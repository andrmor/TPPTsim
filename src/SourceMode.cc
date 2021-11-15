#include "SourceMode.hh"
#include "SessionManager.hh"
#include "G4ParticleGun.hh"
#include "DefinedParticles.hh"
#include "TimeGenerator.hh"
#include "Hist1D.hh"
#include "out.hh"
#include "jstools.hh"
#include "FromFileSource.hh"

#include "G4RandomTools.hh"
#include "G4NistManager.hh"
#include "G4Electron.hh"

#define _USE_MATH_DEFINES
#include <cmath>

SourceModeBase * SourceModeFactory::makeSourceModeInstance(const json11::Json & json)
{
    out("Reading source json");
    std::string Type;
    jstools::readString(json, "Type", Type);

    SourceModeBase * sc = nullptr;

    if      (Type == "PointSource")           sc = new PointSource(json);
    else if (Type == "LineSource")            sc = new LineSource(json);
    else if (Type == "PencilBeam")            sc = new PencilBeam(json);
    else if (Type == "MaterialLimitedSource") sc = new MaterialLimitedSource(json);
    else if (Type == "NaturalLysoSource")     sc = new NaturalLysoSource(json);
    else if (Type == "BlurredPointSource")    sc = new BlurredPointSource(json);
    else if (Type == "FromFileSource")        sc = new FromFileSource(json);
    else
    {
        out("Unknown source type!");
        exit(10);
    }

    return sc;
}

// ---

SourceModeBase::SourceModeBase(ParticleBase * particle, TimeGeneratorBase * timeGenerator) :
    Particle(particle), TimeGenerator(timeGenerator)
{
    ParticleGun = new G4ParticleGun(1);
    //Warning: particle definition can be set only later when physics list is initialized. see initialize() method called by SessionManager
    // the rest of initialization (constructor with json) can also be done there
}

SourceModeBase::~SourceModeBase()
{
    delete ParticleGun;   ParticleGun   = nullptr;
    delete TimeGenerator; TimeGenerator = nullptr;
    delete Particle;      Particle      = nullptr;
}

void SourceModeBase::initialize()
{
    if (Particle)
    {
        ParticleGun->SetParticleDefinition(Particle->getParticleDefinition());

        GammaPair * pair = dynamic_cast<GammaPair*>(Particle);
        bGeneratePair = (bool)pair;

        Isotope * iso = dynamic_cast<Isotope*>(Particle);
        if (iso) bIsotropicDirection = false; // save time by not generating direction

        ParticleGun->SetParticleEnergy(Particle->Energy); // to be changed later if there will be spectra to be sampled from
    }
    customPostInit();
}

void SourceModeBase::GeneratePrimaries(G4Event * anEvent)
{
    ParticleGun->SetParticleTime(TimeGenerator->generateTime());

    if (bIsotropicDirection) Direction = generateDirectionIsotropic(); //else it is fixed
    ParticleGun->SetParticleMomentumDirection(Direction);

    ParticleGun->GeneratePrimaryVertex(anEvent);

    if (bGeneratePair) generateSecondGamma(anEvent);
}

void SourceModeBase::writeToJson(json11::Json::object & json) const
{
    json["Type"] = getTypeName();

    if (Particle)
    {
        json11::Json::object js;
            Particle->writeToJson(js);
        json["Particle"] = js;
    }

    if (TimeGenerator)
    {
        json11::Json::object js;
            TimeGenerator->writeToJson(js);
        json["TimeGenerator"] = js;
    }

    doWriteToJson(json);
}

void SourceModeBase::readFromJson(const json11::Json & json)
{
    if (jstools::contains(json, "Particle"))
    {
        json11::Json::object js;
        jstools::readObject(json, "Particle", js);
        Particle = ParticleFactory::makeParticleInstance(js);
    }

    if (jstools::contains(json, "TimeGenerator"))
    {
        json11::Json::object js;
        jstools::readObject(json, "TimeGenerator", js);
        TimeGenerator = TimeGeneratorFactory::makeTimeGeneratorInstance(js);
    }

    //do not call particular loader, using this function in the constructor!
}

void SourceModeBase::generateSecondGamma(G4Event * anEvent)
{
    ParticleGun->SetParticleMomentumDirection(-Direction);

    ParticleGun->GeneratePrimaryVertex(anEvent);
}

G4ThreeVector SourceModeBase::generateDirectionIsotropic()
{
    //Sphere function of CERN ROOT

    double a = 0, b = 0, r2 = 1.0;
    while (r2 > 0.25)
    {
        a  = G4UniformRand() - 0.5;
        b  = G4UniformRand() - 0.5;
        r2 = a*a + b*b;
    }
    double scale = 8.0 * sqrt(0.25 - r2);

    G4ThreeVector v;
    v[0] = a * scale;
    v[1] = b * scale;
    v[2] = ( -1.0 + 8.0 * r2 );
    return v;
}

// ---

PointSource::PointSource(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin) :
    SourceModeBase(particle, timeGenerator), Origin(origin)
{
    ParticleGun->SetParticlePosition(Origin);
}

PointSource::PointSource(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);

    ParticleGun->SetParticlePosition(Origin);
}

void PointSource::doWriteToJson(json11::Json::object & json) const
{
    json["OriginX"] = Origin.x();
    json["OriginY"] = Origin.y();
    json["OriginZ"] = Origin.z();
}

void PointSource::doReadFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "OriginX", Origin[0]);
    jstools::readDouble(json, "OriginY", Origin[1]);
    jstools::readDouble(json, "OriginZ", Origin[2]);
}

// ---

PencilBeam::PencilBeam(ParticleBase * particle, TimeGeneratorBase * timeGenerator,
                       const G4ThreeVector & origin, const G4ThreeVector & direction,
                       int numParticles, ProfileBase * spread) :
    SourceModeBase(particle, timeGenerator), Origin(origin), NumParticles(numParticles), Profile(spread)
{
    Direction = direction;
    update();
}

PencilBeam::PencilBeam(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);
    update();
}

PencilBeam::~PencilBeam()
{
    delete Profile;
}

void PencilBeam::GeneratePrimaries(G4Event * anEvent)
{
    for (int iP = 0; iP < NumParticles; iP++)
    {
        if (Profile)
        {
            G4ThreeVector pos(Origin);
            Profile->generateOffset(pos);
            //out(pos);
            ParticleGun->SetParticlePosition(pos);
        }
        SourceModeBase::GeneratePrimaries(anEvent);
    }
}

void PencilBeam::update()
{
    if (Profile) Profile->setDirection(Direction);
    ParticleGun->SetParticlePosition(Origin);
    bIsotropicDirection = false;
}

void PencilBeam::doWriteToJson(json11::Json::object & json) const
{
    json["OriginX"] = Origin.x();
    json["OriginY"] = Origin.y();
    json["OriginZ"] = Origin.z();

    json["DirectionX"] = Direction.x();
    json["DirectionY"] = Direction.y();
    json["DirectionZ"] = Direction.z();

    json["NumParticles"] = NumParticles;

    json11::Json::object js;
    if (Profile) Profile->writeToJson(js);
    json["Profile"] = js;
}

void PencilBeam::doReadFromJson(const json11::Json &json)
{
    jstools::readDouble(json, "OriginX", Origin[0]);
    jstools::readDouble(json, "OriginY", Origin[1]);
    jstools::readDouble(json, "OriginZ", Origin[2]);

    jstools::readDouble(json, "DirectionX", Direction[0]);
    jstools::readDouble(json, "DirectionY", Direction[1]);
    jstools::readDouble(json, "DirectionZ", Direction[2]);

    //NumParticles = 1;
    jstools::readInt(json, "NumParticles", NumParticles);

    // !!!*** needs refactoring!
    delete Profile; Profile = nullptr;
    json11::Json::object js;
    jstools::readObject(json, "Profile", js);
    if (!js.empty())
    {
        std::string Type;
        jstools::readString(js, "Type", Type);

        if      (Type == "Uniform") Profile = new UniformProfile(js);
        else if (Type == "Gauss")   Profile = new GaussProfile(js);
        else
        {
            out("Unknown profile type for the beam source!");
            exit(1);
        }
    }
}

// ---

#include "G4Navigator.hh"
MaterialLimitedSource::MaterialLimitedSource(ParticleBase * particle,
                                             TimeGeneratorBase * timeGenerator,
                                             const G4ThreeVector & origin, const G4ThreeVector &boundingBoxFullSize,
                                             const G4String & material,
                                             G4String fileName_EmissionPositions) :
    SourceModeBase(particle, timeGenerator),
    Origin(origin), BoundingBox(boundingBoxFullSize),
    Material(material),
    FileName(fileName_EmissionPositions)
{
    init();
}

MaterialLimitedSource::MaterialLimitedSource(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);
    init();
}

void MaterialLimitedSource::init()
{
    if (!FileName.empty())
    {
        Stream = new std::ofstream();
        Stream->open(FileName);
        if (!Stream->is_open() || Stream->fail() || Stream->bad())
        {
            out("Cannot open file to store emission positions!");
            delete Stream; Stream = nullptr;
        }
        else out("\nSaving source emission positions to file", FileName);
    }
}

MaterialLimitedSource::~MaterialLimitedSource()
{
    if (Stream) Stream->close();
    delete Stream;    Stream    = nullptr;
    delete Navigator; Navigator = nullptr;
}

void MaterialLimitedSource::customPostInit()
{
    Navigator = new G4Navigator();
    SessionManager & SM = SessionManager::getInstance();
    Navigator->SetWorldVolume(SM.physWorld);

    G4NistManager * man = G4NistManager::Instance();
    SourceMat = man->FindMaterial(Material);
}

void MaterialLimitedSource::GeneratePrimaries(G4Event *anEvent)
{
    G4ThreeVector pos;
    int attempts = 0;
    while (true)
    {
        for (int i = 0 ; i < 3; i++)
            pos[i] = Origin[i] + (-0.5 + G4UniformRand()) * BoundingBox[i];

        G4VPhysicalVolume * vol = Navigator->LocateGlobalPointAndSetup(pos);
        if (vol && vol->GetLogicalVolume())
            if (vol->GetLogicalVolume()->GetMaterial() == SourceMat)
            {
                ParticleGun->SetParticlePosition(pos);
                if (Stream)
                    *Stream << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
                break;
            }

        if (attempts++ > 10000)
        {
            out("Made 10000 attempts to generate a position within the source material, but failed!");
            exit(33);
        }
    }

    SourceModeBase::GeneratePrimaries(anEvent);
}

void MaterialLimitedSource::doWriteToJson(json11::Json::object &json) const
{
    json["OriginX"] = Origin.x();
    json["OriginY"] = Origin.y();
    json["OriginZ"] = Origin.z();

    json["BoundingBoxX"] = BoundingBox.x();
    json["BoundingBoxY"] = BoundingBox.y();
    json["BoundingBoxZ"] = BoundingBox.z();

    json["Material"] = Material;
    json["FileName"] = FileName;
}

void MaterialLimitedSource::doReadFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "OriginX", Origin[0]);
    jstools::readDouble(json, "OriginY", Origin[1]);
    jstools::readDouble(json, "OriginZ", Origin[2]);

    jstools::readDouble(json, "BoundingBoxX", BoundingBox[0]);
    jstools::readDouble(json, "BoundingBoxY", BoundingBox[1]);
    jstools::readDouble(json, "BoundingBoxZ", BoundingBox[2]);

    jstools::readString(json, "Material", Material);
    jstools::readString(json, "FileName", FileName);
}

// ---

NaturalLysoSource::NaturalLysoSource(double timeFrom, double timeTo) :
    SourceModeBase(new Hf176exc, new UniformTime(timeFrom, timeTo)), TimeFrom(timeFrom), TimeTo(timeTo)
{
    init();
}

NaturalLysoSource::NaturalLysoSource(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);
    init();
}

void NaturalLysoSource::init()
{
    bIsotropicDirection = false;

    SessionManager & SM = SessionManager::getInstance();
    ScintMaxRadius = sqrt( 0.25*(SM.ScintSizeX*SM.ScintSizeX + SM.ScintSizeY*SM.ScintSizeY + SM.ScintSizeZ*SM.ScintSizeZ) );
    //out("---Max scint radius:", ScintMaxRadius);

    ElectronSpectrum = {{0,3881}, {6,3888}, {12,3969}, {18,4009}, {24,3979}, {30,3994}, {36,3860}, {42,3828}, {48,3753}, {54,3698}, {60,3810}, {66,3706}, {72,3735}, {78,3628}, {84,3569}, {90,3641}, {96,3623}, {102,3618}, {108,3593}, {114,3471}, {120,3547}, {126,3406}, {132,3472}, {138,3279}, {144,3365}, {150,3362}, {156,3299}, {162,3142}, {168,3176}, {174,3070}, {180,3070}, {186,3038}, {192,3094}, {198,2947}, {204,2866}, {210,2872}, {216,2897}, {222,2779}, {228,2669}, {234,2700}, {240,2556}, {246,2498}, {252,2529}, {258,2562}, {264,2400}, {270,2218}, {276,2264}, {282,2197}, {288,2155}, {294,2007}, {300,2084}, {306,1946}, {312,1864}, {318,1775}, {324,1784}, {330,1794}, {336,1702}, {342,1625}, {348,1509}, {354,1516}, {360,1476}, {366,1343}, {372,1332}, {378,1236}, {384,1214}, {390,1219}, {396,1122}, {402,1070}, {408,1025}, {414,917}, {420,913}, {426,832}, {432,770}, {438,758}, {444,695}, {450,626}, {456,629}, {462,567}, {468,523}, {474,472}, {480,399}, {486,389}, {492,340}, {498,285}, {504,250}, {510,245}, {516,231}, {522,179}, {528,131}, {534,145}, {540,114}, {546,77}, {552,67}, {558,57}, {564,28}, {570,16}, {576,13}, {582,5}, {588,2}, {594,0}};
    Hist1D hist(100, 0, 600.0);
    for (const auto & pair : ElectronSpectrum)
        hist.fill(pair.first+0.001, pair.second);

    Sampler = new Hist1DSampler(hist, 12345);
    /*
    hist.report();
    Hist1D tmpHist(100,0,600.0);
    for (int i=0; i<200000; i++) tmpHist.fill(Sampler->getRandom());
    tmpHist.save("/home/andr/WORK/TPPT/testOutput_ED.txt");
    exit(-1);
    */
}

NaturalLysoSource::~NaturalLysoSource()
{
    delete Sampler;
    delete Navigator;
}

void NaturalLysoSource::GeneratePrimaries(G4Event *anEvent)
{
    SessionManager & SM = SessionManager::getInstance();
    const int numScint = SM.ScintRecords.size();
    if (numScint == 0) return;

    const int iScint = G4UniformRand() * numScint; // check: flat() excludes 1 or not
    //out("---Randomly selected scint#", iScint);

    G4ThreeVector pos;
    int iCopy = -1;
    do
    {
        for (int i = 0; i < 3; i++)
            pos[i] = SM.ScintRecords[iScint].CenterPos[i] + ScintMaxRadius * ( -1.0 + 2.0 * G4UniformRand() );

        G4VPhysicalVolume * vol = Navigator->LocateGlobalPointAndSetup(pos);
        iCopy = vol->GetCopyNo();
    }
    while (iCopy != iScint);
    ParticleGun->SetParticlePosition(pos);

    SourceModeBase::GeneratePrimaries(anEvent);   // this will generate Hf176[596.820]

    // have to leave the gun properties ready to fire the next event!
    G4ParticleDefinition * tmpPD     = ParticleGun->GetParticleDefinition();
    double                 tmpEnergy = ParticleGun->GetParticleEnergy();
    const G4ThreeVector    tmpMdir   = ParticleGun->GetParticleMomentumDirection();

    ParticleGun->SetParticleDefinition(G4Electron::Definition());
    ParticleGun->SetParticleEnergy(Sampler->getRandom()*keV);
    ParticleGun->SetParticleMomentumDirection(generateDirectionIsotropic());
    ParticleGun->GeneratePrimaryVertex(anEvent);

    //restoring properties
    ParticleGun->SetParticleMomentumDirection(tmpMdir);
    ParticleGun->SetParticleEnergy(tmpEnergy);
    ParticleGun->SetParticleDefinition(tmpPD);
}

void NaturalLysoSource::doWriteToJson(json11::Json::object & json) const
{
    json["TimeFrom"] = TimeFrom;
    json["TimeTo"]   = TimeTo;
}

void NaturalLysoSource::doReadFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "TimeFrom", TimeFrom);
    jstools::readDouble(json, "TimeTo",   TimeTo);
}

void NaturalLysoSource::customPostInit()
{
    Navigator = new G4Navigator();
    SessionManager & SM = SessionManager::getInstance();
    Navigator->SetWorldVolume(SM.physWorld);
}

BlurredPointSource::BlurredPointSource(ParticleBase *particle, TimeGeneratorBase *timeGenerator, const G4ThreeVector &origin, G4String fileName) :
    PointSource(particle, timeGenerator, origin), FileName(fileName)
{
    init();
}

BlurredPointSource::BlurredPointSource(const json11::Json & json) :
    PointSource(json)
{
    doReadFromJson(json);
    init();
}

void BlurredPointSource::init()
{
    Hist1D dist(21, -10, 10);
    std::ifstream * inStream = new std::ifstream(FileName);
    std::vector<std::pair<double,double>> Input;
    std::string line;

    while (!inStream->eof())
    {
        getline(*inStream, line);
        //out(line);
        std::stringstream ss(line);
        double position, probability;
        ss >> position >> probability;
        Input.push_back({position,probability});
    }
    for (const auto & pair : Input)
        dist.fill(pair.first+0.001, pair.second);

    Sampler = new Hist1DSampler(dist, 12345);
}

BlurredPointSource::~BlurredPointSource()
{
    delete Sampler;
}

void BlurredPointSource::GeneratePrimaries(G4Event *anEvent)
{
    G4ThreeVector vec;
    for (int i=0; i<3; i++)
    {
        vec[i] = Sampler->getRandom() + Origin[i];
    }
    //out(vec);
    ParticleGun->SetParticlePosition(vec);
    SourceModeBase::GeneratePrimaries(anEvent);
}

void BlurredPointSource::doWriteToJson(json11::Json::object &json) const
{
    PointSource::doWriteToJson(json);
    json["FileName"] = FileName;
}

void BlurredPointSource::doReadFromJson(const json11::Json &json)
{
    jstools::readString(json, "FileName", FileName);
}

// ---

LineSource::LineSource(ParticleBase *particle, TimeGeneratorBase *timeGenerator, const G4ThreeVector &startPoint, const G4ThreeVector &endPoint) :
    SourceModeBase(particle, timeGenerator), StartPoint(startPoint), EndPoint(endPoint) {}

LineSource::LineSource(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);
}

void LineSource::GeneratePrimaries(G4Event * anEvent)
{
    double rand = G4UniformRand();
    G4ThreeVector pos;
    for (int i=0; i<3; i++)
        pos[i] = StartPoint[i] + (EndPoint[i] - StartPoint[i]) * rand;
    //out(pos);
    ParticleGun->SetParticlePosition(pos);
    SourceModeBase::GeneratePrimaries(anEvent);
}

void LineSource::doWriteToJson(json11::Json::object & json) const
{
    json["StartPointX"] = StartPoint.x();
    json["StartPointY"] = StartPoint.y();
    json["StartPointZ"] = StartPoint.z();

    json["EndPointX"] = EndPoint.x();
    json["EndPointY"] = EndPoint.y();
    json["EndPointZ"] = EndPoint.z();
}

void LineSource::doReadFromJson(const json11::Json &json)
{
    jstools::readDouble(json, "StartPointX", StartPoint[0]);
    jstools::readDouble(json, "StartPointY", StartPoint[1]);
    jstools::readDouble(json, "StartPointZ", StartPoint[2]);

    jstools::readDouble(json, "EndPointX", EndPoint[0]);
    jstools::readDouble(json, "EndPointY", EndPoint[1]);
    jstools::readDouble(json, "EndPointZ", EndPoint[2]);
}

// ---

void ProfileBase::setDirection(const G4ThreeVector & dir)
{
    Direction = dir.unit();
}

/*
void ProfileBase::addRotation(double Degrees)
{
    Angle = Degrees * M_PI / 180.0;
}
*/

void ProfileBase::writeToJson(json11::Json::object & json) const
{
    json["Type"] = getTypeName();
    doWriteToJson(json);
}

UniformProfile::UniformProfile(const json11::Json & json) : ProfileBase()
{
    readFromJson(json);
}

void UniformProfile::generateOffset(G4ThreeVector & pos) const
{
    double offX = -0.5 * DX + DX * G4UniformRand();
    double offY = -0.5 * DY + DY * G4UniformRand();

    G4ThreeVector posLoc(offX, offY, 0);
    posLoc.rotateUz(Direction);

    pos += posLoc;
}

void UniformProfile::doWriteToJson(json11::Json::object &json) const
{
    json["DX"] = DX;
    json["DY"] = DY;
}

void UniformProfile::readFromJson(const json11::Json &json)
{
    jstools::readDouble(json, "DX", DX);
    jstools::readDouble(json, "DY", DY);
}

GaussProfile::GaussProfile(const json11::Json & json) : ProfileBase()
{
    readFromJson(json);
}

void GaussProfile::generateOffset(G4ThreeVector & pos) const
{
    double offX = G4RandGauss::shoot(0, SigmaX);
    double offY = G4RandGauss::shoot(0, SigmaY);

    G4ThreeVector posLoc(offX, offY, 0);
    posLoc.rotateUz(Direction);

    pos += posLoc;
}

void GaussProfile::doWriteToJson(json11::Json::object &json) const
{
    json["SigmaX"] = SigmaX;
    json["SigmaY"] = SigmaY;
}

void GaussProfile::readFromJson(const json11::Json &json)
{
    jstools::readDouble(json, "SigmaX", SigmaX);
    jstools::readDouble(json, "SigmaY", SigmaY);
}
