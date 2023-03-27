#include "SourceMode.hh"
#include "SessionManager.hh"
#include "G4ParticleGun.hh"
#include "DefinedParticles.hh"
#include "TimeGenerator.hh"
#include "Hist1D.hh"
#include "out.hh"
#include "jstools.hh"
#include "SourceParticleListFile.hh"
#include "SourcePesHistogramFiles.hh"
#include "SourceAnnihilHistFile.hh"

#include "G4RandomTools.hh"
#include "G4NistManager.hh"
#include "G4Electron.hh"

#define _USE_MATH_DEFINES
#include <cmath>

SourceModeBase * SourceModeFactory::makeSourceInstance(const json11::Json & json)
{
    out("Reading source json");
    std::string Type;
    jstools::readString(json, "Type", Type);

    SourceModeBase * sc = nullptr;

    if      (Type == "SourcePoint")              sc = new SourcePoint(json);
    else if (Type == "SourcePointBlurred")       sc = new SourcePointBlurred(json);
    else if (Type == "SourceNa22Point")          sc = new SourceNa22Point(json);
    else if (Type == "SourceLine")               sc = new SourceLine(json);
    else if (Type == "SourceCylinder")           sc = new SourceCylinder(json);
    else if (Type == "SourceBeam")               sc = new SourceBeam(json);
    else if (Type == "SourceMultiBeam")          sc = new SourceMultiBeam(json);
    else if (Type == "SourceMaterialLimited")    sc = new SourceMaterialLimited(json);
    else if (Type == "SourceLysoNatural")        sc = new SourceLysoNatural(json);
    else if (Type == "SourceParticleListFile")   sc = new SourceParticleListFile(json);
    else if (Type == "SourcePesHistogramFiles")  sc = new SourcePesHistogramFiles(json);
    else if (Type == "SourceAnnihilHistFile")    sc = new SourceAnnihilHistFile(json);
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

void SourceModeBase::setParticleEnergy(double energy)
{
    ParticleGun->SetParticleEnergy(energy);
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

SourcePoint::SourcePoint(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin) :
    SourceModeBase(particle, timeGenerator), Origin(origin)
{
    ParticleGun->SetParticlePosition(Origin);
}

SourcePoint::SourcePoint(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);

    ParticleGun->SetParticlePosition(Origin);
}

void SourcePoint::doWriteToJson(json11::Json::object & json) const
{
    json["OriginX"] = Origin.x();
    json["OriginY"] = Origin.y();
    json["OriginZ"] = Origin.z();
}

void SourcePoint::doReadFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "OriginX", Origin[0]);
    jstools::readDouble(json, "OriginY", Origin[1]);
    jstools::readDouble(json, "OriginZ", Origin[2]);
}

// ---

SourceBeam::SourceBeam(ParticleBase * particle, TimeGeneratorBase * timeGenerator,
                       const G4ThreeVector & origin, const G4ThreeVector & direction,
                       int numParticles, ProfileBase * spread) :
    SourceModeBase(particle, timeGenerator), Origin(origin), NumParticles(numParticles), Profile(spread)
{
    Direction = direction;
    update();
}

SourceBeam::SourceBeam(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);
    update();
}

SourceBeam::~SourceBeam()
{
    delete Profile;
}

void SourceBeam::GeneratePrimaries(G4Event * anEvent)
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

void SourceBeam::update()
{
    if (Profile) Profile->setDirection(Direction);
    ParticleGun->SetParticlePosition(Origin);
    bIsotropicDirection = false;
}

void SourceBeam::doWriteToJson(json11::Json::object & json) const
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

void SourceBeam::doReadFromJson(const json11::Json &json)
{
    jstools::readDouble(json, "OriginX", Origin[0]);
    jstools::readDouble(json, "OriginY", Origin[1]);
    jstools::readDouble(json, "OriginZ", Origin[2]);

    jstools::readDouble(json, "DirectionX", Direction[0]);
    jstools::readDouble(json, "DirectionY", Direction[1]);
    jstools::readDouble(json, "DirectionZ", Direction[2]);

    jstools::readInt(json, "NumParticles", NumParticles);

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

void BeamRecord::writeToJson(json11::Json::object & json) const
{
    json["Energy"]     = Energy;
    json["XIsoCenter"] = XIsoCenter;
    json["ZIsoCenter"] = ZIsoCenter;
    json["TimeStart"]  = TimeStart;
    json["TimeSpan"]   = TimeSpan;
    json["StatWeight"] = StatWeight;
}

void BeamRecord::readFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "Energy",     Energy);
    jstools::readDouble(json, "XIsoCenter", XIsoCenter);
    jstools::readDouble(json, "ZIsoCenter", ZIsoCenter);
    jstools::readDouble(json, "TimeStart",  TimeStart);
    jstools::readDouble(json, "TimeSpan",   TimeSpan);
    jstools::readDouble(json, "StatWeight", StatWeight);
}

void BeamRecord::print() const
{
    out("NominalEnergy", Energy/MeV, "MeV; Xiso", XIsoCenter/mm, "mm; Ziso", ZIsoCenter/mm, "mm; Time0", TimeStart/ns, "ns; Tspan", TimeSpan/ns, "ns; StatWeight", StatWeight);
}

SourceMultiBeam::SourceMultiBeam(ParticleBase * particle, const std::vector<BeamRecord> & beams, double totalParticles) :
    SourceModeBase(particle, nullptr), Beams(beams), NumParticles(totalParticles)
{
    loadCalibration();
    calculateParticlesPerStatWeightUnit();
}

SourceMultiBeam::SourceMultiBeam(ParticleBase * particle, const std::string & beamletFileName, double totalParticles) :
    SourceModeBase(particle, nullptr), NumParticles(totalParticles)
{
    loadCalibration();
    loadBeamletData(beamletFileName);
    calculateParticlesPerStatWeightUnit();

    //out(Beams.size());
    //for (const auto & r : Beams) r.print();
}

SourceMultiBeam::SourceMultiBeam(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    readFromJson(json);
    SourceMultiBeam::doReadFromJson(json);
    loadCalibration();
    calculateParticlesPerStatWeightUnit();

    //out(Beams.size());
    //for (const auto & r : Beams) r.print();
}

void SourceMultiBeam::loadCalibration()
{
    std::ifstream inStream(CalibrationFileName); // the file is copied automatically by cmake to the run directory
    if (!inStream.is_open())
    {
        out("Cannot open file with calibration:\n", CalibrationFileName);
        exit(1);
    }

    Calibration.clear();
    double NomEnergy, TrueEnergy, SpotSigma; //[MeV, MeV, mm]
    for (std::string line; std::getline(inStream, line); )
    {
        out(">>>",line);
        if (line.empty()) continue; //allow empty lines
        if (line[0] == '#') continue; //allow comments

        std::stringstream ss(line);  // units in the file are MeV and mbarns
        ss >> NomEnergy >> TrueEnergy >> SpotSigma;
        if (ss.fail())
        {
            out("Unexpected format of a line in the calibration file:\nExpect lines with NominalEnergy[MeV] TrueEnergy[MeV] SpotSigma[mm]");
            exit(3);
        }
        //out(Energy, WaterRange);
        Calibration.push_back({NomEnergy*MeV, TrueEnergy*MeV, SpotSigma*mm});
    }
}

void SourceMultiBeam::loadBeamletData(const std::string & beamletDataFile)
{
    std::ifstream inStream(beamletDataFile);
    if (!inStream.is_open())
    {
        out("Cannot open file with beamlet data:\n", beamletDataFile);
        exit(2);
    }

    Beams.clear();
    double NominalEnergy, XIsoCenter, ZIsoCenter, TimeStart, TimeSpan, StatWeight; // [MeV, mm, mm, ns, ns, number]
    for (std::string line; std::getline(inStream, line); )
    {
        //out(">>>",line);
        if (line.empty()) continue; //allow empty lines
        if (line[0] == '#') continue; //allow comments

        std::stringstream ss(line);  // units in the file are MeV and mbarns
        ss >> NominalEnergy >> XIsoCenter >> ZIsoCenter >> TimeStart >> TimeSpan >> StatWeight;
        if (ss.fail())
        {
            out("Unexpected format of a line in the calibration file:\nExpect lines with NominalEnergy[MeV] TrueEnergy[MeV] SpotSigma[mm]");
            exit(3);
        }
        //out(Energy, WaterRange);
        Beams.push_back( {NominalEnergy*MeV, XIsoCenter*mm, ZIsoCenter*mm, TimeStart*ns, TimeSpan*ns, StatWeight} );
    }
}

void SourceMultiBeam::calculateParticlesPerStatWeightUnit()
{
    double sum = 0;
    for (const BeamRecord & br : Beams) sum += br.StatWeight;

    if (sum <= 0)
    {
        out("Sum stat weight of the beamlets should be a positive number!");
        exit(11);
    }

    ParticlesPerStatWeightUnit = NumParticles / sum;
    out("----->Sum of stat weights:", sum, "  Particles:", NumParticles, "  Particles per stat weight unit:", ParticlesPerStatWeightUnit);
}

void SourceMultiBeam::GeneratePrimaries(G4Event * anEvent)  // optimize! do interpolation once for the same energy! !!!***
{
    if (iRecord >= Beams.size())
    {
        out("All defined beamlets were already generated");
        return;
    }

    if ( iParticle >= round(ParticlesPerStatWeightUnit * Beams[iRecord].StatWeight) )
    {
        iRecord++;
        iParticle = 0;
        return GeneratePrimaries(anEvent);
    }

    //out("Generating particle for beam #", iRecord, " and particle #", iParticle);
    const BeamRecord & Beam = Beams[iRecord];

    const double & nominalEnergy = Beam.Energy;
    double iUpper = 1;
    for (; iUpper < Calibration.size(); iUpper++)
        if (Calibration[iUpper][0] > nominalEnergy) break;
    if (iUpper == Calibration.size())
    {
        out("Energy is outise calibration range");
        exit(5);
    }
    const double interpolationfactor = (nominalEnergy - Calibration[iUpper-1][0]) / ( Calibration[iUpper][0] - Calibration[iUpper-1][0] );
    //out("  Nominal energy:", nominalEnergy, " i.f.:", interpolationfactor);

    //energy
    const double trueEnergy = SessionManager::interpolate(Calibration[iUpper-1][1], Calibration[iUpper][1], interpolationfactor);
    ParticleGun->SetParticleEnergy(trueEnergy);
    //out("  Energy:", trueEnergy, " MeV");

    //time
    const double time = Beam.TimeStart + Beam.TimeSpan * G4UniformRand();
    ParticleGun->SetParticleTime(time);
    //out("  Time:", time, "ns");

    //direction
    const double sigma = SessionManager::interpolate(Calibration[iUpper-1][2], Calibration[iUpper][2], interpolationfactor);
    //out("  Sigma:", sigma);
    const double X = G4RandGauss::shoot(Beam.XIsoCenter, sigma);
    const double Z = G4RandGauss::shoot(Beam.ZIsoCenter, sigma);
    //out(Beam.XIsoCenter, Beam.PositionSigma, X); out(Beam.ZIsoCenter, Beam.PositionSigma, Z); exit(111);
    G4ThreeVector dir = G4ThreeVector(X, 0, Z) - Origin;
    dir = dir.unit();
    ParticleGun->SetParticleMomentumDirection(dir);
    //out("  Direction:", dir);

    //position
    const double fraction = (Origin[1] - StartBeamFromY) / Origin[1];
    const G4ThreeVector xyz = {X * fraction, StartBeamFromY, Z * fraction};
    ParticleGun->SetParticlePosition(xyz);
    //out("  Position:", xyz, "  [mm]");

    ParticleGun->GeneratePrimaryVertex(anEvent);
    iParticle++;
}

std::vector<std::pair<double, double>> SourceMultiBeam::getTimeWindows(double delayAfter, double marginBefore) const
{
    std::vector<std::pair<double, double>> wins;

    for (size_t iB = 0; iB < Beams.size() - 1; iB++)
    {
        const double from     = Beams[iB].TimeStart + Beams[iB].TimeSpan + delayAfter;
        const double duration = Beams[iB+1].TimeStart - marginBefore;
        if (duration > 0) wins.push_back({from, duration});
    }

    return wins;
}

void SourceMultiBeam::doWriteToJson(json11::Json::object & json) const
{
    json11::Json::array ar;
    for (const auto & b : Beams)
    {
        json11::Json::object el;
        b.writeToJson(el);
        ar.push_back(el);
    }
    json["Beams"] = ar;
    json["NumParticles"] = NumParticles;
}

void SourceMultiBeam::doReadFromJson(const json11::Json & json)
{
    Beams.clear();
    json11::Json::array ar;
    jstools::readArray(json, "Beams", ar);
    for (size_t i = 0; i < ar.size(); i++)
    {
        json11::Json::object el = ar[i].object_items();
        BeamRecord rec;
        rec.readFromJson(el);
        Beams.push_back(rec);
    }
    jstools::readDouble(json, "NumParticles", NumParticles);
}

// ---

#include "G4Navigator.hh"
SourceMaterialLimited::SourceMaterialLimited(ParticleBase * particle,
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

SourceMaterialLimited::SourceMaterialLimited(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);
    init();
}

void SourceMaterialLimited::init()
{
    if (!FileName.empty())
    {
        Stream = new std::ofstream();
        SessionManager & SM = SessionManager::getInstance();
        Stream->open(SM.WorkingDirectory + '/' + FileName);
        if (!Stream->is_open() || Stream->fail() || Stream->bad())
        {
            out("Cannot open file to store emission positions!");
            delete Stream; Stream = nullptr;
        }
        else out("\nSaving source emission positions to file", FileName);
    }
}

SourceMaterialLimited::~SourceMaterialLimited()
{
    if (Stream) Stream->close();
    delete Stream;    Stream    = nullptr;
    delete Navigator; Navigator = nullptr;
}

void SourceMaterialLimited::customPostInit()
{
    Navigator = new G4Navigator();
    SessionManager & SM = SessionManager::getInstance();
    Navigator->SetWorldVolume(SM.physWorld);

    G4NistManager * man = G4NistManager::Instance();
    SourceMat = man->FindMaterial(Material);
}

void SourceMaterialLimited::GeneratePrimaries(G4Event *anEvent)
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
                    *Stream << pos[0] << " " << pos[1] << " " << pos[2] << '\n';
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

void SourceMaterialLimited::doWriteToJson(json11::Json::object &json) const
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

void SourceMaterialLimited::doReadFromJson(const json11::Json & json)
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

SourceLysoNatural::SourceLysoNatural(double timeFrom, double timeTo) :
    SourceModeBase(new Hf176exc, new UniformTime(timeFrom, timeTo)), TimeFrom(timeFrom), TimeTo(timeTo)
{
    init();
}

SourceLysoNatural::SourceLysoNatural(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);
    init();
}

void SourceLysoNatural::init()
{
    bIsotropicDirection = false;

    SessionManager & SM = SessionManager::getInstance();
    ScintMaxRadius = sqrt( 0.25*(SM.ScintSizeX*SM.ScintSizeX + SM.ScintSizeY*SM.ScintSizeY + SM.ScintSizeZ*SM.ScintSizeZ) );
    //out("---Max scint radius:", ScintMaxRadius);

    ElectronSpectrum = {{0,3881}, {6,3888}, {12,3969}, {18,4009}, {24,3979}, {30,3994}, {36,3860}, {42,3828}, {48,3753}, {54,3698}, {60,3810}, {66,3706}, {72,3735}, {78,3628}, {84,3569}, {90,3641}, {96,3623}, {102,3618}, {108,3593}, {114,3471}, {120,3547}, {126,3406}, {132,3472}, {138,3279}, {144,3365}, {150,3362}, {156,3299}, {162,3142}, {168,3176}, {174,3070}, {180,3070}, {186,3038}, {192,3094}, {198,2947}, {204,2866}, {210,2872}, {216,2897}, {222,2779}, {228,2669}, {234,2700}, {240,2556}, {246,2498}, {252,2529}, {258,2562}, {264,2400}, {270,2218}, {276,2264}, {282,2197}, {288,2155}, {294,2007}, {300,2084}, {306,1946}, {312,1864}, {318,1775}, {324,1784}, {330,1794}, {336,1702}, {342,1625}, {348,1509}, {354,1516}, {360,1476}, {366,1343}, {372,1332}, {378,1236}, {384,1214}, {390,1219}, {396,1122}, {402,1070}, {408,1025}, {414,917}, {420,913}, {426,832}, {432,770}, {438,758}, {444,695}, {450,626}, {456,629}, {462,567}, {468,523}, {474,472}, {480,399}, {486,389}, {492,340}, {498,285}, {504,250}, {510,245}, {516,231}, {522,179}, {528,131}, {534,145}, {540,114}, {546,77}, {552,67}, {558,57}, {564,28}, {570,16}, {576,13}, {582,5}, {588,2}, {594,0}};
    Hist1D hist(100, 0, 600.0);
    for (const auto & pair : ElectronSpectrum)
        hist.fill(pair.first+0.001, pair.second);

    Sampler = new Hist1DSampler(hist);
}

SourceLysoNatural::~SourceLysoNatural()
{
    delete Sampler;
    delete Navigator;
}

void SourceLysoNatural::GeneratePrimaries(G4Event *anEvent)
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

void SourceLysoNatural::doWriteToJson(json11::Json::object & json) const
{
    json["TimeFrom"] = TimeFrom;
    json["TimeTo"]   = TimeTo;
}

void SourceLysoNatural::doReadFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "TimeFrom", TimeFrom);
    jstools::readDouble(json, "TimeTo",   TimeTo);
}

void SourceLysoNatural::customPostInit()
{
    Navigator = new G4Navigator();
    SessionManager & SM = SessionManager::getInstance();
    Navigator->SetWorldVolume(SM.physWorld);
}

SourcePointBlurred::SourcePointBlurred(ParticleBase *particle, TimeGeneratorBase *timeGenerator, const G4ThreeVector &origin, G4String distributionFileName, int numBins, double range) :
    SourcePoint(particle, timeGenerator, origin),
    FileName(distributionFileName), NumBins(numBins), Range(range)
{
    init();
}

SourcePointBlurred::SourcePointBlurred(const json11::Json & json) :
    SourcePoint(json)
{
    doReadFromJson(json);
    init();
}

void SourcePointBlurred::init()
{
    Hist1D dist(NumBins, -Range, Range);
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

    Sampler = new Hist1DSampler(dist);
}

SourcePointBlurred::~SourcePointBlurred()
{
    delete Sampler;
}

void SourcePointBlurred::GeneratePrimaries(G4Event *anEvent)
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

void SourcePointBlurred::doWriteToJson(json11::Json::object &json) const
{
    SourcePoint::doWriteToJson(json);
    json["FileName"] = FileName;
    json["NumBins"]  = NumBins;
    json["Range"]    = Range;
}

void SourcePointBlurred::doReadFromJson(const json11::Json &json)
{
    SourcePoint::doReadFromJson(json);
    jstools::readString(json, "FileName", FileName);
    jstools::readInt(json, "NumBins", NumBins);
    jstools::readDouble(json, "Range", Range);
}

// ---

SourceLine::SourceLine(ParticleBase *particle, TimeGeneratorBase *timeGenerator, const G4ThreeVector &startPoint, const G4ThreeVector &endPoint) :
    SourceModeBase(particle, timeGenerator), StartPoint(startPoint), EndPoint(endPoint) {}

SourceLine::SourceLine(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);
}

void SourceLine::GeneratePrimaries(G4Event * anEvent)
{
    const G4ThreeVector pos = StartPoint + (EndPoint - StartPoint) * G4UniformRand();
    //out(pos);
    ParticleGun->SetParticlePosition(pos);
    SourceModeBase::GeneratePrimaries(anEvent);
}

void SourceLine::doWriteToJson(json11::Json::object & json) const
{
    json["StartPointX"] = StartPoint.x();
    json["StartPointY"] = StartPoint.y();
    json["StartPointZ"] = StartPoint.z();

    json["EndPointX"] = EndPoint.x();
    json["EndPointY"] = EndPoint.y();
    json["EndPointZ"] = EndPoint.z();
}

void SourceLine::doReadFromJson(const json11::Json &json)
{
    jstools::readDouble(json, "StartPointX", StartPoint[0]);
    jstools::readDouble(json, "StartPointY", StartPoint[1]);
    jstools::readDouble(json, "StartPointZ", StartPoint[2]);

    jstools::readDouble(json, "EndPointX", EndPoint[0]);
    jstools::readDouble(json, "EndPointY", EndPoint[1]);
    jstools::readDouble(json, "EndPointZ", EndPoint[2]);
}

// ---

SourceCylinder::SourceCylinder(ParticleBase * particle, TimeGeneratorBase * timeGenerator, double radius, const G4ThreeVector & startPoint, const G4ThreeVector & endPoint, const std::string & fileName) :
    SourceModeBase(particle, timeGenerator), Radius(radius), StartPoint(startPoint), EndPoint(endPoint), FileName(fileName)
{
    init();
}

SourceCylinder::SourceCylinder(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);

    init();
}

SourceCylinder::~SourceCylinder()
{
    if (Stream)
    {
        Stream->close();
        delete Stream;
    }
}

void SourceCylinder::init()
{
    Axis = EndPoint - StartPoint;
    UnitNormal1 = Axis.unit().orthogonal();
    UnitNormal2 = UnitNormal1;
    UnitNormal2.rotate(0.5 * M_PI, Axis);

    if (!FileName.empty())
    {
        Stream = new std::ofstream();
        SessionManager & SM = SessionManager::getInstance();
        Stream->open(SM.WorkingDirectory + '/' + FileName);
        if (!Stream->is_open() || Stream->fail() || Stream->bad())
        {
            out("Cannot open file to store emission positions!");
            delete Stream; Stream = nullptr;
        }
        else out("\nSaving source emission positions to file", FileName);
    }
}

void SourceCylinder::GeneratePrimaries(G4Event * anEvent)
{
    G4ThreeVector pos = StartPoint + Axis * G4UniformRand();

    double x,y;
    do
    {
        x = -1.0 + 2.0 * G4UniformRand();
        y = -1.0 + 2.0 * G4UniformRand();
    }
    while (x*x + y*y > 1.0);

    pos += UnitNormal1 * x * Radius  +  UnitNormal2 * y * Radius;

    if (Stream) *Stream << pos[0] << " " << pos[1] << " " << pos[2] << '\n';

    ParticleGun->SetParticlePosition(pos);
    SourceModeBase::GeneratePrimaries(anEvent);
}

void SourceCylinder::doWriteToJson(json11::Json::object & json) const
{
    json["Radius"] = Radius;
    json["FileName"] = FileName;

    json["StartPointX"] = StartPoint.x();
    json["StartPointY"] = StartPoint.y();
    json["StartPointZ"] = StartPoint.z();

    json["EndPointX"] = EndPoint.x();
    json["EndPointY"] = EndPoint.y();
    json["EndPointZ"] = EndPoint.z();
}

void SourceCylinder::doReadFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "Radius", Radius);
    jstools::readString(json, "FileName", FileName);

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

RoundProfile::RoundProfile(const json11::Json & json)
{
    readFromJson(json);
}

void RoundProfile::generateOffset(G4ThreeVector & pos) const
{
    double offX, offY;
    do
    {
        offX = -0.5 * Diameter + Diameter * G4UniformRand();
        offY = -0.5 * Diameter + Diameter * G4UniformRand();
    }
    while (offX*offX + offY*offY > 0.25*Diameter*Diameter);

    G4ThreeVector posLoc(offX, offY, 0);
    posLoc.rotateUz(Direction);

    pos += posLoc;
}

void RoundProfile::doWriteToJson(json11::Json::object &json) const
{
    json["Diameter"] = Diameter;
}

void RoundProfile::readFromJson(const json11::Json &json)
{
    jstools::readDouble(json, "Diameter", Diameter);
}

// ---

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

// ---

void readDistFromFile(std::vector<std::pair<double,double>> & dist, const std::string & fileName)
{
    std::ifstream inStream(fileName);
    std::string line;

    dist.clear();
    while (!inStream.eof())
    {
        getline(inStream, line);
        out(line);
        std::stringstream ss(line);
        double position, probability;
        ss >> position >> probability;
        dist.push_back( {position,probability} );
    }
}

CustomProfile::CustomProfile(std::string fileDistributionX, std::string fileDistributionY) : ProfileBase()
{
    readDistFromFile(DistX, fileDistributionX);
    readDistFromFile(DistY, fileDistributionY);
    init();
}

CustomProfile::CustomProfile(const json11::Json & json)
{
    readFromJson(json);
    init();
}

CustomProfile::~CustomProfile()
{
    if (logStream)
    {
        logStream->close();
        delete logStream;
    }
}

void CustomProfile::init()
{
    Hist1D X(DistX.size(), DistX.front().first, DistX.back().first); XSampler = new Hist1DSampler(X); // sampler copies hist
    Hist1D Y(DistY.size(), DistY.front().first, DistY.back().first); YSampler = new Hist1DSampler(Y);

    if (DoLogPositions) logStream = new std::ofstream("CustomProfileGeneratedPositions.txt");
}

void CustomProfile::generateOffset(G4ThreeVector & pos) const
{
    double offX = XSampler->getRandom();
    double offY = YSampler->getRandom();

    G4ThreeVector posLoc(offX, offY, 0);
    posLoc.rotateUz(Direction);

    pos += posLoc;

    if (logStream) *logStream << pos[0] << ' ' << pos[1] << ' ' << pos[2] << '\n';
}

void writeDistToJson(const std::vector<std::pair<double,double>> & dist, json11::Json::array & ar)
{
    for (const auto & pair : dist)
    {
        json11::Json::array el;
            el.push_back(pair.first);
            el.push_back(pair.second);
        ar.push_back(el);
    }
}

void CustomProfile::doWriteToJson(json11::Json::object & json) const
{
    json11::Json::array arX; writeDistToJson(DistX, arX); json["DistributionX"] = arX;
    json11::Json::array arY; writeDistToJson(DistY, arY); json["DistributionY"] = arY;
}

void readDistFromJson(std::vector<std::pair<double,double>> & dist, const json11::Json::array & ar)
{
    dist.clear();
    for (size_t i = 0; i < ar.size(); i++)
    {
        json11::Json::array el = ar[i].array_items();
        if (el.size() != 2)
        {
            out("Custom distribution json should have lines of two numbers");
            exit(111);
        }
        const double pos = el[0].number_value();
        const double sw  = el[1].number_value();
        dist.push_back( {pos, sw} );
    }
}

void CustomProfile::readFromJson(const json11::Json & json)
{
    json11::Json::array arX; jstools::readArray(json, "DistributionX", arX); readDistFromJson(DistX, arX);
    json11::Json::array arY; jstools::readArray(json, "DistributionY", arY); readDistFromJson(DistY, arY);
}

// ---

SourceNa22Point::SourceNa22Point(double timeFrom, double timeTo, const G4ThreeVector & origin) :
    SourceModeBase(new Gamma(1.275*MeV), new UniformTime(timeFrom, timeTo)),
    TimeFrom(timeFrom), TimeTo(timeTo), Origin(origin)
{
    ParticleGun->SetParticlePosition(Origin);
}

SourceNa22Point::SourceNa22Point(const json11::Json & json) : SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);

    ParticleGun->SetParticlePosition(Origin);
}

void SourceNa22Point::GeneratePrimaries(G4Event * anEvent)
{
    SourceModeBase::GeneratePrimaries(anEvent);   // this will generate gamma with 1.275 MeV

    if (G4UniformRand() < 0.096) return; // no beta+

    // generating gamma pair
    double tmpEnergy = ParticleGun->GetParticleEnergy();
    ParticleGun->SetParticleEnergy(0.511*MeV);

    Direction = generateDirectionIsotropic();
    ParticleGun->SetParticleMomentumDirection(Direction);
    ParticleGun->GeneratePrimaryVertex(anEvent);

    ParticleGun->SetParticleMomentumDirection(-Direction);
    ParticleGun->GeneratePrimaryVertex(anEvent);

    //restoring properties
    ParticleGun->SetParticleEnergy(tmpEnergy);
}

void SourceNa22Point::doWriteToJson(json11::Json::object & json) const
{
    json["TimeFrom"] = TimeFrom;
    json["TimeTo"]   = TimeTo;

    json["OriginX"] = Origin.x();
    json["OriginY"] = Origin.y();
    json["OriginZ"] = Origin.z();
}

void SourceNa22Point::doReadFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "TimeFrom", TimeFrom);
    jstools::readDouble(json, "TimeTo",   TimeTo);

    jstools::readDouble(json, "OriginX", Origin[0]);
    jstools::readDouble(json, "OriginY", Origin[1]);
    jstools::readDouble(json, "OriginZ", Origin[2]);
}

// ---

SourceMixer::SourceMixer(std::vector<std::pair<SourceModeBase *, double>> sourcesAndStatWeights) :
    SourceModeBase(nullptr, nullptr), SourcesAndWeights(sourcesAndStatWeights)
{
    init();
}

SourceMixer::SourceMixer(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    init();
}

void SourceMixer::initialize()
{
    for (auto & p : SourcesAndWeights)
        p.first->initialize();
}

void SourceMixer::setParticleEnergy(double energy)
{
    for (auto & p : SourcesAndWeights)
        p.first->setParticleEnergy(energy);
}

void SourceMixer::init()
{
    SumWeight = 0;
    for (const auto & p : SourcesAndWeights)
        SumWeight += p.second;


    if (SourcesAndWeights.empty())
    {
        out("SourceMixer: no sources were defined");
        exit(21);
    }
    if (SumWeight <= 0)
    {
        out("SourceMixer: sum of stat weights for all sources should be positive");
        exit(21);
    }
}

void SourceMixer::GeneratePrimaries(G4Event *anEvent)
{
    double rnd = SumWeight * G4UniformRand();

    for (const auto & source : SourcesAndWeights)
    {
        const double & statWeight = source.second;
        if (rnd <= statWeight)
        {
            source.first->GeneratePrimaries(anEvent);
            return;
        }
        rnd -= statWeight;
    }
}

void SourceMixer::doWriteToJson(json11::Json::object & json) const
{
    json11::Json::array ar;

    for (const auto & p : SourcesAndWeights)
    {
        json11::Json::object js;

            json11::Json::object sjs;
            p.first->writeToJson(sjs);
        js["Source"] = sjs;
        js["StatWeight"] = p.second;

        ar.push_back(js);
    }

    json["Sources"] = ar;
}

void SourceMixer::doReadFromJson(const json11::Json &json)
{
    SourcesAndWeights.clear();

    json11::Json::array ar;
    jstools::readArray(json, "Sources", ar);
    for (size_t i = 0; i < ar.size(); i++)
    {
        json11::Json::object js = ar[i].object_items();

        double statWeight;
        jstools::readDouble(js, "StatWeight", statWeight);

        json11::Json::object sjs;
        jstools::readObject(js, "Source", sjs);
        SourceModeBase * source = SourceModeFactory::makeSourceInstance(sjs);

        SourcesAndWeights.push_back({source, statWeight});
    }
}
