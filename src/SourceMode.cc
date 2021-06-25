#include "SourceMode.hh"
#include "SessionManager.hh"
#include "G4ParticleGun.hh"
#include "DefinedParticles.hh"
#include "TimeGenerator.hh"
#include "Hist1D.hh"
#include "out.hh"

#include "G4RandomTools.hh"
#include "G4NistManager.hh"
#include "G4Electron.hh"

#define _USE_MATH_DEFINES
#include <cmath>

SourceModeBase::SourceModeBase(ParticleBase * particle, TimeGeneratorBase * timeGenerator) :
    Particle(particle), TimeGenerator(timeGenerator)
{
    ParticleGun = new G4ParticleGun(1);

    //Warning: particle definition can be set only later when physics list is initialized. see initialize() method called by SessionManager

    if (Particle)
    {
        bIsotropicDirection = !Particle->bSkipDirection; // can be overwriten by the concrete source type!

        GammaPair * pair = dynamic_cast<GammaPair*>(Particle);
        bGeneratePair = pair;
        //if (pair) bAcollinearity = pair->bAcollineraity;

        ParticleGun->SetParticleEnergy(Particle->Energy); // to be changed later if there will be spectra to be sampled from
    }
}

SourceModeBase::~SourceModeBase()
{
    delete ParticleGun;   ParticleGun   = nullptr;
    delete TimeGenerator; TimeGenerator = nullptr;
    delete Particle;      Particle      = nullptr;
}

void SourceModeBase::initialize()
{
    if (Particle) ParticleGun->SetParticleDefinition(Particle->getParticleDefinition());
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

void SourceModeBase::generateSecondGamma(G4Event * anEvent)
{
    /*
    if (bAcollinearity)
    {
        G4ThreeVector v = -Direction;
        constexpr double Sigma = 0.5*deg / 2.35482;
        double angle = G4RandGauss::shoot(0, Sigma);
        v.rotate(angle, v.orthogonal());
        v.rotate(G4UniformRand()*2.0*M_PI, Direction);
        ParticleGun->SetParticleMomentumDirection(v);
    }
    else
    */
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
    SourceModeBase(particle, timeGenerator)
{
    ParticleGun->SetParticlePosition(origin);
}

// ---

PencilBeam::PencilBeam(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin, const G4ThreeVector & direction) :
    SourceModeBase(particle, timeGenerator)
{
    ParticleGun->SetParticlePosition(origin);

    Direction = direction;
    bIsotropicDirection = false;
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

NaturalLysoSource::NaturalLysoSource(double timeFrom, double timeTo) :
    SourceModeBase(new Lu176, new UniformTime(timeFrom, timeTo))
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

void NaturalLysoSource::customPostInit()
{
    Navigator = new G4Navigator();
    SessionManager & SM = SessionManager::getInstance();
    Navigator->SetWorldVolume(SM.physWorld);
}

BlurredPointSource::BlurredPointSource(ParticleBase *particle, TimeGeneratorBase *timeGenerator, const G4ThreeVector &origin, G4String fileName) :
    PointSource(particle, timeGenerator, origin)
{
    Hist1D dist(21, -10, 10);
    std::ifstream * inStream = new std::ifstream(fileName);
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
    Origin = ParticleGun->GetParticlePosition();
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

// ---

LineSource::LineSource(ParticleBase *particle, TimeGeneratorBase *timeGenerator, const G4ThreeVector &startPoint, const G4ThreeVector &endPoint) :
    SourceModeBase(particle, timeGenerator), StartPoint(startPoint), EndPoint(endPoint) {}

void LineSource::GeneratePrimaries(G4Event *anEvent)
{
    double rand = G4UniformRand();
    G4ThreeVector pos;
    for (int i=0; i<3; i++)
        pos[i] = StartPoint[i] + (EndPoint[i] - StartPoint[i]) * rand;
    //out(pos);
    ParticleGun->SetParticlePosition(pos);
    SourceModeBase::GeneratePrimaries(anEvent);
}
