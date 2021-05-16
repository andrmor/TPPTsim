#include "SourceMode.hh"
#include "SessionManager.hh"
#include "G4ParticleGun.hh"
#include "DefinedParticles.hh"
#include "TimeGenerator.hh"
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

    bIsotropicDirection = !particle->bSkipDirection; // can be overwriten by the concrete source type!

    GammaPair * pair = dynamic_cast<GammaPair*>(particle);
    bGeneratePair = pair;
    //if (pair) bAcollinearity = pair->bAcollineraity;

    ParticleGun->SetParticleEnergy(Particle->Energy); // to be changed later if there will be spectra to be sampled from
}

SourceModeBase::~SourceModeBase()
{
    delete ParticleGun;   ParticleGun   = nullptr;
    delete TimeGenerator; TimeGenerator = nullptr;
    delete Particle;      Particle      = nullptr;
}

void SourceModeBase::initialize()
{
    ParticleGun->SetParticleDefinition(Particle->getParticleDefinition());
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
}

void NaturalLysoSource::GeneratePrimaries(G4Event *anEvent)
{
    G4ThreeVector pos(0,0,0);
    ParticleGun->SetParticlePosition(pos);

    SourceModeBase::GeneratePrimaries(anEvent);   // this will generate Hf176[596.820]

    // have to leave the gun properties ready to fire the next event!
    G4ParticleDefinition * tmpPD     = ParticleGun->GetParticleDefinition();
    double                 tmpEnergy = ParticleGun->GetParticleEnergy();
    const G4ThreeVector    tmpMdir   = ParticleGun->GetParticleMomentumDirection();

    ParticleGun->SetParticleDefinition(G4Electron::Definition());
    ParticleGun->SetParticleEnergy(100.0*keV);
    ParticleGun->SetParticleMomentumDirection(generateDirectionIsotropic());
    ParticleGun->GeneratePrimaryVertex(anEvent);

    //restoring properties
    ParticleGun->SetParticleMomentumDirection(tmpMdir);
    ParticleGun->SetParticleEnergy(tmpEnergy);
    ParticleGun->SetParticleDefinition(tmpPD);
}
