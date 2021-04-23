#include "SourceMode.hh"
#include "SessionManager.hh"
#include "G4ParticleGun.hh"
#include "TimeGenerator.hh"
#include "out.hh"

#include "G4RandomTools.hh"
#include "G4NistManager.hh"

#define _USE_MATH_DEFINES
#include <cmath>

SourceModeBase::SourceModeBase(ParticleBase * particle, TimeGeneratorBase * timeGenerator) :
    Particle(particle), TimeGenerator(timeGenerator)
{
    ParticleGun = new G4ParticleGun(1);
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
    //Warning: particle definition can be set only later when physics list is initialized. see initialize() method

    ParticleGun->SetParticlePosition(origin);

    bSkipDirection = particle->bSkipDirection;
    if (bSkipDirection) ParticleGun->SetParticleMomentumDirection({0,0,1.0});

    ParticleGun->SetParticleEnergy(Particle->Energy);

    bGeneratePair = (bool)dynamic_cast<ParticleGammaPair*>(particle);
}

void PointSource::GeneratePrimaries(G4Event * anEvent)
{
    ParticleGun->SetParticleTime(TimeGenerator->generateTime());

    if (bSkipDirection)
        ParticleGun->GeneratePrimaryVertex(anEvent);
    else
    {
        G4ThreeVector v = generateDirectionIsotropic();
        ParticleGun->SetParticleMomentumDirection(v);

        ParticleGun->GeneratePrimaryVertex(anEvent);

        if (bGeneratePair)
        {
            ParticleGun->SetParticleMomentumDirection(-v);//{-x,-y,-z});
            ParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }
}

// ---

PencilBeam::PencilBeam(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin, const G4ThreeVector & direction) :
    SourceModeBase(particle, timeGenerator)
{
    //Warning: particle definition can be set only later when physics list is initialized
    ParticleGun->SetParticlePosition(origin);
    ParticleGun->SetParticleMomentumDirection(direction);
    ParticleGun->SetParticleEnergy(Particle->Energy);
}

void PencilBeam::GeneratePrimaries(G4Event * anEvent)
{
    ParticleGun->SetParticleTime(TimeGenerator->generateTime());
    ParticleGun->GeneratePrimaryVertex(anEvent);
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
    ParticleGun->SetParticleEnergy(Particle->Energy);

    bSkipDirection = particle->bSkipDirection;
    if (bSkipDirection) ParticleGun->SetParticleMomentumDirection({0,0,1.0});

    ParticleGun->SetParticleEnergy(Particle->Energy);

    bGeneratePair = (bool)dynamic_cast<ParticleGammaPair*>(particle);

    if (!FileName.empty())
    {
        Stream = new std::ofstream();
        Stream->open(FileName);
        if (!Stream->is_open())
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
    delete Stream;
    delete Navigator;
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
    ParticleGun->SetParticleTime(TimeGenerator->generateTime());

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

        if (attempts > 1000)
        {
            out("Made 1000 attempts to generate a position within the source material, but failed!");
            exit(33);
        }
    }

    if (bSkipDirection)
        ParticleGun->GeneratePrimaryVertex(anEvent);
    else
    {
        G4ThreeVector v = generateDirectionIsotropic();
        ParticleGun->SetParticleMomentumDirection(v);
        ParticleGun->GeneratePrimaryVertex(anEvent);

        if (bGeneratePair)
        {
            ParticleGun->SetParticleMomentumDirection(-v);
            ParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }
}
