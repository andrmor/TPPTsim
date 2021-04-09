#include "SourceMode.hh"
#include "SessionManager.hh"
#include "G4ParticleGun.hh"
#include "out.hh"

#include "G4RandomTools.hh"
#include "G4NistManager.hh"

#define _USE_MATH_DEFINES
#include <cmath>

SourceModeBase::SourceModeBase(ParticleBase * particle, int numPerEvent) :
    Particle(particle), NumPerEvent(numPerEvent)
{
    ParticleGun = new G4ParticleGun(numPerEvent);
}

SourceModeBase::~SourceModeBase()
{
    delete ParticleGun; ParticleGun = nullptr;
    delete Particle;    Particle    = nullptr;
}

void SourceModeBase::initialize()
{
    ParticleGun->SetParticleDefinition(Particle->getParticleDefinition());

    customPostInit();
}

// ---

PointSource::PointSource(ParticleBase * particle, const G4ThreeVector & origin, int numPerEvent) :
    SourceModeBase(particle, numPerEvent)
{
    //Warning: particle definition can be set only later when physics list is initialized

    ParticleGun->SetParticlePosition(origin);

    bSkipDirection = particle->bSkipDirection;
    if (bSkipDirection) ParticleGun->SetParticleMomentumDirection({0,0,1.0});

    ParticleGun->SetParticleEnergy(Particle->Energy);

    bGeneratePair = (bool)dynamic_cast<ParticleGammaPair*>(particle);
}

void PointSource::GeneratePrimaries(G4Event * anEvent)
{
    if (bSkipDirection)
        ParticleGun->GeneratePrimaryVertex(anEvent);
    else
    {
        double phi   = acos(-1.0 + 2.0 * G4UniformRand());
        double theta = 2.0 * M_PI * G4UniformRand();

        double x = sin(phi) * cos(theta);
        double y = sin(phi) * sin(theta);
        double z = cos(phi);

        ParticleGun->SetParticleMomentumDirection({x,y,z});
        ParticleGun->GeneratePrimaryVertex(anEvent);

        if (bGeneratePair)
        {
            ParticleGun->SetParticleMomentumDirection({-x,-y,-z});
            ParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }
}

// ---

PointSourceUniformTime::PointSourceUniformTime(ParticleBase *particle, const G4ThreeVector &origin, int numPerEvent, double timeWindow) :
    PointSource(particle, origin, numPerEvent), TimeWindow(timeWindow) {}

void PointSourceUniformTime::GeneratePrimaries(G4Event *anEvent)
{
    ParticleGun->SetParticleTime(G4UniformRand() * TimeWindow);
    PointSource::GeneratePrimaries(anEvent);
}

// ---

PencilBeam::PencilBeam(ParticleBase * particle, const G4ThreeVector & origin, const G4ThreeVector & direction, int numPerEvent) :
    SourceModeBase(particle, numPerEvent)
{
    //Warning: particle definition can be set only later when physics list is initialized
    ParticleGun->SetParticlePosition(origin);
    ParticleGun->SetParticleMomentumDirection(direction);
    ParticleGun->SetParticleEnergy(Particle->Energy);
}

void PencilBeam::GeneratePrimaries(G4Event * anEvent)
{
    ParticleGun->GeneratePrimaryVertex(anEvent);
}

// ---

#include "G4Navigator.hh"
MaterialLimitedSource::MaterialLimitedSource(ParticleBase *particle, const G4ThreeVector &origin, const G4ThreeVector &boundingBoxFullSize, const G4String &material, int numPerEvent) :
    SourceModeBase(particle, numPerEvent), Origin(origin), BoundingBox(boundingBoxFullSize), Material(material)
{
    ParticleGun->SetParticleEnergy(Particle->Energy);

    bSkipDirection = particle->bSkipDirection;
    if (bSkipDirection) ParticleGun->SetParticleMomentumDirection({0,0,1.0});

    ParticleGun->SetParticleEnergy(Particle->Energy);

    bGeneratePair = (bool)dynamic_cast<ParticleGammaPair*>(particle);
}

MaterialLimitedSource::~MaterialLimitedSource()
{
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
        double phi   = acos(-1.0 + 2.0 * G4UniformRand());
        double theta = 2.0 * M_PI * G4UniformRand();

        double x = sin(phi) * cos(theta);
        double y = sin(phi) * sin(theta);
        double z = cos(phi);

        ParticleGun->SetParticleMomentumDirection({x,y,z});
        ParticleGun->GeneratePrimaryVertex(anEvent);

        if (bGeneratePair)
        {
            ParticleGun->SetParticleMomentumDirection({-x,-y,-z});
            ParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }
}

