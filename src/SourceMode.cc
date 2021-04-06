#include "SourceMode.hh"
#include "SessionManager.hh"
#include "G4ParticleGun.hh"
#include "out.hh"

#include "G4RandomTools.hh"

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
