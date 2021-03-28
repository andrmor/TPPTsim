#include "PrimaryGeneratorAction.hh"
#include "SessionManager.hh"
#include "SimMode.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4RandomTools.hh"

#define _USE_MATH_DEFINES
#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
    SessionManager & SM = SessionManager::getInstance();
    fParticleGun = new G4ParticleGun(SM.NumParticlesPerEvent);

    G4ParticleDefinition* particleDefinition = nullptr;
    double Energy = 0;

    switch (SM.SimMode->SourceMode)
    {
        default:;

        case SourceModeEnum::GammaPair :
        {
            particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
            Energy = 511.0*keV;
        }
        break;

        /*
        case SourceModeEnum::C10 :
        {
            particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
            Energy = 0;
            fParticleGun->SetParticleMomentumDirection({0,0,1.0});
        }
        break;
        */

        /*more here*/
    }

    fParticleGun->SetParticleDefinition(particleDefinition);
    fParticleGun->SetParticleEnergy(Energy);
    fParticleGun->SetParticlePosition({0,0,0});
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event * anEvent)
{
    SessionManager & SM = SessionManager::getInstance();

    if (SM.SimMode->SourceMode == SourceModeEnum::GammaPair)
    {
        double phi = acos(-1.0 + 2.0 * G4UniformRand());
        double theta = 2.0 * M_PI * G4UniformRand();
        double x, y, z;

        x = sin(phi) * cos(theta);
        y = sin(phi) * sin(theta);
        z = cos(phi);

        fParticleGun->SetParticleMomentumDirection({x,y,z});
        fParticleGun->GeneratePrimaryVertex(anEvent);
        fParticleGun->SetParticleMomentumDirection({-x,-y,-z});
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else
    {
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}
