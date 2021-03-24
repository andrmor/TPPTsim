#include "PrimaryGeneratorAction.hh"
#include "SessionManager.hh"

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
    fParticleGun = new G4ParticleGun(1);

    SessionManager & SM = SessionManager::getInstance();

    G4ParticleDefinition* particleDefinition = nullptr;
    double Energy = 0;

    switch (SM.SourceMode)
    {
        case SessionManager::GammaPair :
        {
            particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
            Energy = 511*keV;
        }
        break;

        case SessionManager::C10 :
        {
            particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
            Energy = 0;
            fParticleGun->SetParticleMomentumDirection({0,0,1.0});
        }
        break;
    }

    fParticleGun->SetParticleDefinition(particleDefinition);
    fParticleGun->SetParticleEnergy(Energy);
    fParticleGun->SetParticlePosition({0,0,0});
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    SessionManager & SM = SessionManager::getInstance();

    if(SM.SourceMode == SessionManager::GammaPair)
    {
        //Isotropic momentum direction
        for (int i=0; i < SM.NumParticles; i++)
        {
            double phi = acos(-1 + 2 * G4UniformRand());
            double theta = 2 * M_PI * G4UniformRand();
            double x, y, z;

            x = sin(phi) * cos(theta);
            y = sin(phi) * sin(theta);
            z = cos(phi);

            fParticleGun->SetParticleMomentumDirection({x,y,z});
            fParticleGun->GeneratePrimaryVertex(anEvent);
            fParticleGun->SetParticleMomentumDirection({-x,-y,-z});
            fParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }
    else
    {
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}
