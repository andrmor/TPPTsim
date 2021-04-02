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

#include "G4IonTable.hh"
PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
    SessionManager & SM = SessionManager::getInstance();

    ParticleGun = new G4ParticleGun(SM.NumParticlesPerEvent);
    SM.ParticleGun = ParticleGun;

    //the rest is defined in SessionManager::configureSource() -> has to be there since IonTable is recreated during runManager->Initialize();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete ParticleGun; ParticleGun = nullptr;

    //paranoic
    SessionManager & SM = SessionManager::getInstance();
    SM.ParticleGun = nullptr;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event * anEvent)
{
    SessionManager & SM = SessionManager::getInstance();

    if (SM.SourceMode == SourceModeEnum::GammaPair)
    {
        double phi = acos(-1.0 + 2.0 * G4UniformRand());
        double theta = 2.0 * M_PI * G4UniformRand();
        double x, y, z;

        x = sin(phi) * cos(theta);
        y = sin(phi) * sin(theta);
        z = cos(phi);

        ParticleGun->SetParticleMomentumDirection({x,y,z});
        ParticleGun->GeneratePrimaryVertex(anEvent);
        ParticleGun->SetParticleMomentumDirection({-x,-y,-z});
        ParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else
    {
        ParticleGun->GeneratePrimaryVertex(anEvent);
    }
}
