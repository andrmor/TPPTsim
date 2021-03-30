#include "SessionManager.hh"
#include "SimMode.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "out.hh"

#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4VisManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4StepLimiterPhysics.hh"
#include "QGSP_BIC_HP.hh"
#include "G4UImanager.hh"
#include "G4RandomTools.hh"
#include "G4String.hh"
#include "G4ios.hh"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

SessionManager & SessionManager::getInstance()
{
    static SessionManager instance; // Guaranteed to be destroyed, instantiated on first use.
    return instance;
}

SessionManager::SessionManager()
{
    runManager = new G4RunManager;
}

#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
void SessionManager::startSession(int argc, char ** argv)
{
    out("\n\n---------");

    DetectorConstruction * theDetector = new DetectorConstruction();
    runManager->SetUserInitialization(theDetector);

    G4VModularPhysicsList* physicsList = new QGSP_BIC_HP;
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    physicsList->SetDefaultCutValue(0.1*mm);  // Margarida, think about defining Geant4's "regions" - Phantom and the Detector, and using different cut-offs
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserAction(new PrimaryGeneratorAction);

    G4UserSteppingAction * StAct = SimMode->getSteppingAction();
    if (StAct) runManager->SetUserAction(StAct);

    runManager->Initialize();

    configureRandomGenerator();
    configureSource(); //has to be here: after initialize()
    if (SimMode->bNeedGui)    configureGUI(argc, argv);
    if (SimMode->bNeedOutput) configureOutput();
    configureVerbosity();
}

SessionManager::~SessionManager() {}

void SessionManager::configureGUI(int argc, char ** argv)
{
    ui         = new G4UIExecutive(argc, argv);
    visManager = new G4VisExecutive("Quiet");

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/hits/verbose 2");
    UImanager->ApplyCommand("/tracking/verbose 2");
    UImanager->ApplyCommand("/control/saveHistory");

    if (SimMode->DetetctorMode == DetectorModeEnum::WithDetector) scanMaterials();
}

void SessionManager::scanMaterials()
{
    out("-->Scanning materials...");

    std::vector<G4LogicalVolume*> * lvs = G4LogicalVolumeStore::GetInstance();
    for (G4LogicalVolume * lv : *lvs)
    {
        G4Material * mat = lv->GetMaterial();
        out(lv->GetName(), mat->GetName(), mat->GetChemicalFormula(), mat->GetDensity());

        if (mat->GetName() == "G4_Al") lv->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 0, 1.0)));
    }

    out("<--Material scan completed");
}

void SessionManager::startGUI()
{
    visManager->Initialize();
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
}

void SessionManager::configureOutput()
{
    outStream = new std::ofstream();
    outStream->open(WorkingDirectory + "/" + FileName);
    if (!outStream->is_open())
    {
        out("Cannot open file to store output data!");
        delete outStream; outStream = nullptr;
    }
    else out("\nSaving output to file", FileName);
}

void SessionManager::configureRandomGenerator()
{
    randGen = new CLHEP::RanecuEngine();
    randGen->setSeed(Seed);
    G4Random::setTheEngine(randGen);
}

#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
void SessionManager::configureSource()
{
    G4ParticleDefinition * particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
    out(particleDefinition);

    double        Energy    = 0;
    G4ThreeVector Position  = {0, 0, 0};
    G4ThreeVector Direction = {0, 0, 1.0};

    if      (SimMode->SourceMode == SourceModeEnum::Geantino)
    {
        //tests here
        Position  = {0, 0, 100.0};
        Direction = {1.0, 0, 0};
        Energy = 1.0;
    }
    else if (SimMode->SourceMode == SourceModeEnum::GammaPair)
    {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        Energy = 511.0*keV;
    }
    else // assuming one of the PE isotopes
    {
        G4IonTable * ions = G4IonTable::GetIonTable();
        switch (SimMode->SourceMode)
        {
            case SourceModeEnum::C10 : particleDefinition = ions->GetIon(6, 10, 0); break;
            case SourceModeEnum::C11 : particleDefinition = ions->GetIon(6, 11, 0); break;
            case SourceModeEnum::O15 : particleDefinition = ions->GetIon(8, 15, 0); break;
            case SourceModeEnum::N13 : particleDefinition = ions->GetIon(7, 13, 0); break;
            default:;
        };
    };

    ParticleGun->SetParticleDefinition(particleDefinition);
    ParticleGun->SetParticlePosition(Position);
    ParticleGun->SetParticleMomentumDirection(Direction);
    ParticleGun->SetParticleEnergy(Energy);
}

void SessionManager::configureVerbosity()
{
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    if (bVerbose)
    {
        UImanager->ApplyCommand("/hits/verbose 2");
        UImanager->ApplyCommand("/tracking/verbose 2");
        UImanager->ApplyCommand("/control/saveHistory");
    }
    else
    {
        UImanager->ApplyCommand("/hits/verbose 0");
        UImanager->ApplyCommand("/tracking/verbose 0");
        UImanager->ApplyCommand("/control/verbose 0");
        UImanager->ApplyCommand("/run/verbose 0");
    }
    //UImanager->ApplyCommand("/run/initialize");
}

void SessionManager::endSession()
{
    if (outStream) outStream->close();
    delete outStream;

    delete visManager;
    delete runManager;
    delete ui;
}
