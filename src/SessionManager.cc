#include "SessionManager.hh"
#include "SimulationMode.hh"
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
void SessionManager::startSession(int argc, char** argv)
{
    out("\n\n---------");

    DetectorConstruction * theDetector = new DetectorConstruction();
    runManager->SetUserInitialization(theDetector);

    G4VModularPhysicsList* physicsList = new QGSP_BIC_HP;
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    physicsList->SetDefaultCutValue(0.1*mm);  // Margarida, think about defining Geant4's "Regions" - Phantom and the Detector, and using different cut-offs
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserAction(new PrimaryGeneratorAction);

    G4UserSteppingAction * StAct = SimulationMode->getSteppingAction();
    if (StAct) runManager->SetUserAction(StAct);

    //G4UImanager* UImanager = G4UImanager::GetUIpointer();
    //UImanager->ApplyCommand("/control/verbose 0");
    //UImanager->ApplyCommand("/run/verbose 0");
    //UImanager->ApplyCommand("/run/setCut 0.1 mm");   // !!!
    //UImanager->ApplyCommand("/run/initialize");

    runManager->Initialize();

    configureRandomGenerator();
    if (SimulationMode->bNeedGui)    configureGUI(argc, argv);
    if (SimulationMode->bNeedOutput) configureOutput();
    configureVerbosity();
}

SessionManager::~SessionManager() {}

void SessionManager::configureGUI(int argc, char** argv)
{
    ui         = new G4UIExecutive(argc, argv);
    visManager = new G4VisExecutive("Quiet");

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/hits/verbose 2");
    UImanager->ApplyCommand("/tracking/verbose 2");
    UImanager->ApplyCommand("/control/saveHistory");
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
    outStream->open(FileName);
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
