#include "SessionManager.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "out.hh"

#include <sstream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "QGSP_BIC_HP.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc, char** argv)
{
    SessionManager& SM = SessionManager::getInstance();

    // ------ START of user inits ------

    //We should change the run mode here:
    SM.runMode = SessionManager::GUI; //GUI, ScintPosTest, Main, ShowEvent
    int iShowEvent = 9338;

    SM.Seed = 0;

    // ------ END of user inits ------

    G4UIExecutive * ui         = nullptr;
    G4RunManager  * runManager = new G4RunManager;
    G4VisManager  * visManager = nullptr;

    DetectorConstruction * theDetector = new DetectorConstruction();
    runManager->SetUserInitialization(theDetector);

    G4VModularPhysicsList* physicsList = new QGSP_BIC_HP;
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new ActionInitialization());

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/run/initialize");
    UImanager->ApplyCommand("/control/verbose 0");
    UImanager->ApplyCommand("/run/verbose 0");

    UImanager->ApplyCommand("/run/setCut 0.1 mm");   // !!!

    UImanager->ApplyCommand("/run/initialize");

    if (SM.needGui())
    {
        ui         = new G4UIExecutive(argc, argv);
        visManager = new G4VisExecutive("Quiet");
    }

    SM.startSession();

    switch (SM.runMode)
    {
    case SessionManager::GUI:
        SM.setupGUI(UImanager, ui, visManager);
        break;

    case SessionManager::ShowEvent :
        SM.runSimulation(iShowEvent);
        SM.setupGUI(UImanager, ui, visManager);
        break;

    case SessionManager::ScintPosTest:
        SM.NumParticles = 10000;
        SM.bScintPositionTestMode = true;

        UImanager->ApplyCommand("/run/beamOn");
        break;

    case SessionManager::Main :
        SM.runSimulation(10000);
        break;
    }

    SM.endSession();

    delete visManager;
    delete runManager;
    delete ui;
}
