#include "SessionManager.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include <sstream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "QGSP_BIC_HP.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include <QDebug>

int main(int argc, char** argv)
{
    SessionManager& SM = SessionManager::getInstance();

    //SM.bGuiMode         = false; // if true, G4 visualization will start. use beamOn to generate particles

    //We should change the run mode here:
    SM.runMode = SessionManager::ShowEvent;

    SM.Seed = 0;

    G4UIExecutive * ui         = nullptr;
    G4RunManager  * runManager = new G4RunManager;
    G4VisManager* visManager = 0;

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


    if (SM.runMode == SessionManager::GUI || SM.runMode == SessionManager::ShowEvent)
    {
        ui = new G4UIExecutive(argc, argv);
        visManager = new G4VisExecutive("Quiet");

    }

    UImanager->ApplyCommand("/run/setCut 0.1 mm");
    UImanager->ApplyCommand("/run/initialize");

    SM.startSession();


    switch (SM.runMode)
    {
        case SessionManager::GUI:
        {
            SM.runGUI(UImanager, ui, visManager);
        };
        break;

        case SessionManager::ScintPosTest:
        {
            SM.bScintPositionTestMode = true;
            SM.Hits = 0;
            SM.Errors = 0;

            UImanager->ApplyCommand("/run/beamOn");
            std::cout.flush();
            qDebug() << "Total scintillator hits:" << SM.Hits << " Errors:" << SM.Errors;
        };
        break;

        case SessionManager::Main :
        {
            SM.runSimulation(10000);
        };
        break;

        case SessionManager::ShowEvent :
        {
            SM.runSimulation(5263);
            SM.runGUI(UImanager, ui, visManager);
        };
        break;

    }

    qDebug() << "runForrestrun";
    delete visManager;
    delete runManager;
    delete ui;
    qDebug() << "runForrestrun";
    SM.endSession();
}
