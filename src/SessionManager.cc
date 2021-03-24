#include "SessionManager.hh"

#include <iostream>
#include <sstream>
#include <fstream>

#include "G4ParticleDefinition.hh"
#include "G4UImanager.hh"
#include "G4RandomTools.hh"
#include "G4String.hh"

SessionManager & SessionManager::getInstance()
{
    static SessionManager instance; // Guaranteed to be destroyed, instantiated on first use.
    return instance;
}

SessionManager::SessionManager()
{

}

SessionManager::~SessionManager()
{
    delete outStream;
}

void SessionManager::startSession()
{

    outStream = new std::ofstream();
    outStream->open(WorkingDirectory + "/" + BaseFileName);

    if (!outStream->is_open())
    {
        std::cerr << "Cannot open file to store data" << std::endl;
        exit(1);
    }

    randGen = new CLHEP::RanecuEngine();
    randGen->setSeed(Seed);
    G4Random::setTheEngine(randGen);
}

void SessionManager::endSession()
{
    outStream->close();
}

#include <QDebug>
void SessionManager::runSimulation(int NumRuns)
{
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (!bScintPositionTestMode)
    {
        int NumScint = NumScintX * NumScintY * NumRows * NumSegments;
        ScintData.clear();
        std::vector<int> Hits;
        std::stringstream ss;

        for(int i=0; i<NumScint;i++) ScintData.push_back({0,0,0});


        for(int irun=0; irun<NumRuns;irun++)
        {
            UImanager->ApplyCommand("/run/beamOn");
            Hits.clear();

            for(int i=0; i<NumScint;i++)
            {
                if (ScintData[i][1] > 0.500)
                {
                    Hits.push_back(i);
                }
            }

            if (Hits.size() == 2)
            {
                for(int i=0; i<2; i++)
                {
                    int iScint = Hits[i];
                    qDebug()<<iScint<<ScintData[iScint][0]<<ScintData[iScint][1]<<" xyz: "<<ScintPositions[iScint][0]<<ScintPositions[iScint][1]<<ScintPositions[iScint][3]<<irun;
                    *outStream << ScintPositions[iScint][0] / mm << " " << ScintPositions[iScint][1] / mm << " " << ScintPositions[iScint][3] / mm << " " << ScintData[iScint][0] / ns << " " << ScintData[iScint][1] / MeV << std::endl;

                }
                qDebug()<<"---";
            }


            for(int i=0; i<NumScint;i++) ScintData[i] = {0,0,0};
        }

    }
    else
    {
        Hits = 0;
        Errors = 0;

        UImanager->ApplyCommand("/run/beamOn");
        std::cout.flush();
        qDebug() << "Total scintillator hits:" << Hits << " Errors:" << Errors;
        bScintPositionTestMode = false;
    }


}

#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
 void SessionManager::runGUI(G4UImanager *UImanager, G4UIExecutive *ui, G4VisManager* visManager)
{
    UImanager->ApplyCommand("/hits/verbose 2");
    UImanager->ApplyCommand("/tracking/verbose 2");
    UImanager->ApplyCommand("/control/saveHistory");
    visManager->Initialize();
    UImanager->ApplyCommand("/control/execute vis.mac");
    UImanager->ApplyCommand("/hits/verbose 2");
    UImanager->ApplyCommand("/tracking/verbose 2");
    UImanager->ApplyCommand("/control/saveHistory");
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
}

