#include "SessionManager.hh"
#include "out.hh"

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

SessionManager::SessionManager(){}

SessionManager::~SessionManager()
{
    delete outStream;
}

#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
void SessionManager::setupGUI(G4UImanager * UImanager, G4UIExecutive * ui, G4VisManager * visManager)
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

void SessionManager::startSession()
{
    out("\n\n---------");
    outStream = new std::ofstream();
    outStream->open(WorkingDirectory + "/" + BaseFileName);

    if (!outStream->is_open())
    {
        out("Cannot open file to store output data, not saving to file!");
        delete outStream; outStream = nullptr;
    }

    randGen = new CLHEP::RanecuEngine();
    randGen->setSeed(Seed);
    G4Random::setTheEngine(randGen);
}

void SessionManager::endSession()
{
    std::cout.flush();
    G4cout.flush();

    if (runMode == ScintPosTest)
    {
        if (Hits > 1) SumDelta /= Hits;
        out("\nTotal scintillator hits:", Hits, "Max delta:", MaxDelta, " Average delta:", SumDelta);
    }

    if (outStream) outStream->close();
    else if (runMode == Main)
        out("\nOutput stream was not created, nothing was saved");
}

bool SessionManager::needGui() const
{
    return (runMode == GUI || runMode == ShowEvent);
}

void SessionManager::runSimulation(int NumRuns)
{
    const double EnergyThreshold = 0.500*MeV;

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    const int NumScint = NumScintX * NumScintY * NumRows * NumSegments * 2;
    for(int i=0; i < NumScint; i++) ScintData.push_back({0,0,0});
    std::vector<int> hits;

    for (int iRun = 0; iRun < NumRuns; iRun++)
    {
        UImanager->ApplyCommand("/run/beamOn");

        hits.clear();
        for (int iScint = 0; iScint < NumScint; iScint++)
            if (ScintData[iScint][1] > EnergyThreshold)
                hits.push_back(iScint);

        if (hits.size() == 2)
        {
            for (int i = 0; i < 2; i++)
            {
                const int iScint = hits[i];
                const double Time   = ScintData[iScint][0] / ns;
                const double Energy = ScintData[iScint][1] / MeV;
                const G4ThreeVector & Pos   = ScintPositions.at(iScint);
                const double X = Pos[0] / mm;
                const double Y = Pos[1] / mm;
                const double Z = Pos[2] / mm;

                out("Scint#",iScint, Time,"ns ", Energy, "MeV  xyz: (",X,Y,Z,")  Run# ",iRun);

                if (outStream)
                    *outStream << X << " " << Y << " " << Z << " " << Time << " " << Energy << std::endl;

            }
            out("---");
        }

        for (int i = 0; i < NumScint; i++) ScintData[i] = {0,0,0};
    }
}
