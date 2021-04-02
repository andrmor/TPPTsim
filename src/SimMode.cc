#include "SessionManager.hh"
#include "Modes.hh"
#include "SimMode.hh"
#include "SteppingAction.hh"
#include "SensitiveDetectorScint.hh"
#include "out.hh"

#include "G4RunManager.hh"

SimModeGui::SimModeGui()
{
    bNeedGui    = true;
    bNeedOutput = false;

    SessionManager& SM = SessionManager::getInstance();
    SM.NumParticlesPerEvent = 1;
}

void SimModeGui::run()
{
    SessionManager& SM = SessionManager::getInstance();
    SM.startGUI();
}

// ---

SimModeShowEvent::SimModeShowEvent(int EventToShow) : iEvent(EventToShow)
{
    bNeedGui    = true;
    bNeedOutput = false;

    SessionManager& SM = SessionManager::getInstance();
    SM.NumParticlesPerEvent = 1;
}

void SimModeShowEvent::run()
{
    SessionManager& SM = SessionManager::getInstance();
    if (iEvent != 0) SM.runManager->BeamOn(iEvent);
    SimModeGui::run();
}

// ---

SimModeScintPosTest::SimModeScintPosTest()
{
    bNeedGui    = false;
    bNeedOutput = false;

    SessionManager& SM = SessionManager::getInstance();
    SM.NumParticlesPerEvent = 10000;
}

void SimModeScintPosTest::run()
{
    SessionManager& SM = SessionManager::getInstance();

    SM.runManager->BeamOn(1);

    outFlush();
    if (Hits > 1) SumDelta /= Hits;
    out("\n---Test results---\nTotal hits of the scintillators:", Hits, "Max delta:", MaxDelta, " Average delta:", SumDelta, "\n\n");
}

G4UserSteppingAction * SimModeScintPosTest::getSteppingAction()
{
    return new ScintPosTest_SteppingAction;
}

// ---

SimModeSingleEvents::SimModeSingleEvents()
{
    bNeedGui    = false;
    bNeedOutput = true;

    NumEvents   = 10000;

    SessionManager& SM = SessionManager::getInstance();
    SM.NumParticlesPerEvent = 1;
    SM.FileName = "Coincidence-GammaPairs-Test1.txt";
}

void SimModeSingleEvents::run()
{
    SessionManager& SM = SessionManager::getInstance();

    const double EnergyThreshold = 0.500*MeV;

    const int NumScint = SM.countScintillators();
    for(int i=0; i < NumScint; i++) ScintData.push_back({0,0,0});
    std::vector<int> hits;

    for (int iRun = 0; iRun < NumEvents; iRun++)
    {
        SM.runManager->BeamOn(1);

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
                const G4ThreeVector & Pos   = SM.ScintPositions.at(iScint);
                const double X = Pos[0] / mm;
                const double Y = Pos[1] / mm;
                const double Z = Pos[2] / mm;

                out("Scint#",iScint, Time,"ns ", Energy, "MeV  xyz: (",X,Y,Z,")  Run# ",iRun);

                if (SM.outStream)
                    *SM.outStream << X << " " << Y << " " << Z << " " << Time << " " << Energy << std::endl;

            }
            out("---");
        }

        for (int i = 0; i < NumScint; i++) ScintData[i] = {0,0,0};
    }

    outFlush();
    if (!SM.outStream) out("\nOutput stream was not created, nothing was saved");
    else out("Data saved to file:", SM.WorkingDirectory + "/" + SM.FileName);
}

G4VSensitiveDetector * SimModeSingleEvents::getScintDetector()
{
    return new SensitiveDetectorScint_SingleEvents("Scint");
}

// ---

SimModeMultipleEvents::SimModeMultipleEvents()
{
    bNeedGui    = false;
    bNeedOutput = true;

    SessionManager& SM = SessionManager::getInstance();
    SM.NumParticlesPerEvent = 10;
    SM.FileName = "TPPToutput-Test1.txt";
    InitialReserve = 10000;

    bDoCluster     = true;
    MaxTimeDif     = 0.2 * ns;
    double MaxR    = 0.5 * mm;
    MaxR2          = MaxR * MaxR;
}

void SimModeMultipleEvents::run()
{
    SessionManager& SM = SessionManager::getInstance();

    const int NumScint = SM.countScintillators();
    DepositionData.resize(NumScint);
    for(auto & vec : DepositionData) vec.reserve(InitialReserve);


    SM.runManager->BeamOn(1);

    saveData();

    outFlush();
    if (!SM.outStream) out("\nOutput stream was not created, nothing was saved");
    else
    {
        out("\nData saved to file:", SM.WorkingDirectory + "/" + SM.FileName);
        if (bDoCluster) out("Depositions were clustered using",MaxTimeDif,"ns time threshold and maxR of",sqrt(MaxR2),"mm");
    }
}

G4VSensitiveDetector *SimModeMultipleEvents::getScintDetector()
{
    return new SensitiveDetectorScint_MultipleEvents("SD");
}

void SimModeMultipleEvents::saveData()
{
    SessionManager& SM = SessionManager::getInstance();
    const int numScint = SM.countScintillators();

    if (SM.outStream)
    {
        for (int iScint = 0; iScint < numScint; iScint++)
        {
            const G4ThreeVector & sp = SM.ScintPositions[iScint];
            auto & nodes = DepositionData[iScint];

            if (!nodes.empty())
            {
                *SM.outStream << "#" << iScint << " " << sp[0] << " " << sp[1] << " " << sp[2] << std::endl;

                for (const DepositionNodeRecord & n : nodes)
                    *SM.outStream << n.pos[0] << " " << n.pos[1] << " " << n.pos[2] << " " << n.time << " " << n.energy << std::endl;
            }
        }
        for (auto & vec : DepositionData) vec.clear();
    }
}

void DepositionNodeRecord::merge(const DepositionNodeRecord & other)
{
    if (other.energy <= 0) return;

    const double newEnergy = energy + other.energy;
    time = (time * energy  +  other.time * other.energy) / newEnergy;
    energy = newEnergy;
}

bool DepositionNodeRecord::isCluster(const DepositionNodeRecord &other, double maxTimeDelta, double maxR2) const
{
    if ( fabs(time - other.time) > maxTimeDelta ) return false;

    double d2 = 0;
    for (int i = 0; i < 3; i++) d2 += (pos[i] - other.pos[i]) * (pos[i] - other.pos[i]);
    if (d2 > maxR2) return false;

    return true;
}

// ---
