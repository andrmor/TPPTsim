#include "SimulationMode.hh"
#include "out.hh"

#include "G4RunManager.hh"

SimModeGui::SimModeGui(SourceModeEnum sourceMode) :
    SimModeBase(sourceMode)
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void SimModeGui::run()
{
    SessionManager& SM = SessionManager::getInstance();
    SM.startGUI();
}

SimModeShowEvent::SimModeShowEvent(SourceModeEnum sourceMode, int EventIndexToShow) :
    SimModeGui(sourceMode), iEvent(EventIndexToShow)
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void SimModeShowEvent::run()
{
    SessionManager& SM = SessionManager::getInstance();
    if (iEvent != 0) SM.runManager->BeamOn(iEvent);
    SimModeGui::run();
}

SimModeScintPosTest::SimModeScintPosTest(SourceModeEnum sourceMode) :
    SimModeBase(sourceMode)
{
    bNeedGui    = false;
    bNeedOutput = false;
}

void SimModeScintPosTest::run()
{
    SessionManager& SM = SessionManager::getInstance();
    SM.NumParticles = 10000;

    SM.runManager->BeamOn(1);

    outFlush();
    if (SM.Hits > 1) SM.SumDelta /= SM.Hits;
    out("\n---Test results---\nTotal hits of the scintillators:", SM.Hits, "Max delta:", SM.MaxDelta, " Average delta:", SM.SumDelta, "\n\n");
}

#include "SteppingAction.hh"
G4UserSteppingAction * SimModeScintPosTest::getSteppingAction()
{
    return new SteppingAction;
}

SimModeSingleEvents::SimModeSingleEvents(SourceModeEnum sourceMode) :
    SimModeBase(sourceMode)
{
    bNeedGui    = false;
    bNeedOutput = true;
}

void SimModeSingleEvents::run()
{
    SessionManager& SM = SessionManager::getInstance();

    const int    NumPairs = 10000;
    const double EnergyThreshold = 0.500*MeV;

    const int NumScint = SM.NumScintX * SM.NumScintY * SM.NumRows * SM.NumSegments * 2;
    for(int i=0; i < NumScint; i++) SM.ScintData.push_back({0,0,0});
    std::vector<int> hits;

    for (int iRun = 0; iRun < NumPairs; iRun++)
    {
        SM.runManager->BeamOn(1);

        hits.clear();
        for (int iScint = 0; iScint < NumScint; iScint++)
            if (SM.ScintData[iScint][1] > EnergyThreshold)
                hits.push_back(iScint);

        if (hits.size() == 2)
        {
            for (int i = 0; i < 2; i++)
            {
                const int iScint = hits[i];
                const double Time   = SM.ScintData[iScint][0] / ns;
                const double Energy = SM.ScintData[iScint][1] / MeV;
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

        for (int i = 0; i < NumScint; i++) SM.ScintData[i] = {0,0,0};
    }

    outFlush();
    if (!SM.outStream) out("\nOutput stream was not created, nothing was saved");
    else out("Data saved to file:", SM.FileName);
}

#include "SensitiveDetectorScint.hh"
G4VSensitiveDetector * SimModeSingleEvents::getScintDetector()
{
    return new SensitiveDetectorScint("Scint");
}
