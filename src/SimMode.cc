#include "SessionManager.hh"
#include "SimMode.hh"
#include "SteppingAction.hh"
#include "SensitiveDetectorScint.hh"
#include "out.hh"

#include "G4RunManager.hh"

SimModeGui::SimModeGui(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode) :
    SimModeBase(sourceMode, detMode, phantMode)
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void SimModeGui::run()
{
    SessionManager& SM = SessionManager::getInstance();
    SM.startGUI();
}

// ---

SimModeShowEvent::SimModeShowEvent(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode, int EventIndexToShow) :
    SimModeGui(sourceMode, detMode, phantMode), iEvent(EventIndexToShow)
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

// ---

SimModeScintPosTest::SimModeScintPosTest(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode) :
    SimModeBase(sourceMode, detMode, phantMode)
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

SimModeSingleEvents::SimModeSingleEvents(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode) :
    SimModeBase(sourceMode, detMode, phantMode)
{
    bNeedGui    = false;
    bNeedOutput = true;

    SessionManager& SM = SessionManager::getInstance();
    SM.FileName = "Coincidence-GammaPairs-Test1.txt";
}

void SimModeSingleEvents::run()
{
    SessionManager& SM = SessionManager::getInstance();

    const int    NumPairs = 10000;
    const double EnergyThreshold = 0.500*MeV;

    const int NumScint = SM.NumScintX * SM.NumScintY * SM.NumRows * SM.NumSegments * 2;
    for(int i=0; i < NumScint; i++) ScintData.push_back({0,0,0});
    std::vector<int> hits;

    for (int iRun = 0; iRun < NumPairs; iRun++)
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
    return new SingleEvents_SensitiveDetectorScint("Scint");
}