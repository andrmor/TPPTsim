#include "SteppingAction.hh"
#include "SessionManager.hh"
#include "out.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHandle.hh"
#include "G4NavigationHistory.hh"

void SteppingAction::UserSteppingAction(const G4Step *step)
{
    SessionManager & SM = SessionManager::getInstance();

    const G4StepPoint * postP  = step->GetPostStepPoint();
    if (!postP) return;

    if (SM.bScintPositionTestMode)
    {
        if (postP->GetMaterial() == SM.ScintMat)
        {
            SM.Hits++;

            const G4TouchableHandle & touch = postP->GetTouchableHandle(); //to get the physical volumes
            int iScint    = touch->GetVolume(0)->GetCopyNo(); //this volume (scintillator)
            int iAssembly = touch->GetVolume(1)->GetCopyNo(); //container/master of the volume (encapsulation)

            G4ThreeVector globCenterPos = touch->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0,0,0)); //local to global
            double delta = 0;
            for (int i=0; i<3; i++)
            {
                double d = globCenterPos[i] - SM.ScintPositions[iScint][i];
                delta += d*d;
            }
            delta = sqrt(delta);

            SM.SumDelta += delta;
            if (delta > SM.MaxDelta) SM.MaxDelta = delta;

            if (SM.bVerbose)
            {
                out("Index of the scintillator:",iScint, " Index of the assembly:",iAssembly);
                out("Volume center position:", globCenterPos[0], globCenterPos[1], globCenterPos[2]);
                out("   --> from ScintPos:  ",SM.ScintPositions[iScint][0], SM.ScintPositions[iScint][1], SM.ScintPositions[iScint][2]); //it's faster (only once in detector construction)
                out("       --> Delta", delta, "mm");

            }
        }
    }
}
