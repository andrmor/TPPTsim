#include "SteppingAction.hh"
#include "SessionManager.hh"
#include "SimMode.hh"
#include "out.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHandle.hh"
#include "G4NavigationHistory.hh"

void SteppingAction_ScintPosTest::UserSteppingAction(const G4Step * step)
{
    SessionManager & SM = SessionManager::getInstance();
    SimModeScintPosTest * Mode = static_cast<SimModeScintPosTest*>(SM.SimMode);

    const G4StepPoint * postP  = step->GetPostStepPoint();

    if (postP->GetMaterial() == SM.ScintMat)
    {
        Mode->Hits++;

        const G4TouchableHandle & touch = postP->GetTouchableHandle(); //to get the physical volumes
        int iScint    = touch->GetVolume(0)->GetCopyNo(); //this volume (scintillator)
        int iAssembly = touch->GetVolume(1)->GetCopyNo(); //container/master of the volume (encapsulation)

        G4ThreeVector globCenterPos = touch->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0,0,0.5*SM.ScintSizeZ)); //local to global
        double delta = 0;
        for (int i=0; i<3; i++)
        {
            double d = globCenterPos[i] - SM.ScintPositions[iScint][i];
            delta += d*d;
        }
        delta = sqrt(delta);

        Mode->SumDelta += delta;
        if (delta > Mode->MaxDelta) Mode->MaxDelta = delta;

        if (SM.bDebug)
        {
            out("Index of the scintillator:",iScint, " Index of the assembly:",iAssembly);
            out("Volume center position:", globCenterPos[0], globCenterPos[1], globCenterPos[2]);
            out("   --> from ScintPos:  ",SM.ScintPositions[iScint][0], SM.ScintPositions[iScint][1], SM.ScintPositions[iScint][2]); //it's faster (only once in detector construction)
            out("       --> Delta", delta, "mm");

        }
    }
}

// ---

#include "G4VProcess.hh"
#include "G4Track.hh"
void SteppingAction_Tracing::UserSteppingAction(const G4Step *step)
{
    const G4StepPoint * postP  = step->GetPostStepPoint();
    const G4StepPoint * preP   = step->GetPreStepPoint();
    if (!postP) return;

    G4VPhysicalVolume * physVol = preP->GetPhysicalVolume();
    std::string physVolName, logVolName;
    int         copyNum;
    std::string matName;
    G4ThreeVector prePos = preP->GetPosition();

    bool bStart = (step->GetTrack()->GetCurrentStepNumber() == 1);
    if (bStart)
    {
        physVol = preP->GetPhysicalVolume();
        if (physVol)
        {
            physVolName = preP->GetPhysicalVolume()->GetName();
            logVolName  = preP->GetPhysicalVolume()->GetLogicalVolume()->GetName();
            copyNum = preP->GetPhysicalVolume()->GetCopyNo();
            matName = preP->GetMaterial()->GetName();
            out("Created at (", prePos[0], ",",prePos[1],",", prePos[2],")", physVolName, "/", logVolName, "(",copyNum, ")", matName);
        }
        else
            out("Created at (", prePos[0], ",",prePos[1],",", prePos[2],") - seems to be outside the World!" );
    }

    const G4VProcess  * proc = postP->GetProcessDefinedStep();

    std::string procName = "undefined";
    if (proc) procName   = proc->GetProcessName();
    G4ThreeVector pos    = postP->GetPosition();
    physVol = postP->GetPhysicalVolume();

    if (physVol)
    {
        physVolName = postP->GetPhysicalVolume()->GetName();
        logVolName  = postP->GetPhysicalVolume()->GetLogicalVolume()->GetName();
        copyNum = postP->GetPhysicalVolume()->GetCopyNo();
        matName = postP->GetMaterial()->GetName();
        out(procName, " -> (", pos[0], ",",pos[1],",", pos[2],") ", physVolName, "/", logVolName, "(",copyNum, ")", matName);
    }
    else
        out("Exited world at (", pos[0], ",",pos[1],",", pos[2],")" );
}

void SteppingAction_AcollinearityTester::UserSteppingAction(const G4Step *step)
{
    SessionManager & SM = SessionManager::getInstance();
    if (step->GetTrack()->GetParticleDefinition() != SM.GammaPD) return;

    SimModeAcollinTest * Mode = static_cast<SimModeAcollinTest*>(SM.SimMode);
    const G4ThreeVector & v = step->GetPreStepPoint()->GetMomentumDirection();
    out(v);
    //step->GetTrack()->SetTrackStatus(fStopAndKill);

    //if (!step->IsLastStepInVolume()) return;
    //if (step->GetPostStepPoint()->GetPhysicalVolume()) return;
    //exiting World
    //out(step->GetPostStepPoint()->GetPhysicalVolume());
}
