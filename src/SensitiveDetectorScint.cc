#include "SensitiveDetectorScint.hh"
#include "SessionManager.hh"
#include "SimMode.hh"
#include "out.hh"

#include <sstream>

#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

G4bool SensitiveDetectorScint_SingleEvents::ProcessHits(G4Step* step, G4TouchableHistory*)
{  
    const double edep = step->GetTotalEnergyDeposit();
    if (edep == 0) return true;

    SessionManager & SM = SessionManager::getInstance();
    SimModeSingleEvents * Mode = static_cast<SimModeSingleEvents*>(SM.SimMode);

    const G4StepPoint * postP  = step->GetPostStepPoint();

    int iScint = postP->GetPhysicalVolume()->GetCopyNo();

    double & firstTime = Mode->ScintData[iScint][0];
    double & sumEnergy = Mode->ScintData[iScint][1];

    double time = postP->GetGlobalTime();
    if (time < firstTime || firstTime == 0) firstTime = time;
    sumEnergy += edep;

    return true;
}

// ---

G4bool SensitiveDetectorScint_MultipleEvents::ProcessHits(G4Step *step, G4TouchableHistory *)
{
    const double edep = step->GetTotalEnergyDeposit();
    if (edep == 0) return true;

    SessionManager & SM = SessionManager::getInstance();
    SimModeMultipleEvents * Mode = static_cast<SimModeMultipleEvents*>(SM.SimMode);

    const G4StepPoint * postP  = step->GetPostStepPoint();
    int iScint = postP->GetPhysicalVolume()->GetCopyNo();
    std::vector<DepositionNodeRecord> & Nodes = Mode->DepositionData[iScint];

    G4ThreeVector global = postP->GetPosition();
    G4ThreeVector local = postP->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(global);

    //out(iScint, "    ", global[0], global[1], global[2], "->",local[0], local[1], local[2]);

    DepositionNodeRecord newNode(local, postP->GetGlobalTime(), edep);
    if (!Mode->bDoCluster || Nodes.empty() || !Nodes.back().isCluster(newNode, Mode->MaxTimeDif, Mode->MaxR2))
    {
        if (Nodes.size() == Mode->InitialReserve) Mode->saveData();
        Nodes.push_back(newNode);
    }
    else
        Nodes.back().merge(newNode);

    return true;
}
