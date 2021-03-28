#include "SensitiveDetectorScint.hh"
#include "SessionManager.hh"
#include "SimMode.hh"

#include <sstream>

#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

SensitiveDetectorScint::SensitiveDetectorScint(const G4String & name)
    : G4VSensitiveDetector(name) {}

G4bool SensitiveDetectorScint::ProcessHits(G4Step* step, G4TouchableHistory*)
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
