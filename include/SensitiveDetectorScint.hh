#ifndef SensitiveDetectorScint_h
#define SensitiveDetectorScint_h

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

class SensitiveDetectorScint : public G4VSensitiveDetector
{
public:
    SensitiveDetectorScint(const G4String & name);

    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
};

#endif // SensitiveDetectorScint_h
