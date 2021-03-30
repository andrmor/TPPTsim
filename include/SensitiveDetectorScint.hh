#ifndef SensitiveDetectorScint_h
#define SensitiveDetectorScint_h

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

class SingleEvents_SensitiveDetectorScint : public G4VSensitiveDetector
{
public:
    SingleEvents_SensitiveDetectorScint(const G4String & name);

    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
};

// for different modes create separate SDs!

#endif // SensitiveDetectorScint_h
