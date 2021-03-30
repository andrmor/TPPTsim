#ifndef SteppingAction_h
#define SteppingAction_h

#include "G4UserSteppingAction.hh"

class G4Step;

class ScintPosTest_SteppingAction : public G4UserSteppingAction
{
public:
    virtual void UserSteppingAction(const G4Step * step) override;
};

// if other stepping actions are needed, define a separate class for each!

#endif // SteppingAction_h
