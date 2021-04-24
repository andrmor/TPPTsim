#ifndef SteppingAction_h
#define SteppingAction_h

#include "G4UserSteppingAction.hh"

class G4Step;

class SteppingAction_ScintPosTest : public G4UserSteppingAction
{
public:
    virtual void UserSteppingAction(const G4Step * step) override;
};

// ---

class SteppingAction_Tracing : public G4UserSteppingAction
{
public:
    virtual void UserSteppingAction(const G4Step * step) override;
};

// ---

class SteppingAction_AcollinearityTester : public G4UserSteppingAction
{
public:
    virtual void UserSteppingAction(const G4Step * step) override;
};

#endif // SteppingAction_h
