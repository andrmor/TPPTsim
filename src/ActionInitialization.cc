#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"
#include "SessionManager.hh"

ActionInitialization::ActionInitialization()
    : G4VUserActionInitialization() {}

ActionInitialization::~ActionInitialization() {}

void ActionInitialization::Build() const
{
    SessionManager & SM = SessionManager::getInstance();

    SetUserAction(new PrimaryGeneratorAction);
    //SteppingAction is used for tests. The following line can be commented:
    if(SM.runMode == SessionManager::ScintPosTest)
        SetUserAction(new SteppingAction);
}
