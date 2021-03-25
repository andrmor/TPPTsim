#ifndef SimulationMode_h
#define SimulationMode_h

#include "SessionManager.hh"

class G4UserSteppingAction;
class G4VSensitiveDetector;

class SimModeBase
{
public:
    SimModeBase(SourceModeEnum sourceMode) :
        SourceMode(sourceMode) {}
    virtual ~SimModeBase(){}

    virtual void run() {}
    virtual G4UserSteppingAction * getSteppingAction() {return nullptr;}
    virtual G4VSensitiveDetector * getScintDetector()  {return nullptr;}

    SourceModeEnum SourceMode = SourceModeEnum::GammaPair;

    bool bNeedGui    = false;
    bool bNeedOutput = false;
};

class SimModeGui : public SimModeBase
{
public:
    SimModeGui(SourceModeEnum sourceMode);

    void run() override;
};

class SimModeShowEvent : public SimModeGui
{
public:
    SimModeShowEvent(SourceModeEnum sourceMode, int EventIndexToShow);

    void run() override;

    int iEvent = 0;
};

class SimModeScintPosTest : public SimModeBase
{
public:
    SimModeScintPosTest(SourceModeEnum sourceMode);

    void run() override;
    G4UserSteppingAction * getSteppingAction() override;
};

class SimModeSingleEvents : public SimModeBase
{
public:
    SimModeSingleEvents(SourceModeEnum sourceMode);

    void run() override;
    G4VSensitiveDetector * getScintDetector() override;
};



#endif // SimulationMode_h
