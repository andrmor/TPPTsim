#ifndef SimulationMode_h
#define SimulationMode_h

#include "SessionManager.hh"

class G4UserSteppingAction;
class G4VSensitiveDetector;

class SimModeBase
{
public:
    SimModeBase(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode) :
        SourceMode(sourceMode), DetetctorMode(detMode), PhantomMode(phantMode) {}
    virtual ~SimModeBase(){}

    virtual void run() {}
    virtual G4UserSteppingAction * getSteppingAction() {return nullptr;}
    virtual G4VSensitiveDetector * getScintDetector()  {return nullptr;}

    SourceModeEnum   SourceMode    = SourceModeEnum::GammaPair;
    DetectorModeEnum DetetctorMode = DetectorModeEnum::WithDetector;
    PhantomModeEnum  PhantomMode   = PhantomModeEnum::PMMA;

    bool bNeedGui    = false;
    bool bNeedOutput = false;
};

// ---

class SimModeGui : public SimModeBase
{
public:
    SimModeGui(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode);

    void run() override;
};

// ---

class SimModeShowEvent : public SimModeGui
{
public:
    SimModeShowEvent(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode, int EventIndexToShow);

    void run() override;

    int iEvent = 0;
};

// ---
class SimModeScintPosTest : public SimModeBase
{
public:
    SimModeScintPosTest(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode);

    void run() override;

    G4UserSteppingAction * getSteppingAction() override;

    double MaxDelta = 0;
    int    Hits     = 0;
    double SumDelta = 0;
};

// ---

class SimModeSingleEvents : public SimModeBase
{
public:
    SimModeSingleEvents(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode);

    void run() override;

    G4VSensitiveDetector * getScintDetector() override;

    std::vector<G4ThreeVector> ScintData;
};

#endif // SimulationMode_h
