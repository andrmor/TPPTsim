#ifndef SimulationMode_h
#define SimulationMode_h

#include "Modes.hh"

#include <vector>

#include "G4ThreeVector.hh"

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
    DetectorModeEnum DetetctorMode = DetectorModeEnum::OnlyScint;
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
    SimModeShowEvent(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode, int EventToShow);

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

    int NumEvents = 10000;
    std::vector<G4ThreeVector> ScintData;
};

// ---
struct DepositionNodeRecord
{
    DepositionNodeRecord(G4ThreeVector Pos, double Time, double Energy) :
        pos(Pos), time(Time), energy(Energy) {}

    G4ThreeVector pos;
    double        time;
    double        energy;
};
class SimModeMultipleEvents : public SimModeBase
{
public:
    SimModeMultipleEvents(SourceModeEnum sourceMode, DetectorModeEnum detMode, PhantomModeEnum phantMode);

    void run() override;

    G4VSensitiveDetector * getScintDetector() override;

    void saveData(); //might be called multiple times!

    int InitialReserve = 1000; // expected number of nodes per scintillator //later: trigger dump to file when first scint reaches this number of nodes
    std::vector< std::vector<DepositionNodeRecord> > DepositionData;
};

#endif // SimulationMode_h
