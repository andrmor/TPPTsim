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

    void merge(const DepositionNodeRecord & other);
    bool isCluster(const DepositionNodeRecord & other, double maxTimeDelta, double maxR2) const;

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

    bool   bDoCluster     = true; // only considers consequtive nodes!
    double MaxTimeDif     = 0.2;
    double MaxR2          = 0.5 * 0.5;
    size_t InitialReserve = 1000; // trigger dump to file when a scintillator accumulated this number of nodes
    std::vector< std::vector<DepositionNodeRecord> > DepositionData;
};

#endif // SimulationMode_h
