#ifndef SimulationMode_h
#define SimulationMode_h

#include <vector>

#include "G4ThreeVector.hh"

class G4UserSteppingAction;
class G4VSensitiveDetector;

class SimModeBase
{
public:
    virtual ~SimModeBase(){}

    virtual void run() {}
    virtual G4UserSteppingAction * getSteppingAction() {return nullptr;}
    virtual G4VSensitiveDetector * getScintDetector()  {return nullptr;}

    bool bNeedGui    = false;
    bool bNeedOutput = false;
};

// ---

class SimModeGui : public SimModeBase
{
public:
    SimModeGui();

    void run() override;
};

// ---

class SimModeShowEvent : public SimModeGui
{
public:
    SimModeShowEvent(int EventToShow);

    void run() override;

    int iEvent = 0;
};

// ---
class SimModeScintPosTest : public SimModeBase
{
public:
    SimModeScintPosTest();

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
    SimModeSingleEvents();

    void run() override;

    G4VSensitiveDetector * getScintDetector() override;

    int NumEvents = 10000;
    std::vector<G4ThreeVector> ScintData;
};

// ---
struct DepositionNodeRecord
{
    DepositionNodeRecord(double Time, double Energy) :
        time(Time), energy(Energy) {}

    void merge(const DepositionNodeRecord & other);
    bool isCluster(const DepositionNodeRecord & other, double maxTimeDelta) const;

    double        time;
    double        energy;
};
class SimModeMultipleEvents : public SimModeBase
{
public:
    SimModeMultipleEvents(int numRuns, const std::string & FileName, bool bBinary = false);

    void run() override;

    G4VSensitiveDetector * getScintDetector() override;

    void saveData(); //might be called multiple times!

    int    NumRuns        = 1;
    bool   bDoCluster     = true; // only considers consequtive nodes!
    double MaxTimeDif     = 0.2;
    size_t InitialReserve = 1000; // trigger dump to file when a scintillator accumulated this number of nodes
    std::vector< std::vector<DepositionNodeRecord> > DepositionData;
};

// ---

class SimModeTracing : public SimModeGui
{
public:
    SimModeTracing();

    void run() override;

    G4UserSteppingAction * getSteppingAction() override;
};

// ---

class SimModeAcollinTest : public SimModeBase
{
public:
    SimModeAcollinTest(int numEvents, const std::string &fileName);

    void run() override;

    G4UserSteppingAction * getSteppingAction() override;

protected:
    int NumEvents = 1;

};

#endif // SimulationMode_h
