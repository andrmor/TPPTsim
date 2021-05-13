#ifndef SimulationMode_h
#define SimulationMode_h

#include <vector>

#include "G4ThreeVector.hh"

class G4UserSteppingAction;
class G4VSensitiveDetector;
class Hist1D;

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
    SimModeMultipleEvents(int numEvents, const std::string & FileName, bool bBinary = false);

    void run() override;

    G4VSensitiveDetector * getScintDetector() override;

    void saveData(); //might be called multiple times!

    int    NumEvents        = 1;
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
struct DirAndEnergy
{
    DirAndEnergy(G4ThreeVector Dir, double Energy) : dir(Dir), energy(Energy) {}
    DirAndEnergy(){}

    G4ThreeVector dir;
    double        energy;
};
class SimModeAcollinTest : public SimModeBase
{
public:
    SimModeAcollinTest(int numRuns, double range, int numBins, const std::string & fileName); //range: degrees from 180.0

    void run() override;

    G4UserSteppingAction * getSteppingAction() override;

    void addDirection(const G4ThreeVector & v, int parentID, double energy);

protected:
    int NumRuns = 1;
    std::string FileName  = "dummy.txt";
    int From = 0;

    Hist1D * Hist = nullptr;
    std::vector<DirAndEnergy> Gammas;

    int ParentTrackId = -1;
    int numNotTherm = 0;
};

// ---

class SimModeAnnihilTest : public SimModeBase
{
public:
    SimModeAnnihilTest(int numEvents, double range, int numBins, const std::string & fileName);
    ~SimModeAnnihilTest();

    void run() override;

    G4UserSteppingAction * getSteppingAction() override;

    void addPosition(double x);

protected:
    int         NumEvents = 1;
    std::string FileName  = "dummy.txt";

    Hist1D    * Hist      = nullptr;
};

#endif // SimulationMode_h
