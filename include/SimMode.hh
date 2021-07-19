#ifndef SimulationMode_h
#define SimulationMode_h

#include "json11.hh"

#include <vector>

#include "G4ThreeVector.hh"

class G4UserSteppingAction;
class G4VSensitiveDetector;
class Hist1D;

class SimModeBase
{
public:
    virtual ~SimModeBase(){}

    bool bNeedGui    = false;
    bool bNeedOutput = false;

    virtual void run() {}
    virtual G4UserSteppingAction * getSteppingAction() {return nullptr;}
    virtual G4VSensitiveDetector * getScintDetector()  {return nullptr;}

    virtual void onEventStarted() {}

    virtual std::string getTypeName() const = 0;
    void writeToJson(json11::Json::object & json) const;

protected:
    virtual void doWriteToJson(json11::Json::object & json) const = 0;
};

// ---

class SimModeGui : public SimModeBase
{
public:
    SimModeGui();

    void run() override;

    std::string getTypeName() const override {return "SimModeGui";}
protected:
    void doWriteToJson(json11::Json::object &) const override {};
};

// ---

class SimModeShowEvent : public SimModeGui
{
public:
    SimModeShowEvent(int EventToShow);

    void run() override;

    std::string getTypeName() const override {return "SimModeShowEvent";}
protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int iEvent = 0;
};

// ---
class SimModeScintPosTest : public SimModeBase
{
public:
    SimModeScintPosTest();

    void run() override;

    G4UserSteppingAction * getSteppingAction() override;

    std::string getTypeName() const override {return "SimModeScintPosTest";}
protected:
    void doWriteToJson(json11::Json::object &) const override {};

public:
    double MaxDelta = 0;
    int    Hits     = 0;
    double SumDelta = 0;
};

// ---

class SimModeSingleEvents : public SimModeBase
{
public:
    SimModeSingleEvents(int numEvents);

    void run() override;

    G4VSensitiveDetector * getScintDetector() override;

    std::string getTypeName() const override {return "SimModeSingleEvents";}
protected:
    void doWriteToJson(json11::Json::object & json) const override;

public:
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
    SimModeMultipleEvents(int numEvents, const std::string & fileName, bool binary = false);

    void run() override;

    G4VSensitiveDetector * getScintDetector() override;

    void saveData(); // can be called multiple times!

    std::string getTypeName() const override {return "SimModeMultipleEvents";}
protected:
    void doWriteToJson(json11::Json::object & json) const override;

public:
    int         NumEvents   = 1;
    std::string FileName;
    bool        bBinary;

    bool        bDoCluster  = true; // only considers consequtive nodes!
    double      MaxTimeDif  = 0.2;
    size_t      MaxCapacity = 1000; // trigger dump to file when a scintillator accumulated this number of nodes

    std::vector<std::vector<DepositionNodeRecord>> DepositionData;
};

// ---

class SimModeTracing : public SimModeGui
{
public:
    SimModeTracing();

    void run() override;

    G4UserSteppingAction * getSteppingAction() override;

    std::string getTypeName() const override {return "SimModeTracing";}
protected:
    void doWriteToJson(json11::Json::object &) const override {};
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

    std::string getTypeName() const override {return "SimModeAcollinTest";}
protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int NumRuns = 1;
    int From    = 0;
    int NumBins = 100;
    std::string FileName  = "dummy.txt";

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

    std::string getTypeName() const override {return "SimModeAnnihilTest";}
protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int         NumEvents = 1;
    double      Range     = 4.0;
    int         NumBins   = 100;
    std::string FileName  = "dummy.txt";

    Hist1D    * Hist      = nullptr;
};

// ---

class SimModeNatRadTest : public SimModeBase
{
public:
    SimModeNatRadTest(int numEvents, int numBins, const std::string & fileName);
    ~SimModeNatRadTest();

    void run() override;

    G4UserSteppingAction * getSteppingAction() override;

    void addEnergy(int iScint, double energy);

    std::string getTypeName() const override {return "SimModeNatRadTest";}
protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int         NumEvents = 1;
    int         NumBins   = 100;
    std::string FileName  = "dummy.txt";

    Hist1D    * Hist      = nullptr;
    std::vector<double> Deposition;
};

// ---

class SimModeFirstStage : public SimModeBase
{
public:
    SimModeFirstStage(int numEvents, const std::string & fileName, bool bBinary);

    void run() override;

    void saveParticle(const G4String & particle, double energy_keV, double *PosDir, double time);

    void onEventStarted() override;

    std::string getTypeName() const override {return "SimModeFirstStage";}
protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int         NumEvents    = 1;
    std::string FileName     = "dummy.txt";
    int         CurrentEvent = 0;
};

#endif // SimulationMode_h
