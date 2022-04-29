#ifndef SimulationMode_h
#define SimulationMode_h

#include "json11.hh"

#include <vector>

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

class G4UserSteppingAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4VSensitiveDetector;
class Hist1D;

class SimModeBase;

class SimModeFactory
{
public:
    static SimModeBase * makeSimModeInstance(const json11::Json & json);
};

class SimModeBase
{
public:
    virtual ~SimModeBase(){}

    bool bNeedGui    = false;
    bool bNeedOutput = false;

    virtual void preInit() {}  // triggered before Geant4 configuration process!

    virtual void run() {}
    virtual G4UserSteppingAction * getSteppingAction() {return nullptr;}
    virtual G4UserStackingAction * getStackingAction() {return nullptr;}
    virtual G4UserTrackingAction * getTrackingAction() {return nullptr;}
    virtual G4VSensitiveDetector * getScintDetector()  {return nullptr;}

    virtual void onEventStarted() {}

    virtual std::string getTypeName() const = 0;
    void writeToJson(json11::Json::object & json) const;
    virtual void readFromJson(const json11::Json & /*json*/) {}

protected:
    virtual void doWriteToJson(json11::Json::object & /*json*/) const {};
};

// ---

class SimModeGui : public SimModeBase
{
public:
    SimModeGui();

    void run() override;
    std::string getTypeName() const override {return "SimModeGui";}
};

// ---

class DoseExtractorMode : public SimModeBase
{
public:
    DoseExtractorMode(int numEvents, std::array<double,3> binSize, std::array<int,3> numBins, std::array<double,3> origin, std::string fileName);

    void fill(double energy, const G4ThreeVector & pos, double density);

    void run() override;
    std::string getTypeName() const override {return "DoseExtractorMode";}
    G4UserSteppingAction * getSteppingAction() override;
    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int                   NumEvents;
    std::array<double, 3> BinSize; // mm
    std::array<int,    3> NumBins;
    std::array<double, 3> Origin;  // center coordinates of the frame

    double VoxelVolume;

    std::vector<std::vector<std::vector<double>>> Dose;

    bool getVoxel(const G4ThreeVector & pos, std::array<int, 3> & index);

    void saveArray();
    void writeBinningToJson(json11::Json::object & json) const;
};

// ---

class SimModeShowEvent : public SimModeGui
{
public:
    SimModeShowEvent(int EventToShow);

    void run() override;
    std::string getTypeName() const override {return "SimModeShowEvent";}
    void readFromJson(const json11::Json & json) override;

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
    void readFromJson(const json11::Json & json) override;

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
    SimModeMultipleEvents(int numEvents, const std::string & fileName, bool binary,
                          size_t maxCapacity = 10000, bool doCluster  = true, double maxTimeDif  = 0.01*ns);

    void run() override;
    G4VSensitiveDetector * getScintDetector() override;
    void saveData(); // can be called multiple times!
    std::string getTypeName() const override {return "SimModeMultipleEvents";}
    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

public:
    //note that the constructor default-configures the folowing settings:
    int         NumEvents;
    size_t      MaxCapacity;  // trigger dump to file when a scintillator accumulated this number of nodes
    bool        bDoCluster;   // only considers consequtive nodes!
    double      MaxTimeDif;

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
    ~SimModeAcollinTest();

    void addDirection(const G4ThreeVector & v, int parentID, double energy);

    void run() override;
    G4UserSteppingAction * getSteppingAction() override;
    std::string getTypeName() const override {return "SimModeAcollinTest";}
    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void init();

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
    SimModeAnnihilTest(int numEvents, double timeStart, const std::string & fileName, bool binary);

    void saveRecord(const G4ThreeVector & pos, double timeInSeconds);

    void run() override;
    G4UserSteppingAction * getSteppingAction() override;
    std::string getTypeName() const override {return "SimModeAnnihilTest";}
    void readFromJson(const json11::Json & json) override;

    int    NumEvents = 1;
    double TimeStart = 0; // in seconds, ignore all annihilations before that time

protected:
    void doWriteToJson(json11::Json::object & json) const override;
};

// ---

class SimModeNatRadTest : public SimModeBase
{
public:
    SimModeNatRadTest(int numEvents, int numBins, const std::string & fileName);
    ~SimModeNatRadTest();

    void addEnergy(int iScint, double energy);

    void run() override;
    G4UserSteppingAction * getSteppingAction() override;
    std::string getTypeName() const override {return "SimModeNatRadTest";}
    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void init();

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

    void saveParticle(const G4String & particle, double energy_keV, double * PosDir, double time);

    void run() override;
    void onEventStarted() override;
    std::string getTypeName() const override {return "SimModeFirstStage";}
    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int         NumEvents    = 1;
    std::string FileName     = "dummy.txt";
    int         CurrentEvent = 0;
};

// ---

class PesAnalyzerMode : public SimModeBase
{
public:
    PesAnalyzerMode(int numEvents, const std::string & fileName);

    void onNewPrimary();
    void registerTarget(const G4String & target, std::vector<G4String> products);
    void onDecay(const G4String & isotope, std::vector<G4String> products);

    void run() override;
    std::string getTypeName() const override {return "PesAnalyzerMode";}
    //void readFromJson(const json11::Json & json) override;

    G4UserSteppingAction * getSteppingAction() override;

protected:
    //void doWriteToJson(json11::Json::object & json) const override;

    int         NumEvents    = 1;

    bool bPositron = false;
    std::pair<G4String, std::vector<G4String>> Target;               // { TargetIsotopeForPrimary,  {Products}  }
    std::vector<std::pair<G4String, std::vector<G4String>>> Decays;  // {  { DecayingIsotope, {Product} }   }

    std::map<std::string, int> Statistics;
    int numPrim = 0;

};

// ---
struct DepoStatRec
{
    int    iScint;
    double energy;
    double time;
};
class DepoStatMode : public SimModeBase
{
public:
    DepoStatMode(int numEvents, const std::string & fileName);

    std::string getTypeName() const override {return "DepoStatMode";}
    //void readFromJson(const json11::Json & json) override;

    G4UserSteppingAction * getSteppingAction() override;

    void run() override;
    void onEventStarted() override;

    void addRecord(int iScint, double depo, double time);

protected:
    //void doWriteToJson(json11::Json::object & json) const override;

    int NumEvents    = 1;
    std::vector<DepoStatRec> EventRecord;

    void processEventData();

    std::vector<double> Ranges = {0.05, 0.1};

    int num0 = 0;

    int num1 = 0;
    std::vector<int> Single_In;

    int num2 = 0;
    int Two_SameAssembly = 0;
    std::vector<int>    Two_Same_first_In;
    std::vector<int>    Two_Same_second_In;
    std::vector<double> Two_Same_AverageDist_First_In;
    std::vector<double> Two_Same_AverageDist_Second_In;
    std::vector<int>    Two_Same_sum_In;
    std::vector<double> Two_Same_AverageDist_Sum_In;
    std::vector<int>    Two_Same_FirstSmaller_In;
    std::vector<double> Two_Same_AvRatioFirstSecond_In;
    std::vector<int>    Two_Dif_first_In;
    std::vector<int>    Two_Dif_second_In;

    int num3 = 0;

private:
    void initContainers();
    void reportInt(const std::vector<int> & vec, int scaleBy);
    void reportAvDist(const std::vector<double> & vec, const std::vector<int> & scaleVec);
    void reportRatios(const std::vector<int> & vecStat, const std::vector<double> & vecRatio, const std::vector<int> & scaleVec);
    void checkDepoIn(double depo, std::vector<int> & vec);
    void incrementDistance(double depo, const G4ThreeVector & v1, const G4ThreeVector & v2, std::vector<double> & vec);
    void fillRatios(double depo1, double depo2, std::vector<int> & vecStat, std::vector<double> & vecRatio);
};

#endif // SimulationMode_h
