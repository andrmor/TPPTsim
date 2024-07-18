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
class Hist1DRegular;
class G4Material;

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
    virtual void onEventEnded() {}

    virtual std::string getTypeName() const = 0;
    void writeToJson(json11::Json::object & json) const;
    virtual void readFromJson(const json11::Json & /*json*/) {}

protected:
    virtual void doWriteToJson(json11::Json::object & /*json*/) const {}
};

// ---

class ModeGui : public SimModeBase
{
public:
    ModeGui();

    void run() override;
    std::string getTypeName() const override {return "ModeGui";}
};

// ---

class ModeDummy : public SimModeBase
{
public:
    ModeDummy(int numEvents);

    void run() override;
    std::string getTypeName() const override {return "ModeDummy";}

    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int NumEvents;
};

// ---

class ModeDoseExtractor : public SimModeBase
{
public:
    ModeDoseExtractor(int numEvents, std::array<double,3> binSize, std::array<int,3> numBins, std::array<double,3> origin,
                      std::string fileName, bool EnergyDepositionOption = false);

    void fill(double energy, const G4ThreeVector & pos, double density);

    void run() override;
    std::string getTypeName() const override {return "ModeDoseExtractor";}
    G4UserSteppingAction * getSteppingAction() override;
    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int                   NumEvents;
    std::array<double, 3> BinSize; // mm
    std::array<int,    3> NumBins;
    std::array<double, 3> Origin;  // center coordinates of the frame
    std::string           FileName;
    bool                  EnergyDepositionOption = false; // collect energy depsosition (do not divide by voxel mass)

    double VoxelVolume;  // mm3

    std::vector<std::vector<std::vector<double>>> Dose;

    void init();

    bool getVoxel(const G4ThreeVector & pos, std::array<int, 3> & index);

    void saveArray();
    void writeBinningToJson(json11::Json::object & json) const;
};

// ---

class ModeEnergyCalibration : public SimModeBase
{
public:
    ModeEnergyCalibration(int numEvents, double binSize, const std::string & fileName);

    void run() override;
    std::string getTypeName() const override {return "ModeEnergyCalibration";}
    G4UserSteppingAction * getSteppingAction() override;
    void readFromJson(const json11::Json & json) override;

    void fill(double depo, double yPos);

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    void   init();
    void   saveRangeData();
    void   calibrate();
    double extractRange(); // returns zero if failed
    void   loadMdaData(const std::string & fileName);
    void   loadSavedRanges(const std::string & fileName);

    const double RecordedRange = 350.0*mm;
    const double RangeLevel    = 0.925; // according to Falk, between 95% abnd 90%

    const double SpotSigmaFactor = 0.93; // Falk

    int         NumEvents;
    double      BinSize;
    std::string FileName;

    //runtime
    std::vector<double> Deposition;
    std::vector<std::pair<double,double>> Ranges; // pairs of (energy_MeV,range_mm)

    std::vector<std::array<double,3>> MdaData; // Energy, Range, Sigma

};

// ---

class ModeShowEvent : public ModeGui
{
public:
    ModeShowEvent(int EventToShow);

    void run() override;
    std::string getTypeName() const override {return "ModeShowEvent";}
    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int iEvent = 0;
};

// ---
class ModeTestScintPositions : public SimModeBase
{
public:
    ModeTestScintPositions();

    void run() override;
    G4UserSteppingAction * getSteppingAction() override;
    std::string getTypeName() const override {return "ModeTestScintPositions";}

    double MaxDelta = 0;
    int    Hits     = 0;
    double SumDelta = 0;
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
class ModeDepositionScint : public SimModeBase
{
public:
    ModeDepositionScint(int numEvents, const std::string & fileName, bool binary,
                          size_t maxCapacity = 10000, bool doCluster  = true, double maxTimeDif  = 0.01*ns);

    void run() override;
    G4VSensitiveDetector * getScintDetector() override;
    void saveData(); // can be called multiple times!
    std::string getTypeName() const override {return "ModeDepositionScint";}
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

class ModeTracing : public ModeGui
{
public:
    ModeTracing();

    void run() override;
    G4UserSteppingAction * getSteppingAction() override;
    std::string getTypeName() const override {return "ModeTracing";}
};

// ---
struct DirAndEnergy
{
    DirAndEnergy(G4ThreeVector Dir, double Energy) : dir(Dir), energy(Energy) {}
    DirAndEnergy(){}

    G4ThreeVector dir;
    double        energy;
};
class ModeTestAcollinearity : public SimModeBase
{
public:
    ModeTestAcollinearity(int numRuns, double range, int numBins, const std::string & fileName); //range: degrees from 180.0
    ~ModeTestAcollinearity();

    void addDirection(const G4ThreeVector & v, int parentID, double energy);

    void run() override;
    G4UserSteppingAction * getSteppingAction() override;
    std::string getTypeName() const override {return "ModeTestAcollinearity";}
    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void init();

    int NumRuns = 1;
    int From    = 0;
    int NumBins = 100;
    std::string FileName  = "dummy.txt";

    Hist1DRegular * Hist = nullptr;
    std::vector<DirAndEnergy> Gammas;

    int ParentTrackId = -1;
    int numNotTherm = 0;
};

// ---

class ModeTestAnnihilations : public SimModeBase
{
public:
    ModeTestAnnihilations(int numEvents, double timeStart, const std::string & fileName, bool binary);

    void saveRecord(const G4ThreeVector & pos, double timeInSeconds);

    void run() override;
    G4UserSteppingAction * getSteppingAction() override;
    std::string getTypeName() const override {return "ModeTestAnnihilations";}
    void readFromJson(const json11::Json & json) override;

    int    NumEvents = 1;
    double TimeStart = 0; // in seconds, ignore all annihilations before that time

protected:
    void doWriteToJson(json11::Json::object & json) const override;
};

// ---

class ModeTestLysoNatRad : public SimModeBase
{
public:
    ModeTestLysoNatRad(int numEvents, int numBins, const std::string & fileName);
    ~ModeTestLysoNatRad();

    void addEnergy(int iScint, double energy);

    void run() override;
    G4UserSteppingAction * getSteppingAction() override;
    std::string getTypeName() const override {return "ModeTestLysoNatRad";}
    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void init();

    int         NumEvents = 1;
    int         NumBins   = 100;
    std::string FileName  = "dummy.txt";

    Hist1DRegular    * Hist      = nullptr;
    std::vector<double> Deposition;
};

// ---

class ModeParticleLogger : public SimModeBase
{
public:
    ModeParticleLogger(int numEvents, const std::string & fileName, bool bBinary);

    void saveParticle(const G4String & particle, double energy_keV, double * PosDir, double time);

    void run() override;
    void onEventStarted() override;
    std::string getTypeName() const override {return "ModeParticleLogger";}
    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int         NumEvents    = 1;
    std::string FileName     = "dummy.txt";
    int         CurrentEvent = 0;
};

// ---

struct DepoStatRec
{
    int    iScint;
    double energy;
    double time;
};
class ModeTestDepositionStat : public SimModeBase
{
public:
    ModeTestDepositionStat(int numEvents, double thresholdMeV = 0.01, std::vector<double> Ranges = {0.05, 0.1}); // energy windows are 0.511 +- Range[i]

    std::string getTypeName() const override {return "ModeTestDepositionStat";}
    void readFromJson(const json11::Json & json) override;

    G4UserSteppingAction * getSteppingAction() override;

    void run() override;
    void onEventStarted() override;

    void addRecord(int iScint, double depo, double time);

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    int NumEvents    = 1;
    std::vector<DepoStatRec> EventRecord;

    void processEventData();

    double Threshold = 0.01; // in MeV


    std::vector<double> Ranges = {0.05, 0.1}; // WARNING! it is assumed an exclusive (<0.5) window: if one depo node is inside, all other nodes should be outside

    int num0 = 0;

    int num1 = 0;
    std::vector<int> Single_In;

    int num2 = 0;
    int Two_SameAssembly = 0;
    std::vector<int>    Two_Same_First_In;
    std::vector<int>    Two_Same_Second_In;
    std::vector<double> Two_Same_AverageDist_First_In;
    std::vector<double> Two_Same_AverageDist_Second_In;
    std::vector<int>    Two_Same_Sum_In;
    std::vector<double> Two_Same_AverageDist_Sum_In;
    std::vector<Hist1DRegular*>Two_Same_HistDist;
    std::vector<int>    Two_Same_FirstSmaller_In;
    std::vector<Hist1DRegular*>Two_Same_HistFirstOverSum;
    std::vector<int>    Two_Dif_First_In;
    std::vector<int>    Two_Dif_Second_In;
    std::vector<double> Two_Dif_AverageDist_First_In;
    std::vector<double> Two_Dif_AverageDist_Second_In;

    //alternative
    std::vector<int>    NoGroup_In_2;
    std::vector<int>    Assembly_In_2;
    std::vector<int>    Global_In_2;

    int num3 = 0;
    std::vector<int>    NoGroup_In_3;
    std::vector<int>    Assembly_In_3;
    std::vector<int>    Global_In_3;

    int num4 = 0;
    std::vector<int>    NoGroup_In_4;
    std::vector<int>    Assembly_In_4;
    std::vector<int>    Global_In_4;

    int num5plus = 0;
    std::vector<int>    NoGroup_In_5plus;
    std::vector<int>    Assembly_In_5plus;
    std::vector<int>    Global_In_5plus;


private:
    void initContainers();
    void reportInt(const std::vector<int> & vec, int scaleBy);
    void reportAvDist(const std::vector<double> & vec, const std::vector<int> & scaleVec);
    void saveDistHist(const std::vector<Hist1DRegular *> & HistDist);
    void reportRatios(const std::vector<int> & vecStat, const std::vector<int> & scaleVec, std::vector<Hist1DRegular*> & vecHist);
    void fillDepoIn(double depo, std::vector<int> & vec);
    void fillDistHist(double depoSum, const G4ThreeVector & pos1, const G4ThreeVector & pos2, std::vector<Hist1DRegular*> & vecHist);
    void incrementDistance(double depo, const G4ThreeVector & v1, const G4ThreeVector & v2, std::vector<double> & vec);
    void fillRatios(double depo1, double depo2, std::vector<int> & vecStat, std::vector<Hist1DRegular*> & vecHist);
    void fillInByGrouping(std::vector<int> & NoGroup_In, std::vector<int> & Assembly_In, std::vector<int> & Global_In);
    void groupRecords(std::vector<DepoStatRec> & records);
    void reportResults();
};

// ---

class ModePositronTimeLogger : public SimModeBase
{
public:
    ModePositronTimeLogger(int numEvents, double timeFrom, double duration, int numBins, double maxOffsetFromIso, std::string fileName);

    void fillTime(double time, const G4ThreeVector & pos);

    void run() override;

    std::string getTypeName() const override {return "ModePositronTimeLogger";}
    void readFromJson(const json11::Json & json) override;

    G4UserStackingAction * getStackingAction() override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    void init();

    int             NumEvents;
    double          TimeFrom;
    double          Duration;
    int             NumBins;
    double          MaxOffset;
    std::string     FileName;
    Hist1DRegular * Hist      = nullptr;
};

// ---

class ModeTesterForAntonio : public SimModeBase
{
public:
    ModeTesterForAntonio();

    void run() override;

    std::string getTypeName() const override {return "ModeTesterForAntonio";}
    //void readFromJson(const json11::Json & json) override;

    G4UserSteppingAction * getSteppingAction() override;

protected:
    //void doWriteToJson(json11::Json::object & json) const override;

};

// ---

class ModeRadHard : public SimModeBase
{
public:
    ModeRadHard(int numEvents);

    void run() override;
    std::string getTypeName() const override {return "ModeRadHard";}
    void readFromJson(const json11::Json & json) override;
    G4UserSteppingAction * getSteppingAction() override;

    int NumEvents = 1;

    G4Material * MatLYSO = nullptr;
    G4Material * MatSiPM = nullptr;

    Hist1DRegular * HistNeutronEn_LYSO = nullptr;
    Hist1DRegular * HistNeutronEn_SiPM = nullptr;
    size_t NumNeutrons_LYSO = 0;
    size_t NumNeutrons_SiPM = 0;
    double Deposition_LYSO = 0;
    double Deposition_SiPM = 0;

    // deposition histogram binning
    const int    NumBins = 200;
    const double From = 0;    // MeV
    const double To   = 50.0; // MeV

protected:
    void doWriteToJson(json11::Json::object & json) const override;

};

// ---

class ModeScintDepoLogger : public SimModeBase
{
public:
    ModeScintDepoLogger(int numEvents, double maxTime, const std::string & fileName);

    void run() override;

    std::string getTypeName() const override {return "ModeScintDepoLogger";}
    //void readFromJson(const json11::Json & json) override;

    void onEventStarted() override;
    void onEventEnded() override;

    G4UserSteppingAction * getSteppingAction() override;

    void addDepo(int iScint, double depo, double time);

protected:
    //void doWriteToJson(json11::Json::object & json) const override;

    int NumEvents = 0;
    double MaxTime = 0.1e9;
    std::vector<double> Depo;
    bool NoDepoThisEvent = true;

};

#endif // SimulationMode_h
