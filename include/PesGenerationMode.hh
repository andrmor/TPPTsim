#ifndef pesgenerationmode_h
#define pesgenerationmode_h

#include "SimMode.hh"

#include <array>
#include <vector>
#include <map>
#include <string>

class PesGenRecord
{
public:
    PesGenRecord(int targetZ, int targetA, std::string pes) : TargetZ(targetZ), TargetA(targetA), PES(pes) {}
    PesGenRecord(){}

    int TargetZ;      // 6;
    int TargetA;      // 12;

    std::string PES;  // "C11";

    std::vector<std::pair<double,double>> CrossSection; // [MeV] [millibarn]
    double DecayTime = -1.0; // decay time (not half-time!) in seconds; Needed only in the DirectMode

    //runtime
    double NumberDensity;
    double getCrossSection(double energy) const;

    //only for direct mode
    std::vector<std::vector<std::vector<double>>> * ProbArray = nullptr;
    double TimeFactor = -1;
};

// ---

class G4Track;

class PesGenerationMode : public SimModeBase
{
public:
    PesGenerationMode(int numEvents, const std::string & outputFileName, bool binaryOutput);  // MC mode
    PesGenerationMode(int numEvents,
                      std::array<double,3> binSize, std::array<int,3> numBins, std::array<double,3> origin,
                      const std::vector<std::pair<double,double>> & timeWindows);             // Direct mode

    std::string getTypeName() const override {return "PesGenerationMode";}
    G4UserStackingAction * getStackingAction() override;

    void preInit() override;
    void run() override;
    void onEventStarted() override;

    bool modelTrigger(const G4Track * track);

    void saveRecord(const std::string & Pes, double X, double Y, double Z, double Time) const;

    std::vector<PesGenRecord>     BaseRecords;
    std::map<std::string, double> DecayTimes;

    std::vector<std::vector<PesGenRecord>> MaterialRecords; // [indexInMatTable] [Records]

    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    void loadCrossSections(const std::string & fileName);
    void loadLifeTimes(const std::string & fileName);
    void exploreMaterials();
    void updateMatRecords(int iMat, int Z, int A, double IsotopeNumberDensity);
    void updateDecayTimes();

private:
    int    NumEvents;

    bool   bDirectMode = false;
    int    CurrentEvent = 0;

    double LastEnergy;
    double LastTrackLength;
    G4ThreeVector LastPosition;
    int    LastMaterial;

    // MC mode specific
    double SaveDirection[3];
    double SaveEnergy = 0;
    std::vector<double> ProbVec;

    // direct mode specific
    std::array<double, 3> BinSize; // mm
    std::array<int,    3> NumBins;
    std::array<double, 3> Origin;  // center coordinates of the frame
    std::vector<std::pair<double,double>> TimeWindows; // observation windows, [from, to] in seconds

    void commonConstructor();
    bool doTriggerMC(const G4Track * track); // return status (true = kill) is now ignored, proton is traced to the end of track
    void doTriggerDirect(const G4Track * track);

    bool getVoxel(const G4ThreeVector & pos, int * index);
    void addPath(const G4ThreeVector & posFrom, const G4ThreeVector & posTo, std::vector<std::tuple<int,int,int, double>> & path);
    void addPathA(const G4ThreeVector & posFrom, const G4ThreeVector & posTo, std::vector<std::tuple<int, int, int, double> > & path);
    bool isValidVoxel(int * coords) const;
    void initProbArrays();
    void saveArrays();
    void clearTimeFactors();
    void calculateTimeFactor(PesGenRecord & r, double globalTime);
};

#endif // pesgenerationmode_h
