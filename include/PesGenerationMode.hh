#ifndef pesgenerationmode_h
#define pesgenerationmode_h

#include "SimMode.hh"

#include <vector>
#include <string>

class PesGenRecord
{
public:
    PesGenRecord(int targetZ, int targetA, std::string pes) : TargetZ(targetZ), TargetA(targetA), PES(pes) {}
    PesGenRecord(){}

    int         TargetZ; // 6;
    int         TargetA; // 12;

    std::string PES;     // "C11";

    std::vector<std::pair<double,double>> CrossSection; // [MeV] [millibarn]

    //runtime
    double NumberDensity;
    double getCrossSection(double energy) const;

    //only for PES_DIRECT mode!
    std::vector<std::vector<std::vector<double>>> * ProbArray = nullptr;
};

// ---

class G4Track;

class PesGenerationMode : public SimModeBase
{
public:
    PesGenerationMode(int numEvents, const std::string & outputFileName, bool binaryOutput);

    std::string getTypeName() const override {return "PesGenerationMode";}
    G4UserStackingAction * getStackingAction() override;

    void preInit() override;
    void run() override;
    void onEventStarted() override;

    bool modelTrigger(const G4Track * track);

    void saveRecord(const std::string & Pes, double X, double Y, double Z, double Time) const;

    std::vector<PesGenRecord> BaseRecords;

    std::vector<std::vector<PesGenRecord>> MaterialRecords; // [indexInMatTable] [Records]

    //only for PES_DIRECT mode!
    double DX = 1.0; // mm
    double DY = 1.0; // mm
    double DZ = 1.0; // mm
    int Xfrom = -100;
    int Xto   =  100;
    int Yfrom = -100;
    int Yto   =  100;
    int Zfrom = -100;
    int Zto   =  100;
    void saveArrays();

protected:
    void loadCrossSections(const std::string & fileName);
    void exploreMaterials();
    void updateMatRecords(int iMat, int Z, int A, double IsotopeNumberDensity);

private:
    int    NumEvents;
    int    CurrentEvent = 0;
    double SaveDir[3];
    double SaveEnergy = 0;

    double LastEnergy;
    double LastTrackLength;
    G4ThreeVector LastPosition;
    int    LastMaterial;

    std::vector<double> ProbVec;

    bool triggerMC(const G4Track * track);
    bool triggerDirect(const G4Track * track);

    bool getVoxel(const G4ThreeVector & pos, int * index);
    void addPath(const G4ThreeVector & posFrom, const G4ThreeVector & posTo, std::vector<std::tuple<int,int,int, double>> & path);
};

#endif // pesgenerationmode_h
