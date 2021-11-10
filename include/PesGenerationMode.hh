#ifndef pesgenerationmode_h
#define pesgenerationmode_h

#include "SimMode.hh"

#include <vector>
#include <string>

class PesGenRecord
{
public:
    PesGenRecord(int targetZ, int targetA, std::string pes) : TargetZ(targetZ), TargetA(targetA), PES(pes) {}

    int         TargetZ; // 6;
    int         TargetA; // 12;

    std::string PES;     // "C11";

    std::vector<std::pair<double,double>> CrossSection; // [MeV] [millibarn]

    //runtime
    double NumberDensity;
    double getCrossSection(double energy) const;
};

// ---

class PesGenerationMode : public SimModeBase
{
public:
    PesGenerationMode(int numEvents, const std::string & outputFileName, bool binaryOutput);

    std::string getTypeName() const override {return "PesGenerationMode";}
    G4UserStackingAction * getStackingAction() override;

    void preInit() override;
    void run() override;
    void onEventStarted() override;

    void saveRecord(const std::string & Pes, double X, double Y, double Z, double Time) const;

    std::vector<PesGenRecord> BaseRecords;

    std::vector<std::vector<PesGenRecord>> MaterialRecords; // [indexInMatTable] [Records]

protected:
    void exploreMaterials();
    void updateMatRecords(int iMat, int Z, int A, double IsotopeNumberDensity);

private:
    int    NumEvents;
    int    CurrentEvent = 0;
    double SaveDir[3];
    double SaveEnergy = 0;
};

#endif // pesgenerationmode_h
