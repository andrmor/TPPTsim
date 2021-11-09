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
    std::vector<std::pair<double,double>> MFP; // [MeV] [mm]
};

// ---

class PesGenerationMode : public SimModeBase
{
public:
    PesGenerationMode();

    std::string getTypeName() const override {return "PesGenerationMode";}

    void preInit() override;

    void run() override;

    std::vector<PesGenRecord> BaseRecords;

    std::vector<std::vector<PesGenRecord>> MaterialRecords; // [indexInMatTable] [Records]

protected:
    void exploreMaterials();
    void updateMatRecords(int iMat, int Z, int A, double IsotopeNumberDensity);
};

#endif // pesgenerationmode_h
