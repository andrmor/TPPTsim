#ifndef peshistogramsource_h
#define peshistogramsource_h

#include "SourceMode.hh"
#include "GeantParticleGenerator.hh"

#include <string>
#include <vector>

struct PesDataRecord
{
    PesDataRecord(std::string isotope, const std::vector<std::string> & spatialFiles) :
        Isotope(isotope), SpatialFiles(spatialFiles) {}

    std::string Isotope;
    std::vector<std::string> SpatialFiles;
};

class SourcePesHistogramFiles : public SourceModeBase
{
public:
    SourcePesHistogramFiles(const std::string & directory, double multiplier, bool generateUniformOverBin);  // if not uniform, use bin center
    SourcePesHistogramFiles(const json11::Json & json);

    void GeneratePrimaries(G4Event * anEvent) override;

    double CountEvents() override;

    std::string getTypeName() const override {return "SourcePesHistogramFiles";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    void init();

    void checkInputData();

    std::string Directory;
    double      Multiplier;
    bool        GenerateUniformOverBin = false; // false = generate always in the bin center

    std::vector<PesDataRecord> IsotopeBase;
    size_t CurrentRecord = 0;
    size_t CurrentFile = 0;

    std::ifstream * inStream = nullptr;

    GeantParticleGenerator ParticleGenerator;

    static constexpr double TimeSpan = 10e9; // in ns
};

#endif // peshistogramsource_h
