#ifndef gammapairfromannihilhist_h
#define gammapairfromannihilhist_h

#include "SourceMode.hh"
#include "jstools.hh"

#include <string>

class G4ParticleDefinition;

class GammaPairFromAnnihilHist : public SourceModeBase
{
public:
    GammaPairFromAnnihilHist(const std::string & histogramFileName, double activityMultiplier, bool generateUniformOverBin);  // if not uniform, use bin center
    GammaPairFromAnnihilHist(const json11::Json & json);

    void GeneratePrimaries(G4Event * anEvent) override;

    double CountEvents() override;

    std::string getTypeName() const override {return "GammaPairFromAnnihilHist";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    void init();

    std::string HistogramFileName;
    double      ActivityMultiplier;
    bool        GenerateUniformOverBin = false; // false = generate always in the bin center

    static constexpr double TimeSpan = 10e9; // time span of uniform generator, in ns

    // run-time
    size_t CurrentIy = 0;
    BinningParameters Binning;
    std::vector<std::vector<std::vector<double>>> ActivityData;
};

#endif //gammapairfromannihilhist_h
