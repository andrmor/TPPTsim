#ifndef activitygenerationmode_h
#define activitygenerationmode_h

#include "PesGenerationMode.hh"
#include "SimMode.hh"
#include "PesGenRecord.hh"

#include <array>
#include <vector>
#include <string>

class G4Track;

class ActivityGenerationMode : public PesGenerationMode
{
public:
    ActivityGenerationMode(int numEvents,
                           std::array<double,3> binSize, std::array<int,3> numBins, std::array<double,3> origin,
                           const std::vector<std::pair<double,double>> & acquisitionFromTos,
                           const std::string & fileName);

    std::string getTypeName() const override {return "ActivityGenerationMode";}

    void run() override;
    void onEventStarted() override {} // nothing in this sub-mode

    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

private:
    std::vector<std::pair<double,double>> TimeWindows;
    std::string FileName;

    // output data
    std::vector<std::vector<std::vector<double>>> Activity; // total number of decays in the defined time window

    void doTriggerDirect(const G4Track * track);

    double calculateTimeFactor(double t0, double decayTime); // potentially bottleneck -> find a way to use a LUT
    void   saveData();
    void   initActivityArray();
};

#endif // activitygenerationmode_h
