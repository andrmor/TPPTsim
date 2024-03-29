#ifndef activitygenerationmode_h
#define activitygenerationmode_h

#include "ModePesGenerator_Prob.hh"
#include "SimMode.hh"
#include "PesGenRecord.hh"

#include <array>
#include <vector>
#include <string>

class G4Track;

// tracking primary protons and using probability approach to generate activity (positron range is not considered!)
class ModeActivityGenerator : public ModePesGenerator_Prob
{
public:
    ModeActivityGenerator(int numEvents,
                           std::array<double,3> binSize, std::array<int,3> numBins, std::array<double,3> origin,
                           const std::vector<std::pair<double,double>> & acquisitionFromDurationPairs,
                           const std::string & fileName);
                           // assumes that acquisitionFromDurationPairs do not overlap!

    std::string getTypeName() const override {return "ModeActivityGenerator";}

    void run() override;
    void onEventStarted() override {} // nothing in this sub-mode

    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

private:
    std::string FileName;

    // output data
    std::vector<std::vector<std::vector<double>>> Activity; // total number of decays in the defined time windows

    bool doTrigger(const G4Track * track) override;

    void   saveData();
    void   initActivityArray();
};

#endif // activitygenerationmode_h
