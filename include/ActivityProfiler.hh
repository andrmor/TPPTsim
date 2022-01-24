#ifndef activityprofilermode_h
#define activityprofilermode_h

#include "SimMode.hh"

#include <vector>
#include <string>
#include <array>

// TODO -> move to a seprate file
struct BeamTimeWindow
{
    BeamTimeWindow(double from, double to) : From(from), To(to) {}
    BeamTimeWindow(double from, double to, double weight) : From(from), To(to), Weight(weight) {}

    double From   = 0;     // in seconds
    double To     = 1e30;  // in seconds
    double Weight = 1.0;
};
struct ScanTimeWindow
{
    ScanTimeWindow(double from, double to) : From(from), To(to) {}

    double From   = 0;     // in seconds
    double To     = 1e30;  // in seconds
};

struct IsotopeDataRecord
{
    IsotopeDataRecord(std::string isotope, double halfLife, const std::vector<std::string> & spatialFiles) :
        Isotope(isotope), HalfLife(halfLife), SpatialFiles(spatialFiles) {}

    std::string Isotope;
    double      HalfLife;
    std::vector<std::string> SpatialFiles;
};

struct SpatialParameters
{
    std::array<double, 3> BinSize;
    std::array<int,    3> NumBins;
    std::array<double, 3> Origin;

    bool operator!=(const SpatialParameters & other) const;

    void read(const std::string & fileName);
    void report();
};

class ActivityProfilerMode : public SimModeBase
{
public:
    ActivityProfilerMode(const std::vector<BeamTimeWindow> & beamStructure, const std::vector<ScanTimeWindow> & scanWindows,
                         const std::string & dataDirectory, const std::string & baseFileName);

    void run() override;

    std::string getTypeName() const override {return "ActivityProfilerMode";}

protected:
    const std::vector<BeamTimeWindow> BeamTimeWindows;
    const std::vector<ScanTimeWindow> ScanTimeWindows;
    const std::string DataDirectory;
    const std::string BaseFileName;

    std::vector<IsotopeDataRecord> IsotopeBase;
    double SumBeamFactors;

    SpatialParameters Mapping;

    int numTimeRuns = 10000;

    void   init3DArray(std::vector<std::vector<std::vector<double>>> & ar);
    void   init2DArray(std::vector<std::vector<double>> & ar);
    double calculateTimeFactor(double tau);
    double sampleGenerationTime();

    void   save1D(std::vector<double> & ar, const std::string & fileName);
    void   save2D(std::vector<std::vector<double>> & ar, const std::string & fileName);
};

#endif // activityprofilermode_h
