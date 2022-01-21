#include "ActivityProfiler.hh"
#include "out.hh"
#include "jstools.hh"
#include "SessionManager.hh"

#include <iostream>
#include <sstream>
#include <fstream>

ActivityProfilerMode::ActivityProfilerMode(const std::vector<BeamTimeWindow> & beamStructure, const std::vector<ScanTimeWindow> & scanWindows,
                                           const std::string & dataDirectory, const std::string & baseFileName) :
    SimModeBase(),
    BeamTimeWindows(beamStructure), ScanTimeWindows(scanWindows),
    DataDirectory(dataDirectory), BaseFileName(baseFileName)
{
    IsotopeBase.push_back({ "C11", 1220.41518, {"6_12_C11.txt", "8_16_C11.txt"} });
    IsotopeBase.push_back({ "O15", 122.26643,  {"8_16_O15.txt"} });
    IsotopeBase.push_back({ "N13", 597.9,      {"8_16_N13.txt"} });

    SumBeamFactors = 0;
    for (const BeamTimeWindow & win : BeamTimeWindows)
        SumBeamFactors += (win.To - win.From) * win.Weight;
}

void ActivityProfilerMode::run()
{
    // TODO -> check time windows (to > from , weigt > 0 etc)

    // TODO -> check isotopoes are uniques in the database

    std::vector<double> yArr;

    bool bOnStart = true;
    for (const IsotopeDataRecord & iso : IsotopeBase)
    {
        const size_t numRec = iso.SpatialFiles.size();
        out("Isotope", iso.Isotope, " -> ", numRec, "record(s)");
        if (numRec == 0) continue;

        if (bOnStart) // first record of the first isotope -> configure Mapping
        {
            out("OnStart");
            Mapping.read(DataDirectory + '/' + iso.SpatialFiles[0]);
            Mapping.report();
            initArray(DATA);
            bOnStart = false;
        }

        initArray(isoDATA);
        yArr.resize(Mapping.NumBins[1], 0);

        double tauHalf = iso.HalfLife;
        out("Half-time:", tauHalf, "s");

        const double tau = tauHalf / log(2);
        const double timeFactor = calculateTimeFactor(tau);

        for (size_t iR = 0; iR < numRec; iR++)
        {
            std::string fn = DataDirectory + '/' + iso.SpatialFiles[iR];
            out(fn);

            SpatialParameters thisMapping;
            thisMapping.read(fn);
            if (thisMapping != Mapping)
            {
                out("Mismatch in data arrays");
                exit(3);
            }

            std::ifstream in(fn);
            std::string line;
            std::getline(in, line); // first line is the json


            for (int ix = 0; ix < Mapping.NumBins[0]; ix++)
                for (int iy = 0; iy < Mapping.NumBins[1]; iy++)
                {
                    std::getline(in, line);
                    std::string::size_type sz;     // alias of size_t

                    for (int iz = 0; iz < Mapping.NumBins[2]; iz++)
                    {
                        double val = std::stod(line, &sz);
                        line = line.substr(sz);

                        val *= timeFactor;

                        DATA[ix][iy][iz] += val;
                        isoDATA[ix][iy][iz] += val;
                        yArr[iy] += val;
                    }
            }
        }

        for (int iy = 0; iy < Mapping.NumBins[1]; iy++) out(yArr[iy]);

        /*
        hist.NewHist("iso", numY,  OriginY, OriginY + BinSizeY * numY)
                for (var ix = 0; ix < numX; ix++)
                for (var iy = 0; iy < numY; iy++)
        {
            var y = OriginY + BinSizeY * iy
                    for (var iz = 0; iz < numZ; iz++)
            {
                hist.Fill("iso", y, isoDATA[ix][iy][iz])
            }
        }

        hist.Smear("iso", SpatialResolution)
                hist.Draw("iso", "hist")
                grwin.AddToBasket("y_" + isotope)


                core.print("------")
    */
    }
}

bool SpatialParameters::operator!=(const SpatialParameters & other) const
{
    for (int i = 0; i < 3; i++)
    {
        if (BinSize[i] != other.BinSize[i]) return true;
        if (NumBins[i] != other.NumBins[i]) return true;
        if (Origin[i]  != other.Origin[i])  return true;
    }
    return false;
}

void SpatialParameters::read(const std::string & fileName)
{
    if (!SessionManager::isFileExist(fileName))
    {
        out("File", fileName, "does not exist or cannot be open!");
        exit(1);
    }

    std::ifstream in(fileName);
    std::string line;
    std::getline(in, line);
    line.erase(0, 1);

    std::string err;
    json11::Json json = json11::Json::parse(line, err);
    if (!err.empty())
    {
        out(err);
        exit(2);
    }

    json11::Json::array bsar;
    jstools::readArray(json, "BinSize", bsar);
    if (bsar.size() < 3) {out("Cannot read BinSizes"); exit(3);}
    for (int i=0; i<3; i++) BinSize[i] = bsar[i].number_value();

    json11::Json::array nbar;
    jstools::readArray(json, "NumBins", nbar);
    if (nbar.size() < 3) {out("Cannot read NumBins"); exit(3);}
    for (int i=0; i<3; i++) NumBins[i] = nbar[i].int_value();

    json11::Json::array orar;
    jstools::readArray(json, "Origin", orar);
    if (orar.size() < 3) {out("Cannot read Origins"); exit(3);}
    for (int i=0; i<3; i++) Origin[i] = orar[i].number_value();
}

void SpatialParameters::report()
{
    out(
          "BinSizes: (", BinSize[0],',',BinSize[1],',',BinSize[2],") ",
          "NumBins:  (", NumBins[0],',',NumBins[1],',',NumBins[2],") ",
          "Origin:   (", Origin[0], ',',Origin[1], ',',Origin[2], ") "
        );
}

void ActivityProfilerMode::initArray(std::vector<std::vector<std::vector<double>>> & ar)
{
    ar.resize(Mapping.NumBins[0]);
    for (int ix = 0; ix < Mapping.NumBins[0]; ix++)
    {
        ar[ix].resize(Mapping.NumBins[1]);
        for (int iy = 0; iy < Mapping.NumBins[1]; iy++)
        {
            ar[ix][iy].resize(Mapping.NumBins[2]);
            for (int iz = 0; iz < Mapping.NumBins[2]; iz++)
                ar[ix][iy][iz] = 0;
        }
    }
}

double ActivityProfilerMode::calculateTimeFactor(double tau)
{
    double timeFactor = 0;
    for (int iTime = 0; iTime < numTimeRuns; iTime++)
    {
        const double t0 = sampleGenerationTime();
        for (size_t i = 0; i < ScanTimeWindows.size(); i++)
        {
            double timeTo = ScanTimeWindows[i].To - t0;
            if (timeTo  <= 0) continue;
            double timeFrom = ScanTimeWindows[i].From - t0;
            if (timeFrom < 0) timeFrom = 0;

            double from = exp(-timeFrom/tau);
            double to   = exp(-timeTo/tau);
            double delta = from - to;
            timeFactor += delta;
        }
    }
    timeFactor /= numTimeRuns;

    return timeFactor;
}

#include "G4RandomTools.hh"
double ActivityProfilerMode::sampleGenerationTime()
{
    double rnd = G4UniformRand() * SumBeamFactors;

    for (const BeamTimeWindow & win : BeamTimeWindows)
    {
        double val = (win.To - win.From) * win.Weight;
        if (rnd < val)
            return win.From + rnd / win.Weight;

        rnd -= val;
    }

    out("Error in PES time generator");
    exit(3);
}
