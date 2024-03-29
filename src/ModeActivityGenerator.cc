#include "ModeActivityGenerator.hh"
#include "SessionManager.hh"
#include "out.hh"
#include "jstools.hh"

#include "G4MTRunManager.hh"

ModeActivityGenerator::ModeActivityGenerator(int numEvents,
                                               std::array<double, 3> binSize, std::array<int, 3> numBins, std::array<double, 3> origin,
                                               const std::vector<std::pair<double,double>> & acquisitionFromDurationPairs,
                                               const std::string & fileName) :
    ModePesGenerator_Prob(numEvents, binSize, numBins, origin, acquisitionFromDurationPairs),
    FileName(fileName)
{
    initActivityArray();
}

void ModeActivityGenerator::initActivityArray()
{
    // X
    Activity.resize(NumBins[0]);
    for (std::vector<std::vector<double>> & ary : Activity)
    {
        // Y
        ary.resize(NumBins[1]);
        for (std::vector<double> & arz : ary)
        {
            // Z
            arz.resize(NumBins[2], 0);
        }
    }
}

void ModeActivityGenerator::run()
{
    ModePesGenerator_MC::run(); // it is ModePesGenerator_MC indeed!
    saveData();
}

void ModeActivityGenerator::readFromJson(const json11::Json & json)
{
    ModePesGenerator_Prob::readFromJson(json);
    initActivityArray();

    jstools::readString(json, "FileName", FileName);
}

void ModeActivityGenerator::doWriteToJson(json11::Json::object & json) const
{
    ModePesGenerator_Prob::doWriteToJson(json);

    json["FileName"] = FileName;
}

bool ModeActivityGenerator::doTrigger(const G4Track *track)
{
    const std::vector<PesGenRecord> & Records = MaterialRecords[LastMaterial];
    if (Records.empty()) return false;

    std::vector<std::tuple<int, int, int, double>> Path;
    addPathA(LastPosition, track->GetPosition(), Path);
    //out("Path length:", Path.size());
    if (Path.empty()) return false;

    const double meanEnergy = 0.5 * (track->GetKineticEnergy() + LastEnergy);
    for (const PesGenRecord & r : Records)
    {
        const double cs = r.getCrossSection(meanEnergy);
        const double DProbByMM = 1e-25 * cs * r.NumberDensity; // millibarn = 0.001e-28m2 -> 0.001e-22mm2 -> 1e-25 mm2
        const double timeFractionInWindows = calculateAcqusitionTimeFactor(track->GetGlobalTime(), r.DecayTime);

        for (size_t i = 0; i < Path.size(); i++)
        {
            const double numPES = std::get<3>(Path[i]) * DProbByMM;
            const double decays = numPES * timeFractionInWindows;
            //out(decays);

            Activity[std::get<0>(Path[i])][std::get<1>(Path[i])][std::get<2>(Path[i])] += decays;
        }
    }

    return false;
}

void ModeActivityGenerator::saveData()
{
    SessionManager& SM = SessionManager::getInstance();
    std::string fullFileName = SM.WorkingDirectory + '/' + FileName;

    // save resulting 3D array of activity
    std::ofstream stream;
    stream.open(fullFileName);
    if (!stream.is_open())
    {
        out("Cannot open file to store output data!");
        outFlush();
        exit(1);
    }
    else out("\nSaving array data to file", fullFileName);

    json11::Json::object json;
    ModePesGenerator_Prob::doWriteToJson(json);
    json11::Json aa(json);
    std::string str = '#' + aa.dump();
    stream << str << '\n';

    for (const std::vector<std::vector<double>> & ary : Activity)
    {
        for (const std::vector<double> & arz : ary)
        {
            for (const double & val : arz) stream << val << ' ';
            stream << '\n';
        }
    }

    stream.close();
}
