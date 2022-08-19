#include "ActivityGenerationMode.hh"
#include "SessionManager.hh"
#include "out.hh"
#include "jstools.hh"

#include "G4MTRunManager.hh"

ActivityGenerationMode::ActivityGenerationMode(int numEvents,
                                               std::array<double, 3> binSize, std::array<int, 3> numBins, std::array<double, 3> origin,
                                               const std::vector<std::pair<double,double>> & acquisitionFromTos,
                                               const std::string & fileName) :
    PesGenerationMode(numEvents, binSize, numBins, origin),
    TimeWindows(acquisitionFromTos),
    FileName(fileName)
{
    initActivityArray();
}

void ActivityGenerationMode::initActivityArray()
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
            arz.resize(NumBins[2]);
        }
    }
}

void ActivityGenerationMode::run()
{
    SessionManager& SM = SessionManager::getInstance();
    exploreMaterials();

    // this sub-mode is just to debug!
    if (bNeedGui)
    {
        SM.startGUI();
        return;
    }

    SM.runManager->BeamOn(NumEvents);

    saveData();
}

void ActivityGenerationMode::readFromJson(const json11::Json & json)
{
    PesGenerationMode::readFromJson(json);
    initActivityArray();

    TimeWindows.clear();
    json11::Json::array ar;
    jstools::readArray(json, "TimeWindows", ar);
    for (size_t i = 0; i < ar.size(); i++)
    {
        json11::Json::array el = ar[i].array_items();
        TimeWindows.push_back( {el[0].number_value(), el[1].number_value()} );
    }

    jstools::readString(json, "FileName", FileName);
}

void ActivityGenerationMode::doWriteToJson(json11::Json::object & json) const
{
    PesGenerationMode::doWriteToJson(json);

    json11::Json::array ar;
    for (const auto & p : TimeWindows)
    {
        json11::Json::array el;
            el.push_back(p.first);
            el.push_back(p.second);
        ar.push_back(el);
    }
    json["TimeWindows"] = ar;
    json["FileName"] = FileName;
}

void ActivityGenerationMode::doTriggerDirect(const G4Track *track)
{
    const std::vector<PesGenRecord> & Records = MaterialRecords[LastMaterial];
    if (Records.empty()) return;

    std::vector<std::tuple<int, int, int, double>> Path;
    addPathA(LastPosition, track->GetPosition(), Path);
    //out("Path length:", Path.size());
    if (Path.empty()) return;

    const double meanEnergy = 0.5 * (track->GetKineticEnergy() + LastEnergy);
    for (const PesGenRecord & r : Records)
    {
        const double cs = r.getCrossSection(meanEnergy);
        const double DProbByMM = 1e-25 * cs * r.NumberDensity; // millibarn = 0.001e-28m2 -> 0.001e-22mm2 -> 1e-25 mm2

        for (size_t i = 0; i < Path.size(); i++)
        {
            const double numPES = std::get<3>(Path[i]) * DProbByMM;
            const double timeFractionInWindows = calculateTimeFactor(track->GetGlobalTime()/s, r.DecayTime);
            const double decays = numPES * timeFractionInWindows;
            //out(decays);

            Activity[std::get<0>(Path[i])][std::get<1>(Path[i])][std::get<2>(Path[i])] += decays;
        }
    }
}

double ActivityGenerationMode::calculateTimeFactor(double t0, double decayTime)
{
    double timeFactor = 0;

    for (size_t i = 0; i < TimeWindows.size(); i++)
    {
        double timeTo = TimeWindows[i].second - t0;
        if (timeTo  <= 0) continue;
        double timeFrom = TimeWindows[i].first - t0;
        if (timeFrom < 0) timeFrom = 0;

        const double from = exp(-timeFrom/decayTime);
        const double to   = exp(-timeTo/decayTime);
        const double delta = from - to;
        timeFactor += delta;
    }

    //out("Time factor:",timeFactor);
    return timeFactor;
}

void ActivityGenerationMode::saveData()
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
    PesGenerationMode::doWriteToJson(json);
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