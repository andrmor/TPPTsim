#include "SessionManager.hh"
#include "SimMode.hh"
#include "SteppingAction.hh"
#include "SensitiveDetectorScint.hh"
#include "ModePesGenerator_MC.hh"
#include "ModeActivityGenerator.hh"
#include "ModePesGenerator_Prob.hh"
#include "ModeAnnihilationLogger.hh"
#include "Hist1D.hh"
#include "out.hh"
#include "jstools.hh"

#include "G4String.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <sstream>
#include <fstream>

SimModeBase * SimModeFactory::makeSimModeInstance(const json11::Json & json)
{
    out("Reading simulation mode json");
    std::string Type;
    jstools::readString(json, "Type", Type);

    SimModeBase * sm = nullptr;

      // misc mode
    if      (Type == "ModeGui")                  sm = new ModeGui();
    else if (Type == "ModeShowEvent")            sm = new ModeShowEvent(0);
    else if (Type == "ModeParticleLogger")       sm = new ModeParticleLogger(0, "dummy.txt", false);
    else if (Type == "ModeTracing")              sm = new ModeTracing();
    else if (Type == "ModeEnergyCalibration")    sm = new ModeEnergyCalibration(0, 1, "dummy.txt");
      // main functionality
    else if (Type == "ModeDoseExtractor")        sm = new ModeDoseExtractor(0, {1,1,1}, {1,1,1}, {0,0,0}, "dummy.txt");
    else if (Type == "ModeDepositionScint")      sm = new ModeDepositionScint(0, "dummy.txt", false);
    else if (Type == "ModePesGenerator_MC")      sm = new ModePesGenerator_MC(0, "dummy.txt", false);
    else if (Type == "ModePesGenerator_Prob")    sm = new ModePesGenerator_Prob(0, {1,1,1}, {1,1,1}, {0,0,0}, {{0, 1e20}});
    else if (Type == "ModeActivityGenerator")    sm = new ModeActivityGenerator(0, {1,1,1}, {1,1,1}, {0,0,0}, {{0, 1e20}}, "dummy.txt");
    else if (Type == "ModeAnnihilationLogger")   sm = new ModeAnnihilationLogger(0, {1,1,1}, {1,1,1}, {0,0,0}, "dummy.txt");
      // tests
    else if (Type == "ModePositronTimeLogger")   sm = new ModePositronTimeLogger(0, 0, 1, 1, 1.0, "dummy.txt");
    else if (Type == "ModeTestScintPositions")   sm = new ModeTestScintPositions();
    else if (Type == "ModeTestAcollinearity")    sm = new ModeTestAcollinearity(0, 0, 1, "dummy.txt");
    else if (Type == "ModeTestAnnihilations")    sm = new ModeTestAnnihilations(0, 0, "dummy.txt", false);
    else if (Type == "ModeTestLysoNatRad")       sm = new ModeTestLysoNatRad(0, 0, "dummy.txt");
    else if (Type == "ModeTestDepositionStat")   sm = new ModeTestDepositionStat(0, 0.01, {0.05, 0.1});
    else
    {
        out("Unknown simulation mode type!");
        exit(50);
    }

    sm->readFromJson(json);

    return sm;
}

// ---

void SimModeBase::writeToJson(json11::Json::object & json) const
{
    json["Type"] = getTypeName();
    doWriteToJson(json);
}

// ---

ModeGui::ModeGui()
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void ModeGui::run()
{
    SessionManager& SM = SessionManager::getInstance();
    SM.startGUI();
}

// ---

ModeDoseExtractor::ModeDoseExtractor(int numEvents, std::array<double, 3> binSize, std::array<int, 3> numBins, std::array<double, 3> origin,
                                     std::string fileName, bool EnergyDepositionOption) :
    NumEvents(numEvents), BinSize(binSize), NumBins(numBins), Origin(origin),
    FileName(fileName), EnergyDepositionOption(EnergyDepositionOption)
{
    init();
}

void ModeDoseExtractor::init()
{
    bNeedGui    = false;
    bNeedOutput = true;

    SessionManager & SM = SessionManager::getInstance();
    SM.FileName = FileName;
    SM.bBinOutput = false;

    Dose.resize(NumBins[0]);
    for (std::vector<std::vector<double>> & ary : Dose)
    {
        ary.resize(NumBins[1]);
        for (std::vector<double> & arz : ary)
            arz = std::vector<double>(NumBins[2], 0);
    }

    VoxelVolume = 1.0;
    for (int i = 0; i < 3; i++) VoxelVolume *= BinSize[i]; // mm3
}

void ModeDoseExtractor::run()
{
    SessionManager & SM = SessionManager::getInstance();
    SM.runManager->BeamOn(NumEvents);
    saveArray();
}

G4UserSteppingAction * ModeDoseExtractor::getSteppingAction()
{
    return new SteppingAction_Dose();
}

void ModeDoseExtractor::readFromJson(const json11::Json & json)
{
    SessionManager & SM = SessionManager::getInstance();
    jstools::readInt(json, "NumEvents", NumEvents);
    jstools::readString(json, "FileName", SM.FileName);
    jstools::readBool(json, "EnergyDepositionOption", EnergyDepositionOption);

    //BinSize
    {
        json11::Json::array jar;
        jstools::readArray(json, "BinSize", jar);
        for (int i = 0; i < 3; i++) BinSize[i] = jar[i].number_value();
    }
    // NumBins
    {
        json11::Json::array jar;
        jstools::readArray(json, "NumBins", jar);
        for (int i = 0; i < 3; i++) NumBins[i] = jar[i].int_value();
    }
    // Origin
    {
        json11::Json::array jar;
        jstools::readArray(json, "Origin", jar);
        for (int i = 0; i < 3; i++) Origin[i] = jar[i].number_value();
    }

    init();
}

void ModeDoseExtractor::writeBinningToJson(json11::Json::object & json) const
{
    //BinSize
    {
        json11::Json::array jar;
        for (int i = 0; i < 3; i++) jar.push_back(BinSize[i]);
        json["BinSize"] = jar;
    }
    // NumBins
    {
        json11::Json::array jar;
        for (int i = 0; i < 3; i++) jar.push_back(NumBins[i]);
        json["NumBins"] = jar;
    }
    // Origin
    {
        json11::Json::array jar;
        for (int i = 0; i < 3; i++) jar.push_back(Origin[i]);
        json["Origin"] = jar;
    }
}

void ModeDoseExtractor::doWriteToJson(json11::Json::object & json) const
{
    SessionManager & SM = SessionManager::getInstance();
    json["NumEvents"] = NumEvents;
    json["FileName"] = SM.FileName;
    json["EnergyDepositionOption"] = EnergyDepositionOption;

    writeBinningToJson(json);
}

void ModeDoseExtractor::fill(double energy, const G4ThreeVector & pos, double density)
{
    std::array<int,3> index; // on stack, fast
    const bool ok = getVoxel(pos, index);
    if (!ok) return;

    if (EnergyDepositionOption)
        Dose[index[0]][index[1]][index[2]] += energy;
    else
    {
        const double densityKgPerMM3 = density / (kg/mm3);
        if (densityKgPerMM3 > 1e-10) // ignore vacuum
        {
            const double deltaDose = (energy/joule) / ( densityKgPerMM3 * VoxelVolume );
            Dose[index[0]][index[1]][index[2]] += deltaDose;
        }
    }
}

bool ModeDoseExtractor::getVoxel(const G4ThreeVector & pos, std::array<int,3> & index)
{
    for (int i = 0; i < 3; i++)
    {
        index[i] = floor( (pos[i] - Origin[i]) / BinSize[i] );
        if ( index[i] < 0 || index[i] >= NumBins[i]) return false;
    }
    return true;
}

void ModeDoseExtractor::saveArray()
{
    SessionManager & SM = SessionManager::getInstance();

    json11::Json::object json;
    writeBinningToJson(json);
    json11::Json aa(json);
    std::string str = '#' + aa.dump();
    *SM.outStream << str << '\n';

    for (std::vector<std::vector<double>> & ary : Dose)
    {
        for (std::vector<double> & arz : ary)
        {
            for (double d : arz)
                *SM.outStream << d << ' ';
            *SM.outStream << '\n';
        }
    }
}

// ---

ModeEnergyCalibration::ModeEnergyCalibration(int numEvents, double binSize, const std::string & fileName) :
    SimModeBase(), NumEvents(numEvents), BinSize(binSize), FileName(fileName)
{
    init();
}

void ModeEnergyCalibration::init()
{
    if (BinSize <= 0)
    {
        out("\n\n\n------\nError: BinSize has to be positive!");
        exit(13);
    }

    size_t size = RecordedRange / BinSize + 1;
    Deposition.resize(size);

    loadMdaData("EnergyRangeSigma.txt");
    if (MdaData.size() != 94)
    {
        out(MdaData.size());
        out("Not expected size of MDA data file!");
        exit(4);
    }

    loadSavedRanges("/home/andr/WORK/TPPT/MultiBeam/EnergyCalibration/EnergyRangeWater-1mm-50k.txt");

    saveRangeData();
    calibrate();
    exit(111);
}

#include "SourceMode.hh"
void ModeEnergyCalibration::run()
{
    SessionManager & SM = SessionManager::getInstance();

    for (double energy = 70.0*MeV; energy < 225.0*MeV; energy += 1.0*MeV)
    {
        out("Energy: ", energy);
        SM.SourceMode->setParticleEnergy(energy);

        bool failFlag = false;
        double range;
        do
        {
            std::fill(Deposition.begin(), Deposition.end(), 0);
            SM.runManager->BeamOn(NumEvents);
            range = extractRange();
            failFlag = (range == 0);
        }
        while (failFlag);

        Ranges.push_back( {energy, range} );
        out("Range:", Ranges.back().second);
    }

    saveRangeData();
    calibrate();
}

double ModeEnergyCalibration::extractRange()  // returns zero if failed
{
    size_t iMax = 0;
    double max = 0;
    for (size_t i = 0; i < Deposition.size(); i++)
        if (Deposition[i] > max)
        {
            max = Deposition[i];
            iMax = i;
        }

    int numCrosses = 0;
    int iRange = 0;
    for (size_t i = iMax + 1; i < Deposition.size(); i++)
    {
        if (Deposition[i-1] > max * RangeLevel && Deposition[i] <= max * RangeLevel)
        {
            numCrosses++;
            if (numCrosses > 1)
            {
                out("Multiple cross of the threshold line.");
                return 0;
            }
            iRange = i;
        }
    }

    // interpolating in the range (iRange-1, iRange]
    const double iDelta = (Deposition[iRange-1] - max*RangeLevel) / (Deposition[iRange-1] - Deposition[iRange]);
    const double range = BinSize * (iDelta + iRange-1);
    //out("iRange", iRange, "Lev:",max*RangeLevel, " Depos:", Deposition[iRange-1], Deposition[iRange], " iDelta", iDelta, " Result:", range);
    return range;
}

void ModeEnergyCalibration::loadMdaData(const std::string & fileName)
{
    std::ifstream inStream(fileName);
    if (!inStream.is_open())
    {
        out("Cannot open file with MDA data (energy / water_range / spot_sigma):\n", fileName);
        exit(1);
    }

    MdaData.clear();
    for (std::string line; std::getline(inStream, line); )
    {
        //out(">>>",line);
        if (line.empty()) continue; //allow empty lines
        if (line[0] == '#') continue; //allow comments

        std::stringstream ss(line);  // units in the file are MeV and mbarns
        double Energy, WaterRange, Sigma; //[MeV, cm, mm]
        ss >> Energy >> WaterRange >> Sigma;
        if (ss.fail())
        {
            out("Unexpected format of a line in the file with MDA data");
            exit(3);
        }
        //out(Energy, WaterRange, Sigma);
        MdaData.push_back({Energy*MeV, WaterRange*cm, Sigma*mm});
    }
}

void ModeEnergyCalibration::loadSavedRanges(const std::string & fileName)
{
    std::ifstream inStream(fileName);
    if (!inStream.is_open())
    {
        out("Cannot open file with simulated range data (energy / water_range):\n", fileName);
        exit(1);
    }

    Ranges.clear();
    for (std::string line; std::getline(inStream, line); )
    {
        //out(">>>",line);
        if (line.empty()) continue; //allow empty lines

        std::stringstream ss(line);  // units in the file are MeV and mbarns
        double Energy, WaterRange; //[MeV, mm]
        ss >> Energy >> WaterRange;
        if (ss.fail())
        {
            out("Unexpected format of a line in the file with simulated range data");
            exit(3);
        }
        //out(Energy, WaterRange);
        Ranges.push_back({Energy*MeV, WaterRange*mm});
    }
    //out(Ranges.size());
}

G4UserSteppingAction * ModeEnergyCalibration::getSteppingAction()
{
    return new SteppingAction_EnCal();
}

void ModeEnergyCalibration::readFromJson(const json11::Json & json)
{
    jstools::readInt   (json, "NumEvents", NumEvents);
    jstools::readDouble(json, "BinSize",   BinSize);
    jstools::readString(json, "FileName",  FileName);

    init();
}

void ModeEnergyCalibration::fill(double depo, double yPos)
{
    const size_t iBin = -yPos / BinSize;
    if (iBin >= Deposition.size()) return;

    Deposition[iBin] += depo;
}

void ModeEnergyCalibration::doWriteToJson(json11::Json::object & json) const
{
    json["NumEvents"] = NumEvents;
    json["BinSize"]   = BinSize;
    json["FileName"]  = FileName;
}

void ModeEnergyCalibration::saveRangeData()
{
    SessionManager & SM = SessionManager::getInstance();

    std::string fullFileName = SM.WorkingDirectory + '/' + "Ranges-" + FileName;

    std::ofstream outStream;
    outStream.open(fullFileName);
    if (!outStream.is_open() || outStream.fail() || outStream.bad())
    {
        out("Cannot open file:", fullFileName);
        exit(2);
    }

    out("Saving ranges to file", fullFileName);
    out("\n\n\nEnergy(MeV) --> Range(mm)");

    for (const auto & p : Ranges)
    {
        out(p.first, " --> ", p.second);
        outStream << p.first << ' ' << p.second << '\n';
    }

    outStream.close();
}

void ModeEnergyCalibration::calibrate()
{
    SessionManager & SM = SessionManager::getInstance();

    std::string fullFileName = SM.WorkingDirectory + '/' + FileName;

    std::ofstream outStream;
    outStream.open(fullFileName);
    if (!outStream.is_open() || outStream.fail() || outStream.bad())
    {
        out("Cannot open file:", fullFileName);
        exit(2);
    }

    out("Saving calibration data to file", fullFileName);
    out("\n\n\nNominal_energy(MeV), True_energy(MeV), SpotSigma(mm)");

    size_t lastRangeRec = 1;
    for (const auto & ar : MdaData)
    {
        const double & nominalEnergy = ar[0];
        const double & nominalRange  = ar[1];
        const double & nominalSigma  = ar[2];

        double trueEnergy;
        //real energy from range
        for (size_t i = lastRangeRec; i < Ranges.size(); i++)  // Energy(first) / Range(second)
        {
            const double & thisRange     = Ranges[i].second;
            const double & previousRange = Ranges[i-1].second;
            if (thisRange < nominalRange) continue;

            lastRangeRec = i;
            const double fraction = (nominalRange - previousRange) / (thisRange - previousRange);
            trueEnergy = SessionManager::interpolate(Ranges[i-1].first, Ranges[i].first, fraction); //a + fraction * (b−a)
            break;
        }

        //sigma
        const double sigma = nominalSigma * SpotSigmaFactor;

        outStream << nominalEnergy << ' ' << trueEnergy << ' ' << sigma << '\n';
    }

    outStream.close();
}

// ---

ModeShowEvent::ModeShowEvent(int EventToShow) : iEvent(EventToShow)
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void ModeShowEvent::run()
{
    SessionManager& SM = SessionManager::getInstance();
    if (iEvent != 0) SM.runManager->BeamOn(iEvent);
    ModeGui::run();
}

void ModeShowEvent::readFromJson(const json11::Json & json)
{
    jstools::readInt(json, "iEvent", iEvent);
}

void ModeShowEvent::doWriteToJson(json11::Json::object & json) const
{
    json["iEvent"] = iEvent;
}

// ---

ModeTestScintPositions::ModeTestScintPositions()
{
    bNeedGui    = false;
    bNeedOutput = false;
}

void ModeTestScintPositions::run()
{
    SessionManager& SM = SessionManager::getInstance();

    SM.runManager->BeamOn(10000);

    outFlush();
    if (Hits > 1) SumDelta /= Hits;
    out("\n---Test results---\nTotal hits of the scintillators:", Hits, "Max delta:", MaxDelta, " Average delta:", SumDelta, "\n\n");
}

G4UserSteppingAction * ModeTestScintPositions::getSteppingAction()
{
    return new SteppingAction_ScintPosTest;
}

// ---

ModeDepositionScint::ModeDepositionScint(int numEvents, const std::string & fileName, bool binary,
                                             size_t maxCapacity, bool doCluster, double maxTimeDif) :
    NumEvents(numEvents), MaxCapacity(maxCapacity), bDoCluster(doCluster), MaxTimeDif(maxTimeDif)
{
    bNeedGui    = false;
    bNeedOutput = true;

    SessionManager & SM = SessionManager::getInstance();
    SM.bBinOutput  = binary;
    SM.FileName    = fileName;
}

void ModeDepositionScint::run()
{
    SessionManager& SM = SessionManager::getInstance();

    const int NumScint = SM.countScintillators();
    DepositionData.resize(NumScint);
    for(auto & vec : DepositionData) vec.reserve(MaxCapacity);

    SM.runManager->BeamOn(NumEvents);

    saveData();

    outFlush();
    if (!SM.outStream) out("\nOutput stream was not created, nothing was saved");
    else
    {
        out("\nData saved to file:", SM.WorkingDirectory + "/" + SM.FileName);
        if (bDoCluster) out("Depositions were clustered using",MaxTimeDif,"ns time threshold");
    }
}

G4VSensitiveDetector *ModeDepositionScint::getScintDetector()
{
    return new SensitiveDetectorScint_MultipleEvents("SD");
}

void ModeDepositionScint::saveData()
{
    SessionManager& SM = SessionManager::getInstance();
    const int numScint = SM.countScintillators();

    if(SM.bBinOutput)
    {
        for (int iScint = 0; iScint < numScint; iScint++)
        {
            auto & nodes = DepositionData[iScint];

            if (!nodes.empty())
            {
                *SM.outStream << char(0xee);
                SM.outStream->write((char*)&iScint, sizeof(int));

                for (size_t iNodes = 0; iNodes < nodes.size(); iNodes++)
                {
                    *SM.outStream << char(0xff);
                    SM.outStream->write((char*)&DepositionData[iScint][iNodes].time,   sizeof(double));
                    SM.outStream->write((char*)&DepositionData[iScint][iNodes].energy, sizeof(double));
                }
            }
        }
    }
    else
    {
        for (int iScint = 0; iScint < numScint; iScint++)
        {
            auto & nodes = DepositionData[iScint];

            if (!nodes.empty())
            {
                *SM.outStream << "# " << iScint << std::endl;

                for (const DepositionNodeRecord & n : nodes)
                    *SM.outStream << n.time << " " << n.energy << std::endl;
            }
        }
        for (auto & vec : DepositionData) vec.clear();
    }
}

void ModeDepositionScint::readFromJson(const json11::Json &json)
{
    SessionManager & SM = SessionManager::getInstance();

    jstools::readInt   (json, "NumEvents",   NumEvents);
    jstools::readString(json, "FileName",    SM.FileName);
    jstools::readBool  (json, "bBinary",     SM.bBinOutput);

    jstools::readBool  (json, "bDoCluster",  bDoCluster);
    jstools::readDouble(json, "MaxTimeDif",  MaxTimeDif);
    int iMaxCapacity = 10000;
    jstools::readInt   (json, "MaxCapacity", iMaxCapacity);
    MaxCapacity = iMaxCapacity;
}

void ModeDepositionScint::doWriteToJson(json11::Json::object &json) const
{
    SessionManager & SM = SessionManager::getInstance();

    json["NumEvents"]   = NumEvents;
    json["FileName"]    = SM.FileName;
    json["bBinary"]     = SM.bBinOutput;

    json["bDoCluster"]  = bDoCluster;
    json["MaxTimeDif"]  = MaxTimeDif;
    json["MaxCapacity"] = (int)MaxCapacity;
}

void DepositionNodeRecord::merge(const DepositionNodeRecord & other)
{
    if (other.energy <= 0) return;

    const double newEnergy = energy + other.energy;
    time = (time * energy  +  other.time * other.energy) / newEnergy;
    energy = newEnergy;
}

bool DepositionNodeRecord::isCluster(const DepositionNodeRecord &other, double maxTimeDelta) const
{
    if ( fabs(time - other.time) > maxTimeDelta ) return false;
    return true;
}

// ---

ModeTracing::ModeTracing()
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void ModeTracing::run()
{
    ModeGui::run();
}

G4UserSteppingAction * ModeTracing::getSteppingAction()
{
    return new SteppingAction_Tracing;
}

// ---

ModeTestAcollinearity::ModeTestAcollinearity(int numRuns, double range, int numBins, const std::string & fileName) :
    NumRuns(numRuns), NumBins(numBins), FileName(fileName)
{
    From = 180.0 - range;
    init();
}

ModeTestAcollinearity::~ModeTestAcollinearity()
{
    delete Hist;
}

void ModeTestAcollinearity::init()
{
    delete Hist; Hist = new Hist1DRegular(NumBins, From, 180.0);
}

void ModeTestAcollinearity::run()
{
    SessionManager& SM = SessionManager::getInstance();

    for (int iRun = 0; iRun < NumRuns; iRun++)
    {
        ParentTrackId = -1;
        //out("Run #", iRun);

        SM.runManager->BeamOn(1);

        if (Gammas.size() >= 2)
        {
            double angle = Gammas[0].dir.angle(Gammas[1].dir) / deg - 1e-10; // -1e-10 to get 180 deg in the previous bin
            Hist->fill(angle);

            if (angle < From)
            {
                if (Gammas[0].energy > 0.511 || Gammas[0].energy < 0.51 ||
                    Gammas[1].energy > 0.511 || Gammas[1].energy < 0.51)
                {
                    numNotTherm++;
                    //out("Not thermalized positron");
                    //out("Undeflow angle:", angle);
                    //out("Dir/Energies:");
                    //for (const DirAndEnergy & de : Gammas) std::cout << de.dir << " " << de.energy << std::endl;
                }
                else
                {
                    out("strange event!");
                    out("Undeflow angle:", angle);
                    out("Dir/Energies:");
                    for (const DirAndEnergy & de : Gammas) std::cout << de.dir << " " << de.energy << std::endl;
                }
            }
        }
        else
        {
            //out("Unexpected: number of gammas is", Gammas.size());
        }
        Gammas.clear();
    }

    outFlush();
    out("\nDistribution of inter-gamma angles (from", From,"to 180 deg):");
    Hist->report();
    out("NotThermalized:", numNotTherm);
    Hist->save(SM.WorkingDirectory + '/' + FileName);
}

G4UserSteppingAction *ModeTestAcollinearity::getSteppingAction()
{
    return new SteppingAction_AcollinearityTester;
}

void ModeTestAcollinearity::addDirection(const G4ThreeVector & v, int parentID, double energy)
{
    // the first gamma to track will be annihilation one, some of the following ones can be secondary ones!
    //out("ParentId:", parentID, "Energy:", energy);
    if (ParentTrackId == -1) ParentTrackId = parentID;
    else
    {
        if (parentID != ParentTrackId) return;
    }

    Gammas.push_back( DirAndEnergy(v, energy) );
}

void ModeTestAcollinearity::readFromJson(const json11::Json &json)
{
    jstools::readInt   (json, "NumRuns",  NumRuns);
    jstools::readInt   (json, "From",     From);
    jstools::readInt   (json, "NumBins",  NumBins);
    jstools::readString(json, "FileName", FileName);

    init();
}

void ModeTestAcollinearity::doWriteToJson(json11::Json::object &json) const
{
    json["NumRuns"]  = NumRuns;
    json["From"]     = From;
    json["NumBins"]  = NumBins;
    json["FileName"] = FileName;
}

// ---

ModeTestAnnihilations::ModeTestAnnihilations(int numEvents, double timeStart, const std::string & fileName, bool binary) :
    NumEvents(numEvents), TimeStart(timeStart)
{
    bNeedOutput = true;

    SessionManager & SM = SessionManager::getInstance();
    SM.FileName = fileName;
    SM.bBinOutput = binary;
}

void ModeTestAnnihilations::saveRecord(const G4ThreeVector & pos, double timeInSeconds)
{
    SessionManager & SM = SessionManager::getInstance();
    *SM.outStream << pos[0] << ' ' << pos[1] << ' ' << pos[2] << ' ' << timeInSeconds << '\n';
}

G4UserSteppingAction * ModeTestAnnihilations::getSteppingAction()
{
    return new SteppingAction_AnnihilationTester;
}

void ModeTestAnnihilations::run()
{
    SessionManager & SM = SessionManager::getInstance();
    SM.runManager->BeamOn(NumEvents);
}

void ModeTestAnnihilations::readFromJson(const json11::Json &json)
{
    SessionManager & SM = SessionManager::getInstance();

    jstools::readInt   (json, "NumEvents", NumEvents);
    jstools::readDouble(json, "TimeStart", TimeStart);
    jstools::readString(json, "FileName",  SM.FileName);
    jstools::readBool  (json, "Binary",    SM.bBinOutput);
}

void ModeTestAnnihilations::doWriteToJson(json11::Json::object &json) const
{
    SessionManager & SM = SessionManager::getInstance();

    json["NumEvents"] = NumEvents;
    json["TimeStart"] = TimeStart;
    json["FileName"]  = SM.FileName;
    json["Binary"]    = SM.bBinOutput;
}

// ---

ModeTestLysoNatRad::ModeTestLysoNatRad(int numEvents, int numBins, const std::string & fileName) :
    NumEvents(numEvents), NumBins(numBins), FileName(fileName)
{
    init();

    /*
    // This code allows to test a standalone scintillator (will be two on the opposite sides of the ring unless detectorConstructor code is modified!):
    //bNeedGui = true;
    SM.NumScintX = 1;
    SM.NumScintY = 1;
    SM.NumSegments = 1;
    SM.NumRows     = 1;
    SM.ScintSizeX = 57.4;
    SM.ScintSizeY = 57.4;
    SM.ScintSizeZ = 10.0;
    SM.EncapsSizeX = SM.EncapsSizeY = SM.ScintSizeX + 0.2;
    SM.EncapsSizeZ = SM.ScintSizeZ + 0.2;
    */
}

void ModeTestLysoNatRad::init()
{
    delete Hist; Hist = new Hist1DRegular(NumBins, 0, 1.3);
}

ModeTestLysoNatRad::~ModeTestLysoNatRad()
{
    delete Hist;
}

void ModeTestLysoNatRad::run()
{
    SessionManager & SM = SessionManager::getInstance();

    const int numScint = SM.countScintillators();
    Deposition.resize(numScint);

    for (int iEv = 0; iEv < NumEvents; iEv++)
    {
        for (int iSc = 0; iSc < numScint; iSc++)
            Deposition[iSc] = 0;

        SM.runManager->BeamOn(1);

        for (int iSc = 0; iSc < numScint; iSc++)
            if (Deposition[iSc] != 0)
                Hist->fill(Deposition[iSc]);
    }

    outFlush();
    out("\nDistribution of deposited energies[MeV] per scintillator per event:");
    Hist->report();
    Hist->save(SM.WorkingDirectory + '/' + FileName);
}

G4UserSteppingAction * ModeTestLysoNatRad::getSteppingAction()
{
    return new SteppingAction_NatRadTester;
}

void ModeTestLysoNatRad::addEnergy(int iScint, double energy)
{
    Deposition[iScint] += energy;
}

void ModeTestLysoNatRad::readFromJson(const json11::Json &json)
{
    jstools::readInt   (json, "NumEvents", NumEvents);
    jstools::readInt   (json, "NumBins",   NumBins);
    jstools::readString(json, "FileName",  FileName);

    init();
}

void ModeTestLysoNatRad::doWriteToJson(json11::Json::object &json) const
{
    json["NumEvents"] = NumEvents;
    json["NumBins"]   = NumBins;
    json["FileName"]  = FileName;
}

// ---

ModeParticleLogger::ModeParticleLogger(int numEvents, const std::string & fileName, bool bBinary) :
    NumEvents(numEvents)
{
    bNeedOutput = true;

    SessionManager & SM = SessionManager::getInstance();
    SM.FileName   = fileName;
    SM.bBinOutput = bBinary;
}

void ModeParticleLogger::run()
{
    CurrentEvent = 0;
    SessionManager & SM = SessionManager::getInstance();
    SM.runManager->BeamOn(NumEvents);
}

void ModeParticleLogger::saveParticle(const G4String & particle, double energy_keV, double * PosDir, double time)
{
    SessionManager & SM = SessionManager::getInstance();

    if (SM.bBinOutput)
    {
        *SM.outStream << char(0xFF);
        *SM.outStream << particle << char(0x00);
        SM.outStream->write((char*)&energy_keV,  sizeof(double));
        SM.outStream->write((char*)PosDir, 6*sizeof(double));
        SM.outStream->write((char*)&time,    sizeof(double));
    }
    else
    {
        std::stringstream ss;
        ss << particle << ' ';
        ss << energy_keV << ' ';
        ss << PosDir[0] << ' ' << PosDir[1] << ' ' << PosDir[2] << ' ';     //position
        ss << PosDir[3] << ' ' << PosDir[4] << ' ' << PosDir[5] << ' ';     //direction
        ss << time;

        *SM.outStream << ss.rdbuf() << '\n';
    }

    //out("->",particle, energy, "(",PosDir[0],PosDir[1],PosDir[2],")", "(",PosDir[3],PosDir[4],PosDir[5],")",time);
}

void ModeParticleLogger::onEventStarted()
{
    SessionManager & SM = SessionManager::getInstance();
    if (SM.bBinOutput)
    {
        *SM.outStream << char(0xEE);
        SM.outStream->write((char*)&CurrentEvent, sizeof(int));
    }
    else
        *SM.outStream << '#' << CurrentEvent << '\n';

    CurrentEvent++;
}

void ModeParticleLogger::readFromJson(const json11::Json &json)
{
    SessionManager & SM = SessionManager::getInstance();

    jstools::readInt   (json, "NumEvents", NumEvents);
    jstools::readString(json, "FileName",  SM.FileName);
    jstools::readBool  (json, "bBinary",   SM.bBinOutput);
}

void ModeParticleLogger::doWriteToJson(json11::Json::object &json) const
{
    SessionManager & SM = SessionManager::getInstance();

    json["NumEvents"] = NumEvents;
    json["FileName"]  = SM.FileName;
    json["bBinary"]   = SM.bBinOutput;
}

// ---

ModeTestDepositionStat::ModeTestDepositionStat(int numEvents, double thresholdMeV, std::vector<double> ranges) :
    SimModeBase(), NumEvents(numEvents), Threshold(thresholdMeV), Ranges(ranges) {}

void ModeTestDepositionStat::run()
{
    initContainers();

    SessionManager & SM = SessionManager::getInstance();
    SM.runManager->BeamOn(NumEvents);
    processEventData(); // for the last event!

    reportResults();
}

void ModeTestDepositionStat::doWriteToJson(json11::Json::object & json) const
{
    json["NumEvents"] = NumEvents;
    json["Threshold"] = Threshold;

    json11::Json::array jar;
    for (size_t i = 0; i < Ranges.size(); i++) jar.push_back(Ranges[i]);
    json["Ranges"] = jar;
}

void ModeTestDepositionStat::readFromJson(const json11::Json & json)
{
    jstools::readInt   (json, "NumEvents", NumEvents);
    jstools::readDouble(json, "Threshold", Threshold);

    json11::Json::array jar;
    jstools::readArray(json, "Ranges", jar);
    for (size_t i = 0; i < jar.size(); i++) Ranges[i] = jar[i].number_value();
}

void ModeTestDepositionStat::initContainers()
{
    size_t num = Ranges.size();

    Single_In.resize(num, 0);
    Two_Same_First_In.resize(num, 0);
    Two_Same_Second_In.resize(num, 0);
    Two_Same_AverageDist_First_In.resize(num, 0);
    Two_Same_AverageDist_Second_In.resize(num, 0);
    Two_Same_Sum_In.resize(num, 0);
    Two_Same_AverageDist_Sum_In.resize(num, 0);
    Two_Same_HistDist.resize(num, nullptr);
    for (size_t i = 0; i < num; i++)
        Two_Same_HistDist[i] = new Hist1DRegular(100, 0, 30.0);
    Two_Same_FirstSmaller_In.resize(num, 0);
    Two_Same_HistFirstOverSum.resize(num, nullptr);
    for (size_t i = 0; i < num; i++)
        Two_Same_HistFirstOverSum[i] = new Hist1DRegular(100, 0, 1.0);
    Two_Dif_First_In.resize(num, 0);
    Two_Dif_Second_In.resize(num, 0);
    Two_Dif_AverageDist_First_In.resize(num, 0);
    Two_Dif_AverageDist_Second_In.resize(num, 0);

    NoGroup_In_2.resize(num, 0);
    NoGroup_In_3.resize(num, 0);
    NoGroup_In_4.resize(num, 0);
    NoGroup_In_5plus.resize(num, 0);
    Assembly_In_2.resize(num, 0);
    Assembly_In_3.resize(num, 0);
    Assembly_In_4.resize(num, 0);
    Assembly_In_5plus.resize(num, 0);
    Global_In_2.resize(num, 0);
    Global_In_3.resize(num, 0);
    Global_In_4.resize(num, 0);
    Global_In_5plus.resize(num, 0);
}

void ModeTestDepositionStat::reportInt(const std::vector<int> & vec, int scaleBy)
{
    for (size_t i = 0; i < Ranges.size(); i++)
    {
        double factor = Ranges[i];
        out("  -->", 100.0 * vec[i] / scaleBy, "% within +-", std::to_string(factor*100.0), "% of 511 keV");
    }
}

void ModeTestDepositionStat::reportAvDist(const std::vector<double> & vec, const std::vector<int> & scaleVec)
{
    for (size_t i = 0; i < Ranges.size(); i++)
    {
        double factor = Ranges[i];
        out("  -->", vec[i] / scaleVec[i], "mm for +-", std::to_string(factor*100.0), "% of 511 keV");
    }
}

void ModeTestDepositionStat::saveDistHist(const std::vector<Hist1DRegular*> & HistDist)
{
    const SessionManager & SM = SessionManager::getInstance();
    for (size_t i=0; i<Ranges.size(); i++)
    {
        double factor = Ranges[i];
        HistDist[i]->save(SM.WorkingDirectory + "/histDistTwoSumInside-" + std::to_string(factor*100.0) + ".txt");
    }
}

void ModeTestDepositionStat::reportRatios(const std::vector<int> & vecStat, const std::vector<int> & scaleVec, std::vector<Hist1DRegular *> & vecHist)
{
    const SessionManager & SM = SessionManager::getInstance();
    for (size_t i = 0; i < Ranges.size(); i++)
    {
        double factor = Ranges[i];
        out("  -->Fraction with first depo smaller:", (double)vecStat[i] / scaleVec[i],
            " for +-", std::to_string(factor*100.0), "% of 511 keV");

        vecHist[i]->save(SM.WorkingDirectory + "/histFirstOverSum-" + std::to_string(factor*100.0) + ".txt");
    }
}

void ModeTestDepositionStat::onEventStarted()
{
    processEventData(); //previous event!

    EventRecord.clear();
}

G4UserSteppingAction * ModeTestDepositionStat::getSteppingAction()
{
    return new SteppingAction_DepoStatMode();
}

void ModeTestDepositionStat::addRecord(int iScint, double depo, double time)
{
    //out("-->Scint:", iScint, " MeV:", depo, "  ns:", time);

    for (DepoStatRec & r : EventRecord)
    {
        if (r.iScint == iScint)
        {
            r.energy += depo;
            if (time < r.time) r.time = time;
            return;
        }
    }
    EventRecord.emplace_back(DepoStatRec{iScint, depo, time});
}

void ModeTestDepositionStat::fillDepoIn(double depo, std::vector<int> & vec)
{
    for (size_t i = 0; i < Ranges.size(); i++)
    {
        if ( depo > (0.511 * (1.0 - Ranges[i])) &&
             depo < (0.511 * (1.0 + Ranges[i]))
           ) vec[i]++;
    }
}

void ModeTestDepositionStat::fillDistHist(double depoSum, const G4ThreeVector & pos1, const G4ThreeVector & pos2, std::vector<Hist1DRegular *> & vecHist)
{
    for (size_t i = 0; i < Ranges.size(); i++)
    {
        if ( depoSum > (0.511 * (1.0 - Ranges[i])) &&
             depoSum < (0.511 * (1.0 + Ranges[i]))
           )
        {
            const G4ThreeVector d = pos1 - pos2;
            vecHist[i]->fill(d.getR());
        }
    }
}

void ModeTestDepositionStat::incrementDistance(double depo, const G4ThreeVector & v1, const G4ThreeVector & v2, std::vector<double> & vec)
{
    for (size_t i = 0; i < Ranges.size(); i++)
    {
        if ( depo > (0.511 * (1.0 - Ranges[i])) &&
             depo < (0.511 * (1.0 + Ranges[i]))
           )
        {
            const G4ThreeVector d = v1 - v2;
            vec[i] += d.getR();
        }
    }
}

void ModeTestDepositionStat::fillRatios(double depo1, double depo2, std::vector<int> & vecStat, std::vector<Hist1DRegular *> & vecHist)
{
    const double sumDepo = depo1 + depo2;
    for (size_t i = 0; i < Ranges.size(); i++)
    {
        if ( sumDepo > (0.511 * (1.0 - Ranges[i])) &&
             sumDepo < (0.511 * (1.0 + Ranges[i]))
           )
        {
            if (depo1 < depo2) vecStat[i]++;
            vecHist[i]->fill(depo1/sumDepo);
        }
    }
}

void ModeTestDepositionStat::processEventData()
{
    const SessionManager & SM = SessionManager::getInstance();
    //for (const DepoStatRec & r : EventRecord)
    //    out("Scint:", r.iScint, " MeV:", r.energy, "  ns:", r.time);

    // applying threshold
    if (!EventRecord.empty())
    {
        for (int iThis = EventRecord.size()-1; iThis >= 0; iThis--)
        {
            const DepoStatRec & Rec = EventRecord[iThis];
            if (Rec.energy < Threshold)
                EventRecord.erase(EventRecord.begin() + iThis);
        }
    }

    if      (EventRecord.empty())     num0++;
    else if (EventRecord.size() == 1)
    {
        num1++;
        fillDepoIn(EventRecord.front().energy, Single_In);
    }
    else if (EventRecord.size() == 2)
    {
        num2++;

        const DepoStatRec & rFirst  = EventRecord.front();
        const DepoStatRec & rSecond = EventRecord.back();

        bool bSameAssembly = ( SM.ScintRecords[rFirst.iScint].AssemblyNumber == SM.ScintRecords[rSecond.iScint].AssemblyNumber );
        if (bSameAssembly)
        {
            Two_SameAssembly++;
            const double depo1 = rFirst.energy;
            fillDepoIn (depo1, Two_Same_First_In);
            const double depo2 = rSecond.energy;
            fillDepoIn (depo2, Two_Same_Second_In);

            incrementDistance(depo1, SM.ScintRecords[rFirst.iScint].FacePos, SM.ScintRecords[rSecond.iScint].FacePos, Two_Same_AverageDist_First_In);
            incrementDistance(depo2, SM.ScintRecords[rFirst.iScint].FacePos, SM.ScintRecords[rSecond.iScint].FacePos, Two_Same_AverageDist_Second_In);

            const double sumdepo = depo1 + depo2;
            fillDepoIn(sumdepo, Two_Same_Sum_In);
            fillDistHist(sumdepo, SM.ScintRecords[rFirst.iScint].FacePos, SM.ScintRecords[rSecond.iScint].FacePos, Two_Same_HistDist);

            incrementDistance(sumdepo, SM.ScintRecords[rFirst.iScint].FacePos, SM.ScintRecords[rSecond.iScint].FacePos, Two_Same_AverageDist_Sum_In);

            fillRatios(depo1, depo2, Two_Same_FirstSmaller_In, Two_Same_HistFirstOverSum);
        }
        else
        {
            const double depo1 = rFirst.energy;
            fillDepoIn(depo1, Two_Dif_First_In);
            const double depo2 = rSecond.energy;
            fillDepoIn(depo2, Two_Dif_Second_In);

            incrementDistance(depo1, SM.ScintRecords[rFirst.iScint].FacePos, SM.ScintRecords[rSecond.iScint].FacePos, Two_Dif_AverageDist_First_In);
            incrementDistance(depo2, SM.ScintRecords[rFirst.iScint].FacePos, SM.ScintRecords[rSecond.iScint].FacePos, Two_Dif_AverageDist_Second_In);
        }


        fillInByGrouping(NoGroup_In_2, Assembly_In_2, Global_In_2);
    }
    else if (EventRecord.size() == 3)
    {
        num3++;
        fillInByGrouping(NoGroup_In_3, Assembly_In_3, Global_In_3);
    }
    else if (EventRecord.size() == 4)
    {
        num4++;
        fillInByGrouping(NoGroup_In_4, Assembly_In_4, Global_In_4);
    }
    else
    {
        num5plus++;
        fillInByGrouping(NoGroup_In_5plus, Assembly_In_5plus, Global_In_5plus);
    }
}

void ModeTestDepositionStat::fillInByGrouping(std::vector<int> & NoGroup_In, std::vector<int> & Assembly_In, std::vector<int> & Global_In)
{
    for (DepoStatRec & rec : EventRecord)
    {
        const double depo = rec.energy;
        fillDepoIn(depo, NoGroup_In);
    }

    std::vector<DepoStatRec> GrouppedRecord = EventRecord;
    groupRecords(GrouppedRecord);
    for (DepoStatRec & rec : GrouppedRecord)
    {
        const double depo = rec.energy;
        fillDepoIn(depo, Assembly_In);
    }

    double sumDepo = 0;
    for (const DepoStatRec & rec : GrouppedRecord) sumDepo += rec.energy;
    fillDepoIn(sumDepo, Global_In);
}

void ModeTestDepositionStat::groupRecords(std::vector<DepoStatRec> & records)
{
    if (records.size() < 2)  return;

    const SessionManager & SM = SessionManager::getInstance();

    for (int iThis = records.size()-1; iThis > 0; iThis--)
    {
        DepoStatRec & ThisRec = records[iThis];
        for (int iCheck = 0; iCheck < iThis; iCheck++)
        {
            DepoStatRec & CheckRec = records[iCheck];
            if (SM.ScintRecords[ThisRec.iScint].AssemblyNumber == SM.ScintRecords[CheckRec.iScint].AssemblyNumber)
            {
                CheckRec.energy += ThisRec.energy;
                records.erase(records.begin() + iThis);
                break;
            }
        }
    }
}

void ModeTestDepositionStat::reportResults()
{
    // results
    out("\n\n------------------");
    out("Total number of gammas generated:", NumEvents);
    out("Threshold:", Threshold, '\n');

    int remains = NumEvents;

    out("No interaction:",       100.0 * num0/NumEvents, "%"); remains -= num0;

    num1--; // first dummy event disregarded
    out("Single scintillator: ", 100.0 * num1/NumEvents, "%"); remains -= num1;
    reportInt(Single_In, NumEvents);

    out("Two scintillators: ",   100.0 * num2/NumEvents, "%"); remains -= num2;
    out("  Two within the same assembly:", 100.0*Two_SameAssembly/NumEvents, "%");
    out("    First:");
    reportInt(Two_Same_First_In,  NumEvents);
    out("    Second:");
    reportInt(Two_Same_Second_In, NumEvents);
    out("      Average distance between scints");
    out("        when first in:");
    reportAvDist(Two_Same_AverageDist_First_In,  Two_Same_First_In);
    out("        when second in:");
    reportAvDist(Two_Same_AverageDist_Second_In, Two_Same_Second_In);
    out("  Sum of two depositions:");
    reportInt(Two_Same_Sum_In, NumEvents);
    out("      Average distance between scints when sum energy is in window:");
    reportAvDist(Two_Same_AverageDist_Sum_In, Two_Same_Sum_In);
    saveDistHist(Two_Same_HistDist);
    reportRatios(Two_Same_FirstSmaller_In, Two_Same_Sum_In, Two_Same_HistFirstOverSum);
    out("  Two within different assemblies:");
    out("    First:");
    reportInt(Two_Dif_First_In, NumEvents);
    out("    Second:");
    reportInt(Two_Dif_Second_In, NumEvents);
    out("      Average distance between scints:");
    out("        when first in:");
    reportAvDist(Two_Dif_AverageDist_First_In,  Two_Dif_First_In);
    out("        when second in:");
    reportAvDist(Two_Dif_AverageDist_Second_In, Two_Dif_Second_In);

    out("Three scintillators: ", 100.0 * num3/NumEvents, "%"); remains -= num3;
    out("Four  scintillators: ", 100.0 * num4/NumEvents, "%"); remains -= num4;
    out("Five+ scintillators: ", 100.0 * num5plus/NumEvents, "%"); remains -= num5plus;

    //out("Remainer:", remains);

    out("------\n");

    out("If no grouping is made, \"good\" events are:");
    for (size_t i = 0; i < Ranges.size(); i++)
    {
        double factor = Ranges[i];
        out("Window of +-", std::to_string(factor*100.0), "% of 511 keV:");

        double tot1 = 100.0 / NumEvents * Single_In[i];
        double tot2 = 100.0 / NumEvents * NoGroup_In_2[i];
        double tot3 = 100.0 / NumEvents * NoGroup_In_3[i];
        double tot4 = 100.0 / NumEvents * NoGroup_In_4[i];
        double totP = 100.0 / NumEvents * NoGroup_In_5plus[i];
        double tot = tot1 + tot2 + tot3 + tot4 + totP;
        out("  Total", tot, "     1:", tot1, "  2:", tot2, "  3:", tot3, "  4:", tot4, "  5+:", totP);
    }
    out("------\n");
    out("If grouping in assemblies is made, \"good\" events are:");
    for (size_t i = 0; i < Ranges.size(); i++)
    {
        double factor = Ranges[i];
        out("Window of +-", std::to_string(factor*100.0), "% of 511 keV:");

        double tot1 = 100.0 / NumEvents * Single_In[i];
        double tot2 = 100.0 / NumEvents * Assembly_In_2[i];
        double tot3 = 100.0 / NumEvents * Assembly_In_3[i];
        double tot4 = 100.0 / NumEvents * Assembly_In_4[i];
        double totP = 100.0 / NumEvents * Assembly_In_5plus[i];
        double tot = tot1 + tot2 + tot3 + tot4 + totP;
        out("  Total", tot, "     1:", tot1, "  2:", tot2, "  3:", tot3, "  4:", tot4, "  5+:", totP);
    }
    out("------\n");
    out("If global grouping is made, \"good\" events are:");
    for (size_t i = 0; i < Ranges.size(); i++)
    {
        double factor = Ranges[i];
        out("Window of +-", std::to_string(factor*100.0), "% of 511 keV:");

        double tot1 = 100.0 / NumEvents * Single_In[i];
        double tot2 = 100.0 / NumEvents * Global_In_2[i];
        double tot3 = 100.0 / NumEvents * Global_In_3[i];
        double tot4 = 100.0 / NumEvents * Global_In_4[i];
        double totP = 100.0 / NumEvents * Global_In_5plus[i];
        double tot = tot1 + tot2 + tot3 + tot4 + totP;
        out("  Total", tot, "     1:", tot1, "  2:", tot2, "  3:", tot3, "  4:", tot4, "  5+:", totP);
    }

    //out( NumEvents , num0 , num1 , num2 , num3 , num4 , num5plus , remains);
}

// ---

ModePositronTimeLogger::ModePositronTimeLogger(int numEvents, double timeFrom, double duration, int numBins, double maxOffsetFromIso, std::string fileName) :
    NumEvents(numEvents), TimeFrom(timeFrom), Duration(duration), NumBins(numBins), MaxOffset(maxOffsetFromIso), FileName(fileName)
{
    init();
}

void ModePositronTimeLogger::fillTime(double time, const G4ThreeVector & pos)
{
    //out(pos[0],pos[1],pos[2]);
    if (fabs(pos[0]) < MaxOffset && fabs(pos[1]) < MaxOffset && fabs(pos[2]) < MaxOffset)
    {
        Hist->fill(time);
        //out("-->Accept");
    }
}

void ModePositronTimeLogger::run()
{
    SessionManager & SM = SessionManager::getInstance();
    SM.runManager->BeamOn(NumEvents);
    Hist->report();
    Hist->save(SM.WorkingDirectory + '/' + FileName);
}

void ModePositronTimeLogger::readFromJson(const json11::Json &json)
{
    jstools::readInt   (json, "NumEvents", NumEvents);
    jstools::readDouble(json, "TimeFrom",  TimeFrom);
    jstools::readDouble(json, "Duration",  Duration);
    jstools::readInt   (json, "NumBins",   NumBins);
    jstools::readDouble(json, "MaxOffset", MaxOffset);
    jstools::readString(json, "FileName",  FileName);

    init();
}

#include "StackingAction.hh"
G4UserStackingAction * ModePositronTimeLogger::getStackingAction()
{
    return new PositronTimeLoggerStackingAction();
}

void ModePositronTimeLogger::doWriteToJson(json11::Json::object &json) const
{
    json["NumEvents"] = NumEvents;
    json["TimeFrom"]  = TimeFrom;
    json["Duration"]  = Duration;
    json["NumBins"]   = NumBins;
    json["MaxOffset"] = MaxOffset;
    json["FileName"]  = FileName;
}

void ModePositronTimeLogger::init()
{
    delete Hist; Hist = new Hist1DRegular(NumBins, TimeFrom, TimeFrom + Duration);
}

// ---

ModeTesterForAntonio::ModeTesterForAntonio() :
    SimModeBase() {}

void ModeTesterForAntonio::run()
{
    SessionManager & SM = SessionManager::getInstance();

    if (!SM.PhantomLogical)
    {
        out("Phantom is not defined!");
        exit(111);
    }

    SM.runManager->BeamOn(1);

    const G4MaterialTable * theMaterialTable = G4Material::GetMaterialTable();
    const size_t numMat = theMaterialTable->size();
    out("Defined", numMat, "materials");

    for (size_t iMat = 0; iMat < numMat; iMat++)
    {
        G4Material * mat = (*theMaterialTable)[iMat];
        out("Setting material to:", mat->GetName());

        SM.PhantomLogical->SetMaterial(mat);

        SM.runManager->Initialize();

        SM.runManager->BeamOn(1);
    }
}

G4UserSteppingAction * ModeTesterForAntonio::getSteppingAction()
{
    return new SteppingAction_TesterForAntonio();
}

// ---------------------

ModeRadHard::ModeRadHard(int numEvents) :
    NumEvents(numEvents) {}

void ModeRadHard::run()
{
    SessionManager & SM = SessionManager::getInstance();

    HistNeutronEn_LYSO = new Hist1DRegular(200, 0, 50);
    HistNeutronEn_SiPM = new Hist1DRegular(200, 0, 50);

    const G4MaterialTable * theMaterialTable = G4Material::GetMaterialTable();
    const size_t numMat = theMaterialTable->size();
    out("Defined", numMat, "materials");
    for (size_t iMat = 0; iMat < numMat; iMat++)
    {
        G4Material * mat = (*theMaterialTable)[iMat];
        if (mat->GetName() == "SIPM")   MatSiPM = mat;
        if (mat->GetName() == "LYSOCe") MatLYSO = mat;
    }
    if (!MatLYSO || !MatSiPM)
    {
        out("LYSO or SiPM material not defined!");
        exit(222);
    }

    out("Lyso:", MatLYSO, MatLYSO->GetName(), MatSiPM, MatSiPM->GetName());

    SM.runManager->BeamOn(NumEvents);

    size_t numScint = SM.countScintillators();
    //size_t numSiPM  = SM.NumRows * SM.NumSegments * 2.0;

    out("\n", "\n", "\n");
    out("--Deposition--");
    out("  In LYSO", Deposition_LYSO, "MeV", "NumScint:", numScint, "  Depo per one:", Deposition_LYSO / numScint);
    out("  In SiPM", Deposition_SiPM, "MeV", "NumSiPMs:", numScint, "  Depo per one:", Deposition_SiPM / numScint);
    out("--Neutrons--");
    out("  Number in LYSO", NumNeutrons_LYSO, "Per one:", NumNeutrons_LYSO / numScint);
    out("  Number in SiPM", NumNeutrons_SiPM, "Per one:", NumNeutrons_SiPM / numScint);
    out("  Energies for LYSO:");
    HistNeutronEn_LYSO->report();
    out("  Energies for SiPM:");
    HistNeutronEn_SiPM->report();
}

G4UserSteppingAction * ModeRadHard::getSteppingAction()
{
    return new SteppingAction_RadHard();
}
