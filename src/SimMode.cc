#include "SessionManager.hh"
#include "SimMode.hh"
#include "SteppingAction.hh"
#include "SensitiveDetectorScint.hh"
#include "PesGenerationMode.hh"
#include "Hist1D.hh"
#include "out.hh"
#include "jstools.hh"

#include "G4String.hh"
#include "G4RunManager.hh"

#include <iostream>
#include <sstream>
#include <fstream>

SimModeBase * SimModeFactory::makeSimModeInstance(const json11::Json & json)
{
    out("Reading simulation mode json");
    std::string Type;
    jstools::readString(json, "Type", Type);

    SimModeBase * sm = nullptr;

    if      (Type == "SimModeGui")            sm = new SimModeGui();
    else if (Type == "SimModeShowEvent")      sm = new SimModeShowEvent(0);
    else if (Type == "DoseExtractorMode")     sm = new DoseExtractorMode(0, {1,1,1}, {1,1,1}, {0,0,0}, "dummy.txt");
    else if (Type == "SimModeScintPosTest")   sm = new SimModeScintPosTest();
    else if (Type == "SimModeSingleEvents")   sm = new SimModeSingleEvents(0);
    else if (Type == "SimModeMultipleEvents") sm = new SimModeMultipleEvents(0, "dummy.txt", false);
    else if (Type == "SimModeTracing")        sm = new SimModeTracing();
    else if (Type == "SimModeAcollinTest")    sm = new SimModeAcollinTest(0, 0, 1, "dummy.txt");
    else if (Type == "SimModeAnnihilTest")    sm = new SimModeAnnihilTest(0, 0, "dummy.txt", false);
    else if (Type == "SimModeNatRadTest")     sm = new SimModeNatRadTest(0, 0, "dummy.txt");
    else if (Type == "SimModeFirstStage")     sm = new SimModeFirstStage(0, "dummy.txt", false);
    else if (Type == "PesGenerationMode")     sm = new PesGenerationMode(0, "dummy.txt", false);
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

SimModeGui::SimModeGui()
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void SimModeGui::run()
{
    SessionManager& SM = SessionManager::getInstance();
    SM.startGUI();
}

// ---

DoseExtractorMode::DoseExtractorMode(int numEvents, std::array<double, 3> binSize, std::array<int, 3> numBins, std::array<double, 3> origin, std::string fileName) :
    NumEvents(numEvents), BinSize(binSize), NumBins(numBins), Origin(origin)
{
    bNeedGui    = false;
    bNeedOutput = true;

    SessionManager & SM = SessionManager::getInstance();
    SM.FileName = fileName;
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

void DoseExtractorMode::run()
{
    SessionManager & SM = SessionManager::getInstance();
    SM.runManager->BeamOn(NumEvents);
    saveArray();
}

G4UserSteppingAction * DoseExtractorMode::getSteppingAction()
{
    return new SteppingAction_Dose();
}

void DoseExtractorMode::readFromJson(const json11::Json & json)
{
    SessionManager & SM = SessionManager::getInstance();
    jstools::readInt(json, "NumEvents", NumEvents);
    jstools::readString(json, "FileName", SM.FileName);

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
}

void DoseExtractorMode::writeBinningToJson(json11::Json::object & json) const
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

void DoseExtractorMode::doWriteToJson(json11::Json::object & json) const
{
    SessionManager & SM = SessionManager::getInstance();
    json["NumEvents"] = NumEvents;
    json["FileName"] = SM.FileName;

    writeBinningToJson(json);
}

void DoseExtractorMode::fill(double energy, const G4ThreeVector & pos, double density)
{
    const double densityKgPerMM3 = density/kg*mm3;
    if (densityKgPerMM3 < 1e-10) return; // vacuum

    std::array<int,3> index; // on stack, fast
    const bool ok = getVoxel(pos, index);
    if (!ok) return;

    const double deltaDose = (energy/joule) / ( densityKgPerMM3 * VoxelVolume );
    Dose[index[0]][index[1]][index[2]] += deltaDose;
}

bool DoseExtractorMode::getVoxel(const G4ThreeVector & pos, std::array<int,3> & index)
{
    for (int i = 0; i < 3; i++)
    {
        index[i] = floor( (pos[i] - Origin[i]) / BinSize[i] );
        if ( index[i] < 0 || index[i] >= NumBins[i]) return false;
    }
    return true;
}

void DoseExtractorMode::saveArray()
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

SimModeShowEvent::SimModeShowEvent(int EventToShow) : iEvent(EventToShow)
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void SimModeShowEvent::run()
{
    SessionManager& SM = SessionManager::getInstance();
    if (iEvent != 0) SM.runManager->BeamOn(iEvent);
    SimModeGui::run();
}

void SimModeShowEvent::readFromJson(const json11::Json & json)
{
    jstools::readInt(json, "iEvent", iEvent);
}

void SimModeShowEvent::doWriteToJson(json11::Json::object & json) const
{
    json["iEvent"] = iEvent;
}

// ---

SimModeScintPosTest::SimModeScintPosTest()
{
    bNeedGui    = false;
    bNeedOutput = false;
}

void SimModeScintPosTest::run()
{
    SessionManager& SM = SessionManager::getInstance();

    SM.runManager->BeamOn(10000);

    outFlush();
    if (Hits > 1) SumDelta /= Hits;
    out("\n---Test results---\nTotal hits of the scintillators:", Hits, "Max delta:", MaxDelta, " Average delta:", SumDelta, "\n\n");
}

G4UserSteppingAction * SimModeScintPosTest::getSteppingAction()
{
    return new SteppingAction_ScintPosTest;
}

// ---

SimModeSingleEvents::SimModeSingleEvents(int numEvents) : NumEvents(numEvents)
{
    bNeedGui    = false;
    bNeedOutput = true;

    SessionManager& SM = SessionManager::getInstance();
    SM.FileName = "Coincidence-GammaPairs-Test1.txt";
}

void SimModeSingleEvents::run()
{
    SessionManager& SM = SessionManager::getInstance();

    const double EnergyThreshold = 0.500*MeV;

    const int NumScint = SM.countScintillators();
    for(int i=0; i < NumScint; i++) ScintData.push_back({0,0,0});
    std::vector<int> hits;

    for (int iRun = 0; iRun < NumEvents; iRun++)
    {
        SM.runManager->BeamOn(1);

        hits.clear();
        for (int iScint = 0; iScint < NumScint; iScint++)
            if (ScintData[iScint][1] > EnergyThreshold)
                hits.push_back(iScint);

        if (hits.size() == 2)
        {
            for (int i = 0; i < 2; i++)
            {
                const int iScint = hits[i];
                const double Time   = ScintData[iScint][0] / ns;
                const double Energy = ScintData[iScint][1] / MeV;
                const G4ThreeVector & Pos   = SM.ScintRecords.at(iScint).FacePos;
                const double X = Pos[0] / mm;
                const double Y = Pos[1] / mm;
                const double Z = Pos[2] / mm;

                out("Scint#",iScint, Time,"ns ", Energy, "MeV  xyz: (",X,Y,Z,")  Run# ",iRun);

                if (SM.outStream)
                    *SM.outStream << X << " " << Y << " " << Z << " " << Time << " " << Energy << std::endl;

            }
            out("---");
        }

        for (int i = 0; i < NumScint; i++) ScintData[i] = {0,0,0};
    }

    outFlush();
    if (!SM.outStream) out("\nOutput stream was not created, nothing was saved");
    else out("Data saved to file:", SM.WorkingDirectory + "/" + SM.FileName);
}

G4VSensitiveDetector * SimModeSingleEvents::getScintDetector()
{
    return new SensitiveDetectorScint_SingleEvents("Scint");
}

void SimModeSingleEvents::readFromJson(const json11::Json & json)
{
    jstools::readInt(json, "NumEvents", NumEvents);
}

void SimModeSingleEvents::doWriteToJson(json11::Json::object & json) const
{
    json["NumEvents"] = NumEvents;
}

// ---

SimModeMultipleEvents::SimModeMultipleEvents(int numEvents, const std::string & fileName, bool binary,
                                             size_t maxCapacity, bool doCluster, double maxTimeDif) :
    NumEvents(numEvents), MaxCapacity(maxCapacity), bDoCluster(doCluster), MaxTimeDif(maxTimeDif)
{
    bNeedGui    = false;
    bNeedOutput = true;

    SessionManager & SM = SessionManager::getInstance();
    SM.bBinOutput  = binary;
    SM.FileName    = fileName;
}

void SimModeMultipleEvents::run()
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

G4VSensitiveDetector *SimModeMultipleEvents::getScintDetector()
{
    return new SensitiveDetectorScint_MultipleEvents("SD");
}

void SimModeMultipleEvents::saveData()
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

void SimModeMultipleEvents::readFromJson(const json11::Json &json)
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

void SimModeMultipleEvents::doWriteToJson(json11::Json::object &json) const
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

SimModeTracing::SimModeTracing()
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void SimModeTracing::run()
{
    SimModeGui::run();
}

G4UserSteppingAction * SimModeTracing::getSteppingAction()
{
    return new SteppingAction_Tracing;
}

// ---

SimModeAcollinTest::SimModeAcollinTest(int numRuns, double range, int numBins, const std::string & fileName) :
    NumRuns(numRuns), NumBins(numBins), FileName(fileName)
{
    From = 180.0 - range;
    init();
}

SimModeAcollinTest::~SimModeAcollinTest()
{
    delete Hist;
}

void SimModeAcollinTest::init()
{
    delete Hist; Hist = new Hist1D(NumBins, From, 180.0);
}

void SimModeAcollinTest::run()
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

G4UserSteppingAction *SimModeAcollinTest::getSteppingAction()
{
    return new SteppingAction_AcollinearityTester;
}

void SimModeAcollinTest::addDirection(const G4ThreeVector & v, int parentID, double energy)
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

void SimModeAcollinTest::readFromJson(const json11::Json &json)
{
    jstools::readInt   (json, "NumRuns",  NumRuns);
    jstools::readInt   (json, "From",     From);
    jstools::readInt   (json, "NumBins",  NumBins);
    jstools::readString(json, "FileName", FileName);

    init();
}

void SimModeAcollinTest::doWriteToJson(json11::Json::object &json) const
{
    json["NumRuns"]  = NumRuns;
    json["From"]     = From;
    json["NumBins"]  = NumBins;
    json["FileName"] = FileName;
}

// ---

SimModeAnnihilTest::SimModeAnnihilTest(int numEvents, double timeStart, const std::string & fileName, bool binary) :
    NumEvents(numEvents), TimeStart(timeStart)
{
    bNeedOutput = true;

    SessionManager & SM = SessionManager::getInstance();
    SM.FileName = fileName;
    SM.bBinOutput = binary;
}

void SimModeAnnihilTest::saveRecord(const G4ThreeVector & pos, double timeInSeconds)
{
    SessionManager & SM = SessionManager::getInstance();
    *SM.outStream << pos[0] << ' ' << pos[1] << ' ' << pos[2] << ' ' << timeInSeconds << '\n';
}

G4UserSteppingAction * SimModeAnnihilTest::getSteppingAction()
{
    return new SteppingAction_AnnihilationTester;
}

void SimModeAnnihilTest::run()
{
    SessionManager & SM = SessionManager::getInstance();
    SM.runManager->BeamOn(NumEvents);
}

void SimModeAnnihilTest::readFromJson(const json11::Json &json)
{
    SessionManager & SM = SessionManager::getInstance();

    jstools::readInt   (json, "NumEvents", NumEvents);
    jstools::readDouble(json, "TimeStart", TimeStart);
    jstools::readString(json, "FileName",  SM.FileName);
    jstools::readBool  (json, "Binary",    SM.bBinOutput);
}

void SimModeAnnihilTest::doWriteToJson(json11::Json::object &json) const
{
    SessionManager & SM = SessionManager::getInstance();

    json["NumEvents"] = NumEvents;
    json["TimeStart"] = TimeStart;
    json["FileName"]  = SM.FileName;
    json["Binary"]    = SM.bBinOutput;
}

// ---

SimModeNatRadTest::SimModeNatRadTest(int numEvents, int numBins, const std::string & fileName) :
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

void SimModeNatRadTest::init()
{
    delete Hist; Hist = new Hist1D(NumBins, 0, 1.3);
}

SimModeNatRadTest::~SimModeNatRadTest()
{
    delete Hist;
}

void SimModeNatRadTest::run()
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

G4UserSteppingAction * SimModeNatRadTest::getSteppingAction()
{
    return new SteppingAction_NatRadTester;
}

void SimModeNatRadTest::addEnergy(int iScint, double energy)
{
    Deposition[iScint] += energy;
}

void SimModeNatRadTest::readFromJson(const json11::Json &json)
{
    jstools::readInt   (json, "NumEvents", NumEvents);
    jstools::readInt   (json, "NumBins",   NumBins);
    jstools::readString(json, "FileName",  FileName);

    init();
}

void SimModeNatRadTest::doWriteToJson(json11::Json::object &json) const
{
    json["NumEvents"] = NumEvents;
    json["NumBins"]   = NumBins;
    json["FileName"]  = FileName;
}

// ---

SimModeFirstStage::SimModeFirstStage(int numEvents, const std::string & fileName, bool bBinary) :
    NumEvents(numEvents)
{
    bNeedOutput = true;

    SessionManager & SM = SessionManager::getInstance();
    SM.FileName   = fileName;
    SM.bBinOutput = bBinary;
}

void SimModeFirstStage::run()
{
    CurrentEvent = 0;
    SessionManager & SM = SessionManager::getInstance();
    SM.runManager->BeamOn(NumEvents);
}

void SimModeFirstStage::saveParticle(const G4String & particle, double energy_keV, double * PosDir, double time)
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

void SimModeFirstStage::onEventStarted()
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

void SimModeFirstStage::readFromJson(const json11::Json &json)
{
    SessionManager & SM = SessionManager::getInstance();

    jstools::readInt   (json, "NumEvents", NumEvents);
    jstools::readString(json, "FileName",  SM.FileName);
    jstools::readBool  (json, "bBinary",   SM.bBinOutput);
}

void SimModeFirstStage::doWriteToJson(json11::Json::object &json) const
{
    SessionManager & SM = SessionManager::getInstance();

    json["NumEvents"] = NumEvents;
    json["FileName"]  = SM.FileName;
    json["bBinary"]   = SM.bBinOutput;
}

// ---

PesAnalyzerMode::PesAnalyzerMode(int numEvents, const std::string & fileName) :
    NumEvents(numEvents)
{
    SessionManager & SM = SessionManager::getInstance();

    SM.FileName = fileName;
    bNeedOutput = true;
}

void PesAnalyzerMode::onNewPrimary()
{
    numPrim++;

    if (bPositron)
    {
        std::string products;
        for (const G4String & p : Target.second)
        {
            if (p == "gamma") continue;
            products += p + " ";
        }

        std::string key = Target.first ;

        if (Decays.size() == 1)
        {
            std::string secondaries;
            for (const G4String & p : Decays.front().second)
            {
                if (p == "gamma") continue;
                secondaries += p + " ";
            }
            key += " -> (" + products + ")";
        }
        else
        {
            //out("=> ",Target.first, " -> (" + products + ")");
            key += "->(" + products + ")";
            for (const auto & pair : Decays)
            {
                std::string secondaries;
                for (const G4String & p : pair.second)
                {
                    if (p == "gamma") continue;
                    secondaries += p + " ";
                }
                //out("  " + pair.first, " -> (" + secondaries + ")");
                key += " + " + pair.first + "->(" + secondaries + ")";
            }
        }

        auto it = Statistics.find(key);
        if (it == Statistics.end()) Statistics[key] = 1;
        else Statistics[key]++;


    }

    Target = {"None", {}};
    Decays.clear();
    bPositron = false;
}

void PesAnalyzerMode::registerTarget(const G4String & target, std::vector<G4String> products)
{
    std::sort(products.begin(), products.end());
    Target = {target, products};
}

void PesAnalyzerMode::onDecay(const G4String & isotope, std::vector<G4String> products)
{
    std::sort(products.begin(), products.end());
    Decays.push_back({isotope, products});

    for (const G4String & p : products)
        if (p == "e+") bPositron = true;
}

void PesAnalyzerMode::run()
{
    SessionManager & SM = SessionManager::getInstance();

    SM.runManager->BeamOn(NumEvents);
    onNewPrimary(); // to trigger update for the last primary

    for (auto it : Statistics)
    {
        out(it.first, it.second);
    }
    out("Num primaries:", --numPrim);
}

G4UserSteppingAction * PesAnalyzerMode::getSteppingAction()
{
    return new SteppingAction_PesAnalyzer();
}

// ---

DepoStatMode::DepoStatMode(int numEvents, const std::string & fileName) :
    SimModeBase(), NumEvents(numEvents)
{
    SessionManager & SM = SessionManager::getInstance();

    SM.FileName = fileName;
    bNeedOutput = true;

    //SM.EvNumberInterval = 1;
}

void DepoStatMode::run()
{
    SessionManager & SM = SessionManager::getInstance();
    SM.runManager->BeamOn(NumEvents);

    processEventData(); // for the last event!

    // results
    out("\n\n------------------");
    int remains = NumEvents;
    out("No interaction:", 100.0*numNothing/NumEvents, "%");
        remains -= numNothing;
    out("Single scintillator: ", 100.0*numSingle/NumEvents, "%");
    out("  --> Within 5% / 10% of 511 keV:", 100.0*numSingle_5percent/NumEvents, "% /",100.0*numSingle_10percent/NumEvents,"%");
        remains -= numSingle;
    out("Two scintillators: ", 100.0*numTwo/NumEvents, "%", "  Same assembly:", 100.0*numTwo_SameAssembly/NumEvents, "%");
        remains -= numTwo;
    out("  Two within the same assembly:");
    out("  --> First  within 5% / 10% of 511 keV:", 100.0*numTwo_Same_first_5percent/NumEvents, "% /",100.0*numTwo_Same_first_10percent/NumEvents,"%");
    out("  --> Second within 5% / 10% of 511 keV:", 100.0*numTwo_Same_second_5percent/NumEvents, "% /",100.0*numTwo_Same_second_10percent/NumEvents,"%",
        "average distance between scints:", averageScintDistSecond5/numTwo_Same_second_5percent, "and", averageScintDistSecond10/numTwo_Same_second_10percent, "mm");
    out("  --> Sum within 5% / 10% of 511 keV:", 100.0*numTwo_Same_sum_5percent/NumEvents, "% /",100.0*numTwo_Same_sum_10percent/NumEvents,"%",
        "average distance between scints:", averageScintDistSum5/numTwo_Same_sum_5percent, "and", averageScintDistSum10/numTwo_Same_sum_10percent, "mm");
    out("     --> For 5%,  fraction when first deposition is smaller:", (double)two_FirstSmaller5 /numTwo_Same_sum_5percent,  "Average first/second ratio:",two_AvRatioFirstSecond5/numTwo_Same_sum_5percent);
    out("     --> For 10%, fraction when first deposition is smaller:", (double)two_FirstSmaller10/numTwo_Same_sum_10percent, "Average first/second ratio:",two_AvRatioFirstSecond5/numTwo_Same_sum_10percent);
    out("  Two within different assemblies:");
    out("  --> First  within 5% / 10% of 511 keV:", 100.0*numTwo_Dif_first_5percent/NumEvents, "% /",100.0*numTwo_Dif_first_10percent/NumEvents,"%");
    out("  --> Second within 5% / 10% of 511 keV:", 100.0*numTwo_Dif_second_5percent/NumEvents, "% /",100.0*numTwo_Dif_second_10percent/NumEvents,"%");

    //"  Within 5% / 10% of 511 keV:", 100.0*numSingle_5percent/NumEvents, "% /",100.0*numSingle_10percent/NumEvents,"%");


    out("Unclassified:", 100.0*remains/NumEvents, "%");
}

void DepoStatMode::onEventStarted()
{
    processEventData(); //previous event!

    EventRecord.clear();
}

G4UserSteppingAction * DepoStatMode::getSteppingAction()
{
    return new SteppingAction_DepoStatMode();
}

void DepoStatMode::addRecord(int iScint, double depo, double time)
{
    //out("-->Scint:", iScint, " MeV:", depo, "  ns:", time);

    for (DepoStatRec & r : EventRecord)
    {
        if (iScint == r.iScint)
        {
            r.energy += depo;
            if (time < r.time) r.time = time;
            return;
        }
    }
    EventRecord.emplace_back(DepoStatRec{iScint, depo, time});
}

void DepoStatMode::processEventData()
{
    SessionManager & SM = SessionManager::getInstance();
    //for (const DepoStatRec & r : EventRecord)
    //    out("Scint:", r.iScint, " MeV:", r.energy, "  ns:", r.time);

    if      (EventRecord.empty())     numNothing++;
    else if (EventRecord.size() == 1)
    {
        numSingle++;
        const double depo = EventRecord.front().energy;
        if (depo > 0.511*0.95 && depo < 0.511*1.05) numSingle_5percent++;
        if (depo > 0.511*0.90 && depo < 0.511*1.10) numSingle_10percent++;
    }
    else if (EventRecord.size() == 2)
    {
        numTwo++;

        const DepoStatRec & rFirst  = EventRecord.front();
        const DepoStatRec & rSecond = EventRecord.back();

        bool bSameAssembly = ( SM.ScintRecords[rFirst.iScint].AssemblyNumber == SM.ScintRecords[rSecond.iScint].AssemblyNumber );
        if (bSameAssembly)
        {
            numTwo_SameAssembly++;
            const double depo1 = rFirst.energy;
            if (depo1 > 0.511*0.95 && depo1 < 0.511*1.05) numTwo_Same_first_5percent++;
            if (depo1 > 0.511*0.90 && depo1 < 0.511*1.10) numTwo_Same_first_10percent++;

            const double depo2 = rSecond.energy;
            if (depo2 > 0.511*0.95 && depo2 < 0.511*1.05)
            {
                numTwo_Same_second_5percent++;
                G4ThreeVector d = SM.ScintRecords[rFirst.iScint].FacePos - SM.ScintRecords[rSecond.iScint].FacePos;
                averageScintDistSecond5 += d.getR();
            }
            if (depo2 > 0.511*0.90 && depo2 < 0.511*1.10)
            {
                numTwo_Same_second_10percent++;
                G4ThreeVector d = SM.ScintRecords[rFirst.iScint].FacePos - SM.ScintRecords[rSecond.iScint].FacePos;
                averageScintDistSecond10 += d.getR();
            }

            const double sumdepo = depo1 + depo2;
            if (sumdepo > 0.511*0.95 && sumdepo < 0.511*1.05)
            {
                numTwo_Same_sum_5percent++;
                G4ThreeVector d = SM.ScintRecords[rFirst.iScint].FacePos - SM.ScintRecords[rSecond.iScint].FacePos;
                averageScintDistSum5 += d.getR();

                if (depo1 < depo2) two_FirstSmaller5++;
                two_AvRatioFirstSecond5 += depo1/depo2;
            }
            if (sumdepo > 0.511*0.90 && sumdepo < 0.511*1.10)
            {
                numTwo_Same_sum_10percent++;
                G4ThreeVector d = SM.ScintRecords[rFirst.iScint].FacePos - SM.ScintRecords[rSecond.iScint].FacePos;
                averageScintDistSum10 += d.getR();

                if (depo1 < depo2) two_FirstSmaller10++;
                two_AvRatioFirstSecond10 += depo1/depo2;
            }
        }
        else
        {
            const double depo1 = rFirst.energy;
            if (depo1 > 0.511*0.95 && depo1 < 0.511*1.05) numTwo_Dif_first_5percent++;
            if (depo1 > 0.511*0.90 && depo1 < 0.511*1.10) numTwo_Dif_first_10percent++;
            const double depo2 = rSecond.energy;
            if (depo2 > 0.511*0.95 && depo2 < 0.511*1.05) numTwo_Dif_second_5percent++;
            if (depo2 > 0.511*0.90 && depo2 < 0.511*1.10) numTwo_Dif_second_10percent++;
        }
    }
    else if (EventRecord.size() == 3)
    {

    }
}
