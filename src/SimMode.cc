#include "SessionManager.hh"
#include "SimMode.hh"
#include "SteppingAction.hh"
#include "SensitiveDetectorScint.hh"
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
    else if (Type == "SimModeScintPosTest")   sm = new SimModeScintPosTest();
    else if (Type == "SimModeSingleEvents")   sm = new SimModeSingleEvents(0);
    else if (Type == "SimModeMultipleEvents") sm = new SimModeMultipleEvents(0, "dummy.txt", false);
    else if (Type == "SimModeTracing")        sm = new SimModeTracing();
    else if (Type == "SimModeAcollinTest")    sm = new SimModeAcollinTest(0, 0, 1, "dummy.txt");
    else if (Type == "SimModeAnnihilTest")    sm = new SimModeAnnihilTest(0, 0, 1, "dummy.txt");
    else if (Type == "SimModeNatRadTest")     sm = new SimModeNatRadTest(0, 0, "dummy.txt");
    else if (Type == "SimModeFirstStage")     sm = new SimModeFirstStage(0, "dummy.txt", false);
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

SimModeAnnihilTest::SimModeAnnihilTest(int numEvents, double range, int numBins, const std::string & fileName) :
    NumEvents(numEvents), Range(range), NumBins(numBins), FileName(fileName)
{
    init();
}

void SimModeAnnihilTest::init()
{
    delete Hist; Hist = new Hist1D(NumBins, -Range, Range);
}

SimModeAnnihilTest::~SimModeAnnihilTest()
{
    delete Hist;
}

G4UserSteppingAction * SimModeAnnihilTest::getSteppingAction()
{
    return new SteppingAction_AnnihilationTester;
}

void SimModeAnnihilTest::run()
{
    SessionManager & SM = SessionManager::getInstance();
    SM.runManager->BeamOn(NumEvents);

    outFlush();
    out("\nDistribution of annihilation positions:");
    Hist->report();
    Hist->save(SM.WorkingDirectory + '/' + FileName);
}

void SimModeAnnihilTest::addPosition(double x)
{
    Hist->fill(x);
}

void SimModeAnnihilTest::readFromJson(const json11::Json &json)
{
    jstools::readInt   (json, "NumEvents", NumEvents);
    jstools::readDouble(json, "Range",     Range);
    jstools::readInt   (json, "NumBins",   NumBins);
    jstools::readString(json, "FileName",  FileName);

    init();
}

void SimModeAnnihilTest::doWriteToJson(json11::Json::object &json) const
{
    json["NumEvents"] = NumEvents;
    json["Range"]     = Range;
    json["NumBins"]   = NumBins;
    json["FileName"]  = FileName;
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

        *SM.outStream << ss.rdbuf() << std::endl;
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
        *SM.outStream << '#' << CurrentEvent << std::endl;

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
