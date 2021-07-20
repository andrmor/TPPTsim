#include "SessionManager.hh"
#include "PhantomMode.hh"
#include "SourceMode.hh"
#include "SimMode.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "ScintRecord.hh"
#include "out.hh"

#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4VisManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4StepLimiterPhysics.hh"
#include "QGSP_BIC_HP.hh"
#include "G4UImanager.hh"
#include "G4RandomTools.hh"
#include "G4String.hh"
#include "G4FastSimulationPhysics.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Gamma.hh"
#include "G4RegionStore.hh"
#include "G4AutoDelete.hh"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

SessionManager & SessionManager::getInstance()
{
    static SessionManager instance; // Guaranteed to be destroyed, instantiated on first use.
    return instance;
}

SessionManager::SessionManager()
{
    runManager = new G4RunManager;
}

void SessionManager::startSession(int argc, char ** argv)
{
    out("\n\n---------");

    if (!isDirExists(WorkingDirectory))
    {
        out("Provided working directory does not exist:\n", WorkingDirectory);
        exit(1);
    }
    if (!SourceMode)
    {
        out("Source mode not provided!");
        exit(2);
    }
    if (!SimMode)
    {
        out("Simulation mode not provided!");
        exit(3);
    }
    if (!PhantomMode)
    {
        out("Phantom mode not provided!");
        exit(4);
    }

    DetectorConstruction * theDetector = new DetectorConstruction();
    runManager->SetUserInitialization(theDetector);

    G4VModularPhysicsList* physicsList = new QGSP_BIC_HP;
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    physicsList->SetDefaultCutValue(0.1*mm);  // see createPhantomRegion and createScintRegion for specific cuts!
    runManager->SetUserInitialization(physicsList);
    if (bSimAcollinearity || bKillNeutrinos) createFastSimulationPhysics(physicsList);

    runManager->SetUserAction(new PrimaryGeneratorAction); // SourceMode cannot be directly inherited from G4VUserPrimaryGeneratorAction due to initialization order

    G4UserSteppingAction * StAct = SimMode->getSteppingAction();
    if (StAct) runManager->SetUserAction(StAct);

    runManager->SetUserAction(new EventAction);

    // ---
    runManager->Initialize();
    // ---

    GammaPD = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    configureRandomGenerator();
    initializeSource();
    if (SimMode->bNeedGui)    configureGUI(argc, argv);
    if (SimMode->bNeedOutput) configureOutput();
    configureVerbosity();

    saveConfig(WorkingDirectory + "/Sim-lastRunConfig.json");
}

SessionManager::~SessionManager() {}

void SessionManager::configureGUI(int argc, char ** argv)
{
    ui         = new G4UIExecutive(argc, argv);
    visManager = new G4VisExecutive("Quiet");

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/hits/verbose 2");
    UImanager->ApplyCommand("/tracking/verbose 2");
    UImanager->ApplyCommand("/control/saveHistory");

    if ( detectorContains(DetComp::GDML) ) scanMaterials();
}

void SessionManager::scanMaterials()
{
    out("-->Scanning materials...");

    std::vector<G4LogicalVolume*> * lvs = G4LogicalVolumeStore::GetInstance();
    for (G4LogicalVolume * lv : *lvs)
    {
        G4Material * mat = lv->GetMaterial();
        out(lv->GetName(), mat->GetName(), mat->GetChemicalFormula(), mat->GetDensity());

        if (mat->GetName() == "G4_Al") lv->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 0, 1.0)));
    }

    out("<--Material scan completed");
}

void SessionManager::createFastSimulationPhysics(G4VModularPhysicsList * physicsList)
{
    G4FastSimulationPhysics * fsm = new G4FastSimulationPhysics();
    //https://indico.cern.ch/event/789510/contributions/3297180/attachments/1817759/2973421/G4Tutorial_fastSim_vFin.pdf
    // see registerAcollinearGammaModel() and registerParticleKillerModel()

    if (bSimAcollinearity) fsm->ActivateFastSimulation("gamma");
    if (bKillNeutrinos)    fsm->ActivateFastSimulation("nu_e");
    physicsList->RegisterPhysics(fsm);
}

#include "AcollinearGammaModel.hh"
void SessionManager::registerAcollinearGammaModel(G4Region * region)
{
    AcollinearGammaModel * mod = new AcollinearGammaModel("AcollinearGammas", region);
    G4AutoDelete::Register(mod);
}

#include "ParticleKiller.hh"
void SessionManager::registerParticleKillerModel(G4Region *region)
{
    ParticleKillerModel * mod = new ParticleKillerModel("ParticleKiller", region);
    G4AutoDelete::Register(mod);
}

void SessionManager::createPhantomRegion(G4LogicalVolume * logVolPhantom)
{
    regPhantom = new G4Region("Phantom");
    regPhantom->AddRootLogicalVolume(logVolPhantom);

    if (bSimAcollinearity) registerAcollinearGammaModel(regPhantom);
    if (bKillNeutrinos)    registerParticleKillerModel(regPhantom);

    G4ProductionCuts * cuts = new G4ProductionCuts();
    cuts->SetProductionCut(CutPhantomGamma,    G4ProductionCuts::GetIndex("gamma"));
    cuts->SetProductionCut(CutPhantomElectron, G4ProductionCuts::GetIndex("e-"));
    cuts->SetProductionCut(CutPhantomPositron, G4ProductionCuts::GetIndex("e+"));
    regPhantom->SetProductionCuts(cuts);
}

void SessionManager::createScintillatorRegion(G4LogicalVolume * logVolScint)
{
    regScint = new G4Region("Scintillators");
    regScint->AddRootLogicalVolume(logVolScint);

    G4ProductionCuts * cuts = new G4ProductionCuts();
    cuts->SetProductionCut(CutScintGamma,    G4ProductionCuts::GetIndex("gamma"));
    cuts->SetProductionCut(CutScintElectron, G4ProductionCuts::GetIndex("e-"));
    cuts->SetProductionCut(CutScintPositron, G4ProductionCuts::GetIndex("e+"));
    regScint->SetProductionCuts(cuts);
}

int SessionManager::countScintillators() const
{
    return NumScintX * NumScintY * NumRows * NumSegments * 2.0;
}

int SessionManager::getNumberNatRadEvents(double timeFromInNs, double timeToInNs) const
{
    double volume_cm3 = ScintSizeX/cm * ScintSizeY/cm * ScintSizeZ/cm;
    out("---Volume of one scintillator =", volume_cm3, "cm3");
    volume_cm3 *= countScintillators();
    out("---Total LYSO volume = ", volume_cm3, "cm3");
    double decaysPerSecond = activityLYSO * volume_cm3;
    out("---Decays per second = ", decaysPerSecond);

    const int numEvents = decaysPerSecond * 1e-9 * (timeToInNs - timeFromInNs);
    out("---Decays for time range from", timeFromInNs, "ns to", timeToInNs, "ns =", numEvents);
    return numEvents;
}

bool SessionManager::detectorContains(DetComp component) const
{
    return std::count(DetectorComposition.begin(), DetectorComposition.end(), component); // pre-c++20 ugly version of "contains"
}

void SessionManager::saveScintillatorTable(const std::string & fileName)
{
    std::ofstream stream;
    stream.open(fileName);
    if (!stream.is_open() || stream.fail() || stream.bad())
    {
        std::cout  << "Cannot open file to store scintillator data: " << fileName << std::endl;
        return;
    }

    for (ScintRecord & rec : ScintRecords) rec.write(stream);

    stream.close();
}

void SessionManager::startGUI()
{
    visManager->Initialize();
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
}

void SessionManager::configureOutput()
{
    outStream = new std::ofstream();

    std::string fullFileName = WorkingDirectory + "/" + FileName;
    if (bBinOutput) outStream->open(fullFileName, std::ios::out | std::ios::binary);
    else outStream->open(fullFileName);

    if (!outStream->is_open() || outStream->fail() || outStream->bad())
    {
        out("Cannot open file to store output data!");
        outFlush();
        exit(1);
        delete outStream; outStream = nullptr;
    }
    else out("\nSaving output to file", fullFileName);
}

void SessionManager::configureRandomGenerator()
{
    randGen = new CLHEP::RanecuEngine();
    randGen->setSeed(Seed);
    G4Random::setTheEngine(randGen);
}

void SessionManager::initializeSource()
{
    SourceMode->initialize();
}

void SessionManager::configureVerbosity()
{
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    if (bG4Verbose)
    {
        UImanager->ApplyCommand("/hits/verbose 2");
        UImanager->ApplyCommand("/tracking/verbose 2");
        UImanager->ApplyCommand("/control/saveHistory");
    }
    else
    {
        UImanager->ApplyCommand("/hits/verbose 0");
        UImanager->ApplyCommand("/tracking/verbose 0");
        UImanager->ApplyCommand("/control/verbose 0");
        UImanager->ApplyCommand("/run/verbose 0");
    }
    //UImanager->ApplyCommand("/run/initialize");
}

void SessionManager::endSession()
{
    if (outStream) outStream->close();
    delete outStream;

    delete visManager;
    delete runManager;
    delete ui;
}

#include <sys/types.h>
#include <sys/stat.h>

int SessionManager::isDirExists(const std::string & dirName)
{
    struct stat info;

    if (stat(dirName.data(), &info) != 0) return false;
    else if (info.st_mode & S_IFDIR)      return true;
    else                                  return false;
}

#include "json11.hh"
void SessionManager::saveConfig(const std::string & fileName) const
{
    json11::Json::object json;

    json["Seed"] = Seed;

    json["bSimAcollinearity"] = bSimAcollinearity;
    json["bKillNeutrinos"]    = bKillNeutrinos;

    json11::Json::object jsCuts;
        jsCuts["CutPhantomGamma"]    = CutPhantomGamma;
        jsCuts["CutPhantomElectron"] = CutPhantomElectron;
        jsCuts["CutPhantomPositron"] = CutPhantomPositron;
        jsCuts["CutScintGamma"]      = CutScintGamma;
        jsCuts["CutScintElectron"]   = CutScintElectron;
        jsCuts["CutScintPositron"]   = CutScintPositron;
    json["Cuts"] = jsCuts;

    json["WorkingDirectory"] = WorkingDirectory;

    json["bG4Verbose"] = bG4Verbose;
    json["bDebug"]     = bDebug;

    json["bShowEventNumber"] = bShowEventNumber;
    json["EvNumberInterval"] = EvNumberInterval;

    //Phantom mode
    {
        json11::Json::object js;
        PhantomMode->writeToJson(js);
        json["PhantomMode"] = js;
    }

    // Detector composition
    {
        json11::Json::array ar;
        for (const DetComp & el : DetectorComposition) ar.push_back(static_cast<int>(el));
        json["DetectorComposition"] = ar;
    }

    // Source
    {
        json11::Json::object js;
        SourceMode->writeToJson(js);
        json["SourceMode"] = js;
    }

    // Simulation mode
    {
        json11::Json::object js;
        SimMode->writeToJson(js);
        json["SimMode"] = js;
    }

    std::string json_str = json11::Json(json).dump();
    std::ofstream confStream;
    confStream.open(fileName);
    if (confStream.is_open())
        confStream << json_str << std::endl;
    confStream.close();
}

void assertKey(const json11::Json & json, const std::string & key)
{
    if (json.object_items().count(key) == 0)
    {
        out("Config json does not contain required key:", key);
        exit(1);
    }
}

void SessionManager::loadConfig(const std::string & fileName)
{
    out("Reading config file:", fileName);

    std::ifstream in(fileName);
    std::stringstream sstr;
    sstr << in.rdbuf();
    std::string cs = sstr.str();

    std::string err;
    json11::Json json = json11::Json::parse(cs, err);
    if (!err.empty())
    {
        out(err);
        exit(1);
    }

    assertKey(json, "Seed");
    Seed = json["Seed"].int_value();
    out("Seed:", Seed);

    assertKey(json, "bSimAcollinearity");
    bSimAcollinearity = json["bSimAcollinearity"].bool_value();
    out("bSimAcollinearity:", bSimAcollinearity);

    assertKey(json, "bKillNeutrinos");
    bKillNeutrinos = json["bKillNeutrinos"].bool_value();
    out("bKillNeutrinos:", bKillNeutrinos);

    assertKey(json, "Cuts");
    json11::Json::object jsCuts = json["Cuts"].object_items();
    {
        assertKey(json, "CutPhantomGamma");
        CutPhantomGamma = json["CutPhantomGamma"].number_value();
        out("CutPhantomGamma:", CutPhantomGamma);
        assertKey(json, "CutPhantomElectron");
        CutPhantomElectron = json["CutPhantomElectron"].number_value();
        out("CutPhantomElectron:", CutPhantomElectron);
        assertKey(json, "aaa");
        CutPhantomPositron = json["CutPhantomPositron"].number_value();
        out("CutPhantomPositron:", CutPhantomPositron);
        assertKey(json, "CutPhantomPositron");
        CutScintGamma = json["CutScintGamma"].number_value();
        out("CutScintGamma:", CutScintGamma);
        assertKey(json, "CutScintGamma");
        CutScintElectron = json["CutScintElectron"].number_value();
        out("CutScintElectron:", CutScintElectron);
        assertKey(json, "CutScintPositron");
        CutScintPositron = json["CutScintPositron"].number_value();
        out("CutScintPositron:", CutScintPositron);
    }

    assertKey(json, "WorkingDirectory");
    WorkingDirectory = json["WorkingDirectory"].string_value();
    out("WorkingDirectory:", WorkingDirectory);

    assertKey(json, "bG4Verbose");
    bG4Verbose = json["bG4Verbose"].bool_value();
    out("bG4Verbose:", bG4Verbose);

    assertKey(json, "bDebug");
    bDebug = json["bDebug"].bool_value();
    out("bDebug:", bDebug);

    assertKey(json, "bShowEventNumber");
    bShowEventNumber = json["bShowEventNumber"].bool_value();
    out("bShowEventNumber:", bShowEventNumber);

    assertKey(json, "EvNumberInterval");
    EvNumberInterval = json["EvNumberInterval"].int_value();
    out("EvNumberInterval:", EvNumberInterval);

    //Phantom mode
    {
        assertKey(json, "PhantomMode");
        json11::Json::object js = json["PhantomMode"].object_items();
        //PhantomMode->writeToJson(js);
    }

    // Detector composition
    {
        assertKey(json, "DetectorComposition");
        json11::Json::array ar = json["DetectorComposition"].array_items();
        DetectorComposition.clear();

        out("Detector composition items:");
        for (size_t i = 0; i < ar.size(); i++)
        {
            // TODO: solve possible problem: index is larger than the defined enum! !!!***
            const json11::Json & arEl = ar[i];
            int ci = arEl.int_value();
            out("-->", ci);
            DetComp el = static_cast<DetComp>(ci);
            DetectorComposition.push_back(el);
        }
        out("Detector composition items end.");
    }

    // Source
    {
        assertKey(json, "SourceMode");
        json11::Json::object js = json["SourceMode"].object_items();
        //SourceMode->writeToJson(js);
    }

    // Simulation mode
    {
        assertKey(json, "SimMode");
        json11::Json::object js = json["SimMode"].object_items();
        //SimMode->writeToJson(js);
    }

    out("Load success!");
}
