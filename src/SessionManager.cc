#include "SessionManager.hh"
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
#include "AcollinearGammaModel.hh"
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

    if (!SourceMode)
    {
        out("Source mode not provided!");
        exit(1);
    }
    if (!SimMode)
    {
        out("Simulation mode not provided!");
        exit(2);
    }
    if (!PhantomMode)
    {
        out("Phantom mode not provided!");
        exit(3);
    }

    DetectorConstruction * theDetector = new DetectorConstruction();
    runManager->SetUserInitialization(theDetector);

    G4VModularPhysicsList* physicsList = new QGSP_BIC_HP;
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    physicsList->SetDefaultCutValue(0.1*mm);  // see createPhntomRegion and createScintRegion for specific cuts!
    if (bSimAcollinearity)
    {
        G4FastSimulationPhysics * fsm = new G4FastSimulationPhysics();
        // see registerAcollinearGammaModel()
        //https://indico.cern.ch/event/789510/contributions/3297180/attachments/1817759/2973421/G4Tutorial_fastSim_vFin.pdf
        fsm->ActivateFastSimulation("gamma");
        physicsList->RegisterPhysics(fsm);
    }
    runManager->SetUserInitialization(physicsList);

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

void SessionManager::registerAcollinearGammaModel(G4Region * region)
{
    AcollinearGammaModel * mod = new AcollinearGammaModel("AcollinearGammas", region);
    G4AutoDelete::Register(mod);
}

void SessionManager::createPhantomRegion(G4LogicalVolume * logVolPhantom)
{
    regPhantom = new G4Region("Phantom");
    regPhantom->AddRootLogicalVolume(logVolPhantom);

    if (bSimAcollinearity) registerAcollinearGammaModel(regPhantom);

    G4ProductionCuts * cuts = new G4ProductionCuts();
    cuts->SetProductionCut(0.5*mm, G4ProductionCuts::GetIndex("gamma"));
    cuts->SetProductionCut(1.0*mm, G4ProductionCuts::GetIndex("e-"));
    regPhantom->SetProductionCuts(cuts);
}

void SessionManager::createScintillatorRegion(G4LogicalVolume * logVolScint)
{
    regScint = new G4Region("Scintillators");
    regScint->AddRootLogicalVolume(logVolScint);

    G4ProductionCuts * cuts = new G4ProductionCuts();
    cuts->SetProductionCut(0.1*mm, G4ProductionCuts::GetIndex("gamma"));
    cuts->SetProductionCut(0.1*mm, G4ProductionCuts::GetIndex("e-"));
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
