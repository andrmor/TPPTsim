#include "SessionManager.hh"
#include "SourceMode.hh"
#include "SourceParticleListFile.hh"
#include "TimeGenerator.hh"
#include "DefinedParticles.hh"
#include "SimMode.hh"
#include "PhantomMode.hh"
#include "DicomPhantom.hh"
#include "DetComp.hh"
#include "out.hh"
#include "ModePesGenerator_MC.hh"
#include "ModeActivityGenerator.hh"
#include "ActivityProfiler.hh"
#include "ModeAnnihilationLogger.hh"
#include "SourcePesHistogramFiles.hh"
#include "SourceAnnihilHistFile.hh"
#include "SourcePositronium.hh"
#include "Positronium.hh"

#include <string>
#include <chrono>

int main(int argc, char** argv)
{
    SessionManager & SM = SessionManager::getInstance(); // side effect: outputs Geant4 version

    std::string filename;//  = "/home/andr/WORK/tmp/SimConfigBox.json"; //here you can directly provide the config file name
    // WARNING: the filename can be overriden with a command line arguments, e.g. sim -f /path/filename.json

    SM.parseRunArguments(argc, argv, filename); // checks for override of the random generator seed and config file name

    if (!filename.empty())
    {
        SM.loadConfig(filename);
    }
    else
    {
     // --- START of user init ---

        // General settings
        SM.Seed = 100;                // WARNING: the seed can be overriden with a command line argument, e.g. sim -s 123456
        SM.SimAcollinearity   = true;  // only for the phantom region!
        SM.KillNeutrinos      = true;
        SM.UseStepLimiter     = true; SM.PhantomStepLimit = 0.25*mm;
        //SM.UseStepLimiter     = true; SM.PhantomStepLimit = 0.325*mm;

        SM.CutPhantomGamma    = 10.0*mm; SM.CutPhantomElectron = 10.0*mm; SM.CutPhantomPositron = 0.1*mm;
        SM.CutScintGamma      =  0.1*mm; SM.CutScintElectron   = 0.1*mm;  SM.CutScintPositron   = 0.1*mm;

        SM.Verbose          = false;
        SM.ShowEventNumber  = true; SM.EvNumberInterval = 10000;

        SM.WorkingDirectory  = "/home/andr/WORK/TPPT/tmp1";

        // Phantom
        SM.Phantom = new PhantomNone;
        //SM.Phantom = new PhantomCylinder(25.4, 100.0, "G4_PLEXIGLASS");
        //SM.Phantom = new PhantomCylinder(25.4, 100.0, "G4_POLYETHYLENE");
        //SM.Phantom = new PhantomMarekWater();
        //SM.Phantom = new PhantomCylinder(6.35, 20.0, "G4_Cu");
        //SM.Phantom = new PhantomCylinder(6.35, 20.0, EMaterial::Ni400);

        // Detector components
        SM.DetectorComposition.add(DetComp::Scintillators);
        //SM.DetectorComposition.add(DetComp::DetComp::SIPM);
        //SM.DetectorComposition.add({DetComp::Base, DetComp::ClosedStructure, DetComp::SIPM, DetComp::PCB, DetComp::CopperStructure, DetComp::CoolingAssemblies});
        //SM.DetectorComposition.add({DetComp::Base, DetComp::SIPM, DetComp::PCB, DetComp::CopperStructure, DetComp::CoolingAssemblies});
        //SM.DetectorComposition.add(DetComp::ParticleLogger); // required only for ModeParticleLogger

        // Source
        //SM.SourceMode = new SourcePoint(new Positronium(1.0, "gammas.txt"), new UniformTime(0, 0.1*s), {0,0,0});
        SM.SourceMode = new SourceCylinder(new Positronium(1.0, "gammas.txt"), new UniformTime(0, 0.1*s),  10.0, {0,0,-50.0}, {0,0,50.0}, "positions.txt");
        //SM.SourceMode = new SourceCylinder(new Positronium(1.0, true, "gammas.txt"), new UniformTime(0, 0.1*s),  10.0, {0,0,-50.0}, {0,0,50.0}, "positions.txt");

        //SM.SourceMode = new SourcePositronium(1.0, new UniformTime(0, 0.1*s), "GeneratedGammaData.txt"); // file name is optional: if not defined, data are not saved!
        //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -55.0*mm}, {0,0,1.0}, 1, new RoundProfile(20.0*mm));
        //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -55.0*mm}, {0,0,1.0}, 1, new CustomRaidalProfile("/home/andr/WORK/TPPT/FlashBeamProfile.txt", true) );
        //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -55.0*mm}, {0,0,1.0}, 1, new CustomRaidalProfile("/home/andr/WORK/TPPT/FlashBeamProfile.txt", false) );
        //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -55.0*mm}, {0,0,1.0}, 1, new CustomRaidalProfile("/home/andr/WORK/TPPT/Flat.txt", true) );
        //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -15.0*mm}, {0,0,1.0}, 1, new RoundProfile(20.0*mm));

        // Simulation mode
        SM.SimMode          = new ModeGui(); //SM.SimMode          = new ModeTracing();
        //SM.SimMode          = new ModeDummy(1e6);
        //SM.SimMode          = new ModeRadHard(1e6);
        //SM.SimMode          = new ModePesGenerator_Prob(1e5, {110.0, 110.0, 110.0}, {1, 1, 1}, {-55, -55, -55}, { {0.1*s, 1.0*s} });
        //SM.SimMode          = new ModeActivityGenerator(1e6, {1000.0, 1.0, 1.0}, {1, 101, 100}, {-500, -50.5, -50.0}, { {0.1*s, 155.0*s} }, "Activity.dat");
        //SM.SimMode          = new ModeDoseExtractor(1e5, {0.5, 0.5, 0.5}, {201, 201, 200}, {-50.25, -50.25, -50}, "Dose.txt");

        //SM.SimMode          = new ModeDoseExtractor(1e5, {0.5, 0.5, 0.5}, {200, 200, 200}, {-50, -50, -50}, "Dose.txt", false);
        //SM.SimMode          = new ModeDoseExtractor(1e5, {0.25, 0.25, 0.25}, {201, 201, 400}, {-25.125, -25.125, -50}, "Dose.txt", false);
        //SM.SimMode          = new ModeDoseExtractor(1e6, {0.5, 0.5, 0.5}, {101, 101, 200}, {-25.25, -25.25, -50}, "Dose.txt", false);
        //SM.SimMode          = new ModeDoseExtractor(1e5, {100.0, 100.0, 0.25}, {1, 1, 410}, {-50.0, -50.0, -50.0}, "Dose.txt", false);
        //SM.SimMode          = new ModeDoseExtractor(1e5, {100.0, 100.0, 0.325}, {1, 1, 350}, {-50.0, -50.0, -52.1}, "Dose.txt", false);
        //SM.SimMode          = new ModeDoseExtractor(1e5, {100.0, 100.0, 0.25}, {1, 1, 410}, {-50.0, -50.0, -10.0}, "Dose.txt", false);


    // --- END of user init ---
    }

    SM.startSession(); // warning: automatically saves config in working directory as SimConfig.json
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    SM.SimMode->run();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    out("Run time", (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count())*1e-6, "s");
    SM.endSession();
}

// Examples of phantoms
//    SM.Phantom = new PhantomNone;
//    SM.Phantom = new PhantomCylinder(200.0*mm, 200.0*mm, EMaterial::PMMA);
//    SM.Phantom = new PhantomBox(90.0, 300.0, 90.0, EMaterial::PE);
//    SM.Phantom = new PhantomDerenzo(200.0, 100.0, {1.8, 2.0, 2.2, 2.5, 3.0, 6.0}, 10.0, 5.0, 60.0);
//    SM.Phantom = new PhantomDICOM("/home/andr/WORK/TPPT/DicomPhant", "headCT_", 84, 252, 8, 155.0, {0,0,0}, {180.0*deg,0,0});
//    SM.Phantom = new PhantomDICOM("/home/andr/WORK/TPPT/DicomPhant", "headCT_", 150, 210, 2, 155.0, {0,0,0}, {180.0*deg,0,0});
// "Custom" phantom
//    SM.Phantom = new PhantomEnergyCalibration;
//    SM.Phantom = new PhantomEspana();
//    SM.Phantom = new PhantomBauerGel();
//    SM.Phantom = new PhantomParam;

// Examples of sources
//    SM.SourceMode       = new SourcePoint(new O15, new ConstantTime(0), {20.0, 20.0, -20.0});
//    SM.SourceMode       = new SourcePoint(new GammaPair, new ExponentialTime(0, 2.034*60*s), {1.2, 2.3, 2});
//    SM.SourceMode       = new SourceMultiBeam(new Proton(), "/home/andr/WORK/TPPT/MultiBeam/BeamletData.txt", 10); // NomEnergy[MeV] XIso[mm] ZIso[mm] Time0[ns] TimeSpan[ns] StatWeight
//    SM.SourceMode       = new SourceBeam(new Geantino(), new ConstantTime(0), {0*mm, 0*mm, 3500.0*mm}, {0,0,-1.0});
//    SM.SourceMode       = new SourcePointBlurred(new GammaPair, new ExponentialTime(0, 2.034*60*s), {0, 0, 0}, "/data/margarida/Data/AnnihilTest.txt", 21, 10.0);
//    SM.SourceMode       = new SourceNa22Point(0,1.0*s, {0, 0, 0});
//    SM.SourceMode       = new SourceLine(new GammaPair, new UniformTime(0, 1e10), {20.0, 20.0, -20.0}, {20.0, 20.0, 20.0});
//    SM.SourceMode       = new SourceBeam(new Proton(116.0*MeV), new UniformTime(0, 372*s), {0, -150.0, 0}, {0,1.0,0}, 1, new UniformProfile(70.0*mm, 70.0*mm));
//    SM.SourceMode       = new SourceMaterialLimited(new O15, new ConstantTime(0), {0, 0, 0}, {200.0,200.0,200.0}, "G4_WATER", "der.txt");
//    SM.SourceMode       = new SourceLysoNatural(0, 1e-5*s);
//    SM.SourceMode       = new SourceParticleListFile("/home/andr/WORK/TPPT/FirstStage.bin", true);
//    SM.SourceMode       = new SourceCylinder(new GammaPair, new UniformTime(0, 500.0*s), 0.5*330, {0,0,-0.5*105}, {0,0,0.5*105});//, "testPos.txt" );
//    SM.SourceMode       = new SourcePesHistogramFiles("/home/andr/WORK/tmp", 1000, true);
//    SM.SourceMode       = new SourceAnnihilHistFile("/home/andr/WORK/TPPT/PESGen/test100.txt", 1, true);
//    SM.SourceMode       = new SourceBeam(new Proton(160.0*MeV), new UniformTime(0, 238.0*s), {0, -150.0, 0}, {0,1.0,0}, 1);
//    SM.SourceMode       = new SourceBeam(new Proton(157.43*MeV), new UniformTime(0, 238.0*s), {0, -150.0, 0}, {0,1.0,0}, 100);   // PE-E3
//    SM.SourceMode       = new SourceBeam(new Proton(106.82*MeV), new UniformTime(0, 254.0*s), {0, -200.0, 0}, {0,1.0,0}, 100);   // PE-E1
//    SM.SourceMode       = new SourceBeam(new Proton(125.67*MeV), new UniformTime(0, 207.0*s), {0, -200.0, 0}, {0,1.0,0}, 100);   // Gel-E2
//    SM.SourceMode       = new SourceBeam(new Proton(176.75*MeV), new UniformTime(0, 203.0*s), {0, -200.0, 0}, {0,1.0,0}, 100);   // Gel-E4
//    SM.SourceMode       = new SourceBeam(new Proton(125.67*MeV), new UniformTime(0, 194.0*s), {0, -200.0, 0}, {0,1.0,0}, 100);   // PMMA-E2

// Examples of simulation modes
//    SM.SimMode          = new ModeTestScintPositions();
//    SM.SimMode          = new ModeParticleLogger(20, "TestParticleSaver.txt", false);
//    SM.SimMode          = new ModeTracing();
//    SM.SimMode          = new ModeDoseExtractor(1e5, {1,1,1}, {121,120,121}, {-60.5, -60, -60.5}, "DoseEspana.txt", false);
//    SM.SimMode          = new ModeDepositionScint(SM.SourceMode->CountEvents(), "FromProb-1000m.txt", false);
//    SM.SimMode          = new ModePesGenerator_MC(SM.SourceMode->CountEvents(), "Pes.dat", false);
//    SM.SimMode          = new ModePesGenerator_Prob(1e5, {1.0, 1.0, 1.0}, {201, 201, 201}, {-100.5, -100, -100.5}, { {0, 1e20*s} });
//    SM.SimMode          = new ModeAnnihilationLogger(SM.SourceMode->CountEvents(), {1.0, 1.0, 1.0}, {201, 201, 201}, {-100.5, -100, -100.5}, "test.txt");
//    SM.SimMode          = new ModeShowEvent(119);
//    SM.SimMode          = new ModeTestAcollinearity(10000, 2.0, 100, "AcolTest.txt");
//    SM.SimMode          = new ModeTestAnnihilations(SM.SourceMode->CountEvents(), 0, "Annihil.txt", false);
//    SM.SimMode          = new ModeDepositionScint(1000, "SimOutput.txt", false);
//    SM.SimMode          = new ModeDepositionScint(SM.SourceMode->CountEvents(), "SimOutput.txt", false); // if using FromFileSource to use all events in the file
//    SM.SimMode          = new ModeDepositionScint(SM.getNumberNatRadEvents(0, 1e-5*s), "SimOutput.bin", true);
//    SM.SimMode          = new ModeParticleLogger(1e3, "FirstStage.bin", true);
//    SM.SimMode          = new ModeTestLysoNatRad(1000000, 500, "natRadEnergyDistr.txt");
//    SM.SimMode          = new ModeTestDepositionStat(1e6, 0.01, {0.05, 0.1});
//    SM.SimMode          = new ModeActivityGenerator(SM.SourceMode->CountEvents(), {1.0, 1.0, 1.0}, {201, 201, 201}, {-100.5, -100, -100.5},  { {0, 1e10*s} }, "multiNew.dat");
//    SM.SimMode           = new ModeEnergyCalibration(50000, 1, "Convertion.txt");

// Obsolete detector components
//    SM.DetectorComposition.add(DetComp::GDML); SM.GdmlFileName = "detector.gdml";
//    SM.DetectorComposition.add(DetComp::Nozzle);

// Automatic name generation using seed
//    std::string filename = SM.generateName("base", "txt")   // e.g. SM.generateName("base", "txt") when seed is 100 --> "base_100.txt"
