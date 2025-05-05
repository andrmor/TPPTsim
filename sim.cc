#include "SessionManager.hh"
#include "SourceMode.hh"
#include "SourceParticleListFile.hh"
#include "TimeGenerator.hh"
#include "DefinedParticles.hh"
#include "SimMode.hh"
#include "BeamCollimatorMode.hh"
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
        SM.Seed = 1000;                // WARNING: the seed can be overriden with a command line argument, e.g. sim -s 123456
        SM.SimAcollinearity   = true;  // only for the phantom region!
        SM.KillNeutrinos      = true;
        //SM.UseStepLimiter     = true; SM.PhantomStepLimit = 0.25*mm;
        //SM.UseStepLimiter     = true; SM.PhantomStepLimit = 0.325*mm;

        SM.CutPhantomGamma    = 10.0*mm; SM.CutPhantomElectron = 10.0*mm; SM.CutPhantomPositron = 0.1*mm;
        SM.CutScintGamma      =  0.1*mm; SM.CutScintElectron   = 0.1*mm;  SM.CutScintPositron   = 0.1*mm;

        SM.Verbose          = false;
        SM.ShowEventNumber  = true; SM.EvNumberInterval = 10000;
        //SM.ShowEventNumber  = true; SM.EvNumberInterval = 10;

        //SM.WorkingDirectory  = "/home/andr/WORK/TPPT_summer2024/tmp2/Scanner4heads3mm";
        //SM.WorkingDirectory  = "/media/andr/HDD/work";
        //SM.WorkingDirectory  = "/home/andr/WORK/TPPT_summer2024/johnTest/tpptsim_tmp";
        //SM.WorkingDirectory  = "/home/andr/WORK/UsProposalSim";
        SM.WorkingDirectory  = "/home/andr/WORK/Na22sourcePosi";

        // Beam collimator (optional)
        //SM.BeamCollimator = new BeamCollimatorMarek(BeamCollimatorMarek::Holes369, {0,0,-75*mm}, 90*deg);

        // Phantom
        //SM.Phantom = new PhantomNone;
        SM.Phantom = new PhantomCylinder(25.4, 100.0, "G4_PLEXIGLASS");
        //SM.Phantom = new PhantomBox(2*1.5*25.4, 100.0, 2*1.5*25.4, "G4_PLEXIGLASS");
        //SM.Phantom = new PhantomCylinder(30.0, 50.0, "G4_PLEXIGLASS");
        //SM.Phantom = new PhantomCylinder(40.0, 50.0, "G4_PLEXIGLASS");
        //SM.Phantom = new PhantomCylinder(25.4, 100.0, "G4_POLYETHYLENE");
        //SM.Phantom = new PhantomMarekCompartments();
        //SM.Phantom = new PhantomCylinder(6.35, 20.0, "G4_Cu");
        //SM.Phantom = new PhantomCylinder(6.35, 20.0, EMaterial::Ni400);
        //SM.Phantom = new PhantomBeamDerenzoAndr3inv();
        //SM.Phantom = new PhantomSphere(5.0*mm, 4.0*mm, {0,0,0}, "G4_PLEXIGLASS");

        // Detector components
        //SM.DetectorComposition.add(DetComp::Scintillators);
        //SM.DetectorComposition.add({DetComp::Base, DetComp::ClosedStructure, DetComp::SIPM, DetComp::PCB, DetComp::CopperStructure, DetComp::CoolingAssemblies});
        //SM.DetectorComposition.add({DetComp::Base, DetComp::SIPM, DetComp::PCB, DetComp::CopperStructure, DetComp::CoolingAssemblies});
        //SM.DetectorComposition.add(DetComp::ParticleLogger); // required only for ModeParticleLogger
        // special scintillator arrangements
        //SM.DetectorComposition.add(DetComp::MiniPET);
        //SM.DetectorComposition.add(DetComp::MicroPET);
        //SM.DetectorComposition.add(DetComp::DoiPET);
        //SM.DetectorComposition.add(DetComp::FlatPanelPET);
    //    SM.DetectorComposition.add(DetComp::TungstenCubes2);
        //SM.DetectorComposition.add(DetComp::PLoggerMicroPET);

        // Source
        //SM.SourceMode = new SourceMultiBeam(new Proton(), {BeamRecord{72.5*MeV, 0,0, 0,0.1*s, 1.0}}, 1e6); // Energy XIsoCenter ZIsoCenter TimeStart TimeSpan StatWeight
        //SM.SourceMode = new SourceMultiBeam(new Proton(), {BeamRecord{89.6*MeV, 0,0, 0,0.1*s, 1.0}}, 1e6); // Energy XIsoCenter ZIsoCenter TimeStart TimeSpan StatWeight
        //SM.SourceMode = new SourceMultiBeam(new Proton(), {BeamRecord{103.8*MeV, 0,0, 0,0.1*s, 1.0}}, 1e6); // Energy XIsoCenter ZIsoCenter TimeStart TimeSpan StatWeight
        //SM.SourceMode = new SourceMultiBeam(new Proton(), {BeamRecord{103.8*MeV, 0,0, 0,30*s, 1.0}}, 1e6); // Energy XIsoCenter ZIsoCenter TimeStart TimeSpan StatWeight

        SourcePositronium * source = new SourcePositronium(0.5, new ConstantTime(0));
        source->setNa22Origin(); // adds deexcitation gamma of Na22
        source->setParaLifetime(0*ns);
        source->setOrtoLifetime(138.6*ns);
        SM.SourceMode = source;

        //SM.SourceMode = new SourcePoint(new Na22_decay(), new ConstantTime(0), {0*mm, 0*mm, 0*mm});
        //SM.SourceMode = new SourceBeam(new Geantino(), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -55.0*mm}, {0,0,1.0});
        //SM.SourceMode = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -55.0*mm}, {0,0,1.0});
        //SM.SourceMode = new SourceCylinder(new GammaPair(), new UniformTime(0, 10*s), 0.5*24*mm, {0,0,-25*mm}, {0,0,25*mm});
        //SM.SourceMode = new SourcePoint(new Positronium(1.0, "gammas.txt"), new UniformTime(0, 0.1*s), {0,0,0});
        //SM.SourceMode = new SourceCylinder(new Positronium(1.0, true, "gammas.txt"), new UniformTime(0, 0.1*s),  10.0, {0,0,-50.0}, {0,0,50.0}, "positions.txt");
        //SM.SourceMode = new SourceCylinder(new Positronium(1.0, true, "gammas.txt"), new UniformTime(0, 0.1*s),  10.0, {0,0,-50.0}, {0,0,50.0}, "positions.txt");
        //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -55.0*mm}, {0,0,1.0}, 1, new RoundProfile(20.0*mm));
        //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -15.0*mm}, {0,0,1.0}, 1, new RoundProfile(20.0*mm));
        //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -55.0*mm}, {0,0,1.0}, 1, new CustomRaidalProfile("/home/andr/WORK/TPPT/Flat.txt", true) );
        //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -55.0*mm}, {0,0,1.0}, 1, new CustomRaidalProfile("/home/andr/WORK/TPPT/FlashBeamProfile.txt", true) );
       //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -55.0*mm}, {0,0,1.0}, 1, new CustomRaidalProfile("/home/andr/WORK/TPPT/FlashBeamProfile.txt", false) );
       //SM.SourceMode       = new SourceBeam(new Proton(75.8*MeV), new UniformTime(0, 0.1*s), {0*mm, 0*mm, -110.0*mm}, {0,0,1.0}, 1, new CustomRaidalProfile("/home/andr/WORK/TPPT/FlashBeamProfile.txt", false) );
        //SM.SourceMode = new SourceAnnihilHistFile("/home/andr/WORK/TPPT_summer2024/tmp/Activity1e6.dat", 3.5e4, true);
        //SM.SourceMode = new SourceAnnihilHistFile("/home/andr/WORK/TPPT_summer2024/tmp/Activity1e6_ph2.dat", 3.5e4, true);
        //SM.SourceMode = new SourceAnnihilHistFile("/home/andr/WORK/TPPT_summer2024/tmp/Activity1e6_ph2.dat", 3.5e4, true, {50.0*mm, 50.0*mm, 25.0*mm});
        //SM.SourceMode = new SourceAnnihilHistFile("/media/andr/HDD/work/Activity123vert_1e7.dat", 3.5e3, true, {50.0*mm, 50.0*mm, 0.0*mm});
        //SM.SourceMode = new SourceAnnihilHistFile("/media/andr/HDD/work/Activity123hor_1e7.dat", 3.5e3, true, {50.0*mm, 50.0*mm, 0.0*mm});
        //SM.SourceMode = new SourceAnnihilHistFile("/media/andr/HDD/work/Activity1e6_ph3inv.dat", 3.5e4, true, {50.0*mm, 50.0*mm, 0.0*mm});
//        SM.SourceMode = new SourceCylinder(new GammaPair(), new UniformTime(0, 1000*s), 0.5*mm, {50.0,50.0,-50.0}, {50.0,50.0,50.0});
//        SM.SourceMode = new SourceAnnihilHistFile("/home/andr/WORK/TPPT_summer2024/tmp/Activity1e6_ph2inv_a.dat", 3.5e4, true);
        //SM.SourceMode = new SourceAnnihilHistFile("/home/andr/WORK/TPPT_summer2024/Sphere/ActivitySphereRad5i4_1e7.dat", 3.5e3, true, {50.0*mm, 50.0*mm, 0.0*mm});
//       SM.SourceMode = new SourceAnnihilHistFile("/media/andr/HDD/work/Activity1e6_MarekCross.dat", 3.5e4, true);
    //    SM.SourceMode = new SourceAnnihilHistFile("/media/andr/HDD/work/Activity1e6_Marek3holes.dat", 3.5e4, true);
    //     SM.SourceMode = new SourceAnnihilHistFile("/media/andr/HDD/work/Activity1e6_ph3inv.dat", 3.5e4, true);
    //    SM.SourceMode = new SourceAnnihilHistFile("/media/andr/HDD/work/Activity1e6_Micro3holes.dat", 3.5e4, true);
   //     SM.SourceMode = new SourceAnnihilHistFile("/media/andr/HDD/work/Activity1e6_Ring.dat", 3.5e4, true);
   //     SM.SourceMode = new SourceParticleListFileMicroPET("/media/andr/HDD/work/PartLog.txt", false);
        //SM.SourceMode = new SourceParticleListFile("/media/andr/HDD/work/PartLog.txt", false);
        //SM.SourceMode = new SourceParticleListFileMicroPET("/media/andr/HDD/work/PartLog_400.txt", false, 30*mm);


        // Simulation mode
       //SM.SimMode          = new ModeGui();
        SM.SimMode          = new SourceTester(100000, 100,0,1000*ns, "time.dat");
        //SM.SimMode = new ModeScintDepoLogger(1e8, 1e10*ns, "DepoTester_100keV_2cubes.txt");
        //SM.SimMode = new ModeScintDepoLogger(1e7, 1e10*ns, "DepoTester_Rotated.txt");
        //SM.SimMode = new ModeScintDepoLogger(1e8, 0.1*s, "Depo_PartLog_400_2Cubes.txt");
        //SM.SimMode          = new ModeTracing();
        //SM.SimMode          = new ModeDummy(1e5);
        //SM.SimMode          = new ModeRadHard(1e6);
        //SM.SimMode          = new ModeActivityGenerator(1e6, {1000.0, 1.0, 1.0}, {1, 101, 100}, {-500, -50.5, -50.0}, { {0.1*s, 155.0*s} }, "Activity.dat");
        //SM.SimMode          = new ModeDoseExtractor(1e5, {0.5, 0.5, 0.5}, {201, 201, 200}, {-50.25, -50.25, -50}, "Dose.txt");
        //SM.SimMode          = new ModeDoseExtractor(1e5, {0.5, 0.5, 0.5}, {201, 201, 200}, {-50.25, -50.25, -35}, "Dose.txt");
        //SM.SimMode          = new ModeDoseExtractor(1e5, {0.5, 0.5, 0.5}, {201, 201, 200}, {-50.25, -50.25, -35}, "DepoE.txt", true);
        //SM.SimMode          = new ModePesGenerator_Prob(1e5, {0.5, 0.5, 0.5}, {201, 201, 200}, {-50.25, -50.25, -35}, { {0.1*s, 1000*s} });

        //SM.SimMode = new ModeActivityGenerator(1e6, {1.0, 1.0, 1.0}, {70, 100, 70}, {-35.0, -50.0, -35.0}, { {0.1*s, 1000.0*s} }, "Activity1e6_72i5.dat");
        //SM.SimMode = new ModeActivityGenerator(1e6, {1.0, 1.0, 1.0}, {70, 100, 70}, {-35.0, -50.0, -35.0}, { {0.1*s, 1000.0*s} }, "Activity1e6_89i6.dat");
        //SM.SimMode = new ModeActivityGenerator(1e6, {1.0, 1.0, 1.0}, {70, 100, 70}, {-35.0, -50.0, -35.0}, { {0.1*s, 1000.0*s} }, "Activity1e6_103i8.dat");

        //SM.SimMode = new ModeActivityGenerator(1e6, {1.0, 1.0, 1.0}, {70, 100, 70}, {-35.0, -50.0, -35.0}, { {30.1*s, 15*60*s} }, "Activity1e6_103i8_30i15.dat");

        //SM.SimMode = new ModeActivityGenerator(1e6, {0.2, 0.2, 1.0}, {200, 200, 40}, {-20.0, -20.0, -20.0}, { {0.1*s, 1000.0*s} }, "Activity1e6_ph2_a.dat");
        //SM.SimMode = new ModeDepositionScint(200, "DepoFromActivityFlash_v2_doi_50_50_25.dat", true);
        //SM.SimMode = new ModeDepositionScint(5e7, "DepoLine_5e7_3mmScanner.dat", true);
//        SM.SimMode = new ModeActivityGenerator(1e6, {0.2, 0.2, 1.0}, {200, 200, 50}, {-20.0, -20.0, -25.0}, { {0.1*s, 1000.0*s} }, "Activity1e6_ph2inv_a.dat");
//        SM.SimMode = new ModeActivityGenerator(1e6, {0.2, 0.2, 1.0}, {200, 200, 50}, {-20.0, -20.0, -25.0}, { {0.1*s, 1000.0*s} }, "Activity1e6_Marek3holes.dat");
   //    SM.SimMode = new ModeActivityGenerator(1e6, {0.2, 0.2, 1.0}, {200, 200, 50}, {-20.0, -20.0, -25.0}, { {0.1*s, 1000.0*s} }, "Activity1e6_ph3inv.dat");
   //    SM.SimMode = new ModeActivityGenerator(1e6, {0.2, 0.2, 1.0}, {200, 200, 50}, {-20.0, -20.0, -25.0}, { {0.1*s, 1000.0*s} }, "Activity1e6_Micro3holes.dat");
   //    SM.SimMode = new ModeActivityGenerator(1e7, {0.2, 0.2, 1.0}, {200, 200, 50}, {-20.0, -20.0, -25.0}, { {0.1*s, 1000.0*s} }, "Activity123hor_1e7.dat");
//     SM.SimMode = new ModeActivityGenerator(1e7, {0.2, 0.2, 0.2}, {50, 50, 50}, {-5.0, -5.0, -5.0}, { {0.1*s, 1000.0*s} }, "ActivitySphereRad5i4_1e7.dat");
//     SM.SimMode = new ModeActivityGenerator(1, {0.2, 0.2, 0.2}, {50, 50, 50}, {-5.0, -5.0, -5.0}, { {0.1*s, 1000.0*s} }, "Test.dat");
//        SM.SimMode = new ModeDepositionScint(5e7, "Depo_Activity_ph2inv_a.dat", true);
//       SM.SimMode = new ModeDepositionScint(200, "Depo_Activity1e6_MarekCross.dat", true);
 //       SM.SimMode = new ModeDepositionScint(200, "Depo_Activity1e6_Marek3holes.dat", true);
    //    SM.SimMode = new ModeDepositionScint(200, "Depo_Activity1e6_ph3inv.dat", true);
    //     SM.SimMode = new ModeDepositionScint(200, "Depo_Activity1e6_Micro3holes.dat", true);
       //    SM.SimMode = new ModeParticleLogger(1e8, "1PartLog_400.txt", false);
      //SM.SimMode = new ModeDepositionScint(200, "Depo_Activity123vert_1e7_to3i5e10_xyz50_50_0.dat", true);
      //SM.SimMode = new ModeDepositionScint(200, "Depo_Activity123hor_1e7_to3i5e10_xyz50_50_0.dat", true);
      //SM.SimMode = new ModeDepositionScint(200, "Depo_Activity1e6_ph3inv_to3i5e10_xyz50_50_0.dat", true);
//      SM.SimMode = new ModeDepositionScint(50, "Depo_ActivitySphereRad5i4_1e7_to3i5e10_xyz50_50_0.dat", true);
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

/*
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
*/
