#include "SessionManager.hh"
#include "SourceMode.hh"
#include "FromFileSource.hh"
#include "TimeGenerator.hh"
#include "DefinedParticles.hh"
#include "SimMode.hh"
#include "PhantomMode.hh"
#include "DicomPhantom.hh"
#include "DetComp.hh"
#include "out.hh"

#include <chrono>

int main(int argc, char** argv)
{
    std::string filename;
    if (argc > 1)
    {
        filename = std::string(argv[1]);
        out("\nLoading config from file:", filename);
    }
    else out("\nNo config file provided as argument, using configuration defined in the main of sim");

    SessionManager & SM = SessionManager::getInstance(); // side effect: outputs Geant4 version

        // warning: automatically saves config in working directory as SimConfig.json
        // beware of possible overright!

        //here you can directly provide the config file name
    //filename = "/home/andr/WORK/TPPT/SimConfig1.json";

    if (!filename.empty())
    {
        SM.loadConfig(filename);
    }
    else
    {
     // --- START of user init ---

        // General settings
        SM.Seed               = 1000;
        SM.SimAcollinearity   = true;  // only for the phantom region!
        SM.KillNeutrinos      = true;

        SM.CutPhantomGamma    = 10.0*mm;
        SM.CutPhantomElectron = 10.0*mm;
        SM.CutPhantomPositron = 0.1 *mm;
        SM.CutScintGamma      = 0.1 *mm;
        SM.CutScintElectron   = 0.1 *mm;
        SM.CutScintPositron   = 0.1 *mm;

        //double timeFrom = 0;
        //double timeTo   = 1e-5*s;  // currently implemented only for the natural rad from LYSO!

        //SM.WorkingDirectory  = "/home/andr/WORK/TPPT";
        SM.WorkingDirectory = "/data/margarida/Data";

        SM.Verbose          = false;
        SM.Debug            = false;
        SM.ShowEventNumber  = true;
        SM.EvNumberInterval = 10000;

        // Phantom
        //SM.PhantomMode      = new PhantomNone;
        SM.PhantomMode      = new PhantomPMMA;
        //SM.PhantomMode      = new PhantomTinyCube;
        //SM.PhantomMode      = new PhantomDerenzo(200.0, 100.0, {1.8, 2.0, 2.2, 2.5, 3.0, 6.0}, 20.0, 10.0, 45.0);
        //SM.PhantomMode      = new PhantomParam;
        //SM.PhantomMode      = new PhantomModeDICOM(155.0, {0,0,50}, "Data.dat", true);

        // Enabled detector components - it is also possible to use .set( {comp1, comp2, ...} )
        SM.DetectorComposition.add(DetComp::Scintillators);
        //SM.DetectorComposition.add(DetComp::Base);
        //SM.DetectorComposition.add(DetComp::ClosedStructure);
        //SM.DetectorComposition.add(DetComp::SIPM);
        //SM.DetectorComposition.add(DetComp::PCB);
        SM.DetectorComposition.add(DetComp::CopperStructure);
        //SM.DetectorComposition.add(DetComp::CoolingAssemblies);
            // Need special care using the following component - might be not cumulative
        //SM.DetectorComposition.add(DetComp::FirstStageMonitor);
        //SM.DetectorComposition.add(DetComp::GDML); SM.GdmlFileName = "detector.gdml";

        // Source
        SM.SourceMode       = new PointSource(new GammaPair, new ExponentialTime(0, 2.034*60*s), {0, 0, SM.GlobalZ0});
        //SM.SourceMode       = new BlurredPointSource(new GammaPair, new ExponentialTime(0, 2.034*60*s), {0, 0, SM.GlobalZ0}, "/data/margarida/Data/AnnihilTest.txt");
        //SM.SourceMode       = new PointSource(new O15, new ConstantTime(0), {0, 0, SM.GlobalZ0});
        //SM.SourceMode       = new LineSource(new O15, new ConstantTime(0), {20.0, 20.0, SM.GlobalZ0-20.0}, {20.0, 20.0, SM.GlobalZ0+20.0});
        //SM.SourceMode       = new PencilBeam(new Gamma(511.0*keV), new ConstantTime(0), {100.0, 160.0, 50.0+SM.GlobalZ0}, {0,0,-1.0});
        //SM.SourceMode       = new PencilBeam(new Geantino, new ConstantTime(0), {0, 0, SM.GlobalZ0+12.3*mm}, {1.0,0,0});
        //SM.SourceMode       = new MaterialLimitedSource(new O15, new ConstantTime(0), {0, 0, SM.GlobalZ0}, {200.0,200.0,200.0}, "G4_WATER", "/home/andr/WORK/TPPT/der.txt");
        //SM.SourceMode       = new MaterialLimitedSource(new GammaPair, new ExponentialTime(0, 2.034*60.0*s), {0, 0, SM.GlobalZ0}, {200.0,200.0,200.0}, "G4_WATER", "/home/andr/WORK/TPPT/der.txt");
        //SM.SourceMode       = new NaturalLysoSource(timeFrom, timeTo);
        //SM.SourceMode       = new FromFileSource("/home/andr/WORK/TPPT/FirstStage.bin", true);

        // Operation mode
        SM.SimMode          = new SimModeGui();
        //SM.SimMode          = new SimModeShowEvent(119);
        //SM.SimMode          = new SimModeScintPosTest();
        //SM.SimMode          = new SimModeTracing();
        //SM.SimMode          = new SimModeAcollinTest(10000, 2.0, 100, "AcolTest.txt");
        //SM.SimMode          = new SimModeAnnihilTest(1e6, 10, 1000, "AnnihilTest.txt");
        //SM.SimMode          = new SimModeNatRadTest(1000000, 500, "natRadEnergyDistr.txt");
        //SM.SimMode          = new SimModeSingleEvents(10000);
        //SM.SimMode          = new SimModeMultipleEvents(1000, "1.txt", false);
        //SM.SimMode          = new SimModeMultipleEvents(10000000, "All10MSim.bin", true);
        //SM.SimMode          = new SimModeMultipleEvents(SM.getNumberNatRadEvents(timeFrom, timeTo), "SimOutput.bin", true);
        //SM.SimMode          = new SimModeFirstStage(1e3, "FirstStage.bin", true);
        //SM.SimMode          = new SimModeMultipleEvents(SM.SourceMode->CountEvents(), "SimOutput.txt", false); // if using FromFileSource to use all events in the file
        //SM.SimMode          = new SimModeDoseExtractor(1e6, 20.0, 100, "Dose.txt");

    // --- END of user init ---
    }

    SM.startSession();
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    SM.SimMode->run();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    out("Run time", (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count())*1e-6, "s");
    SM.endSession();
}
