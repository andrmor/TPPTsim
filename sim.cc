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
#include "PesGenerationMode.hh"
#include "ActivityProfiler.hh"

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
        // beware of possible overwright!

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
        SM.Seed               = 100;
        SM.SimAcollinearity   = true;  // only for the phantom region!
        SM.KillNeutrinos      = true;
        SM.UseStepLimiter     = true; SM.PhantomStepLimt = 1.0*mm;

        SM.CutPhantomGamma    = 10.0*mm;
        SM.CutPhantomElectron = 10.0*mm;
        SM.CutPhantomPositron = 0.1 *mm;
        SM.CutScintGamma      = 0.1 *mm;
        SM.CutScintElectron   = 0.1 *mm;
        SM.CutScintPositron   = 0.1 *mm;

        //double timeFrom = 0;
        //double timeTo   = 1e-5*s;  // currently implemented only for the natural rad from LYSO!

        SM.WorkingDirectory  = "/home/andr/WORK/TPPT";
        //SM.WorkingDirectory = "/data/margarida/Data";

        SM.Verbose          = false;
        SM.Debug            = false;
        SM.ShowEventNumber  = true;
        SM.EvNumberInterval = 10000;

        // Phantom
        //SM.PhantomMode      = new PhantomNone;
        SM.PhantomMode      = new PhantomPMMA;
        //SM.PhantomMode      = new PhantomDerenzo(200.0, 100.0, {1.8, 2.0, 2.2, 2.5, 3.0, 6.0}, 20.0, 10.0, 45.0);
        //SM.PhantomMode      = new PhantomParam();
        //SM.PhantomMode      = new PhantomDICOM("/home/andr/WORK/TPPT/DicomPhant", "headCT_", 84, 252, 8, 155.0, {0,0,0});
        //SM.PhantomMode      = new PhantomDICOM("/home/andr/WORK/TPPT/DicomPhant", "headCT_", 150, 210, 2, 155.0, {0,0,0});
        //SM.PhantomMode      = new PhantomCustomBox(150.0, 200.0, 150.0, PhantomCustomBox::HDPE);
        //SM.PhantomMode      = new PhantomEspana();
        //SM.PhantomMode      = new PhantomCustomBox(90.0, 300.0, 90.0, PhantomCustomBox::PE);
        //SM.PhantomMode      = new PhantomBauerGel();
        //SM.PhantomMode      = new PhantomCustomBox(90.0, 300.0, 90.0, PhantomCustomBox::PMMA);
        //SM.PhantomMode      = new PhantomCustomBox(90.0, 300.0, 90.0, PhantomCustomBox::Brain);
        //SM.PhantomMode      = new PhantomRT;

        // Enabled detector components - it is also possible to use .set( {comp1, comp2, ...} )
        SM.DetectorComposition.add(DetComp::Scintillators);
        SM.DetectorComposition.add(DetComp::Base);
        SM.DetectorComposition.add(DetComp::ClosedStructure);
        SM.DetectorComposition.add(DetComp::SIPM);
        SM.DetectorComposition.add(DetComp::PCB);
        SM.DetectorComposition.add(DetComp::CopperStructure);
        SM.DetectorComposition.add(DetComp::CoolingAssemblies);

            // Need special care using the following component - might be not cumulative
        //SM.DetectorComposition.add(DetComp::FirstStageMonitor);
        //SM.DetectorComposition.add(DetComp::GDML); SM.GdmlFileName = "/home/margarida/Downloads/GDML_version3/mother.gdml";

        // Source
        //SM.SourceMode       = new PointSource(new GammaPair, new ExponentialTime(0, 2.034*60*s), {1.2, 2.3, 2});
        //SM.SourceMode       = new BlurredPointSource(new GammaPair, new ExponentialTime(0, 2.034*60*s), {0, 0, 0}, "/data/margarida/Data/AnnihilTest.txt");
        //SM.SourceMode       = new PointSource(new O15, new ConstantTime(0), {0, 0, 0});
        //SM.SourceMode       = new LineSource(new O15, new ConstantTime(0), {20.0, 20.0, -20.0}, {20.0, 20.0, 20.0});
        //SM.SourceMode       = new PencilBeam(new Proton(116.0*MeV), new UniformTime(0, 372*s), {0, -150.0, 0}, {0,1.0,0}, 1, new UniformProfile(70.0*mm, 70.0*mm));
        //SM.SourceMode       = new PencilBeam(new Proton(160.0*MeV), new UniformTime(0, 238.0*s), {0, -150.0, 0}, {0,1.0,0}, 1);
        SM.SourceMode       = new PencilBeam(new Geantino, new UniformTime(0, 238.0*s), {0,5,5}, {1,0,0}, 1);
        //SM.SourceMode       = new PencilBeam(new Proton(55.0*MeV), new UniformTime(0, 238.0*s), {0, -150.0, 0}, {0,1.0,0}, 1);
        //SM.SourceMode       = new PencilBeam(new Proton(157.43*MeV), new UniformTime(0, 238.0*s), {0, -150.0, 0}, {0,1.0,0}, 100);   // PE-E3
        //SM.SourceMode       = new PencilBeam(new Proton(106.82*MeV), new UniformTime(0, 254.0*s), {0, -200.0, 0}, {0,1.0,0}, 100);   // PE-E1
        //SM.SourceMode       = new PencilBeam(new Proton(125.67*MeV), new UniformTime(0, 207.0*s), {0, -200.0, 0}, {0,1.0,0}, 100);   // Gel-E2
        //SM.SourceMode       = new PencilBeam(new Proton(176.75*MeV), new UniformTime(0, 203.0*s), {0, -200.0, 0}, {0,1.0,0}, 100);   // Gel-E4
        //SM.SourceMode       = new PencilBeam(new Proton(125.67*MeV), new UniformTime(0, 194.0*s), {0, -200.0, 0}, {0,1.0,0}, 100);   // PMMA-E2
        //SM.SourceMode       = new PencilBeam(new Geantino, new ConstantTime(0), {184.0*mm, 10.0*mm, -100.0*mm}, {0,0,1.0});
        //SM.SourceMode       = new MaterialLimitedSource(new O15, new ConstantTime(0), {0, 0, 0}, {200.0,200.0,200.0}, "G4_WATER", "/home/andr/WORK/TPPT/der.txt");
        //SM.SourceMode       = new NaturalLysoSource(timeFrom, timeTo);
        //SM.SourceMode       = new FromFileSource("/home/andr/WORK/TPPT/FirstStage.bin", true);
        //SM.SourceMode       = new MaterialLimitedSource(new O15, new UniformTime(0, 500.0*s), {0, 20.0, 0}, {40.0, 6.0, 66.0}, "G4_WATER", "/home/andr/WORK/TPPT/source.txt");

        // Simulation mode
        //SM.SimMode          = new SimModeGui();
        //SM.SimMode          = new SimModeShowEvent(119);
        //SM.SimMode          = new DoseExtractorMode(1e5, {1,1,1}, {121,120,121}, {-60.5, -60, -60.5}, "DoseEspana.txt");
        //SM.SimMode          = new SimModeScintPosTest();
        SM.SimMode          = new SimModeTracing();
        //SM.SimMode          = new SimModeAcollinTest(10000, 2.0, 100, "AcolTest.txt");
        //SM.SimMode          = new SimModeAnnihilTest(SM.SourceMode->CountEvents(), 0, "Annihil.txt", false);
        //SM.SimMode          = new SimModeNatRadTest(1000000, 500, "natRadEnergyDistr.txt");
        //SM.SimMode          = new SimModeSingleEvents(10000);
        //SM.SimMode          = new SimModeMultipleEvents(1000, "SimOutput.txt", false);
        //SM.SimMode          = new SimModeMultipleEvents(1e7, "SimOutput.bin", true);
        //SM.SimMode          = new SimModeMultipleEvents(SM.getNumberNatRadEvents(timeFrom, timeTo), "SimOutput.bin", true);
        //SM.SimMode          = new SimModeFirstStage(1e3, "FirstStage.bin", true);
        //SM.SimMode          = new SimModeMultipleEvents(SM.SourceMode->CountEvents(), "SimOutput.txt", false); // if using FromFileSource to use all events in the file
        //SM.SimMode          = new PesGenerationMode(1e6, "Pes.dat", false); // MC PES mode, number of protons = events times last argument in PencilBeam!
        //SM.SimMode          = new PesGenerationMode(1e3, {1.0, 1.0, 1.0}, {91, 200, 91}, {-45.5, -150, -45.5}); // Direct PES mode; number of protons = events times last argument in PencilBeam!
        //SM.SimMode          = new PesGenerationMode(1e3, {1.0, 1.0, 1.0}, {101, 250, 101}, {-50.5, -200, -50.5}); // Direct PES mode; number of protons = events times last argument in PencilBeam!
        //SM.SimMode          = new ActivityProfilerMode({{0,194,1}}, {{417,418}}, "/home/andr/WORK/TPPT", "bench");
        //SM.SimMode          = new PesAnalyzerMode(100000, "aaaaa.txt");

    // --- END of user init ---
    }

    SM.startSession();
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    SM.SimMode->run();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    out("Run time", (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count())*1e-6, "s");
    SM.endSession();
}
