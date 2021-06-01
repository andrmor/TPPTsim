#include "SessionManager.hh"
#include "SourceMode.hh"
#include "TimeGenerator.hh"
#include "DefinedParticles.hh"
#include "SimMode.hh"
#include "PhantomMode.hh"
#include "out.hh"
#include "Hist1D.hh"

#include <fstream>
#include <sstream>

#include <chrono>

int main(int argc, char** argv)
{
    SessionManager & SM = SessionManager::getInstance();

// --- START of user init ---

  // General settings
    SM.Seed              = 0;
    SM.bSimAcollinearity = true;  // only for the phantom region!

    double timeFrom = 0;
    double timeTo   = 1e-5*s;  // currently implemented only for the natural rad from LYSO!

    //SM.WorkingDirectory  = "/home/andr/WORK/TPPT";
    SM.WorkingDirectory = "/data/margarida/Data";

    SM.bG4Verbose        = false;
    SM.bDebug            = false;
    SM.bShowEventNumber  = false;

  // Phantom
    //SM.PhantomMode      = new PhantomNone;
    SM.PhantomMode      = new PhantomPMMA;
    //SM.PhantomMode      = new PhantomTinyCube;
    //SM.PhantomMode      = new PhantomDerenzo(200.0, 100.0, {1.8, 2.0, 2.2, 2.5, 3.0, 6.0}, 20.0, 10.0, 45.0);
    //SM.PhantomMode      = new PhantomParam;

  // Detector
    //SM.DetectorComposition = {};
    //SM.DetectorComposition = {DetComp::Scintillators};
    SM.DetectorComposition = {DetComp::Scintillators, DetComp::GDML};

  // Source
    //SM.SourceMode       = new PointSource(new GammaPair, new ExponentialTime(0, 2.034*60*s), {0, 0, SM.GlobalZ0});
    Hist1D dist(21, -10, 10);
    std::string inputFileName  = "/data/margarida/Data/AnnihilTest.txt";
    std::ifstream * inStream = new std::ifstream(inputFileName);
    std::vector<std::pair<double,double>> Input;
    std::string line;

    while (!inStream->eof())
    {
        getline(*inStream, line);
        //out(line);
        std::stringstream ss(line);
        double position, probability;
        ss >> position >> probability;
        Input.push_back({position,probability});
    }
    for (const auto & pair : Input)
        dist.fill(pair.first+0.001, pair.second);
    SM.SourceMode       = new BlurredPointSource(new GammaPair, new ExponentialTime(0, 2.034*60*s), {0, 0, SM.GlobalZ0}, dist, 12345);
    //SM.SourceMode       = new PointSource(new O15, new ConstantTime(0), {0, 0, SM.GlobalZ0});
    //SM.SourceMode       = new PencilBeam(new GammaPair(511.0*keV, true), new ConstantTime(0), {0, 0, SM.GlobalZ0}, {1.0,0,0});
    //SM.SourceMode       = new PencilBeam(new Geantino, new ConstantTime(0), {70.0, 190.0, SM.GlobalZ0+52.5}, {1.0,-1.1,0});
    //SM.SourceMode       = new MaterialLimitedSource(new O15, new ConstantTime(0), {0, 0, SM.GlobalZ0}, {200.0,200.0,200.0}, "G4_WATER", "/home/andr/WORK/TPPT/der.txt");
    //SM.SourceMode       = new MaterialLimitedSource(new GammaPair, new ExponentialTime(0, 2.034*60.0*s), {0, 0, SM.GlobalZ0}, {200.0,200.0,200.0}, "G4_WATER", "/home/andr/WORK/TPPT/der.txt");
    //SM.SourceMode       = new NaturalLysoSource(timeFrom, timeTo);

  // Operation mode
    SM.SimMode          = new SimModeGui();
    //SM.SimMode          = new SimModeShowEvent(119);
    //SM.SimMode          = new SimModeScintPosTest();
    //SM.SimMode          = new SimModeTracing();
    //SM.SimMode          = new SimModeAcollinTest(10000, 2.0, 100, "AcolTest.txt");
    //SM.SimMode          = new SimModeNatRadTest(1000000, 500, "natRadEnergyDistr.txt");
    //SM.SimMode          = new SimModeSingleEvents();
    //SM.SimMode          = new SimModeMultipleEvents(100, "SimOutput.txt", false);
    //SM.SimMode          = new SimModeMultipleEvents(10e6, "SimOutput.bin", true);
    //SM.SimMode          = new SimModeMultipleEvents(SM.getNumberNatRadEvents(timeFrom, timeTo), "SimOutput.bin", true);

// --- END of user init ---
    SM.startSession(argc, argv);
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    SM.SimMode->run();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    out("Run time", (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count())*1e-6, "s");
    SM.endSession();
}
