#include "SessionManager.hh"
#include "SourceMode.hh"
#include "TimeGenerator.hh"
#include "DefinedParticles.hh"
#include "SimMode.hh"
#include "PhantomMode.hh"
#include "out.hh"

#include <sstream>

int main(int argc, char** argv)
{
    SessionManager & SM = SessionManager::getInstance();

// --- START of user init ---

  // General settings
    SM.Seed              = 0;
    SM.bSimAcollinearity = true;
    SM.WorkingDirectory  = "/home/andr/WORK/TPPT";
    //SM.WorkingDirectory = "/data/margarida/Data";
    SM.bG4Verbose        = false;
    SM.bDebug            = false;

  // Phantom
    //SM.PhantomMode      = new PhantomNone;
    //SM.PhantomMode      = new PhantomPMMA;
    SM.PhantomMode      = new PhantomTinyCube;
    //SM.PhantomMode      = new PhantomDerenzo(200.0, 100.0, {1.8, 2.0, 2.2, 2.5, 3.0, 6.0}, 20.0, 10.0, 45.0);

  // Detector
    SM.DetectorComposition = {};
    //SM.DetectorComposition = {DetComp::Scintillators};
    //SM.DetectorComposition = {DetComp::Scintillators, DetComp::GDML};

  // Source
    //SM.SourceMode       = new PointSource(new GammaPair, new ConstantTime(0), {0, 0, SM.GlobalZ0});
    //SM.SourceMode       = new PointSource(new GammaPair, new ExponentialTime(0, 2.034*60.0*s), {0, 0, SM.GlobalZ0});
    SM.SourceMode       = new PointSource(new O15, new ConstantTime(0), {0, 0, SM.GlobalZ0});
    //SM.SourceMode       = new PencilBeam(new GammaPair(511.0*keV, true), new ConstantTime(0), {0, 0, SM.GlobalZ0}, {1.0,0,0});
    //SM.SourceMode       = new PencilBeam(new Geantino, new ConstantTime(0), {150.0, 150.0, SM.GlobalZ0-100.0}, {0,0,1.0});
    //SM.SourceMode       = new MaterialLimitedSource(new GammaPair, new ConstantTime(0), {0, 0, SM.GlobalZ0}, {200.0,200.0,200.0}, "G4_WATER", "/home/andr/WORK/TPPT/der.txt");

  // Operation mode
    //SM.SimMode          = new SimModeGui();
    //SM.SimMode          = new SimModeShowEvent(100000);
    //SM.SimMode          = new SimModeScintPosTest();
    SM.SimMode          = new SimModeAcollinTest(1, "AcolTest.txt");
    //SM.SimMode          = new SimModeSingleEvents();
    //SM.SimMode          = new SimModeMultipleEvents(100, "SimOutput.txt", false);
    //SM.SimMode          = new SimModeMultipleEvents(1e6, "SimOutput.bin", true);
    //SM.SimMode          = new SimModeTracing();

// --- END of user init ---

    SM.startSession(argc, argv);
    SM.SimMode->run();
    SM.endSession();
}
