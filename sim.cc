#include "SessionManager.hh"
#include "SourceMode.hh"
#include "SimMode.hh"
#include "PhantomMode.hh"
#include "out.hh"

#include <sstream>

int main(int argc, char** argv)
{
    SessionManager & SM = SessionManager::getInstance();

// --- START of user init ---

  // General settings
    SM.Seed             = 0;
    SM.WorkingDirectory = "/home/andr/WORK/TPPT";
    //SM.WorkingDirectory = "/data/margarida/Data";
    SM.bG4Verbose       = false;
    SM.bDebug           = false;

  // Phantom
    SM.PhantomMode      = new PhantomModePMMA;
    //SM.PhantomMode      = new PhantomModeNone;
    //SM.PhantomMode      = new PhantomModeDerenzo(200.0, 200.0, {1.8, 2.0, 2.2, 2.5, 3.0, 6.0}, 20.0, 10.0, 45.0);
    //SM.PhantomMode      = new PhantomModeDerenzo(200.0, 200.0, {4.0}, 90.0, 5.0, 45.0);

  // Detector
    SM.DetectorComposition = {DetComp::Scint};
    //SM.DetectorComposition = {DetComp::Scint, DetComp::GDML};
    //SM.DetectorComposition = {};

  // Source
    SM.SourceMode       = new PointSource(new ParticleGammaPair, {0, 0, SM.GlobalZ0}, 1);
    //SM.SourceMode       = new PointSource(new ParticleC11, {0, 0, SM.GlobalZ0}, 100);
    //SM.SourceMode       = new PointSource(new ParticleN12, {0, 0, SM.GlobalZ0}, 100); // 11ms
    //SM.SourceMode       = new PointSourceUniformTime(new ParticleGammaPair, {0, 0, SM.GlobalZ0}, 1, 1000.0);
    //SM.SourceMode       = new PencilBeam(new ParticleGamma(511.0*keV), {0, 0, SM.GlobalZ0}, {1.0,0,0}, 1);
    //SM.SourceMode       = new PencilBeam(new ParticleGeantino, {150.0, 150.0, SM.GlobalZ0-100.0}, {0,0,1.0}, 1);
    //SM.SourceMode       = new MaterialLimitedSource(new ParticleGammaPair, {0, 0, SM.GlobalZ0}, {200.0,200.0,200.0}, "G4_WATER", "/home/andr/WORK/TPPT/der.txt");
    //SM.SourceMode       = new MaterialLimitedSource(new ParticleGamma, {0, 0, SM.GlobalZ0}, {200.0,200.0,200.0}, "G4_WATER");

  // Operation mode
    //SM.SimMode          = new SimModeGui();
    //SM.SimMode          = new SimModeShowEvent(100000);
    SM.SimMode          = new SimModeScintPosTest();
    //SM.SimMode          = new SimModeSingleEvents();
    //SM.SimMode          = new SimModeMultipleEvents(10000, "SimOutput.bin", true);
    //SM.SimMode          = new SimModeTracing();

// --- END of user init ---

    SM.startSession(argc, argv); // has to be after setting up of the simulation mode
    SM.SimMode->run();
    SM.endSession();
}
