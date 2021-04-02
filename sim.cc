#include "SessionManager.hh"
#include "SourceMode.hh"
#include "SimMode.hh"
#include "out.hh"

#include <sstream>

int main(int argc, char** argv)
{
    SessionManager & SM = SessionManager::getInstance();

// --- START of user init ---

  // General settings
    SM.Seed             = 0;
    SM.WorkingDirectory = "/home/andr/WORK/TPPT";
    SM.bG4Verbose       = false;
    SM.bDebug           = false;

  // Phantom
    SM.PhantomMode      = PhantomModeEnum::PMMA;

  // Detector
    SM.DetetctorMode    = DetectorModeEnum::OnlyScints;
    //SM.DetetctorMode    = DetectorModeEnum::ScintsAndGDML;

  // Source
    //SM.SourceMode       = new PointSource(new ParticleC11, {0, 0, SM.GlobalZ0}, 1);
    //SM.SourceMode       = new PointSource(new ParticleGammaPair, {0, 0, SM.GlobalZ0}, 1);
    //SM.SourceMode       = new PencilBeam(new ParticleGamma(511.0*keV), {0, 0, SM.GlobalZ0}, {1.0,0,0}, 1);
    SM.SourceMode       = new PencilBeam(new ParticleGeantino, {150.0, 150.0, SM.GlobalZ0-100.0}, {0,0,1.0}, 1);

  // Operation mode
    //SM.SimMode          = new SimModeGui();
    //SM.SimMode          = new SimModeShowEvent(9643);
    //SM.SimMode          = new SimModeScintPosTest();
    //SM.SimMode          = new SimModeSingleEvents();
    //SM.SimMode          = new SimModeMultipleEvents();
    SM.SimMode          = new SimModeTracing();

// --- END of user init ---

    SM.startSession(argc, argv); // has to be after setting up of the simulation mode
    SM.SimMode->run();
    SM.endSession();
}
