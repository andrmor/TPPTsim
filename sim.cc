#include "SessionManager.hh"
#include "SimMode.hh"
#include "out.hh"

#include <sstream>

int main(int argc, char** argv)
{
    SessionManager & SM = SessionManager::getInstance();

    // --- START of user init ---

    SM.Seed             = 0;
    SM.WorkingDirectory = "/home/andr/WORK/TPPT";
    SM.bG4Verbose       = false;
    SM.bDebug           = true;

    // Phantom
    SM.PhantomMode      = PhantomModeEnum::PMMA;

    // Detector
    SM.DetetctorMode    = DetectorModeEnum::OnlyScint;
    //SM.DetetctorMode    = DetectorModeEnum::ScintsAndGDML;

    // Source
    SM.SourceMode       = SourceModeEnum::GammaPair;
    //SM.SourceMode       = SourceModeEnum::C11;

    // Operation mode
    //SM.SimMode          = new SimModeGui();
    //SM.SimMode          = new SimModeShowEvent(9643);
    //SM.SimMode          = new SimModeScintPosTest();
    //SM.SimMode          = new SimModeSingleEvents();
    SM.SimMode          = new SimModeMultipleEvents();

    // --- END of user init ---

    SM.startSession(argc, argv); // has to be after setting up of the simulation mode
    SM.SimMode->run();
    SM.endSession();
}
