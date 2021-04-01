#include "SessionManager.hh"
#include "SimMode.hh"
#include "out.hh"

#include <sstream>

int main(int argc, char** argv)
{
    SessionManager & SM = SessionManager::getInstance();

    // --- START of user init ---

    SM.Seed = 0;
    SM.WorkingDirectory = "/home/andr/WORK/TPPT";

    //SM.SimMode = new SimModeMultipleEvents(SourceModeEnum::GammaPair, DetectorModeEnum::OnlyScint, PhantomModeEnum::PMMA);
    SM.SimMode = new SimModeGui(SourceModeEnum::GammaPair, DetectorModeEnum::OnlyScint, PhantomModeEnum::PMMA);
    //SM.SimMode = new SimModeGui(SourceModeEnum::C11, DetectorModeEnum::WithDetector, PhantomModeEnum::PMMA);
    //SM.SimMode = new SimModeGui(SourceModeEnum::GammaPair, DetectorModeEnum::WithDetector, PhantomModeEnum::PMMA);
    //SM.SimMode = new SimModeShowEvent(SourceModeEnum::GammaPair, DetectorModeEnum::OnlyScint, PhantomModeEnum::PMMA, 9643);
    //SM.SimMode = new SimModeScintPosTest(SourceModeEnum::GammaPair, DetectorModeEnum::OnlyScint, PhantomModeEnum::PMMA);
    //SM.SimMode = new SimModeSingleEvents(SourceModeEnum::GammaPair, DetectorModeEnum::OnlyScint, PhantomModeEnum::PMMA);

    // --- END of user init ---

    SM.startSession(argc, argv); // has to be after setting up of the simulation mode
    SM.SimMode->run();
    SM.endSession();
}
