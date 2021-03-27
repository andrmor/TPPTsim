#include "SessionManager.hh"
#include "SimulationMode.hh"
#include "out.hh"

#include <sstream>

int main(int argc, char** argv)
{
    SessionManager& SM = SessionManager::getInstance();

    // --- START of user init ---

    SM.Seed = 0;
    SM.FileName = "/data/margarida/Data/Test.txt";

    SM.SimulationMode = new SimModeGui(SourceModeEnum::GammaPair, DetectorModeEnum::WithDetector, PhantomModeEnum::PMMA);
    //SM.SimulationMode = new SimModeShowEvent(SourceModeEnum::GammaPair, DetectorModeEnum::OnlyScint, PhantomModeEnum::PMMA, 9338);
    //SM.SimulationMode = new SimModeScintPosTest(SourceModeEnum::GammaPair, DetectorModeEnum::OnlyScint, PhantomModeEnum::PMMA);
    //SM.SimulationMode = new SimModeSingleEvents(SourceModeEnum::GammaPair, DetectorModeEnum::OnlyScint, PhantomModeEnum::PMMA);

    // --- END of user init ---

    SM.startSession(argc, argv); // has to be after setting up of the simulation mode

    SM.SimulationMode->run();

    SM.endSession();
}
