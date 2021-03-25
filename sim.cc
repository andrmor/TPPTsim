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

    //SM.SimulationMode = new SimModeGui(SourceModeEnum::GammaPair);
    //SM.SimulationMode = new SimModeShowEvent(SourceModeEnum::GammaPair, 9338);
    //SM.SimulationMode = new SimModeScintPosTest(SourceModeEnum::GammaPair);
    SM.SimulationMode = new SimModeSingleEvents(SourceModeEnum::GammaPair);

    // --- END of user init ---

    SM.startSession(argc, argv); // has to be after setting up of the simulation mode

    SM.SimulationMode->run();

    SM.endSession();
}
