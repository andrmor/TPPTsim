#include "TrackingAction.hh"
#include "SessionManager.hh"
#include "PesGenerationMode.hh"

void PesGeneratorTrackingAction::PreUserTrackingAction(const G4Track * /*track*/)
{
    SessionManager & SM = SessionManager::getInstance();
    PesGenerationMode * Mode = static_cast<PesGenerationMode*>(SM.SimMode);
    Mode->bNewTrackStarted = true;
}
