#include "TrackingAction.hh"
#include "SessionManager.hh"
#include "PesGenerationMode.hh"
#include "AnnihilationLoggerMode.hh"
#include "out.hh"

#include "G4Track.hh"
#include "G4Positron.hh"

void PesGeneratorTrackingAction::PreUserTrackingAction(const G4Track * /*track*/)
{
    SessionManager & SM = SessionManager::getInstance();
    PesGenerationMode * Mode = static_cast<PesGenerationMode*>(SM.SimMode);
    Mode->bNewTrackStarted = true;
}

void AnnihilationLoggerTrackingAction::PostUserTrackingAction(const G4Track * track)
{
    if (track->GetParticleDefinition() != G4Positron::Definition()) return;

    //out(track->GetPosition()[0], track->GetPosition()[1], track->GetPosition()[2]);

    SessionManager & SM = SessionManager::getInstance();
    AnnihilationLoggerMode * mode = static_cast<AnnihilationLoggerMode*>(SM.SimMode);
    mode->fillPosition(track->GetPosition());
}
