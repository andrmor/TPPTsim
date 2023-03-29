#include "StackingAction.hh"
#include "SessionManager.hh"
#include "G4Positron.hh"

#include "G4Track.hh"

G4ClassificationOfNewTrack PesGeneratorStackingAction::ClassifyNewTrack(const G4Track * track)
{
    if (track->GetParentID() == 0) return fUrgent;
    return fKill;
}

G4ClassificationOfNewTrack AnnihilationLoggerStackingAction::ClassifyNewTrack(const G4Track * track)
{
    if (track->GetParentID() == 0) return fUrgent;
    if (track->GetParticleDefinition() == G4Positron::Definition()) return fUrgent;
    return fKill;
}

#include "SimMode.hh"
G4ClassificationOfNewTrack PositronTimeLoggerStackingAction::ClassifyNewTrack(const G4Track * track)
{
    if (track->GetParticleDefinition() != G4Positron::Definition()) return fUrgent;

    SessionManager & SM = SessionManager::getInstance();
    ModePositronTimeLogger * Mode = static_cast<ModePositronTimeLogger*>(SM.SimMode);
    Mode->fillTime(track->GetGlobalTime());
    return fKill;
}
