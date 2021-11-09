#include "StackingAction.hh"

#include "G4Track.hh"

G4ClassificationOfNewTrack PesGeneratorStackingAction::ClassifyNewTrack(const G4Track * track)
{
    if (track->GetParentID() == 0) return fUrgent;
    return fKill;
}
