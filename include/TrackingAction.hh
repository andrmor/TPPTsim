#ifndef TrackingAction_h
#define TrackingAction_h

#include "G4UserTrackingAction.hh"

class TrackingAction : public G4UserTrackingAction
{
public:
    void PreUserTrackingAction(const G4Track * track) override;

protected:
    int PrevID       = -1;
    int PrevParentID = -1;
};

#endif // TrackingAction_h
