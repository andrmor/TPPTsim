#include "TrackingAction.hh"
#include "SessionManager.hh"
#include "out.hh"

#include "G4Track.hh"
#include "G4VProcess.hh"

void TrackingAction::PreUserTrackingAction(const G4Track *track)
{
    SessionManager & SM = SessionManager::getInstance();

    if (track->GetParticleDefinition() != SM.GammaPD) return;

    const G4VProcess * proc = track->GetCreatorProcess();
    if (!proc) return;

    if (proc->GetProcessSubType() != 5) return; // need to confirm: 5 is always annihilation?

    //out("-->", track->GetParticleDefinition()->GetParticleName(), proc->GetProcessName(), proc->GetProcessType(), proc->GetProcessSubType());

    int Id       = track->GetTrackID();
    int parentId = track->GetParentID();

    //out(Id, parentId, PrevID, PrevParentID);

    if (Id == (PrevID - 1) && parentId == PrevParentID)
    {
        //out("Rotating momentum direction for this gamma!");
        G4ThreeVector v = track->GetMomentumDirection();
        G4ThreeVector vCopy(v);
        constexpr double Sigma = 0.5*deg / 2.35482;
        double angle = G4RandGauss::shoot(0, Sigma);
        v.rotate(angle, v.orthogonal());
        v.rotate(G4UniformRand()*2.0*M_PI, vCopy);
        const_cast<G4Track*>(track)->SetMomentumDirection(v);

        PrevID = -1;
        PrevParentID = -1;
    }
    else
    {
        PrevID       = Id;
        PrevParentID = parentId;
    }
}
