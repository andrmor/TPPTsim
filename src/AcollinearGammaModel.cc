#include "AcollinearGammaModel.hh"
#include "SessionManager.hh"
#include "G4Gamma.hh"
#include "G4Track.hh"
#include "G4RandomTools.hh"
#include "out.hh"

G4bool AcollinearGammaModel::IsApplicable(const G4ParticleDefinition & particle)
{
    bool applicable = (&particle == G4Gamma::GammaDefinition());
    out("isApplicable called", particle.GetParticleName(), applicable);
    return applicable;
}

G4bool AcollinearGammaModel::ModelTrigger(const G4FastTrack & track)
{
    bool flag = (track.GetPrimaryTrack()->GetCurrentStepNumber() == 1);
    out("ModelTrigger called", flag);
    return flag;
}

void AcollinearGammaModel::DoIt(const G4FastTrack & fastTrack, G4FastStep &step)
{
    out("DoIt called");

    SessionManager & SM = SessionManager::getInstance();
    const G4Track * track = fastTrack.GetPrimaryTrack();

        int Id       = track->GetTrackID();
        int parentId = track->GetParentID();

        out(Id, parentId, PrevID, PrevParentID);

        if (Id == (PrevID - 1) && parentId == PrevParentID)
        {
            out("Rotating momentum direction for this gamma!");
            G4ThreeVector v = track->GetMomentumDirection();
            G4ThreeVector vCopy(v);
            constexpr double Sigma = 0.5*deg / 2.35482;
            double angle = G4RandGauss::shoot(0, Sigma);
            v.rotate(angle, v.orthogonal());
            v.rotate(G4UniformRand()*2.0*M_PI, vCopy);
            step.ProposePrimaryTrackFinalMomentumDirection (v, false);

            PrevID = -1;
            PrevParentID = -1;
        }
        else
        {
            PrevID       = Id;
            PrevParentID = parentId;
        }
}
