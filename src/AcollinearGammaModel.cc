#include "AcollinearGammaModel.hh"
#include "SessionManager.hh"
#include "G4Gamma.hh"
#include "G4Track.hh"
#include "G4RandomTools.hh"
#include "out.hh"

G4bool AcollinearGammaModel::IsApplicable(const G4ParticleDefinition & particle)
{
    const bool applicable = (&particle == G4Gamma::GammaDefinition());
    //out("Acollin isApplicable called", particle.GetParticleName(), applicable);
    return applicable;
}

G4bool AcollinearGammaModel::ModelTrigger(const G4FastTrack & fastTrack)
{
    const G4Track * track = fastTrack.GetPrimaryTrack();
    if (track->GetCurrentStepNumber() != 1) return false; // not the first step

    const double energy = track->GetKineticEnergy();
    if (energy < 0.510 || energy > 0.512) return false;  // annihilationis -> 0.510999

    const int Id       = track->GetTrackID();
    const int parentId = track->GetParentID();
    //out(Id, parentId, PrevID, PrevParentID);

    if (Id == (PrevID - 1) && parentId == PrevParentID)
    {
        PrevID = -1;
        PrevParentID = -1;
        return true;
    }
    else
    {
        PrevID       = Id;
        PrevParentID = parentId;
        return false;
    }
}

void AcollinearGammaModel::DoIt(const G4FastTrack & fastTrack, G4FastStep & step)
{
    //out("Rotating momentum direction for this gamma!");
    const G4Track * track = fastTrack.GetPrimaryTrack();
    G4ThreeVector v = track->GetMomentumDirection();
    G4ThreeVector vCopy(v);
    constexpr double Sigma = 0.5*deg / 2.35482;
    double angle = G4RandGauss::shoot(0, Sigma);
    v.rotate(angle, v.orthogonal());
    v.rotate(G4UniformRand()*2.0*M_PI, vCopy);
    step.ProposePrimaryTrackFinalMomentumDirection(v, false);

    //later add counter
    //SessionManager & SM = SessionManager::getInstance();
}
