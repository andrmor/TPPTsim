#include "FastPesGenerator.hh"
#include "SessionManager.hh"
#include "G4Proton.hh"
#include "G4Track.hh"
#include "G4RandomTools.hh"
#include "out.hh"
#include "Randomize.hh"

G4bool FastPesGeneratorModel::IsApplicable(const G4ParticleDefinition & particle)
{
    const bool applicable = (&particle == G4Proton::ProtonDefinition());
    out("PES isApplicable called", particle.GetParticleName(), applicable);
    return applicable;
}

G4bool FastPesGeneratorModel::ModelTrigger(const G4FastTrack & fastTrack)
{
    const G4Track * track = fastTrack.GetPrimaryTrack();
    if (track->GetParentID() != 0) return false;
    const int StepNumber = track->GetCurrentStepNumber();
    out("PES call", StepNumber);

    if (StepNumber == 1)
    {
        LastEnergy      = track->GetKineticEnergy();
        LastTrackLength = track->GetTrackLength();
        LastPosition    = track->GetPosition();
        return false;
    }

    double        Energy   = track->GetKineticEnergy();
    double        Length   = track->GetTrackLength();
    G4ThreeVector Position = track->GetPosition();

    if (LastEnergy > Energy)
    {
        double stepLength = Length - LastTrackLength;
        double meanEnergy = 0.5 * (Energy + LastEnergy);
        out("Step", stepLength, "meanE", meanEnergy);

        double trigStep = -Lambda * log(G4UniformRand());
        if (trigStep < stepLength)
        {
            G4ThreeVector TriggerPosition = LastPosition + trigStep/stepLength*(Position - LastPosition);
            out("Triggered! Position:", TriggerPosition);

            // tmp
            sumLength += LastTrackLength + trigStep/stepLength*(Length - LastTrackLength);
            num++;
            out("Estimated MFP (need good statistics):", sumLength/num);

            return true;
        }
    }
    else out("Zero dEnergy");

    LastEnergy      = Energy;
    LastTrackLength = Length;
    LastPosition    = Position;
    return false;
}

void FastPesGeneratorModel::DoIt(const G4FastTrack &, G4FastStep & step)
{
    out("PES doit");
    step.KillPrimaryTrack();
}
