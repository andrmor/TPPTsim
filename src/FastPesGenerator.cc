#include "FastPesGenerator.hh"
#include "SessionManager.hh"
#include "G4Proton.hh"
#include "G4Track.hh"
#include "G4RandomTools.hh"
#include "out.hh"
#include "Randomize.hh"
#include "PesGenerationMode.hh"
#include "SessionManager.hh"

G4bool FastPesGeneratorModel::IsApplicable(const G4ParticleDefinition & particle)
{
    const bool applicable = (&particle == G4Proton::ProtonDefinition());
    //out("PES isApplicable called", particle.GetParticleName(), applicable);
    return applicable;
}

G4bool FastPesGeneratorModel::ModelTrigger(const G4FastTrack & fastTrack)
{
    const G4Track * track = fastTrack.GetPrimaryTrack();
    if (track->GetParentID() != 0) return false;
    const int StepNumber = track->GetCurrentStepNumber();
    //out("PES call", StepNumber);

    if (StepNumber == 1)
    {
        LastEnergy      = track->GetKineticEnergy();
        LastTrackLength = track->GetTrackLength();
        LastPosition    = track->GetPosition();
        LastMaterial    = track->GetMaterial()->GetIndex();
        return false;
    }

    double        Energy   = track->GetKineticEnergy();
    double        Length   = track->GetTrackLength();
    G4ThreeVector Position = track->GetPosition();

    if (LastEnergy > Energy)
    {
        double stepLength = Length - LastTrackLength;
        double meanEnergy = 0.5 * (Energy + LastEnergy);
        //out("Step", stepLength, "MeanEenergy", meanEnergy, " Material index", LastMaterial);

        SessionManager & SM = SessionManager::getInstance();
        const PesGenerationMode * PGM = static_cast<PesGenerationMode*>(SM.SimMode);

        const std::vector<PesGenRecord> & Records = PGM->MaterialRecords[LastMaterial];
        if (Records.empty()) return false;

        //probability is proportional to CS * NumberDensity
        ProbVec.clear(); ProbVec.reserve(Records.size());
        double sumProb = 0;
        for (const PesGenRecord & r : Records)
        {
            const double cs = r.getCrossSection(meanEnergy);
            const double relProb = cs * r.NumberDensity;
            sumProb += relProb;
            ProbVec.push_back(relProb);
        }

        // selecting the reaction
        size_t index = 0;
        double val = sumProb * G4UniformRand();
        for (; index+1 < ProbVec.size(); index++)
        {
            if (val < ProbVec[index]) break;
            val -= ProbVec[index];
        }

        const double mfp = 1e25 / ProbVec[index]; // millibarn = 0.001e-28m2 -> 0.001e-22mm2 -> 1e-25 mm2

        double trigStep = -mfp * log(G4UniformRand());
        if (trigStep < stepLength)
        {
            G4ThreeVector TriggerPosition = LastPosition + trigStep/stepLength*(Position - LastPosition);
            //out("Triggered! Position:", TriggerPosition);
            PGM->saveRecord(Records[index].PES, TriggerPosition[0], TriggerPosition[1], TriggerPosition[2], track->GetGlobalTime());
            return true;
        }
    }
    //else out("Zero dEnergy");

    LastEnergy      = Energy;
    LastTrackLength = Length;
    LastPosition    = Position;
    LastMaterial    = track->GetMaterial()->GetIndex();
    return false;
}

void FastPesGeneratorModel::DoIt(const G4FastTrack &, G4FastStep & step)
{
    //out("PES doit");
    step.KillPrimaryTrack();
}
