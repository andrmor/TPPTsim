#include "TrackingAction.hh"
#include "SessionManager.hh"
#include "ModePesGenerator_MC.hh"
#include "ModeAnnihilationLogger.hh"
#include "out.hh"

#include "G4Track.hh"
#include "G4Positron.hh"
#include "G4VProcess.hh"

void PesGeneratorTrackingAction::PreUserTrackingAction(const G4Track * /*track*/)
{
    SessionManager & SM = SessionManager::getInstance();
    ModePesGenerator_MC * Mode = static_cast<ModePesGenerator_MC*>(SM.SimMode);
    Mode->bNewTrackStarted = true;
}

void AnnihilationLoggerTrackingAction::PostUserTrackingAction(const G4Track * track)
{
    if (track->GetParticleDefinition() != G4Positron::Definition()) return;

    //out(track->GetPosition()[0], track->GetPosition()[1], track->GetPosition()[2]);

    //const G4VProcess * proc = track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep();
    //std::string procName = (proc ? proc->GetProcessName() : "undefined");
    //if (procName != "annihil") out(procName);

    SessionManager & SM = SessionManager::getInstance();
    ModeAnnihilationLogger * mode = static_cast<ModeAnnihilationLogger*>(SM.SimMode);
    mode->fillPosition(track->GetPosition());
}

void SourceTester_TrackingAction::PreUserTrackingAction(const G4Track *track)
{
    if (track->GetParentID() == 0)
    {
        SessionManager & SM = SessionManager::getInstance();
        SourceTester * Mode = static_cast<SourceTester*>(SM.SimMode);

        Mode->registerParticle(track->GetParticleDefinition()->GetParticleName(), track->GetGlobalTime());
    }
    const_cast<G4Track*>(track)->SetTrackStatus(fStopAndKill);
}
