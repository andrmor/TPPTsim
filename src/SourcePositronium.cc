#include "SourcePositronium.hh"
#include "DefinedParticles.hh"

#include "G4ParticleGun.hh"

SourceThreeGammas::SourceThreeGammas(TimeGeneratorBase * timeGenerator) :
    SourceModeBase(new Gamma(0.511*MeV), timeGenerator)
{
    ParticleGun->SetParticlePosition({0,0,0});
}

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TimeGenerator.hh"

void SourceThreeGammas::GeneratePrimaries(G4Event * anEvent)
{
    ParticleGun->SetParticleTime(TimeGenerator->generateTime());

    //ROOT example from  http://physik.uibk.ac.at/hephy/praktikum/tau/TGenPhaseSpace_docu.pdf
    //see also https://root.cern.ch/doc/master/PhaseSpace_8C.html
    //Momentum units GeV/c, Energy units GeV

    TLorentzVector target(0.0, 0.0, 0.0, 0.938);
    TLorentzVector beam(0.0, 0.0, .4, .4);
    TLorentzVector W = beam + target;
    Double_t masses[2] = { 0.938, 0.135};
    TGenPhaseSpace event;
    event.SetDecay(W, 2, masses);

    event.Generate();

    for (size_t iPart = 0; iPart < 2; iPart++)
    {
        TLorentzVector * part = event.GetDecay(iPart);

        ParticleGun->SetParticleEnergy(part->Energy()*GeV);

        const TVector3 dir = part->Vect().Unit();
        ParticleGun->SetParticleMomentumDirection({dir[0], dir[1], dir[2]});

        ParticleGun->GeneratePrimaryVertex(anEvent);
    }
}
