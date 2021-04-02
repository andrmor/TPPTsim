#include "DefinedParticles.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

G4ParticleDefinition * ParticleGeantino::getParticleDefinition() const
{
    G4ParticleTable * ptab = G4ParticleTable::GetParticleTable();
    return ptab->FindParticle("geantino");
}

G4ParticleDefinition * ParticleGamma::getParticleDefinition() const
{
    G4ParticleTable * ptab = G4ParticleTable::GetParticleTable();
    return ptab->FindParticle("gamma");
}

G4ParticleDefinition * ParticleGammaPair::getParticleDefinition() const
{
    G4ParticleTable * ptab = G4ParticleTable::GetParticleTable();
    return ptab->FindParticle("gamma");
}

ParticlePES::ParticlePES(int z, int a) :
    ParticleBase("ion_" + std::to_string(z) + "_" + std::to_string(a), true), Z(z), A(a)
{
    Energy = 0;
}

G4ParticleDefinition * ParticlePES::getParticleDefinition() const
{
    G4IonTable * itab = G4IonTable::GetIonTable();
    return itab->GetIon(Z, A, 0);
}

// ---

G4ParticleDefinition *ParticleProton::getParticleDefinition() const
{
    G4ParticleTable * ptab = G4ParticleTable::GetParticleTable();
    return ptab->FindParticle("proton");
}
