#include "DefinedParticles.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

G4ParticleDefinition * Geantino::getParticleDefinition() const
{
    G4ParticleTable * ptab = G4ParticleTable::GetParticleTable();
    return ptab->FindParticle("geantino");
}

G4ParticleDefinition * Gamma::getParticleDefinition() const
{
    G4ParticleTable * ptab = G4ParticleTable::GetParticleTable();
    return ptab->FindParticle("gamma");
}

G4ParticleDefinition * GammaPair::getParticleDefinition() const
{
    G4ParticleTable * ptab = G4ParticleTable::GetParticleTable();
    return ptab->FindParticle("gamma");
}

Isotope::Isotope(int z, int a) :
    ParticleBase("ion_" + std::to_string(z) + "_" + std::to_string(a), true), Z(z), A(a)
{
    Energy = 0;
}

Isotope::Isotope(int z, int a, double excitationEnergy) :
    ParticleBase("ion_" + std::to_string(z) + "_" + std::to_string(a) + '[' + std::to_string(excitationEnergy) + ']', true),
    Z(z), A(a), ExcitationEnergy(excitationEnergy)
{
    Energy = 0;
}

G4ParticleDefinition * Isotope::getParticleDefinition() const
{
    G4IonTable * itab = G4IonTable::GetIonTable();
    return itab->GetIon(Z, A, ExcitationEnergy);
}

void Isotope::doWriteToJson(json11::Json::object & json) const
{
    json["Z"] = Z;
    json["A"] = A;
    json["ExcitationEnergy"] = ExcitationEnergy;
}

// ---

G4ParticleDefinition *Proton::getParticleDefinition() const
{
    G4ParticleTable * ptab = G4ParticleTable::GetParticleTable();
    return ptab->FindParticle("proton");
}

void ParticleBase::writeToJson(json11::Json::object & json) const
{
    json["Type"] = getTypeName();
    json["Energy"] = Energy;
    doWriteToJson(json);
}
