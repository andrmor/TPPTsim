#include "DefinedParticles.hh"
#include "jstools.hh"
#include "out.hh"
#include "Positronium.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

ParticleBase * ParticleFactory::makeParticleInstance(const json11::Json & json)
{
    out("Reading particle json");
    std::string Type;
    jstools::readString(json, "Type", Type);

    ParticleBase * p = nullptr;

    if      (Type == "Geantino")    p = new Geantino();
    else if (Type == "Gamma")       p = new Gamma();
    else if (Type == "GammaPair")   p = new GammaPair();
    else if (Type == "Isotope")     p = new Isotope(0, 0, 0);
    else if (Type == "Proton")      p = new Proton();
    else if (Type == "Positronium") p = new Positronium(1.0);
    else
    {
        out("Unknown particle type!");
        exit(20);
    }

    p->readFromJson(json);
    return p;
}

// ---

void ParticleBase::writeToJson(json11::Json::object & json) const
{
    json["Type"] = getTypeName();
    json["Energy"] = Energy;
    doWriteToJson(json);
}

void ParticleBase::readFromJson(const json11::Json &json)
{
    jstools::readDouble(json, "Energy", Energy);
    doReadFromJson(json);
}

// ---

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
    ParticleBase("ion_" + std::to_string(z) + "_" + std::to_string(a)), Z(z), A(a)
{
    Energy = 0;
}

Isotope::Isotope(int z, int a, double excitationEnergy) :
    ParticleBase("ion_" + std::to_string(z) + "_" + std::to_string(a) + '[' + std::to_string(excitationEnergy) + ']'),
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

void Isotope::doReadFromJson(const json11::Json & json)
{
    jstools::readInt   (json, "Z", Z);
    jstools::readInt   (json, "A", A);
    jstools::readDouble(json, "ExcitationEnergy", ExcitationEnergy);
}

// ---

G4ParticleDefinition *Proton::getParticleDefinition() const
{
    G4ParticleTable * ptab = G4ParticleTable::GetParticleTable();
    return ptab->FindParticle("proton");
}

// ---

#include "G4ParticleGun.hh"
#include "TimeGenerator.hh"
#include "SourceMode.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4HadDecayGenerator.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
Na22_decay::Na22_decay() :
    Isotope(11,22)
{
    StandardPrimaryGeneration = false;
}

void Na22_decay::customInit()
{
    if (!Generator)
    {
        Generator = new G4HadDecayGenerator();
        G4IonTable * ionTable = G4IonTable::GetIonTable();
        Sodium_ion = ionTable->GetIon(11, 22);
        G4ParticleDefinition * neon_ion = ionTable->GetIon(10, 22, 1274.537*keV);

        Na22Mass = Sodium_ion->GetPDGMass();

        Masses.push_back( G4Positron::Definition()->GetPDGMass() );
        Masses.push_back( G4NeutrinoE::Definition()->GetPDGMass() );
        Masses.push_back( neon_ion->GetPDGMass() );

        //Generator->SetVerboseLevel(3);
    }
}

//#include <QDebug>
void Na22_decay::generatePrimaries(G4ParticleGun *particleGun, TimeGeneratorBase *timeGenerator, G4Event *anEvent)
{
    particleGun->SetParticleTime(timeGenerator->generateTime());

    particleGun->SetParticleDefinition( G4Gamma::Definition() );
    particleGun->SetParticleMomentumDirection(SourceModeBase::generateDirectionIsotropic());
    particleGun->SetParticleEnergy(1274.537*keV);
    particleGun->GeneratePrimaryVertex(anEvent);

    if (G4UniformRand() < 0.096) return; // no beta+

    //qDebug() << "Beta+ decay:";

    //bool ok = Generator->Generate(Na22Mass, Masses, FinalState);
    bool ok = Generator->Generate(Sodium_ion, Masses, FinalState);
    //qDebug() << ok;
    if (!ok)
    {
        out("Decay generation failed!");
        exit(1234);
    }

    //for (size_t i = 0; i < FinalState.size(); i++)
    //    qDebug() << (FinalState[i].e()-Masses[i])/keV << "mom:" << FinalState[i].px() << FinalState[i].py() << FinalState[i].pz();

    particleGun->SetParticleDefinition( G4Positron::Definition() );
    particleGun->SetParticleMomentumDirection(SourceModeBase::generateDirectionIsotropic());
    particleGun->SetParticleEnergy( FinalState[0].e() - Masses[0] );
    particleGun->GeneratePrimaryVertex(anEvent);

}
