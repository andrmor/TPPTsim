#include "SourcePositronium.hh"
#include "DefinedParticles.hh"
#include "out.hh"

#include "G4ParticleGun.hh"
#include "G4LorentzVector.hh"

#include "G4Event.hh"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TimeGenerator.hh"
#include "G4RandomTools.hh"
#include "G4Gamma.hh"

SourceThreeGammas::SourceThreeGammas(TimeGeneratorBase * timeGenerator) :
    SourceModeBase(new Gamma(0.511*MeV), timeGenerator)
{
    //ParticleGun->SetParticlePosition({0,0,0});
}

void SourceThreeGammas::GeneratePrimaries(G4Event * anEvent)
{
    double particle_time = TimeGenerator->generateTime();
    const G4ThreeVector particle_position = {0,0,0};

    //the rest inherited frpom: pModel->GeneratePrimaryVertices( event, particle_time, particle_position);
    //selection of type from PreparePositroniumParametrization() is directly in GetPrimaryVertexFromPositroniumAnnihilation

    if (ModelType == WithPrompt)
        anEvent->AddPrimaryVertex( GetPrimaryVertexFromDeexcitation(particle_time, particle_position) );

    anEvent->AddPrimaryVertex( GetPrimaryVertexFromPositroniumAnnihilation(particle_time, particle_position) );
}

void SourceThreeGammas::customPostInit()
{
    //fParaPs  = new Positronium("pPs", 0.1244*ns, 2);
    //fOrthoPs = new Positronium("oPs", 138.6*ns, 3);
}

G4PrimaryVertex * SourceThreeGammas::GetPrimaryVertexFromDeexcitation(double particle_time, const  G4ThreeVector & particle_position )
{
    G4PrimaryVertex* vertex = new G4PrimaryVertex(particle_position, particle_time);
    vertex->SetPrimary( GetGammaFromDeexcitation() );
    return vertex;
}

G4PrimaryVertex * SourceThreeGammas::GetPrimaryVertexFromPositroniumAnnihilation(G4double particle_time, const G4ThreeVector & particle_position)
{
    bool ThreeGammas;
    if      (fThreeGammaFraction == 1.0) ThreeGammas = true;
    else if (fThreeGammaFraction == 0.0) ThreeGammas = false;
    else ThreeGammas = (G4UniformRand() < fThreeGammaFraction);

    G4double shifted_particle_time = particle_time;// + G4RandExponential::shoot( pInfoPs->GetLifeTime() );  !!!***

    G4PrimaryVertex * vertex = new G4PrimaryVertex(particle_position, shifted_particle_time);
    std::vector<G4PrimaryParticle*> gammas = GetGammasFromPositroniumAnnihilation(ThreeGammas);
    std::for_each( gammas.begin(), gammas.end(), [&](G4PrimaryParticle * gamma) { vertex->SetPrimary(gamma); } );
    return vertex;
}

G4PrimaryParticle * SourceThreeGammas::GetGammaFromDeexcitation()
{
    //inherited from: G4PrimaryParticle* gamma = GetSingleGamma( fPromptGammaEnergy );

    G4PrimaryParticle * gamma = new G4PrimaryParticle( G4Gamma::Definition() );

    G4ThreeVector momentum_direction = GetUniformOnSphere();

    G4LorentzVector lv_gamma( momentum_direction.x(), momentum_direction.y(), momentum_direction.z(), 1.0 );
    lv_gamma *= fPromptGammaEnergy;

    gamma->Set4Momentum( lv_gamma.px(), lv_gamma.py(), lv_gamma.pz(), lv_gamma.e() );
    gamma->SetPolarization( GetPolarization( gamma->GetMomentumDirection() ) );
    //gamma->SetUserInformation( GetPrimaryParticleInformation( gamma, GateEmittedGammaInformation::GammaKind::Prompt ) );

    return gamma;
}

#include "G4DynamicParticle.hh"
std::vector<G4PrimaryParticle*> SourceThreeGammas::GetGammasFromPositroniumAnnihilation(bool threeGammas)
{
    const size_t numGammas = (threeGammas ? 3 : 2);
    std::vector<G4PrimaryParticle*> gammas(numGammas) ;

    G4DecayProducts * decay_products = pInfoPs->GetDecayProducts();

    for (size_t i = 0; i < numGammas; i++)
    {
        G4PrimaryParticle * gamma = new G4PrimaryParticle( G4Gamma::Definition() );

        G4DynamicParticle * dynamic_gamma = (*decay_products)[i];
        G4LorentzVector lv = dynamic_gamma->Get4Momentum();
        gamma->Set4Momentum( lv.px(), lv.py(), lv.pz(), lv.e() );
        gamma->SetPolarization( dynamic_gamma->GetPolarization() );
        //gamma->SetUserInformation( GetPrimaryParticleInformation(  gamma, GateEmittedGammaInformation::GammaKind::Annihilation ) );
        gammas[i] = gamma;
    }
    delete decay_products;

    return gammas;
}




// ----------------------

#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"

Positronium::Positronium(G4String name, G4double life_time, G4int annihilation_gammas_number) :
    fName(name), fLifeTime(life_time), fAnnihilationGammasNumber(annihilation_gammas_number)
{
    G4ParticleDefinition * positronium_def = G4ParticleTable::GetParticleTable()->FindParticle(name);
    if (!positronium_def)
    {
        out("positronium_def for name", name, "not found!");
        exit(111);
    }
    G4DecayTable * positronium_decay_table = positronium_def->GetDecayTable();
    pDecayChannel = positronium_decay_table->GetDecayChannel(0);
}

G4DecayProducts * Positronium::GetDecayProducts()
{
    return pDecayChannel->DecayIt();
}

// ------ misc -------

G4ThreeVector SourceThreeGammas::GetUniformOnSphere() const
{
    //Based on TRandom::Sphere
    G4double a = 0,b = 0, r2 = 1;
    while ( r2 > 0.25 )
    {
        a  = G4UniformRand() - 0.5;
        b  = G4UniformRand() - 0.5;
        r2 =  a*a + b*b;
    }

    G4double scale = 8.0 * sqrt(0.25 - r2);
    return G4ThreeVector( a * scale, b * scale, -1. + 8.0 * r2 );
}

G4ThreeVector SourceThreeGammas::GetPolarization(const G4ThreeVector & momentum) const
{
    G4ThreeVector polarization(0.0,0.0,0.0);

    G4ThreeVector a0,b0,d0;
    d0 = momentum.unit();
    a0 = GetPerpendicularVector( d0 ).unit();
    b0 = d0.cross( a0 ).unit();
    G4double angle_radians = G4UniformRand() * M_PI;
    polarization = std::cos( angle_radians ) * a0 + std::sin( angle_radians ) * b0;
    polarization.unit();
    return polarization;
}

G4ThreeVector SourceThreeGammas::GetPerpendicularVector(const G4ThreeVector & v) const
{
    G4double dx = v.x();
    G4double dy = v.y();
    G4double dz = v.z();

    G4double x = dx < 0.0 ? -dx : dx;
    G4double y = dy < 0.0 ? -dy : dy;
    G4double z = dz < 0.0 ? -dz : dz;

    if (x < y) { return x < z ? G4ThreeVector(-dy,dx,0) : G4ThreeVector(0,-dz,dy); }
    else { return y < z ? G4ThreeVector(dz,0,-dx) : G4ThreeVector(-dy,dx,0); }
}

