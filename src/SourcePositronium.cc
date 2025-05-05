#include "SourcePositronium.hh"
#include "DefinedParticles.hh"
#include "SessionManager.hh"
#include "PositroniumDecayChannel.hh"
#include "TimeGenerator.hh"
#include "jstools.hh"
#include "out.hh"

#include <string>

#include "G4ParticleGun.hh"
#include "G4LorentzVector.hh"
#include "G4Event.hh"
#include "G4RandomTools.hh"
#include "G4Gamma.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"

SourcePositronium::SourcePositronium(double orthoDecayFraction, TimeGeneratorBase * timeGenerator, const std::string & fileName_GeneratedGammas) :
    SourceModeBase(nullptr, timeGenerator), OrthoFraction(orthoDecayFraction), GammaOutputFileName(fileName_GeneratedGammas)
{
    if (!fileName_GeneratedGammas.empty()) setupStream(fileName_GeneratedGammas);
}

SourcePositronium::SourcePositronium(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);
}

void SourcePositronium::setupStream(const std::string & baseFileName)
{
    Stream = new std::ofstream();
    SessionManager & SM = SessionManager::getInstance();
    Stream->open(SM.WorkingDirectory + '/' + baseFileName);
    if (!Stream->is_open() || Stream->fail() || Stream->bad())
    {
        out("Cannot open file to store emission positions!");
        delete Stream; Stream = nullptr;
    }
    else out("\nSaving source emission positions to file", baseFileName);
}

SourcePositronium::~SourcePositronium()
{
    if (Stream) Stream->close();
    delete Stream; Stream = nullptr;
}

void SourcePositronium::GeneratePrimaries(G4Event * anEvent)
{
    double particle_time = TimeGenerator->generateTime();

    //the rest inherited frpom: pModel->GeneratePrimaryVertices( event, particle_time, particle_position);
    //selection of type from PreparePositroniumParametrization() is directly in GetPrimaryVertexFromPositroniumAnnihilation

    if (AddPromptGamma)
        anEvent->AddPrimaryVertex( GetPrimaryVertexFromDeexcitation(particle_time, Position) );

    anEvent->AddPrimaryVertex( GetPrimaryVertexFromPositroniumAnnihilation(particle_time, Position) );
}

void SourcePositronium::setTimeGenerator(TimeGeneratorBase * generator)
{
    TimeGenerator = generator;
}

void SourcePositronium::setNa22Origin()
{
    AddPromptGamma = true;
    PromptGammaEnergy = 1.275*MeV;
}

void SourcePositronium::doWriteToJson(json11::Json::object & json) const
{
    json["OrthoFraction"]       = OrthoFraction;
    json["AddPromptGamma"]      = AddPromptGamma;
    json["PromptGammaEnergy"]   = PromptGammaEnergy;
    json["GammaOutputFileName"] = GammaOutputFileName;
}

void SourcePositronium::doReadFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "OrthoFraction",       OrthoFraction);
    jstools::readBool  (json, "AddPromptGamma",      AddPromptGamma);
    jstools::readDouble(json, "PromptGammaEnergy",   PromptGammaEnergy);
    jstools::readString(json, "GammaOutputFileName", GammaOutputFileName);
}

void SourcePositronium::customPostInit()
{
    OrthoDecayChannel = new PositroniumDecayChannel(true);
    ParaDecayChannel  = new PositroniumDecayChannel(false);
}

G4PrimaryVertex * SourcePositronium::GetPrimaryVertexFromDeexcitation(double particle_time, const  G4ThreeVector & particle_position )
{
    G4PrimaryVertex* vertex = new G4PrimaryVertex(particle_position, particle_time);
    vertex->SetPrimary( GetGammaFromDeexcitation() );
    return vertex;
}

G4PrimaryVertex * SourcePositronium::GetPrimaryVertexFromPositroniumAnnihilation(G4double particle_time, const G4ThreeVector & particle_position)
{
    bool ThreeGammas;
    if      (OrthoFraction == 1.0) ThreeGammas = true;
    else if (OrthoFraction == 0.0) ThreeGammas = false;
    else ThreeGammas = (G4UniformRand() < OrthoFraction);

    double lifetime = (ThreeGammas ? OrthoLifetime : ParaLifetime) / log(2); // from halftime to decay time
    G4double shifted_particle_time = particle_time + G4RandExponential::shoot(lifetime);

    G4PrimaryVertex * vertex = new G4PrimaryVertex(particle_position, shifted_particle_time);
    std::vector<G4PrimaryParticle*> gammas = GetGammasFromPositroniumAnnihilation(ThreeGammas);
    std::for_each( gammas.begin(), gammas.end(), [&](G4PrimaryParticle * gamma) { vertex->SetPrimary(gamma); } );
    return vertex;
}

G4PrimaryParticle * SourcePositronium::GetGammaFromDeexcitation()
{
    //inherited from: G4PrimaryParticle* gamma = GetSingleGamma( fPromptGammaEnergy );

    G4PrimaryParticle * gamma = new G4PrimaryParticle( G4Gamma::Definition() );

    G4ThreeVector momentum_direction = GetUniformOnSphere();

    G4LorentzVector lv_gamma( momentum_direction.x(), momentum_direction.y(), momentum_direction.z(), 1.0 );
    lv_gamma *= PromptGammaEnergy;

    gamma->Set4Momentum( lv_gamma.px(), lv_gamma.py(), lv_gamma.pz(), lv_gamma.e() );
    gamma->SetPolarization( GetPolarization( gamma->GetMomentumDirection() ) );
    //gamma->SetUserInformation( GetPrimaryParticleInformation( gamma, GateEmittedGammaInformation::GammaKind::Prompt ) );

    if (Stream)
        *Stream << lv_gamma.px() << " "
                << lv_gamma.py() << " "
                << lv_gamma.pz() << " "
                << lv_gamma.e()  << ':';

    return gamma;
}

std::vector<G4PrimaryParticle*> SourcePositronium::GetGammasFromPositroniumAnnihilation(bool threeGammas)
{
    const size_t numGammas = (threeGammas ? 3 : 2);
    std::vector<G4PrimaryParticle*> gammas(numGammas) ;

    PositroniumDecayChannel * decayChannel = (threeGammas ? OrthoDecayChannel : ParaDecayChannel);
    G4DecayProducts * decay_products = decayChannel->DecayIt(0);

    for (size_t i = 0; i < numGammas; i++)
    {
        G4PrimaryParticle * gamma = new G4PrimaryParticle( G4Gamma::Definition() );

        G4DynamicParticle * dynamic_gamma = (*decay_products)[i];
        G4LorentzVector lv = dynamic_gamma->Get4Momentum();
        gamma->Set4Momentum( lv.px(), lv.py(), lv.pz(), lv.e() );
        gamma->SetPolarization( dynamic_gamma->GetPolarization() );
        //gamma->SetUserInformation( GetPrimaryParticleInformation(  gamma, GateEmittedGammaInformation::GammaKind::Annihilation ) );
        gammas[i] = gamma;

        if (Stream)
        {
            *Stream << lv.px() << " "
                    << lv.py() << " "
                    << lv.pz() << " "
                    << lv.e();
            if (i != numGammas-1) *Stream << ':';
        }
    }
    delete decay_products;

    if (Stream) *Stream << '\n';

    return gammas;
}

// ------ misc -------

G4ThreeVector SourcePositronium::GetUniformOnSphere() const
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

G4ThreeVector SourcePositronium::GetPolarization(const G4ThreeVector & momentum) const
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

G4ThreeVector SourcePositronium::GetPerpendicularVector(const G4ThreeVector & v) const
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
