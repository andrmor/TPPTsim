#ifndef sourcepositronium_h
#define sourcepositronium_h

#include "SourceMode.hh"

#include <vector>

#include "G4SystemOfUnits.hh"

class Positronium;
class G4PrimaryVertex;
class G4PrimaryParticle;
class PositroniumDecayChannel;

class SourcePositronium : public SourceModeBase
{
public:
    enum DecayModel {Standard, WithPrompt};
    SourcePositronium(double orthoDecayFraction, TimeGeneratorBase * timeGenerator);
    //SourcePositronium(const json11::Json & json);

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "SourcePositronium";}

protected:
    void doWriteToJson(json11::Json::object & json) const override {}
    //void doReadFromJson(const json11::Json & json);

    void customPostInit() override;

    DecayModel ModelType = Standard;
    double fPromptGammaEnergy = 1.275*MeV;
    double fThreeGammaFraction = 0.5;

    const G4double ParaLifetime  = 0.1244*ns;
    const G4double OrthoLifetime = 138.6*ns;

    PositroniumDecayChannel * OrthoDecayChannel = nullptr; // 3 annihilation gammas
    PositroniumDecayChannel * ParaDecayChannel  = nullptr; // 2 annihilation gammas

private:
    G4PrimaryVertex * GetPrimaryVertexFromDeexcitation(double particle_time, const G4ThreeVector & particle_position);
    G4PrimaryVertex * GetPrimaryVertexFromPositroniumAnnihilation(G4double particle_time, const G4ThreeVector & particle_position);
    G4PrimaryParticle * GetGammaFromDeexcitation();
    std::vector<G4PrimaryParticle*> GetGammasFromPositroniumAnnihilation(bool threeGammas);
    G4ThreeVector GetUniformOnSphere() const;
    G4ThreeVector GetPolarization(const G4ThreeVector & momentum) const;
    G4ThreeVector GetPerpendicularVector(const G4ThreeVector & v) const;
};

#endif // sourcepositronium_h
