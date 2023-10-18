#ifndef sourcepositronium_h
#define sourcepositronium_h

#include "SourceMode.hh"

#include <vector>

class Positronium;
class G4PrimaryVertex;
class G4PrimaryParticle;

class SourceThreeGammas : public SourceModeBase
{
public:
    enum DecayModel {Standard, WithPrompt};
    SourceThreeGammas(TimeGeneratorBase * timeGenerator);
    //SourceThreeGammas(const json11::Json & json);

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "SourceThreeGammas";}


protected:
    void doWriteToJson(json11::Json::object & json) const override {}
    //void doReadFromJson(const json11::Json & json);

    void customPostInit() override;

    DecayModel ModelType = Standard;
    double fPromptGammaEnergy = 1.275; // in MeV
    double fThreeGammaFraction = 1.0;

    Positronium * fParaPs  = nullptr;
    Positronium * fOrthoPs = nullptr;
    Positronium * pInfoPs  = nullptr;

private:
    G4PrimaryVertex * GetPrimaryVertexFromDeexcitation(double particle_time, const G4ThreeVector & particle_position);
    G4PrimaryVertex * GetPrimaryVertexFromPositroniumAnnihilation(G4double particle_time, const G4ThreeVector & particle_position);
    G4PrimaryParticle * GetGammaFromDeexcitation();
    std::vector<G4PrimaryParticle*> GetGammasFromPositroniumAnnihilation(bool threeGammas);
    G4ThreeVector GetUniformOnSphere() const;
    G4ThreeVector GetPolarization(const G4ThreeVector &momentum) const;
    G4ThreeVector GetPerpendicularVector(const G4ThreeVector &v) const;
};

class G4DecayProducts;
class G4VDecayChannel;
class Positronium
{
public:
    Positronium(G4String name, G4double life_time, G4int annihilation_gammas_number);
    void SetLifeTime(const G4double & life_time) {fLifeTime = life_time;}
    G4double GetLifeTime() const {return fLifeTime;}
    G4String GetName() const {return fName;}
    G4int GetAnnihilationGammasNumber() const {return fAnnihilationGammasNumber;}
    G4DecayProducts * GetDecayProducts();
private:
    G4String fName = "";
    G4double fLifeTime = 0.0; //[ns]
    G4int fAnnihilationGammasNumber = 0;
    G4VDecayChannel * pDecayChannel = nullptr;
};

enum PositroniumKind {pPs, oPs};
/*
 Decay model descibes decay of positronium.
    * For example for PositroniumKind::pPs we have decays:
    * 1) Standard: pPs -> 2 gamma
    * 2) WithPrompt: pPs* -> 2 gamma + prompt_gamma
    * Only prompt gamma has nod modified time, other gammas always has modified time.
*/

#endif // sourcepositronium_h
