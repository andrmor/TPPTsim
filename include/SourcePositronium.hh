#ifndef sourcepositronium_h
#define sourcepositronium_h

#include "SourceMode.hh"

#include <vector>

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

class G4PrimaryVertex;
class G4PrimaryParticle;
class PositroniumDecayChannel;

class SourcePositronium : public SourceModeBase
{
public:
    // orthoDecayFraction defines fraction (from 0 to 1) of events from ortho-Positronium (3 annihil gammas)
    // if fileName is defined, the generated gammas are saved to the file in the working directory, one event per line with the format:
    // ix1 iy1 iz1 e1[MeV]:ix2 iy2 iz2 e2[MeV] etc
    SourcePositronium(double orthoDecayFraction, TimeGeneratorBase * timeGenerator, const std::string & fileName_GeneratedGammas = "");
    SourcePositronium(const json11::Json & json);
    ~SourcePositronium();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "SourcePositronium";}

    void setTimeGenerator(TimeGeneratorBase * generator);
    G4ThreeVector Position = {0,0,0};

    void setNa22Origin();

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    void customPostInit() override;
    void setupStream(const std::string & baseFileName);

    double OrthoFraction = 0.5;
    bool   AddPromptGamma = false;
    double PromptGammaEnergy = 1.275*MeV;

    std::string GammaOutputFileName;

    const G4double ParaLifetime  = 0.1244*ns;
    const G4double OrthoLifetime = 138.6*ns;

    PositroniumDecayChannel * OrthoDecayChannel = nullptr; // 3 annihilation gammas
    PositroniumDecayChannel * ParaDecayChannel  = nullptr; // 2 annihilation gammas

    std::ofstream * Stream = nullptr;

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
