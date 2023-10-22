#ifndef positronium_h
#define positronium_h

#include "DefinedParticles.hh"
#include "json11.hh"

#include <string>

class SourcePositronium;

class Positronium : public ParticleBase
{
public:
    Positronium(double threeGammaDecayFraction);
    Positronium(double threeGammaDecayFraction, const std::string & fileName_GammaEnergyPositions);
    ~Positronium();

    G4ParticleDefinition * getParticleDefinition() const override;

    std::string getTypeName() const override {return "Positronium";}

    void customInit() override;
    void generatePrimaries(G4ParticleGun * particleGun, TimeGeneratorBase * timeGenerator, G4Event * anEvent) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json) override;

    double ThreeGammaDecayFraction = 1.0;
    std::string LogFileName;

    SourcePositronium * Source = nullptr;
};

#endif // positronium_h
