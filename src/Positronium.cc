#include "Positronium.hh"
#include "G4Geantino.hh"
#include "jstools.hh"
#include "SourcePositronium.hh"

Positronium::Positronium(double threeGammaDecayFraction) :
    ParticleBase("Positronium"), ThreeGammaDecayFraction(threeGammaDecayFraction){}

Positronium::Positronium(double threeGammaDecayFraction, const std::string & fileName_LogGammaMomentumAndEnergy) :
    ParticleBase("Positronium"), ThreeGammaDecayFraction(threeGammaDecayFraction), LogFileName(fileName_LogGammaMomentumAndEnergy){}

Positronium::Positronium(double threeGammaDecayFraction, bool assumeNa22Origin, const std::string &fileName_LogGammaMomentumAndEnergy) :
    ParticleBase("Positronium"), ThreeGammaDecayFraction(threeGammaDecayFraction), Na22origin(assumeNa22Origin), LogFileName(fileName_LogGammaMomentumAndEnergy){}

Positronium::~Positronium()
{
    delete Source; Source = nullptr;
}

G4ParticleDefinition * Positronium::getParticleDefinition() const
{
    return G4Geantino::Definition();
}

void Positronium::customInit()
{
    StandardPrimaryGeneration = false;
    Source = new SourcePositronium(ThreeGammaDecayFraction, nullptr, LogFileName);
    if (Na22origin) Source->setNa22Origin();
    Source->initialize();
}

void Positronium::doWriteToJson(json11::Json::object & json) const
{
    json["ThreeGammaDecayFraction"] = ThreeGammaDecayFraction;
    json["LogFileName"] = LogFileName;
}

void Positronium::doReadFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "ThreeGammaDecayFraction", ThreeGammaDecayFraction);
    jstools::readString(json, "LogFileName", LogFileName);
}

#include "G4ParticleGun.hh"
void Positronium::generatePrimaries(G4ParticleGun * particleGun, TimeGeneratorBase * timeGenerator, G4Event * anEvent)
{
    Source->setTimeGenerator(timeGenerator);
    Source->Position = particleGun->GetParticlePosition();
    Source->GeneratePrimaries(anEvent);
}
