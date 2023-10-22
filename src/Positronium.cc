#include "Positronium.hh"
#include "G4Geantino.hh"
#include "jstools.hh"
#include "SourcePositronium.hh"

Positronium::Positronium(double threeGammaDecayFraction) :
    ParticleBase("Positronium"), ThreeGammaDecayFraction(threeGammaDecayFraction)
{
    StandardPrimaryGeneration = false;
}

Positronium::Positronium(double threeGammaDecayFraction, const std::string & fileName_GammaEnergyPositions) :
    ParticleBase("Positronium"), ThreeGammaDecayFraction(threeGammaDecayFraction), LogFileName(fileName_GammaEnergyPositions)
{
    StandardPrimaryGeneration = false;
}

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
    Source = new SourcePositronium(ThreeGammaDecayFraction, nullptr, LogFileName);
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
