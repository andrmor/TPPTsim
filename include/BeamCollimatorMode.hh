#include "json11.hh"

#include <string>
#include <vector>

#include "G4ThreeVector.hh"

class BeamCollimatorBase;
class G4LogicalVolume;

class BeamCollimatorFactory // when adding a new entry, do not forget to modify the implementation of the factory!
{
public:
    static BeamCollimatorBase * makeBeamCollimatorInstance(const json11::Json & json);
};

class BeamCollimatorBase
{
public:
    virtual ~BeamCollimatorBase(){}

    virtual std::string getTypeName() const = 0;
    virtual void defineBeamCollimator(G4LogicalVolume * /*logicWorld*/) = 0;

    void writeToJson(json11::Json::object & json) const;
    virtual void readFromJson(const json11::Json & /*json*/) {}

protected:
    virtual void doWriteToJson(json11::Json::object & /*json*/) const {}
};

// ---

class BeamCollimatorNone : public BeamCollimatorBase
{
public:
    std::string getTypeName() const override {return "BeamCollimatorNone";}
    void defineBeamCollimator(G4LogicalVolume *) override {}
};

// ---

class BeamCollimatorMarek : public BeamCollimatorBase
{
public:
    enum EOpeningOptions {Blind, Hole3, Holes369, Holes123, Cross, Ring};
    BeamCollimatorMarek(EOpeningOptions openingType, const G4ThreeVector & innerFacePosition, double rotationAngle = 0);

    std::string getTypeName() const override {return "BeamCollimatorMarek";}
    void defineBeamCollimator(G4LogicalVolume * logicWorld) override;

    void readFromJson(const json11::Json & json) override;

protected:
    EOpeningOptions OpeningType = Blind;
    G4ThreeVector   FacePosition = {0,0,-100.0};
    double          Angle = 0;

    void doWriteToJson(json11::Json::object & /*json*/) const override;
};
