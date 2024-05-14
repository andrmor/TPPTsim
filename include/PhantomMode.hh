#ifndef PhantomMode_h
#define PhantomMode_h

#include "MaterialBuilder.hh"
#include "json11.hh"

#include <string>
#include <vector>

class PhantomModeBase;
class G4LogicalVolume;

class PhantomModeFactory // when adding a new mode, do not forget to modify the implementation of the factory!
{
public:
    static PhantomModeBase * makePhantomModeInstance(const json11::Json & json);
};

class PhantomModeBase
{
public:
    virtual ~PhantomModeBase(){}

    virtual std::string getTypeName() const = 0;
    virtual G4LogicalVolume * definePhantom(G4LogicalVolume * /*logicWorld*/) = 0; //returns the phantom's logical volume

    void writeToJson(json11::Json::object & json) const;
    virtual void readFromJson(const json11::Json & /*json*/) {}

protected:
    virtual void doWriteToJson(json11::Json::object & /*json*/) const {}
};

// ---

class PhantomNone : public PhantomModeBase
{
    std::string getTypeName() const override {return "PhantomNone";}
    G4LogicalVolume * definePhantom(G4LogicalVolume *) override {return nullptr;}
};

// ---

class PhantomCylinder : public PhantomModeBase
{
public:
    PhantomCylinder(double diameter, double length, std::string g4_material_name);
    PhantomCylinder(double diameter, double length, EMaterial material);
    PhantomCylinder();

    std::string getTypeName() const override {return "PhantomCylinder";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;

    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    double Diameter = 200.0;
    double Length   = 200.0;

    MaterialBuilder * MatBuilder = nullptr;
};

class PhantomBox : public PhantomModeBase
{
public:
    PhantomBox(double sizeX, double sizeY, double sizeZ, std::string g4_material_name);
    PhantomBox(double sizeX, double sizeY, double sizeZ, EMaterial material);
    PhantomBox();

    std::string getTypeName() const override {return "PhantomBox";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;

    void readFromJson(const json11::Json & json) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    double SizeX = 100.0;
    double SizeY = 100.0;
    double SizeZ = 100.0;

    MaterialBuilder * MatBuilder = nullptr;
};

// ---

class PhantomDerenzo : public PhantomModeBase
{
public:
    PhantomDerenzo(double diameter, double height, const std::vector<double> & holeDiameters, double radialOffset, double margin, double dPhi) :
        Diameter(diameter), Height(height), HoleDiameters(holeDiameters), RadialOffset(radialOffset), Margin(margin), DPhi(dPhi) {}

    std::string getTypeName() const override {return "PhantomDerenzo";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;

    void readFromJson(const json11::Json & json) override;

protected:
    double Diameter     = 200.0;
    double Height       = 200.0;
    std::vector<double> HoleDiameters;
    double RadialOffset = 50.0;
    double Margin       = 50.0;
    double DPhi         = 0;     // rotation angle around vertical axis

    void doWriteToJson(json11::Json::object & json) const override;
};

// ---

class PhantomParam : public PhantomModeBase
{
public:
    std::string getTypeName() const override {return "PhantomParam";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;

    void readFromJson(const json11::Json &) override {}

protected:
    void doWriteToJson(json11::Json::object &) const override {}
};

// ---

class PhantomEnergyCalibration : public PhantomModeBase
{
public:
    std::string getTypeName() const override {return "PhantomEnergyCalibration";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
};

// ---

class PhantomEspana : public PhantomModeBase
{
public:
    std::string getTypeName() const override {return "PhantomEspana";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
};

class PhantomBauerGel : public PhantomModeBase
{
public:
    std::string getTypeName() const override {return "PhantomBauerGel";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
};

class PhantomBauerCa : public PhantomModeBase
{
public:
    std::string getTypeName() const override {return "PhantomBauerCa";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
};

class PhantomRT : public PhantomModeBase
{
public:
    std::string getTypeName() const override {return "PhantomRT";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
};

class PhantomMarekWater : public PhantomModeBase
{
public:
    std::string getTypeName() const override {return "PhantomMarekWater";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
};

class PhantomMarekCompartments : public PhantomModeBase
{
public:
    std::string getTypeName() const override {return "PhantomMarekCompartments";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
};

class PhantomBeamDerenzoAndr : public PhantomModeBase
{
public:
    std::string getTypeName() const override {return "PhantomBeamDerenzoAndr";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
};

class PhantomBeamDerenzoAndr2 : public PhantomModeBase
{
public:
    std::string getTypeName() const override {return "PhantomBeamDerenzoAndr2";}
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;

private:
    void placeArray(G4LogicalVolume * element, double step, double x, double y, G4LogicalVolume * container, const std::string & name);
};


#endif // PhantomMode_h
