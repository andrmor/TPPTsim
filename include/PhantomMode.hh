#ifndef PhantomMode_h
#define PhantomMode_h

#include <vector>

class G4LogicalVolume;

class PhantomModeBase
{
public:
    virtual ~PhantomModeBase(){}

    virtual G4LogicalVolume * definePhantom(G4LogicalVolume * /*logicWorld*/) {return nullptr;}
};

// ---

class PhantomNone : public PhantomModeBase
{
    //nothing to add :)
};

// ---

class PhantomPMMA : public PhantomModeBase
{
public:
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
};

// ---

class PhantomTinyCube : public PhantomModeBase
{
public:
    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
};

// ---

class PhantomDerenzo : public PhantomModeBase
{
public:
    PhantomDerenzo(double diameter, double height, const std::vector<double> & holeDiameters, double radialOffset, double margin, double dPhi) :
        Diameter(diameter), Height(height), HoleDiameters(holeDiameters), RadialOffset(radialOffset), Margin(margin), DPhi(dPhi) {}

    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;

protected:
    double Diameter     = 200.0;
    double Height       = 200.0;
    std::vector<double> HoleDiameters;
    double RadialOffset = 50.0;
    double Margin       = 50.0;
    double DPhi         = 0;

};

#endif // PhantomMode_h
