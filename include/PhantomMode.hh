#ifndef PhantomMode_h
#define PhantomMode_h

#include <vector>

class G4LogicalVolume;

class PhantomModeBase
{
public:
    virtual ~PhantomModeBase(){}

    virtual void definePhantom(G4LogicalVolume * logicWorld) {}
};

// ---

class PhantomModeNone : public PhantomModeBase
{
    //nothing to add :)
};

// ---

class PhantomModePMMA : public PhantomModeBase
{
public:
    void definePhantom(G4LogicalVolume * logicWorld) override;
};

// ---

class PhantomModeDerenzo : public PhantomModeBase
{
public:
    PhantomModeDerenzo(double height, double diameter, const std::vector<double> & holeDiameters, double radialOffset, double margin, double dPhi) :
        Height(height), Diameter(diameter), HoleDiameters(holeDiameters), RadialOffset(radialOffset), Margin(margin), DPhi(dPhi) {}

    void definePhantom(G4LogicalVolume * logicWorld) override;

protected:
    double Height       = 200.0;
    double Diameter     = 200.0;
    std::vector<double> HoleDiameters;
    double RadialOffset = 50.0;
    double Margin       = 50.0;
    double DPhi         = 0;

};

#endif // PhantomMode_h
