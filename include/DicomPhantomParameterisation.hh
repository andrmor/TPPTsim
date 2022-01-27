#ifndef DicomPhantomParameterisation_HH
#define DicomPhantomParameterisation_HH

#include <map>

#include "G4PhantomParameterisation.hh"
#include "G4ThreeVector.hh"

class G4VisAttributes;
class DicomPhantomZSliceHeader;

class DicomPhantomParameterisation : public G4VPVParameterisation //G4PhantomParameterisation
{
public:
    typedef std::map<G4String,G4VisAttributes*> ColourMap_t;

public:
    DicomPhantomParameterisation(std::vector<std::pair<double,double>> coord2D, double zStart,
                                 const std::vector<G4Material *> & materials, G4String colourFile);

    G4Material * ComputeMaterial(const G4int repNo, G4VPhysicalVolume *currentVol, const G4VTouchable *parentTouch = nullptr) override;
    void         ComputeTransformation (const int copyNo, G4VPhysicalVolume * physVol) const override;

    //const ColourMap_t& GetColourMap() const { return fColours; }
    //ColourMap_t& GetColourMap() { return fColours; }

    // ANDR
    void SetVoxelDimensions(double DX, double DY, double DZ) {HalfVoxelX = DX; HalfVoxelY = DY; HalfVoxelZ = DZ;}
    // ----

protected:
    std::vector<std::pair<double,double>> XY;
    double ZStart;
    const std::vector<G4Material*> & Materials;

private:
    void ReadColourData(G4String colourFile);

private:
    ColourMap_t fColours;
    std::map<G4int, G4VisAttributes*> mColours;

    DicomPhantomZSliceHeader* ZSliceHeader;
    //G4int control;

    //ANDR
    double HalfVoxelX;
    double HalfVoxelY;
    double HalfVoxelZ;
    //----

};

#endif
