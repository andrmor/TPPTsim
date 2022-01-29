#ifndef DicomPhantomParameterisation_HH
#define DicomPhantomParameterisation_HH

#include <vector>
#include <map>

#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"

class G4VisAttributes;

class DicomPhantomParameterisation : public G4VPVParameterisation
{
public:
    DicomPhantomParameterisation(const std::vector<std::pair<double,double>> & coord2D, double zStart,
                                 const std::vector<G4Material*> & materials,
                                 const std::map<G4String,G4VisAttributes*> & colourMap);

    G4Material * ComputeMaterial(const G4int repNo, G4VPhysicalVolume *currentVol, const G4VTouchable *parentTouch = nullptr) override;
    void         ComputeTransformation (const int copyNo, G4VPhysicalVolume * physVol) const override;

    void         setVoxelHalfSizeZ(double dz) {HalfVoxelZ = dz;}

protected:
    const std::vector<std::pair<double,double>> & XY;
    const double ZStart;
    const std::vector<G4Material*> & Materials;
    const std::map<G4String, G4VisAttributes*> & ColourMap;

    double HalfVoxelZ;
    int BoxesPerSlice = 1;
    G4VisAttributes * defaultVisAttributes = nullptr;
};

#endif
