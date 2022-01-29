#include "DicomPhantomParameterisation.hh"
#include "out.hh"

#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VVisManager.hh"

DicomPhantomParameterisation::DicomPhantomParameterisation(const std::vector<std::pair<double, double> > & coord2D, double zStart,
                                                           const std::vector<G4Material*> & materials,
                                                           const std::map<G4String,G4VisAttributes*> & colourMap) :
    G4VPVParameterisation(), XY(coord2D), ZStart(zStart),
    Materials(materials), ColourMap(colourMap)
{
    BoxesPerSlice = XY.size();
    defaultVisAttributes = new G4VisAttributes(G4Colour(1.0, 0, 0));
}

G4Material * DicomPhantomParameterisation::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume * physVol, const G4VTouchable *)
{
    G4Material * mat = Materials[copyNo];

    if (G4VVisManager::GetConcreteInstance() && physVol)
    {
        const G4String & matName = mat->GetName();
        auto it = ColourMap.find(matName);
        if (it != ColourMap.end())
            physVol->GetLogicalVolume()->SetVisAttributes(it->second);
        else
            physVol->GetLogicalVolume()->SetVisAttributes(defaultVisAttributes);
    }

    return mat;
}

void DicomPhantomParameterisation::ComputeTransformation(const int copyNo, G4VPhysicalVolume * physVol) const
{
    double Z = ZStart + std::ceil(copyNo/BoxesPerSlice) * 2.0 * HalfVoxelZ;
    int numInSlice = copyNo % BoxesPerSlice;
    const std::pair<double, double> & XYpair = XY[numInSlice];

    G4ThreeVector origin(XYpair.first, XYpair.second, Z);
    physVol->SetTranslation(origin);
    physVol->SetRotation(nullptr);
}
