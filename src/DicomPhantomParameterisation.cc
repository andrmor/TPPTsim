#include "DicomPhantomParameterisation.hh"
#include "out.hh"

//#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VVisManager.hh"

DicomPhantomParameterisation::DicomPhantomParameterisation(std::vector<Voxel> & voxels,
                                                           const std::map<G4String,G4VisAttributes*> & colourMap) :
    G4VPVParameterisation(), Voxels(voxels), ColourMap(colourMap)
{
    defaultVisAttributes = new G4VisAttributes(G4Colour(1.0, 0, 0));
}

DicomPhantomParameterisation::~DicomPhantomParameterisation()
{
    delete defaultVisAttributes;
}

G4Material * DicomPhantomParameterisation::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume * physVol, const G4VTouchable *)
{
    G4Material * mat = Voxels[copyNo].Mat;

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
    const Voxel & vox = Voxels[copyNo];

    physVol->SetTranslation( {vox.X, vox.Y, vox.Z} );
    physVol->SetRotation(nullptr);
}
