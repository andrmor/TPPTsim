#include "DicomPhantomParameterisation.hh"
#include "DicomHandler.hh"
#include "DicomPhantomZSliceHeader.hh"

#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

DicomPhantomParameterisation::DicomPhantomParameterisation(std::vector<std::pair<double,double>> coord2D, double zStart,
                                                           const std::vector<G4Material*> & materials, G4String colourFile)
//    :G4PhantomParameterisation(), XY(coord2D), ZStart(zStart)
    : G4VPVParameterisation(),
      XY(coord2D), ZStart(zStart), Materials(materials)
{
    ReadColourData(colourFile);
//    SetSkipEqualMaterials(false);
}

void DicomPhantomParameterisation::ReadColourData(G4String colourFile)
{
    G4VisAttributes * blankAtt = new G4VisAttributes;
    blankAtt->SetVisibility(false);
    fColours["Default"] = blankAtt;

    std::ifstream fin(colourFile.c_str());
    int nMate;
    G4String mateName;
    G4double cred, cgreen, cblue, copacity;
    fin >> nMate;
    for (int ii = 0; ii < nMate; ii++)
    {
        fin >> mateName;
        if (fin.eof()) break;
        fin >> cred >> cgreen >> cblue >> copacity;
        G4Colour colour( cred, cgreen, cblue, copacity );
        G4VisAttributes * visAtt = new G4VisAttributes(colour);
        visAtt->SetVisibility(true);
        fColours[mateName] = visAtt;
        mColours[ii] = new G4VisAttributes(*visAtt);
    }
}

#include "out.hh"
G4Material * DicomPhantomParameterisation::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume * physVol, const G4VTouchable *)
{
    //return (copyNo % 3 == 0 ? air : air1);

    //G4Material* mate = G4PhantomParameterisation::ComputeMaterial(copyNo, physVol, 0);
    G4Material * mat = Materials[copyNo];

    if (G4VVisManager::GetConcreteInstance() && physVol)
    {
        G4String matName = mat->GetName();
//        std::string::size_type iuu = matName.find("__");
//        if (iuu != std::string::npos) matName = matName.substr(0, iuu);

        auto it = fColours.find(matName);
        if (it != fColours.end())
            physVol->GetLogicalVolume()->SetVisAttributes(it->second);
        else
        {
            bool found = false;
            for (const auto & itr : fColours)
            {
                G4String mat_color = itr.first;
                auto len = mat_color.length();
                if (matName.find(mat_color) == 0 && matName.length() > len && matName[len] == '_')
                {
                    physVol->GetLogicalVolume()->SetVisAttributes( fColours.find(mat_color)->second );
                    found = true;
                }
                if (found) break;
            }
//            if (!found)
//            {
//                G4int matIndex = G4int(GetMaterialIndex(copyNo));
//                static uintmax_t n = 0;
//                if (n++ < 100) G4cout << "Unknown material name " << matName << " for index " << matIndex << G4endl;
//                if (mColours.find(matIndex) != mColours.end())
//                    physVol->GetLogicalVolume()->SetVisAttributes(mColours.find(matIndex)->second);
//                else
//                    physVol->GetLogicalVolume()->SetVisAttributes(fColours.begin()->second);
//            }
        }
    }

    return mat;
}

void DicomPhantomParameterisation::ComputeTransformation(const int copyNo, G4VPhysicalVolume * physVol) const
{
//    G4double HalfVoxelZ = GetVoxelHalfZ();
    int BoxesPerSlice = XY.size();

    double Z = ZStart + std::ceil(copyNo/BoxesPerSlice) * 2.0*HalfVoxelZ;
    int numInSlice = copyNo % BoxesPerSlice;
    const std::pair<double, double> & XYpair = XY[numInSlice];

    G4ThreeVector origin(XYpair.first, XYpair.second, Z);
    physVol->SetTranslation(origin);
    physVol->SetRotation(nullptr);
}
