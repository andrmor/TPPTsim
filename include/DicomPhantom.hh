//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file medical/DICOM/include/DicomDetectorConstruction.hh
/// \brief Definition of the DicomDetectorConstruction class
//

#ifndef DicomPhantom_h
#define DicomPhantom_h

#include "PhantomMode.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DicomPhantomZSliceHeader.hh"
#include "G4ThreeVector.hh"
#include "G4IntersectionSolid.hh"
#include "json11.hh"

#include <vector>
#include <set>
#include <map>

class G4Material;
class G4Box;
class G4Tubs;
class G4LogicalVolume;
class DicomPhantomZSliceMerged;

//*******************************************************
/// Dicom detector construction
///
///      - Start the building of the geometry
///      - Initialisation of materials
///      - Creation of the world
///      - Reading of the DICOM data
///
/// History: 30.11.07  First version
/// \author  P. Arce
//*******************************************************

struct matInfo
{
    G4double fSumdens;
    G4int fNvoxels;
    G4int fId;
};

class PhantomModeDICOM : public PhantomModeBase
{
public:
    PhantomModeDICOM(G4double phantRadius, const std::vector<double> & posContainer, G4String dataFile, bool delFiles);

    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
    std::string       getTypeName() const override {return "PhantomModeDICOM";}
    void              readFromJson(const json11::Json & json) override;

protected:
    void              doWriteToJson(json11::Json::object & json) const override;


    void ConstructPhantom();
    void InitialisationOfMaterials();

    void ReadPhantomData();

    // read the DICOM files describing the phantom
    void ReadVoxelDensities(std::ifstream & fin);

    void ReadPhantomDataFile(const G4String & fname);
    // read one of the DICOM files describing the phantom
    // (usually one per Z slice).
    //  Build a DicomPhantomZSliceHeader for each file

    void MergeZSliceHeaders();
    // merge the slice headers of all the files

    G4Material * BuildMaterialWithChangingDensity(const G4Material * origMate, G4float density, G4String newMateName);
    // build a new material if the density of the voxel is different
    // to the other voxels

    void ConstructPhantomContainer(G4LogicalVolume * logicWorld);

protected:
    G4double            PhantRadius;
    std::vector<double> PosContainer;
    G4String            DicomPath;
    G4String            DataFile;
    bool                DelFiles;
    double              zStart;

    std::vector<std::pair<double,double>> BoxXY;

    int                 BoxesPerSlice;
    size_t            * reducedIDs;

    //std::vector<unsigned long> myMateIDs;
    //size_t            * testIDs;

    G4Material        * fAir = nullptr;

    // World ...
    //G4Tubs            * fWorld_solid;
    //G4LogicalVolume   * fWorld_logic;
    //G4VPhysicalVolume * fWorld_phys;

    G4Tubs            * fContainer_solid;   // AM: not needed here
    G4LogicalVolume   * fContainer_logic;   // AM: not needed here
    G4VPhysicalVolume * fContainer_phys;

    G4int               fNoFiles; // number of DICOM files

    std::vector<G4Material*> fOriginalMaterials;  // list of original materials
    std::vector<G4Material*> fMaterials;
    // list of new materials created to distinguish different density
    //  voxels that have the same original materials

    size_t * fMateIDs; // index of material of each voxel
    //unsigned int* fMateIDs; // index of material of each voxel

    std::map<G4int,G4double> fDensityDiffs;
    // Density difference to distinguish material for each original
    // material (by index)

    std::vector<DicomPhantomZSliceHeader*> fZSliceHeaders;
    // list of z slice header (one per DICOM files)
    DicomPhantomZSliceHeader * fZSliceHeaderMerged;
    // z slice header resulted from merging all z slice headers

    G4int fNVoxelX, fNVoxelY, fNVoxelZ;
    G4double fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ;
    G4double fMinX,fMinY,fMinZ; // minimum extension of voxels (position wall)
    G4double fMaxX,fMaxY,fMaxZ; // maximum extension of voxels (position wall)

    std::map<G4int,G4Material*> thePhantomMaterialsOriginal;
    // map numberOfMaterial to G4Material. They are the list of materials as
    // built from .geom file

    DicomPhantomZSliceMerged * fMergedSlices;

    std::set<G4LogicalVolume*> fScorers;
    G4IntersectionSolid * test;// =new G4IntersectionSolid("bx2CyleafTr", box2leaf, cyLeafTr);
    G4bool fConstructed;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
inline G4int PhantomModeDICOM::GetTotalVoxels() const
{
    return fNVoxelX * fNVoxelY * fNVoxelZ;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
