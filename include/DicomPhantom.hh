#ifndef DicomPhantom_h
#define DicomPhantom_h

// includes code from P. Arce

#include "PhantomMode.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DicomPhantomZSliceHeader.hh"
//#include "G4IntersectionSolid.hh"
#include "json11.hh"

#include <vector>
#include <map>
//#include <set>

class G4Material;
class G4LogicalVolume;
//class DicomPhantomZSliceMerged;

struct matInfo
{
    double fSumdens;
    int    fNvoxels;
    int    fId;
};

class PhantomDICOM : public PhantomModeBase
{
public:
    PhantomDICOM(double phantRadius, const std::vector<double> & posContainer, std::string dataFile, bool delFiles);

    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
    std::string       getTypeName() const override {return "PhantomDICOM";}
    void              readFromJson(const json11::Json & json) override;

protected:
    void              doWriteToJson(json11::Json::object & json) const override;

    // read the DICOM files describing the phantom
    //void ReadVoxelDensities(std::ifstream & fin);   // NOT IN USE

    void ReadPhantomDataFile(const G4String & fname);
    // read one of the DICOM files describing the phantom
    // (usually one per Z slice).
    //  Build a DicomPhantomZSliceHeader for each file

    void MergeZSliceHeaders();
    // merge the slice headers of all the files

    G4Material * BuildMaterialWithChangingDensity(const G4Material * origMate, G4float density, G4String newMateName);
    // build a new material if the density of the voxel is different
    // to the other voxels


protected:
    double              PhantRadius;
    std::vector<double> PosContainer;
    std::string         DicomPath;
    std::string         DataFile;
    bool                DelFiles;
    double              zStart;

    std::vector<std::pair<double,double>> BoxXY;

    int                 BoxesPerSlice;
    size_t            * reducedIDs;

    G4Material        * fAir = nullptr;

    G4LogicalVolume   * fContainer_logic;
    G4VPhysicalVolume * fContainer_phys;

    int                 fNoFiles; // number of DICOM files

    std::vector<G4Material*> fOriginalMaterials;  // list of original materials
    std::vector<G4Material*> fMaterials;
    // list of new materials created to distinguish different density
    //  voxels that have the same original materials

    size_t * fMateIDs; // index of material of each voxel
    //unsigned int* fMateIDs; // index of material of each voxel

    std::map<int, double> fDensityDiffs; // Density difference to distinguish material for each original material (by index)

    std::vector<DicomPhantomZSliceHeader*> fZSliceHeaders; // list of z slice header (one per DICOM files)
    DicomPhantomZSliceHeader * fZSliceHeaderMerged;        // z slice header resulted from merging all z slice headers

    int    fNVoxelX, fNVoxelY, fNVoxelZ;
    double fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ;
    double fMinX,fMinY,fMinZ; // minimum extension of voxels (position wall)
    double fMaxX,fMaxY,fMaxZ; // maximum extension of voxels (position wall)

    std::map<int, G4Material*> thePhantomMaterialsOriginal; // map numberOfMaterial to G4Material. They are the list of materials as built from .geom file

    //DicomPhantomZSliceMerged * fMergedSlices;
    //std::set<G4LogicalVolume*> fScorers;
    //G4IntersectionSolid * test;// =new G4IntersectionSolid("bx2CyleafTr", box2leaf, cyLeafTr);
    //bool fConstructed;

    double densityDiff = -1.0; // former from "DICOM_CHANGE_MATERIAL_DENSITY" --> -1.0 inicates this mechanism is disabled

private:
    void buildMaterials();
    void readPhantomData();
    void computePhantomVoxelization(); // TODO : refactor
    void constructPhantomContainer(G4LogicalVolume * logicWorld);
    void constructPhantom();
};

#endif
