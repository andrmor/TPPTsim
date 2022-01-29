#ifndef DicomPhantom_h
#define DicomPhantom_h

// includes code from P. Arce

#include "PhantomMode.hh"
//#include "G4IntersectionSolid.hh"
#include "globals.hh"
#include "json11.hh"

#include <vector>
#include <map>

class G4Material;
class G4LogicalVolume;
class DicomPhantomZSliceHeader;
class G4VPhysicalVolume;
class G4VisAttributes;

class PhantomDICOM : public PhantomModeBase
{
public:
    PhantomDICOM(std::string dataDir, std::string sliceBaseFileName, int sliceFrom, int sliceTo,
                 int lateralCompression, double containerRadius, const std::vector<double> & posInContainer);

    G4LogicalVolume * definePhantom(G4LogicalVolume * logicWorld) override;
    std::string       getTypeName() const override {return "PhantomDICOM";}
    void              readFromJson(const json11::Json & json) override;

protected:
    void              doWriteToJson(json11::Json::object & json) const override;

    std::string              DataDir;
    std::string              SliceFileBase;
    int                      SliceFrom;
    int                      SliceTo;
    int                      LateralCompression;
    double                   PhantRadius;
    std::vector<double>      PosContainer;

    const std::string        ConvertionFileName = "CT2Density.dat";
    const double             densityDiff = -1.0;      // former was read from env variable "DICOM_CHANGE_MATERIAL_DENSITY" --> -1.0 inicates this mechanism is disabled

    bool                     ContainerInvisible = true;
    std::vector<std::pair<std::string, float>> MatUpperDens;
    std::map<G4String, G4VisAttributes*> ColourMap;                 // in use during tracking!
    std::vector<std::string> SliceFiles;
    double                   zStart;
    std::vector<std::pair<double,double>> BoxXY;                    // in use during tracking!
    int                      BoxesPerSlice;

    G4Material             * AirMat = nullptr;

    G4LogicalVolume        * fContainer_logic;
    G4VPhysicalVolume      * fContainer_phys;

    int                      fNoFiles; // number of DICOM files

    std::vector<G4Material*> fOriginalMaterials;  // list of original materials
    std::vector<G4Material*> fMaterials;    // list of new materials created to distinguish different density voxels that have the same original materials

    std::vector<size_t>      MaterialIDs;        // index of material of each voxel

    std::vector<G4Material*> ReducedMaterials;                     // in use during tracking!

    std::map<int, double>    fDensityDiffs; // Density difference to distinguish material for each original material (by index)

    std::vector<DicomPhantomZSliceHeader*> fZSliceHeaders; // list of z slice headers (one per DICOM files)
    DicomPhantomZSliceHeader * fZSliceHeaderMerged = nullptr;        // z slice header resulted from merging all z slice headers

    int    fNVoxelX, fNVoxelY, fNVoxelZ;
    double fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ;
    double fMinX, fMinY, fMinZ; // minimum extension of voxels (position wall)
    double fMaxX, fMaxY, fMaxZ; // maximum extension of voxels (position wall)

    //std::map<int, G4Material*> thePhantomMaterialsOriginal; // map numberOfMaterial to G4Material. They are the list of materials as built from .geom file
    //DicomPhantomZSliceMerged * fMergedSlices;
    //std::set<G4LogicalVolume*> fScorers;
    //G4IntersectionSolid * test;// =new G4IntersectionSolid("bx2CyleafTr", box2leaf, cyLeafTr);
    //bool fConstructed;

private:
    void buildMaterials();
    void readPhantomData();
    void readPhantomDataFile(const G4String & fname); // read one of the DICOM files describing the phantom (usually one per Z slice) and builds a corresponding DicomPhantomZSliceHeader
    void mergeZSliceHeaders();
    void computePhantomVoxelization();
    void constructPhantomContainer(G4LogicalVolume * logicWorld);
    void constructPhantom();
    void readMaterialFile(const std::string & fileName);
    void readColorMap(const std::string & fileName);
    void generateSliceFileNames();
    void clearTmpAndOptimizeContainers();

    G4Material * buildMaterialWithChangingDensity(const G4Material * origMate, G4float density, G4String newMateName); // build a new material if the density of the voxel is different to the other voxels

    // NOT IN USE
    //void ReadVoxelDensities(std::ifstream & fin); // read the DICOM files describing the phantom
};

/*
struct matInfo
{
    double fSumdens;
    int    fNvoxels;
    int    fId;
};
*/

#endif
