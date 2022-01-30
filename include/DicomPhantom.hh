#ifndef DicomPhantom_h
#define DicomPhantom_h

// includes code from P. Arce

#include "PhantomMode.hh"
#include "Voxel.hh"
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
                 int lateralCompression, double containerRadius, const std::vector<double> & posInWorld);

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
    std::vector<double>      PosInWorld;

    const std::string        ConvertionFileName = "CT2Density.dat";

    bool                     ContainerInvisible = true;
    std::vector<std::pair<std::string, float>> MatUpperDens;
    std::map<G4String, G4VisAttributes*> ColourMap;   // --->PERSISTENT: in use during gui session!
    std::vector<std::string> SliceFiles;
    double                   zStart;

    std::vector<Voxel>       Voxels;      // --->PERSISTENT: in use during tracking!

    G4Material             * AirMat = nullptr;

    G4LogicalVolume        * fContainer_logic;
    G4VPhysicalVolume      * fContainer_phys;

    int                      fNoFiles;            // number of DICOM files

    std::vector<G4Material*> fOriginalMaterials;  // list of original materials, small
    std::vector<size_t>      MaterialIDs;         // index of material of each voxel


    std::vector<DicomPhantomZSliceHeader*> fZSliceHeaders;      // list of z slice headers (one per DICOM files)
    DicomPhantomZSliceHeader * fZSliceHeaderMerged = nullptr;   // z slice header resulted from merging all z slice headers

    int    fNVoxelX, fNVoxelY, fNVoxelZ;
    double fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ;
    double fMinX, fMinY, fMinZ; // minimum extension of voxels (position wall)
    double fMaxX, fMaxY, fMaxZ; // maximum extension of voxels (position wall)

private:
    void buildMaterials();
    void readPhantomData();
    void readPhantomDataFile(const G4String & fname); // read one of the DICOM files describing the phantom (usually one per Z slice) and builds a corresponding DicomPhantomZSliceHeader
    void mergeZSliceHeaders();
    void computePhantomVoxelization();
    void constructPhantomContainer(G4LogicalVolume * logicWorld);
    void constructPhantom();
    void optimizeMemory();
    void readMaterialFile(const std::string & fileName);
    void readColorMap(const std::string & fileName);
    void generateSliceFileNames();

};

#endif
