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

    std::string         DataDir;
    std::string         SliceFileBase;
    int                 SliceFrom;
    int                 SliceTo;
    int                 LateralCompression;
    double              PhantRadius;
    std::vector<double> PosInWorld;

    const std::string   ConvertionFileName = "CT2Density.dat";
    const bool          ContainerInvisible = true;

    std::vector<std::pair<std::string, float>> MatUpperDens;
    std::map<G4String, G4VisAttributes*>       ColourMap;     // --->PERSISTENT: in use during gui session!

    std::vector<std::string> SliceFiles;
    double                   zStart;

    std::vector<Voxel>       Voxels;                          // --->PERSISTENT: in use during tracking!

    std::vector<G4Material*> Materials;
    G4Material             * AirMat = nullptr;                // voxels of this material are not created
    std::vector<size_t>      MaterialIDs;                     // index of the material for each voxel in the dicom

    std::vector<DicomPhantomZSliceHeader*> fZSliceHeaders;    // list of z slice headers (one per DICOM files)
    DicomPhantomZSliceHeader * fZSliceHeaderMerged = nullptr; // z slice header resulted from merging all z slice headers

    int    fNVoxelX, fNVoxelY, fNVoxelZ;
    double fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ;

private:
    void buildMaterials();
    void readPhantomData();
    void readPhantomDataFile(const std::string & fname); // read one of the DICOM files describing the phantom (usually one per Z slice) and builds a corresponding DicomPhantomZSliceHeader
    void mergeZSliceHeaders();
    void prepareParameterizationParameters();
    G4LogicalVolume * makeContainer(G4LogicalVolume * logicWorld);
    void constructPhantom(G4LogicalVolume * PhantomLogical);
    void optimizeMemory();
    void readMaterialFile(const std::string & fileName);
    void readColorMap(const std::string & fileName);
    void generateSliceFileNames();

};

#endif
