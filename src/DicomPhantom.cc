#include "DicomPhantom.hh"
#include "DicomPhantomParameterisation.hh"
#include "DicomPhantomZSliceHeader.hh"
#include "DicomHandler.hh"
#include "DicomPhantomZSliceHeader.hh"
#include "jstools.hh"
#include "out.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4UIcommand.hh"
#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVParameterised.hh"
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

PhantomDICOM::PhantomDICOM(std::string dataDir, std::string sliceBaseFileName, int sliceFrom, int sliceTo,
                           int lateralCompression, double containerRadius, const std::vector<double> & posInWorld) :
                           DataDir(dataDir), SliceFileBase(sliceBaseFileName), SliceFrom(sliceFrom), SliceTo(sliceTo),
                           LateralCompression(lateralCompression),
                           PhantRadius(containerRadius), PosInWorld(posInWorld)
{
    ContainerInvisible = true;
}

G4LogicalVolume * PhantomDICOM::definePhantom(G4LogicalVolume * logicWorld)
{
    // clear old g4dcm files in the phantom directory - if compression changes, they are invalid!
    // we can consider saving compression factor in a file in the same directory not to force-rebuild existent files
    std::string cmd = "rm -f " + DataDir + "/*.g4dcm*";
    system(cmd.data());

    readMaterialFile(DataDir + '/' + "Materials.dat");
    readColorMap(DataDir + '/' + "ColorMap.dat");

    generateSliceFileNames();

    DicomHandler dcmHandler;
    dcmHandler.processFiles(DataDir, ConvertionFileName, LateralCompression, MatUpperDens, SliceFiles);

    buildMaterials();
    readPhantomData();
    mergeZSliceHeaders();
    computePhantomVoxelization();
    constructPhantomContainer(logicWorld);
    constructPhantom();

    optimizeMemory();

    return fContainer_logic;
}

void PhantomDICOM::optimizeMemory()
{
    MaterialIDs.clear(); MaterialIDs.shrink_to_fit();
    out("MaterialIDs", MaterialIDs.capacity());

    for (auto & sh : fZSliceHeaders) delete sh;
    fZSliceHeaders.clear();

    delete fZSliceHeaderMerged; fZSliceHeaderMerged = nullptr;
}

void PhantomDICOM::doWriteToJson(json11::Json::object & json) const
{
    json["PhantRadius"] = PhantRadius;
    // TODO: Hugo!
}

void PhantomDICOM::readFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "PhantRadius", PhantRadius);
    // TODO: Hugo!
}

void PhantomDICOM::buildMaterials()
{
    double z, a;
    std::string name, symbol;

    // Hydrogen
    z = 1; a = 1.008*g/mole;
    G4Element* elH = new G4Element("Hydrogen", "H", z, a);

    // Carbon
    z = 6, a = 12.011*g/mole;
    G4Element* elC = new G4Element("Carbon", "C", z, a);

    // Nitrogen
    z = 7; a = 14.007*g/mole;
    G4Element* elN = new G4Element("Nitrogen", "N", z, a);

    // Oxygen
    z = 8; a = 16.00*g/mole;
    G4Element* elO = new G4Element("Oxygen", "O", z, a);

    // Sodium
    z = 11; a = 22.98977*g/mole;
    G4Element* elNa = new G4Element("Sodium", "Na", z, a);

    // Magnesium
    z = 12; a = 24.305*g/mole;
    G4Element* elMg = new G4Element("Magnesium", "Mg", z, a);

    // Aluminium
    z = 13; a = 26.98154*g/mole;
    G4Element* elAl = new G4Element("Aluminium", "Al", z, a);

    // Phosphorus
    z = 15; a = 30.97376*g/mole;
    G4Element* elP = new G4Element("Phosphorus", "P", z, a);

    // Sulfur
    z = 16; a = 32.06*g/mole;
    G4Element* elS = new G4Element("Sulfur", "S", z, a);

    // Chlorine
    z = 17; a = 35.45*g/mole;
    G4Element* elCl = new G4Element("Chlorine", "Cl", z, a);

    // Argon
    z = 18; a = 39.95*g/mole;
    G4Element* elAr = new G4Element("Argon", "Ar", z, a);

    // Potassium
    z = 19; a = 39.09833*g/mole;
    G4Element* elK = new G4Element("Potassium", "K", z, a);

    // Calcium
    z = 20; a = 40.07844*g/mole;
    G4Element* elCa = new G4Element("Calcium", "Ca", z, a);

    // Iron
    z = 26; a = 55.84522*g/mole;
    G4Element* elFe = new G4Element("Iron", "Fe", z, a);

    // Zinc
    z = 30; a = 65.38222*g/mole;
    G4Element* elZn = new G4Element("Zinc", "Zn", z, a);

    int numberofElements;
    double density;

    // Air
    density = 0.00121*g/cm3;  numberofElements = 3;
    AirMat = new G4Material("Air", density, numberofElements);
    AirMat->AddElement(elN, 75.5*perCent);
    AirMat->AddElement(elO, 23.3*perCent);
    AirMat->AddElement(elAr, 1.3*perCent);

    // Lung
    density = 0.500*g/cm3; numberofElements = 9;
    G4Material* lung = new G4Material("Lung", density, numberofElements);
    lung->AddElement(elH, 10.3*perCent);
    lung->AddElement(elC, 10.5*perCent);
    lung->AddElement(elN, 3.1*perCent);
    lung->AddElement(elO, 74.9*perCent);
    lung->AddElement(elNa,0.2*perCent);
    lung->AddElement(elP, 0.2*perCent);
    lung->AddElement(elS, 0.3*perCent);
    lung->AddElement(elCl,0.3*perCent);
    lung->AddElement(elK, 0.2*perCent);

    // Adipose tissue
    density = 0.95*g/cm3; numberofElements = 7;
    G4Material* adiposeTissue = new G4Material("AdiposeTissue", density, numberofElements);
    adiposeTissue->AddElement(elH, 11.4*perCent);
    adiposeTissue->AddElement(elC, 59.8*perCent);
    adiposeTissue->AddElement(elN, 0.7*perCent);
    adiposeTissue->AddElement(elO, 27.8*perCent);
    adiposeTissue->AddElement(elNa,0.1*perCent);
    adiposeTissue->AddElement(elS, 0.1*perCent);
    adiposeTissue->AddElement(elCl,0.1*perCent);

    // Muscle
    density = 1.05*g/cm3; numberofElements = 9;
    G4Material* muscle = new G4Material("Muscle", density, numberofElements);
    muscle->AddElement(elH, 10.2*perCent);
    muscle->AddElement(elC, 14.3*perCent);
    muscle->AddElement(elN, 3.4*perCent);
    muscle->AddElement(elO, 71.0*perCent);
    muscle->AddElement(elNa,0.1*perCent);
    muscle->AddElement(elP, 0.2*perCent);
    muscle->AddElement(elS, 0.3*perCent);
    muscle->AddElement(elCl,0.1*perCent);
    muscle->AddElement(elK, 0.4*perCent);
    
    // Cartilage
    density = 1.1*g/cm3; numberofElements = 8;
    G4Material* cartilage = new G4Material("Cartilage", density, numberofElements);
    cartilage->AddElement(elH, 9.6*perCent);
    cartilage->AddElement(elC, 9.9*perCent);
    cartilage->AddElement(elN, 2.2*perCent);
    cartilage->AddElement(elO, 74.4*perCent);
    cartilage->AddElement(elNa,0.5*perCent);
    cartilage->AddElement(elP, 2.2*perCent);
    cartilage->AddElement(elS, 0.9*perCent);
    cartilage->AddElement(elCl,0.3*perCent);

    // 2/3 Cartilage, 1/3 Bone
    density = 1.35*g/cm3; numberofElements = 11;
    G4Material* cartilBone = new G4Material("CartilBone", density, numberofElements);
    cartilBone->AddElement(elH, 7.9745*perCent);
    cartilBone->AddElement(elC, 11.411*perCent);
    cartilBone->AddElement(elN, 2.8663*perCent);
    cartilBone->AddElement(elO, 64.4699*perCent);
    cartilBone->AddElement(elNa,0.3333*perCent);
    cartilBone->AddElement(elMg,0.0733*perCent);
    cartilBone->AddElement(elP, 4.9657*perCent);
    cartilBone->AddElement(elS, 0.705*perCent);
    cartilBone->AddElement(elCl,0.2*perCent);
    cartilBone->AddElement(elCa,6.9977*perCent);
    cartilBone->AddElement(elZn, 0.0033*perCent);
    
    // 1/3 Cartilage, 2/3 Bone
    density = 1.6*g/cm3; numberofElements = 11;
    G4Material* boneCartil = new G4Material("BoneCartil", density, numberofElements);
    boneCartil->AddElement(elH, 6.3489*perCent);
    boneCartil->AddElement(elC, 12.922*perCent);
    boneCartil->AddElement(elN, 3.5327*perCent);
    boneCartil->AddElement(elO, 54.5397*perCent);
    boneCartil->AddElement(elNa,0.1667*perCent);
    boneCartil->AddElement(elMg,0.1467*perCent);
    boneCartil->AddElement(elP, 7.7313*perCent);
    boneCartil->AddElement(elS, 0.51*perCent);
    boneCartil->AddElement(elCl,0.1*perCent);
    boneCartil->AddElement(elCa,13.9953*perCent);
    boneCartil->AddElement(elZn,0.0067*perCent);
    
    // Bone (ICRP - NIST)
    density = 1.85*g/cm3; numberofElements = 9;
    G4Material* bone = new G4Material ("Bone", density, numberofElements);
    bone->AddElement(elH, 4.7234*perCent);
    bone->AddElement(elC, 14.433*perCent);
    bone->AddElement(elN, 4.199*perCent);
    bone->AddElement(elO, 44.6096*perCent);
    bone->AddElement(elMg,0.22*perCent);
    bone->AddElement(elP, 10.497*perCent);
    bone->AddElement(elS, 0.315*perCent);
    bone->AddElement(elCa,20.993*perCent);
    bone->AddElement(elZn,0.01*perCent);

    // Denser Bone (ICRP - NIST)
    density = 2.1*g/cm3; numberofElements = 9;
    G4Material* denserBone = new G4Material ("DenserBone", density, numberofElements);
    denserBone->AddElement(elH, 4.7234*perCent);
    denserBone->AddElement(elC, 14.433*perCent);
    denserBone->AddElement(elN, 4.199*perCent);
    denserBone->AddElement(elO, 44.6096*perCent);
    denserBone->AddElement(elMg,0.22*perCent);
    denserBone->AddElement(elP, 10.497*perCent);
    denserBone->AddElement(elS, 0.315*perCent);
    denserBone->AddElement(elCa,20.993*perCent);
    denserBone->AddElement(elZn,0.01*perCent);

    // 1/3 Bone, 2/3 Aluminium
    density = 2.4*g/cm3; numberofElements = 10;
    G4Material* alBone = new G4Material("AlBone", density, numberofElements);
    alBone->AddElement(elH, 1.5745*perCent);
    alBone->AddElement(elC, 4.8110*perCent);
    alBone->AddElement(elN, 1.3997*perCent);
    alBone->AddElement(elO, 14.8699*perCent);
    alBone->AddElement(elMg,0.0733*perCent);
    alBone->AddElement(elAl,66.6666*perCent);
    alBone->AddElement(elP, 3.4990*perCent);
    alBone->AddElement(elS, 0.1050*perCent);
    alBone->AddElement(elCa,6.9977*perCent);
    alBone->AddElement(elZn,0.0033*perCent);

    // Aluminium
    density = 2.7*g/cm3; numberofElements = 1;
    G4Material* aluminum = new G4Material ("Aluminum", density, numberofElements);
    aluminum->AddElement(elAl, 100.0*perCent);

    // Denser Aluminium
    density =2.83*g/cm3; numberofElements = 1;
    G4Material* denserAluminum = new G4Material ("DenserAluminum", density, numberofElements);
    denserAluminum->AddElement(elAl, 100.0*perCent);

    // Iron
    density = 7.87*g/cm3; numberofElements = 1;
    G4Material* iiron = new G4Material ("IIron", density, numberofElements);
    iiron->AddElement(elFe, 100.0*perCent);

    fOriginalMaterials.push_back(AirMat);           // 0.00129 g/cm3
    fOriginalMaterials.push_back(lung);           // 0.50 g/cm3
    fOriginalMaterials.push_back(adiposeTissue);  // 0.95 g/cm3
    fOriginalMaterials.push_back(muscle);         // 1.05 g/cm3
    fOriginalMaterials.push_back(cartilage);      // 1.10 g/cm3
    fOriginalMaterials.push_back(cartilBone);     // 1.35 g/cm3
    fOriginalMaterials.push_back(boneCartil);     // 1.60 g/cm3
    fOriginalMaterials.push_back(bone);           // 1.85 g/cm3
    fOriginalMaterials.push_back(denserBone);     // 2.10 g/cm3
    fOriginalMaterials.push_back(alBone);         // 2.40 g/cm3
    fOriginalMaterials.push_back(aluminum);       // 2.70 g/cm3
    fOriginalMaterials.push_back(denserAluminum); // 2.83 g/cm3
    fOriginalMaterials.push_back(iiron);          // 7.87 g/cm3
}

void PhantomDICOM::readPhantomData()
{
    fNoFiles = SliceFiles.size();
    for (const auto & fname : SliceFiles)
        readPhantomDataFile(DataDir + '/' + fname + ".g4dcm");
}

void PhantomDICOM::readPhantomDataFile(const G4String & fname)
{
    out("Reading phantom data file", fname);
    std::ifstream fin(fname.c_str(), std::ios_base::in);
    if (!fin.is_open())
    {
        out("Cannot open phantom data file", fname);
        exit(10);
    }
    DicomPhantomZSliceHeader* sliceHeader = new DicomPhantomZSliceHeader(fin);
    fZSliceHeaders.push_back(sliceHeader);

    const size_t nVoxels = sliceHeader->GetNoVoxels();
    if (fZSliceHeaders.size() == 1) MaterialIDs.resize(fNoFiles * nVoxels); // if it is the first slice, initiliaze MaterialIDs

    size_t voxelCopyNo = (fZSliceHeaders.size() - 1) * nVoxels;
    size_t matID;
    for (size_t ii = 0; ii < nVoxels; ii++, voxelCopyNo++)
    {
        fin >> matID;
        //out("--->", fOriginalMaterials[matID]->GetName());
        MaterialIDs[voxelCopyNo] = matID;
    }
}

void PhantomDICOM::mergeZSliceHeaders()
{
    fZSliceHeaderMerged = new DicomPhantomZSliceHeader( *fZSliceHeaders[0] );
    for (size_t ii = 1; ii < fZSliceHeaders.size(); ii++)
        *fZSliceHeaderMerged += *fZSliceHeaders[ii];
}

void PhantomDICOM::computePhantomVoxelization()
{
    fNVoxelX = fZSliceHeaderMerged->GetNoVoxelX();
    fNVoxelY = fZSliceHeaderMerged->GetNoVoxelY();
    fNVoxelZ = fZSliceHeaderMerged->GetNoVoxelZ();

    fVoxelHalfDimX = fZSliceHeaderMerged->GetVoxelHalfX();
    fVoxelHalfDimY = fZSliceHeaderMerged->GetVoxelHalfY();
    fVoxelHalfDimZ = fZSliceHeaderMerged->GetVoxelHalfZ();

    zStart = -fNVoxelZ * fVoxelHalfDimZ + fVoxelHalfDimZ;
}

void PhantomDICOM::constructPhantomContainer(G4LogicalVolume * logicWorld)
{
    G4Tubs * solid   = new G4Tubs("PhContTube", 0, PhantRadius * mm, fNVoxelZ * fVoxelHalfDimZ * mm, 0, 360 * deg);
    fContainer_logic = new G4LogicalVolume(solid, AirMat, "PhContL", 0, 0, 0);

    const G4ThreeVector pos(PosInWorld[0]*mm, PosInWorld[1]*mm, PosInWorld[2]*mm);

    fContainer_phys = new G4PVPlacement(nullptr, pos, fContainer_logic, "PhCont", logicWorld, false, 1);

    if (ContainerInvisible) fContainer_logic->SetVisAttributes(G4VisAttributes::Invisible);
    else                    fContainer_logic->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 1.0)));
}

void PhantomDICOM::constructPhantom()
{
    const int voxelsPerSlice = fNVoxelX * fNVoxelY;

    const double radius = PhantRadius - 2.0 * fVoxelHalfDimX;
    const double radius2 = radius * radius;
    double XC = 0.5 * fNVoxelX;
    double YC = 0.5 * fNVoxelY;

    const int totNumVoxels = fNVoxelX * fNVoxelY * fNVoxelZ;

    for (int iVoxel = 0; iVoxel < totNumVoxels; iVoxel++)
    {
        const int iSlice = iVoxel / voxelsPerSlice;
        const int iVoxelIndexInSlice = iVoxel - iSlice * voxelsPerSlice;

        const int ix = iVoxelIndexInSlice / fNVoxelX;
        const int iy = iVoxelIndexInSlice % fNVoxelX;

        const double x = (ix - XC) * 2.0 * fVoxelHalfDimX;
        const double y = (iy - YC) * 2.0 * fVoxelHalfDimY;

        if (x*x + y*y < radius2)
        {
            G4Material * mat = fOriginalMaterials[MaterialIDs[iVoxel]];
            if (mat == AirMat) continue;

            const double z = zStart + 2.0 * iSlice * fVoxelHalfDimX;
            Voxels.push_back({-y, x, z, mat});
        }
    }

    G4Box           * voxel_solid = new G4Box("Voxel", fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ);
    G4LogicalVolume * voxel_logic = new G4LogicalVolume(voxel_solid, AirMat, "VoxelL",  0,0,0);

    DicomPhantomParameterisation * param = new DicomPhantomParameterisation(Voxels, ColourMap);
    new G4PVParameterised("phantom", voxel_logic, fContainer_logic, kUndefined, Voxels.size(), param);
}

void PhantomDICOM::readMaterialFile(const std::string & fileName)
{
    std::ifstream inStream(fileName);
    if (!inStream.is_open())
    {
        out("Cannot open data file with materials", fileName);
        exit(10);
    }

    out("Reading materials and the corresponding upper densities:");
    std::string name;
    float density;
    for (std::string line; std::getline(inStream, line); )
    {
        //out(">>>",line);
        if (line.empty()) continue; //allow empty lines

        std::stringstream ss(line);
        ss >> name >> density;
        if (ss.fail())
        {
            out("Unexpected format in data file with materials:", fileName);
            exit(10);
        }
        out(name, density);
        MatUpperDens.push_back({name, density});
    }
}

void PhantomDICOM::readColorMap(const std::string & fileName)
{
    std::ifstream inStream(fileName);
    if (!inStream.is_open())
    {
        out("Cannot open color map file", fileName);
        exit(10);
    }

    out("Reading color map");
    G4String name;
    double red, green, blue, opacity;
    for (std::string line; std::getline(inStream, line); )
    {
        out(">>>",line);
        if (line.empty()) continue;   //allow empty lines
        if (line[0] == '#') continue; //comment

        std::stringstream ss(line);
        ss >> name >> red >> green >> blue >> opacity;
        if (ss.fail())
        {
            out("Unexpected format in color map data file:", fileName);
            exit(10);
        }
        out(name, red, green, blue, opacity);
        G4Colour colour(red, green, blue, opacity);
        G4VisAttributes * visAtt = new G4VisAttributes(colour);
        visAtt->SetVisibility(opacity != 0);
        ColourMap[name] = visAtt;
    }
}

void PhantomDICOM::generateSliceFileNames()
{
    SliceFiles.clear();

    for (int iSlice = SliceFrom; iSlice < SliceTo; iSlice++)
    {
        std::string numStr = std::to_string(iSlice);
        if      (numStr.size() == 1) numStr = "00" + numStr;
        else if (numStr.size() == 2) numStr = "0" + numStr;
        std::string name = SliceFileBase + numStr;
        SliceFiles.push_back(name);
    }
}
