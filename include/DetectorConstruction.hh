#ifndef DetectorConstruction_H
#define DetectorConstruction_H

#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    G4VPhysicalVolume * Construct();

private:
    G4LogicalVolume * createAssembly(int & iScint, G4RotationMatrix * AssemblyRot, G4ThreeVector AssemblyPos, double Angle, int headNumber);
    void addFSM(G4Material *material);
    void addScintillators();
    void positionAssembly(G4RotationMatrix * rot, G4ThreeVector pos, double angle, int & iScint, int iAssembly, int headNumber);
    void addBase();
    void addSIPM();
    void addPCB();
    void addCopperStructure();

    G4VSolid * solidEncaps = nullptr;
    G4VSolid * solidScint  = nullptr;

    G4Material * EncapsMat = nullptr;

    G4LogicalVolume * logicWorld = nullptr;
    G4LogicalVolume * logicScint = nullptr;
};

#endif //DetectorConstruction_H
