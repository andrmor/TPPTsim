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
    G4LogicalVolume * createAssembly(int & iScint, G4RotationMatrix * AssemblyRot, G4ThreeVector AssemblyPos);
    void positionAssembly(G4RotationMatrix * rot, G4ThreeVector pos, int & iScint, int iAssembly);

    G4VSolid * solidEncaps;
    G4VSolid * solidScint;

    G4Material * EncapsMat;

    G4LogicalVolume * logicWorld;
    G4LogicalVolume * logicScint;
};

#endif //DetectorConstruction_H
