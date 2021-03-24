#include "DetectorConstruction.hh"
#include "SensitiveDetectorScint.hh"
#include "SessionManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4Navigator.hh"
#include "G4SDManager.hh"
#include <QDebug>

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4NistManager * man = G4NistManager::Instance();
    SessionManager & SM = SessionManager::getInstance();

    // Materials
    G4Material * matVacuum = man->FindOrBuildMaterial("G4_Galactic");

    std::vector<G4int> natoms;
    std::vector<G4String> elements;
    elements.push_back("C"); natoms.push_back(5);
    elements.push_back("H"); natoms.push_back(8);
    elements.push_back("O"); natoms.push_back(2);
    G4Material * matPMMA = man->ConstructNewMaterial("PMMA", elements, natoms, 1.18*g/cm3);

    //LYSO:Ce - https://www.crystals.saint-gobain.com/sites/imdf.crystals.com/files/documents/lyso-material-data-sheet.pdf
    natoms.clear();
    elements.clear();
    elements.push_back("Lu"); natoms.push_back(45);
    elements.push_back("Y") ; natoms.push_back(5);
    elements.push_back("Si"); natoms.push_back(25);
    elements.push_back("O") ; natoms.push_back(125);
    elements.push_back("Ce"); natoms.push_back(1); //<1% Cerium - https://authors.library.caltech.edu/8519/1/CHEieeetns05.pdf
    G4Material * matLYSO = man->ConstructNewMaterial("LYSO", elements, natoms, 7.1*g/cm3);

    SM.ScintMat = matLYSO;

    EncapsMat = man->FindOrBuildMaterial("G4_TEFLON");

    // Sensitive Detector
    SensitiveDetectorScint * pSD_Scint = new SensitiveDetectorScint("Scint");
    G4SDManager::GetSDMpointer()->AddNewDetector(pSD_Scint);

    // Geometry
    G4Box             * solidWorld = new G4Box("World", 500.0*mm, 500.0*mm, 500.0*mm);
                        logicWorld = new G4LogicalVolume(solidWorld, matVacuum, "World");
    G4VPhysicalVolume * physWorld  = new G4PVPlacement(nullptr, {0, 0, 0}, logicWorld, "World", nullptr, false, 0);
    logicWorld->SetVisAttributes(G4VisAttributes({0, 1, 0}));
    logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

    G4VSolid          * solidPmma = new G4Tubs("Cyl", 0, 100.0*mm, 100.0*mm, 0, 360.0*deg);
    G4LogicalVolume   * logicPmma = new G4LogicalVolume(solidPmma, matPMMA, "Cyl");
    //G4LogicalVolume   * logicPmma = new G4LogicalVolume(solidPmma, matVacuum, "Cyl");
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, 0}, logicPmma, "Target", logicWorld, false, 0);
    logicPmma->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 1.0)));

    solidScint = new G4Box("Scint", 0.5 * SM.ScintSizeX, 0.5 * SM.ScintSizeY, 0.5 * SM.ScintSizeZ);
    logicScint = new G4LogicalVolume(solidScint, SM.ScintMat, "Scint");
    logicScint->SetSensitiveDetector(pSD_Scint);
    logicScint->SetVisAttributes(G4VisAttributes({1, 0, 0}));

    solidEncaps = new G4Box("Encaps",  0.5 * SM.EncapsSizeX, 0.5 * SM.EncapsSizeY, 0.5 * SM.EncapsSizeZ);

    SM.ScintPositions.clear();

    int iAssembly = 0;
    int iScint    = 0;
    for (int iA = 0; iA < SM.NumSegments; iA++)
    {
        double Radius = 0.5 * (SM.InnerDiam + SM.EncapsSizeZ);
        double Angle  = SM.AngularStep * iA + SM.Angle0;
        double X = Radius * sin(Angle);
        double Y = Radius * cos(Angle);
        G4RotationMatrix * rot = new CLHEP::HepRotation(-Angle, 90.0*deg, 0);

        for (int iZ = 0; iZ < SM.NumRows; iZ++)
        {
            double RowPitch = SM.EncapsSizeY + SM.RowGap;
            double Z = -0.5 * (SM.NumRows - 1) * RowPitch  +  iZ * RowPitch;

            positionAssembly(rot, G4ThreeVector( X,  Y, Z), iScint, iAssembly++);
            positionAssembly(rot, G4ThreeVector(-X, -Y, Z), iScint, iAssembly++);
        }
    }

    return physWorld;
}

G4LogicalVolume * DetectorConstruction::createAssembly(int & iScint, G4RotationMatrix * AssemblyRot, G4ThreeVector AssemblyPos)
{
    SessionManager & SM = SessionManager::getInstance();

    G4LogicalVolume * logicEncaps = new G4LogicalVolume(solidEncaps, EncapsMat, "Encaps");
    logicEncaps->SetVisAttributes(G4VisAttributes({0, 1, 0}));

    for (int ix = 0; ix < SM.NumScintX; ix++)
        for (int iy = 0; iy < SM.NumScintY; iy++)
        {
            double X = -0.5 * (SM.NumScintX - 1) * SM.ScintPitchX  +  SM.ScintPitchX * ix;
            double Y = -0.5 * (SM.NumScintY - 1) * SM.ScintPitchY  +  SM.ScintPitchY * iy;
            G4ThreeVector ScintPos(X, Y, 0);
            new G4PVPlacement(nullptr, ScintPos, logicScint, "TT", logicEncaps, true, iScint++);

            //qDebug() << "Assembly position:"<< AssemblyPos[0] << AssemblyPos[1] << AssemblyPos[2];
            //qDebug() << "Local:" << ScintPos[0] << ScintPos[1] << ScintPos[2];
            G4ThreeVector glob = (*AssemblyRot).inverse()(ScintPos);
            //qDebug() << "Glob before shift:" << glob[0] << glob[1] << glob[2];
            glob[0] += AssemblyPos[0];
            glob[1] += AssemblyPos[1];
            glob[2] += AssemblyPos[2];
            //qDebug() << "Glob:" << glob[0] << glob[1] << glob[2];
            //qDebug() << "----";

            SM.ScintPositions.push_back(glob);
        }

    return logicEncaps;
}

void DetectorConstruction::positionAssembly(G4RotationMatrix * rot, G4ThreeVector pos, int & iScint, int iAssembly)
{
    new G4PVPlacement(rot, pos, createAssembly(iScint, rot, pos), "Encaps"+std::to_string(iAssembly), logicWorld, true, iAssembly);
}

