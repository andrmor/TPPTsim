#include "DetectorConstruction.hh"
#include "DetComp.hh"
#include "SensitiveDetectorScint.hh"
#include "SensitiveDetectorFSM.hh"
#include "SessionManager.hh"
#include "SimMode.hh"
#include "PhantomMode.hh"
#include "out.hh"
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

#include "G4GDMLParser.hh"

#define _USE_MATH_DEFINES
#include <math.h>

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

    //EncapsMat = man->FindOrBuildMaterial("G4_TEFLON");
    EncapsMat = matPMMA;

    // Geometry
    if ( SM.detectorContains(DetComp::GDML) )
    {
        G4GDMLParser parser;
        parser.Read("mother.gdml", false);
        SM.physWorld  = parser.GetWorldVolume();
        logicWorld = SM.physWorld->GetLogicalVolume();
    }
    else
    {
        G4Box * solidWorld   = new G4Box("World", 500.0*mm, 500.0*mm, 500.0*mm);
                logicWorld   = new G4LogicalVolume(solidWorld, matVacuum, "World");
                SM.physWorld = new G4PVPlacement(nullptr, {0, 0, 0}, logicWorld, "World", nullptr, false, 0);
    }
    //logicWorld->SetVisAttributes(G4VisAttributes({0, 1, 0}));
    logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

    G4LogicalVolume * phantom = SM.PhantomMode->definePhantom(logicWorld);
    if (phantom) SM.createPhantomRegion(phantom);

    if (SM.detectorContains(DetComp::Scintillators))     addScintillators();
    if (SM.detectorContains(DetComp::FirstStageMonitor)) addFSM(matVacuum);

    return SM.physWorld;
}

void DetectorConstruction::addFSM(G4Material * material)
{
    SessionManager & SM = SessionManager::getInstance();

    double OuterRadius = 0.5 * SM.InnerDiam - 2.0*mm;
    double InnerRadius = OuterRadius - 1.0*mm;
    double Height      = 5.0 * SM.EncapsSizeX; // Margarida, please calculate the minimum size

    G4VSolid          * solidFSM = new G4Tubs("FSM_Cyl", InnerRadius, OuterRadius, 0.5 * Height, 0, 360.0*deg);
    G4LogicalVolume   * logicFSM = new G4LogicalVolume(solidFSM, material, "FSM");
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, SM.GlobalZ0}, logicFSM, "FSM_PV", logicWorld, false, 0);
    logicFSM->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 1.0)));

    G4VSensitiveDetector * pSD_FSM = new SensitiveDetectorFSM("SD_FSM");
    G4SDManager::GetSDMpointer()->AddNewDetector(pSD_FSM);
    logicFSM->SetSensitiveDetector(pSD_FSM);
}

void DetectorConstruction::addScintillators()
{
    SessionManager & SM = SessionManager::getInstance();

    solidScint = new G4Box("Scint", 0.5 * SM.ScintSizeX, 0.5 * SM.ScintSizeY, 0.5 * SM.ScintSizeZ);
    logicScint = new G4LogicalVolume(solidScint, SM.ScintMat, "Scint"); //SiPM are interfaced at the local negative Z
    logicScint->SetVisAttributes(G4VisAttributes({1, 0, 0}));

    SM.createScintillatorRegion(logicScint);

    /*  just to check that positive Z is towards the phantom:
    G4Box * b   = new G4Box("b", 1.0*mm, 1.0*mm, 1.0*mm);
    G4LogicalVolume * lb  = new G4LogicalVolume(b, matVacuum, "bb");
    new G4PVPlacement(nullptr, {0, 0, 5}, lb, "ss", logicScint, false, 0);
    */

    solidEncaps = new G4Box("Encaps",  0.5 * SM.EncapsSizeX, 0.5 * SM.EncapsSizeY, 0.5 * SM.EncapsSizeZ);

    SM.ScintRecords.clear();

    int iAssembly = 0;
    int iScint    = 0;
    for (int iA = 0; iA < SM.NumSegments; iA++)
    {
        double Radius = 0.5 * (SM.InnerDiam + SM.EncapsSizeZ);
        double Angle  = SM.AngularStep * iA + SM.Angle0;
        double X = Radius * sin(Angle);
        double Y = Radius * cos(Angle);
        G4RotationMatrix * rot  = new CLHEP::HepRotation(-Angle,             90.0*deg, 0);
        G4RotationMatrix * rot1 = new CLHEP::HepRotation(-Angle + 180.0*deg, 90.0*deg, 0);

        for (int iZ = 0; iZ < SM.NumRows; iZ++)
        {
            double RowPitch = SM.EncapsSizeY + SM.RowGap;
            double Z = -0.5 * (SM.NumRows - 1) * RowPitch  +  iZ * RowPitch + SM.GlobalZ0;

            positionAssembly(rot,  G4ThreeVector( X,  Y, Z), Angle,             iScint, iAssembly++, 0);
            positionAssembly(rot1, G4ThreeVector(-X, -Y, Z), Angle + 0.5*M_PI , iScint, iAssembly++, 1);
        }
    }

    // Sensitive Detector
    G4VSensitiveDetector * pSD_Scint = SM.SimMode->getScintDetector();
    if (pSD_Scint)
    {
        G4SDManager::GetSDMpointer()->AddNewDetector(pSD_Scint);
        logicScint->SetSensitiveDetector(pSD_Scint);
    }

    SM.saveScintillatorTable(SM.WorkingDirectory + '/' + "LUT.txt");
}

G4LogicalVolume * DetectorConstruction::createAssembly(int & iScint, G4RotationMatrix * AssemblyRot, G4ThreeVector AssemblyPos, double Angle, int headNumber)
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
            new G4PVPlacement(nullptr, ScintPos, logicScint, "Scint", logicEncaps, true, iScint++);

            ScintRecord rec;

            rec.CenterPos = (*AssemblyRot).inverse()(ScintPos);
            for (int i=0; i<3; i++) rec.CenterPos[i] += AssemblyPos[i];

            ScintPos[2] += 0.5 * SM.ScintSizeZ;
            rec.FacePos   = (*AssemblyRot).inverse()(ScintPos);
            for (int i=0; i<3; i++) rec.FacePos[i]   += AssemblyPos[i];

            rec.Angle = Angle;
            rec.HeadNumber = headNumber;

            SM.ScintRecords.push_back(rec);
        }

    return logicEncaps;
}

void DetectorConstruction::positionAssembly(G4RotationMatrix * rot, G4ThreeVector pos, double angle, int & iScint, int iAssembly, int headNumber)
{
    new G4PVPlacement(rot, pos, createAssembly(iScint, rot, pos, angle, headNumber), "Encaps"+std::to_string(iAssembly), logicWorld, true, iAssembly);
}
