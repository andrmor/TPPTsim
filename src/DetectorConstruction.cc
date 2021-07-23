#include "DetectorConstruction.hh"
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
    if (SM.detectorContains(DetComp::Base))              addBase();
    if (SM.detectorContains(DetComp::SIPM))              addSIPM();
    if (SM.detectorContains(DetComp::PCB))               addPCB();
    if (SM.detectorContains(DetComp::CopperStructure))   addCopperStructure();

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
    logicScint->SetVisAttributes(G4VisAttributes({0, 0, 1}));

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
    logicEncaps->SetVisAttributes(G4VisAttributes({0, 1, 1}));

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

void DetectorConstruction::addBase()
{
    SessionManager & SM = SessionManager::getInstance();
    G4NistManager * man = G4NistManager::Instance();

    //Material: ABS, (C8H8·C4H6·C3H3N)n
    std::vector<G4int> natoms;
    std::vector<G4String> elements;
    elements.push_back("C"); natoms.push_back(15);
    elements.push_back("H"); natoms.push_back(17);
    elements.push_back("N"); natoms.push_back(1);
    G4Material * matABS = man->ConstructNewMaterial("ABS", elements, natoms, 1.07*g/cm3);

    G4Tubs * solidBase   = new G4Tubs("Base", SM.RMin, SM.RMax, SM.BaseHeight * 0.5, SM.Angle0 - 10.5 * deg, 120.0 * deg);
    G4LogicalVolume * logicBase   = new G4LogicalVolume(solidBase, matABS, "Base");
    G4Colour grey(0.5, 0.5, 0.5);
    logicBase ->SetVisAttributes(new G4VisAttributes(grey));

    G4RotationMatrix * rot  = new CLHEP::HepRotation(90*deg, 0, 0);
    G4RotationMatrix * rot1 = new CLHEP::HepRotation(-90.0*deg, 0, 0);

    new G4PVPlacement(rot, {0, 0, 0}, logicBase, "Base_PV", logicWorld, false, 0);
    new G4PVPlacement(rot1, {0, 0, 0}, logicBase, "Base_PV", logicWorld, false, 0);
    new G4PVPlacement(rot, {0, 0, SM.BaseHeight + SM.SystHeight}, logicBase, "Base_PV", logicWorld, false, 0);
    new G4PVPlacement(rot1, {0, 0, SM.BaseHeight + SM.SystHeight}, logicBase, "Base_PV", logicWorld, false, 0);
}

void DetectorConstruction::addSIPM()
{
    SessionManager & SM = SessionManager::getInstance();
    G4NistManager * man = G4NistManager::Instance();

    //Material: Epoxy
    std::vector<G4int> natoms;
    std::vector<G4String> elements;

    elements.push_back("C"); natoms.push_back(21);
    elements.push_back("H") ; natoms.push_back(25);
    elements.push_back("Cl") ; natoms.push_back(1);
    elements.push_back("O") ; natoms.push_back(5);
    G4Material * matEpoxy = man->ConstructNewMaterial("Epoxy", elements, natoms, 1.4*g/cm3);

    G4Box * solidSIPM   = new G4Box("SIPM", 0.5 * 25.8*mm, 0.5 * 25.8*mm, 0.5 * 1.35*mm);
    G4LogicalVolume * logicSIPM   = new G4LogicalVolume(solidSIPM, matEpoxy, "SIPM");
    G4Colour grey(0.5, 0.5, 0.5);
    logicSIPM ->SetVisAttributes(new G4VisAttributes(grey));

    for (int iA = 0; iA < SM.NumSegments; iA++)
    {
        double Radius = 0.5 * (SM.InnerDiam + 2 * SM.EncapsSizeZ + 1.35*mm);
        double Angle  = SM.AngularStep * iA + SM.Angle0;
        double X = Radius * sin(Angle);
        double Y = Radius * cos(Angle);
        G4RotationMatrix * rot  = new CLHEP::HepRotation(-Angle,             90.0*deg, 0);
        G4RotationMatrix * rot1 = new CLHEP::HepRotation(-Angle + 180.0*deg, 90.0*deg, 0);

        for (int iZ = 0; iZ < SM.NumRows; iZ++)
        {
            double RowPitch = SM.EncapsSizeY + SM.RowGap;
            double Z = -0.5 * (SM.NumRows - 1) * RowPitch  +  iZ * RowPitch + SM.GlobalZ0;

            new G4PVPlacement(rot, G4ThreeVector( X,  Y, Z), logicSIPM, "SIPM", logicWorld, false, 0);
            new G4PVPlacement(rot1, G4ThreeVector(-X, -Y, Z), logicSIPM, "SIPM", logicWorld, false, 0);
        }
    }
}

void DetectorConstruction::addPCB()
{
    SessionManager & SM = SessionManager::getInstance();
    G4NistManager * man = G4NistManager::Instance();

    //Material: PBC (FR4 -> Epoxy + Fiber Glass), https://agenda.infn.it/event/14179/contributions/23405/attachments/16712/18966/Geant_FR4_Raffaella_v1.pdf
    std::vector<G4int> natoms;
    std::vector<G4String> elements;
    int ncomponents;

    elements.push_back("C"); natoms.push_back(21);
    elements.push_back("H") ; natoms.push_back(25);
    elements.push_back("Cl") ; natoms.push_back(1);
    elements.push_back("O") ; natoms.push_back(5);
    G4Material * matEpoxy = man->ConstructNewMaterial("Epoxy", elements, natoms, 1.4*g/cm3);

    natoms.clear();
    elements.clear();
    elements.push_back("Si"); natoms.push_back(1);
    elements.push_back("O") ; natoms.push_back(2);
    G4Material * matFiberGlass = man->ConstructNewMaterial("FiberGlass", elements, natoms, 2.329*g/cm3);

    G4Material * matFR4 = new G4Material("FR4", 1.86*g/cm3, ncomponents = 2);
    matFR4 -> AddMaterial(matEpoxy, 0.472);
    matFR4 -> AddMaterial(matFiberGlass, 0.528);

    G4Box * solidPCB1   = new G4Box("PCB1", 0.5 * 28.9*mm, 0.5 * 52.2*mm, 0.5 * 3.175*mm);
    G4LogicalVolume * logicPCB1   = new G4LogicalVolume(solidPCB1, matFR4, "PCB1");
    logicPCB1 ->SetVisAttributes(G4VisAttributes({0, 1, 0}));

    G4Box * solidPCB2   = new G4Box("PCB2", 0.5 * 28.9*mm, 0.5 * 2*mm, 0.5 * 45*mm);
    G4LogicalVolume * logicPCB2   = new G4LogicalVolume(solidPCB2, matFR4, "PCB2");
    logicPCB2 ->SetVisAttributes(G4VisAttributes({0, 1, 0}));

    G4Box * solidPCB3   = new G4Box("PCB3", 0.5 * 28.9*mm, 0.5 * 52.2*mm, 0.5 * 2*mm);
    G4LogicalVolume * logicPCB3   = new G4LogicalVolume(solidPCB3, matFR4, "PCB3");
    logicPCB3 ->SetVisAttributes(G4VisAttributes({0, 0, 1}));

    for (int iA = 0; iA < SM.NumSegments; iA++)
    {
        double Radius1 = 0.5 * (SM.InnerDiam + 2 * SM.EncapsSizeZ + 2 * 1.35*mm + 2 * 0.010*mm + 3.175); //1.35 mm of SIPM + 0.010 mm between the SIPM and the PCB layers
        double Radius2 = 0.5 * (SM.InnerDiam + 2 * SM.EncapsSizeZ + 2 * 1.35*mm + 2 * 0.010*mm + 2 * 3.175 + 5*cm);
        double Radius3 = 0.5 * (SM.InnerDiam + 2 * SM.EncapsSizeZ + 2 * 1.35*mm + 2 * 0.010*mm + 2 * 3.175 + 2 * 5*cm + 2*mm);
        double Angle  = SM.AngularStep * iA + SM.Angle0;
        double X1 = Radius1 * sin(Angle);
        double Y1 = Radius1 * cos(Angle);
        double X2 = Radius2 * sin(Angle);
        double Y2 = Radius2 * cos(Angle);
        double X3 = Radius3 * sin(Angle);
        double Y3 = Radius3 * cos(Angle);
        G4RotationMatrix * rot  = new CLHEP::HepRotation(-Angle,             90.0*deg, 0);
        G4RotationMatrix * rot1 = new CLHEP::HepRotation(-Angle + 180.0*deg, 90.0*deg, 0);

        for (int iZ = 0; iZ < SM.NumRows / 2; iZ++) //Criar uma NumRows = 2 para os elementos de PCB
        {
            double RowPitch = 52.2*mm + SM.RowGap; //+ SM.RowGap; //SM.EncapsSizeY + SM.RowGap;
            double Z = -0.5 * (SM.NumRows / 2 - 1) * RowPitch  +  iZ * RowPitch + SM.GlobalZ0;

            new G4PVPlacement(rot, G4ThreeVector( X1,  Y1, Z), logicPCB1, "PCB1", logicWorld, false, 0);
            new G4PVPlacement(rot1, G4ThreeVector(-X1, -Y1, Z), logicPCB1, "PCB1", logicWorld, false, 0);

            new G4PVPlacement(rot, G4ThreeVector( X3,  Y3, Z), logicPCB3, "PCB3", logicWorld, false, 0);
            new G4PVPlacement(rot1, G4ThreeVector(-X3, -Y3, Z), logicPCB3, "PCB3", logicWorld, false, 0);
        }

        for (int iZ = 0; iZ < SM.NumRows; iZ++)
        {
            double RowPitch = 26.25 * mm; //+ SM.RowGap; //SM.EncapsSizeY + SM.RowGap;
            double Z = -0.5 * (SM.NumRows - 1) * RowPitch  +  iZ * RowPitch + SM.GlobalZ0;

            new G4PVPlacement(rot, G4ThreeVector( X2,  Y2, Z), logicPCB2, "PCB2", logicWorld, false, 0);
            new G4PVPlacement(rot1, G4ThreeVector(-X2, -Y2, Z), logicPCB2, "PCB2", logicWorld, false, 0);
        }
    }
}

void DetectorConstruction::addCopperStructure()
{
    SessionManager & SM = SessionManager::getInstance();
    G4NistManager * man = G4NistManager::Instance();

    //Material: Copper
    G4Material * matCopper = man->FindOrBuildMaterial("G4_Cu");

    //Columns between scintillators:
    G4Trd * solidCopper   = new G4Trd("Copper", 0.5 * 3.10 * mm, 0.5 * 2.40 * mm, 0.5 * 105 * mm, 0.5 * 105 * mm, 0.5 * 1.01 * mm);
    G4LogicalVolume * logicCopper   = new G4LogicalVolume(solidCopper, matCopper, "Copper");
    G4Colour brown(0.45,0.25,0.0);
    logicCopper ->SetVisAttributes(new G4VisAttributes(brown));

    for (int iA = 0; iA < 13; iA++)
    {
        double Radius = 0.5 * (SM.InnerDiam + 2 * SM.EncapsSizeZ + 1*mm);
        double Angle  = SM.AngularStep * iA + SM.Angle0 - 4.5*deg;
        double X1 = Radius * sin(Angle);
        double Y1 = Radius * cos(Angle);
        G4RotationMatrix * rot  = new CLHEP::HepRotation(-Angle,             90.0*deg, 0);
        G4RotationMatrix * rot1 = new CLHEP::HepRotation(-Angle + 180.0*deg, 90.0*deg, 0);

        for (int iZ = 0; iZ < 1; iZ++)
        {
            double RowPitch = 28.9 * mm;
            double Z = -0.5 * (1 - 1) * RowPitch  +  iZ * RowPitch + SM.GlobalZ0;

            new G4PVPlacement(rot, G4ThreeVector( X1,  Y1, Z), logicCopper, "Copper", logicWorld, false, 0);
            new G4PVPlacement(rot1, G4ThreeVector(-X1, -Y1, Z), logicCopper, "Copper", logicWorld, false, 0);
        }
    }

    G4Trd * solidCopper2   = new G4Trd("Copper2", 0.5 * 0.9 * mm, 0.5 * 0.7 * mm, 0.5 * 105 * mm, 0.5 * 105 * mm, 0.5 * 3.175 * mm);
    G4LogicalVolume * logicCopper2   = new G4LogicalVolume(solidCopper2, matCopper, "Copper");
    logicCopper2 ->SetVisAttributes(new G4VisAttributes(brown));

    for (int iA = 0; iA < 13; iA++)
    {
        double Radius = 0.5 * (SM.InnerDiam + 2 * SM.EncapsSizeZ + 2 * 1*mm + 2 * 0.010*mm + 3.175); //1 mm of SIPM + 0.010 mm between the SIPM and the PCB layers
        double Angle  = SM.AngularStep * iA + SM.Angle0 - 4.5 * deg;
        double X1 = Radius * sin(Angle);
        double Y1 = Radius * cos(Angle);
        G4RotationMatrix * rot  = new CLHEP::HepRotation(-Angle,             90.0*deg, 0);
        G4RotationMatrix * rot1 = new CLHEP::HepRotation(-Angle + 180.0*deg, 90.0*deg, 0);

        for (int iZ = 0; iZ < 1; iZ++)
        {
            double RowPitch = 0.35 * mm + 28.9 * mm;
            double Z = -0.5 * (1 - 1) * RowPitch  +  iZ * RowPitch + SM.GlobalZ0;

            new G4PVPlacement(rot, G4ThreeVector( X1,  Y1, Z), logicCopper2, "Copper2", logicWorld, false, 0);
            new G4PVPlacement(rot1, G4ThreeVector(-X1, -Y1, Z), logicCopper2, "Copper2", logicWorld, false, 0);
        }
    }

    //External columns:
    G4Box * solidCopper3   = new G4Box("Copper3", 0.5 * 3.0*mm, 0.5 * 3.0*mm, 0.5 * 105.0*mm);
    G4LogicalVolume * logicCopper3   = new G4LogicalVolume(solidCopper3, matCopper, "Copper3");
    logicCopper3 ->SetVisAttributes(new G4VisAttributes(brown));
    //new G4PVPlacement(0, G4ThreeVector(0, 0, SM.GlobalZ0), logicCopper3, "Copper3", logicWorld, false, 0);

    for (int iA = 0; iA < 13; iA++)
    {
        double Radius = 0.5 * (SM.InnerDiam + 2 * SM.EncapsSizeZ + 2 * 1.35*mm + 2 * 0.010*mm + 2 * 3.175 + 2 * 5*cm + 3.0*mm);
        double Angle  = SM.AngularStep * iA + SM.Angle0 - 4.5 * deg;
        double X1 = Radius * sin(Angle);
        double Y1 = Radius * cos(Angle);
        G4RotationMatrix * rot  = new CLHEP::HepRotation(-Angle,             0*deg, 0);
        G4RotationMatrix * rot1 = new CLHEP::HepRotation(-Angle + 180.0*deg, 0*deg, 0);

        for (int iZ = 0; iZ < 1; iZ++)
        {
            double RowPitch = 0.35 * mm + 28.9 * mm;
            double Z = -0.5 * (1 - 1) * RowPitch  +  iZ * RowPitch + SM.GlobalZ0;

            new G4PVPlacement(rot, G4ThreeVector( X1,  Y1, Z), logicCopper3, "Copper3", logicWorld, false, 0);
            new G4PVPlacement(rot1, G4ThreeVector(-X1, -Y1, Z), logicCopper3, "Copper3", logicWorld, false, 0);
        }
    }

    //Holders:
    G4Box * solidHolder   = new G4Box("Holder", 0.5 * 3.0*mm, 0.5 * 3.0*mm, 0.5 * 10.0*mm);
    G4LogicalVolume * logicHolder   = new G4LogicalVolume(solidHolder, matCopper, "Holder");
    logicHolder ->SetVisAttributes(new G4VisAttributes(brown));

    for (int iA = 0; iA < 13; iA++)
    {
        double Radius = 0.5 * (SM.InnerDiam + 2 * SM.EncapsSizeZ + 2 * 1.35*mm + 2 * 0.010*mm + 2 * 3.175 + 2 * 5*cm + 2 * 3.0*mm + 10.0*mm);
        double Angle  = SM.AngularStep * iA + SM.Angle0 - 4.5 * deg;
        double X1 = Radius * sin(Angle);
        double Y1 = Radius * cos(Angle);
        G4RotationMatrix * rot  = new CLHEP::HepRotation(-Angle,             90.0*deg, 0);
        G4RotationMatrix * rot1 = new CLHEP::HepRotation(-Angle + 180.0*deg, 90.0*deg, 0);

        for (int iZ = 0; iZ < 3; iZ++)
        {
            double RowPitch = 35*mm;
            double Z = -0.5 * (3 - 1) * RowPitch  +  iZ * RowPitch + SM.GlobalZ0;

            new G4PVPlacement(rot, G4ThreeVector( X1,  Y1, Z), logicHolder, "Holder", logicWorld, false, 0);
            new G4PVPlacement(rot1, G4ThreeVector(-X1, -Y1, Z), logicHolder, "Holder", logicWorld, false, 0);
        }
    }
}
