#include "PhantomMode.hh"
#include "SessionManager.hh"
#include "DicomPhantom.hh"
#include "out.hh"
#include "jstools.hh"

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

#include <math.h>

PhantomModeBase * PhantomModeFactory::makePhantomModeInstance(const json11::Json & json)
{
    out("Reading phantom json");
    std::string Type;
    jstools::readString(json, "Type", Type);

    PhantomModeBase * ph = nullptr;

    if      (Type == "PhantomNone")       ph = new PhantomNone();
    else if (Type == "PhantomCylinder")   ph = new PhantomCylinder();
    else if (Type == "PhantomBox")        ph = new PhantomBox();
    else if (Type == "PhantomDerenzo")    ph = new PhantomDerenzo(100.0, 100.0, {}, 0, 0, 0);
    else if (Type == "PhantomParam")      ph = new PhantomParam();
    else if (Type == "PhantomDICOM")      ph = new PhantomDICOM("", "", 0,0, 1, 100.0, {0,0,50.0}, {0,0,0});
    else if (Type == "PhantomEspana")     ph = new PhantomEspana();
    else if (Type == "PhantomBauerGel")   ph = new PhantomBauerGel();
    else if (Type == "PhantomBauerCa")    ph = new PhantomBauerCa();
    else if (Type == "PhantomMarekWater") ph = new PhantomMarekWater();
    else
    {
        out("Unknown phantom type!");
        exit(10);
    }

    if (ph) ph->readFromJson(json);
    return ph;
}

// ---

void PhantomModeBase::writeToJson(json11::Json::object & json) const
{
    json["Type"] = getTypeName();
    doWriteToJson(json);
}

// ---

PhantomCylinder::PhantomCylinder(double diameter, double length, std::string g4_material_name) :
    Diameter(diameter), Length(length), MatBuilder(new MaterialBuilder(g4_material_name)) {}

PhantomCylinder::PhantomCylinder(double diameter, double length, EMaterial material) :
    Diameter(diameter), Length(length), MatBuilder(new MaterialBuilder(material)) {}

PhantomCylinder::PhantomCylinder() : MatBuilder(new MaterialBuilder(EMaterial::Undefined)) {}

G4LogicalVolume * PhantomCylinder::definePhantom(G4LogicalVolume * logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();

    G4Material * mat = MatBuilder->build();
    G4VSolid          * solidPmma = new G4Tubs("Phantom_Cyl", 0, 0.5*Diameter, 0.5*Length, 0, 360.0*deg);
    G4LogicalVolume   * logicPmma = new G4LogicalVolume(solidPmma, mat, "Phantom");
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, SM.GlobalZ0}, logicPmma, "Phantom_PV", logicWorld, false, 0);
    logicPmma->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    return logicPmma;
}

void PhantomCylinder::readFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "Diameter", Diameter);
    jstools::readDouble(json, "Length",   Length);

    MatBuilder->readFromJson(json);
}

void PhantomCylinder::doWriteToJson(json11::Json::object & json) const
{
    json["Diameter"] = Diameter;
    json["Length"]   = Length;

    MatBuilder->writeToJson(json);
}

// ---

PhantomBox::PhantomBox(double sizeX, double sizeY, double sizeZ, std::string g4_material_name) :
    SizeX(sizeX), SizeY(sizeY), SizeZ(sizeZ), MatBuilder(new MaterialBuilder(g4_material_name)) {}

PhantomBox::PhantomBox(double sizeX, double sizeY, double sizeZ, EMaterial material) :
    SizeX(sizeX), SizeY(sizeY), SizeZ(sizeZ), MatBuilder(new MaterialBuilder(material)) {}

PhantomBox::PhantomBox() : MatBuilder(new MaterialBuilder(EMaterial::Undefined)) {}

G4LogicalVolume * PhantomBox::definePhantom(G4LogicalVolume * logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();

    G4Material * mat = MatBuilder->build();
    G4VSolid          * solid = new G4Box("Phantom_Box", 0.5 * SizeX * mm, 0.5 * SizeY * mm, 0.5 * SizeZ * mm);
    G4LogicalVolume   * logic = new G4LogicalVolume(solid, mat, "Phantom");
    new G4PVPlacement(nullptr, {0, 0, SM.GlobalZ0}, logic, "Phantom_PV", logicWorld, false, 0);
    logic->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    return logic;
}

void PhantomBox::readFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "SizeX", SizeX);
    jstools::readDouble(json, "SizeY", SizeY);
    jstools::readDouble(json, "SizeZ", SizeZ);

    MatBuilder->readFromJson(json);
}

void PhantomBox::doWriteToJson(json11::Json::object & json) const
{
    json["SizeX"] = SizeX;
    json["SizeY"] = SizeY;
    json["SizeZ"] = SizeZ;

    MatBuilder->writeToJson(json);
}

// ---

G4LogicalVolume * PhantomDerenzo::definePhantom(G4LogicalVolume *logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();
    G4NistManager * man = G4NistManager::Instance();

    std::vector<G4int> natoms;
    std::vector<G4String> elements;
    elements.push_back("C"); natoms.push_back(5);
    elements.push_back("H"); natoms.push_back(8);
    elements.push_back("O"); natoms.push_back(2);
    G4Material * matBulk = man->ConstructNewMaterial("PMMA_phantom", elements, natoms, 1.18*g/cm3); //man->FindOrBuildMaterial("G4_Galactic");

    G4Material * matEmitter = man->FindOrBuildMaterial("G4_WATER"); //"G4_AIR"

    G4VSolid          * solidPmma = new G4Tubs("Phantom_Cyl", 0, 0.5*Diameter*mm, 0.5*Height*mm, 0, 360.0*deg);
    G4LogicalVolume   * logicPmma = new G4LogicalVolume(solidPmma, matBulk, "Phantom");
    //new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, SM.GlobalZ0}, logicPmma, "Phantom_PV", logicWorld, false, 0);
    new G4PVPlacement(nullptr, {0, 0, SM.GlobalZ0}, logicPmma, "Phantom_PV", logicWorld, false, 0);
    logicPmma->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 1.0)));

    const int numSectors = HoleDiameters.size();
    const double deltaAngle = 2.0 * M_PI / numSectors;

    for (int iSec = 0; iSec < numSectors; iSec++)
    {
        const double holeDiameter = HoleDiameters[iSec];
        double a = deltaAngle * iSec + DPhi * M_PI / 180.0;

        //x and y are in not rotated coordinates, sector is upward: centered at x=0, y+)
        double y = RadialOffset;
        int iRowCounter = 0;
        while (y + 0.5*holeDiameter < 0.5*Diameter - 1.0)
        {
            y = RadialOffset + iRowCounter * 2.0*holeDiameter * cos(M_PI / 6.0);

            G4VSolid        * solidHole = new G4Tubs("sh", 0, 0.5*holeDiameter*mm, 0.5*Height*mm - 1.0*mm, 0, 360.0*deg);
            G4LogicalVolume * logicHole = new G4LogicalVolume(solidHole, matEmitter, "Hole" + std::to_string(iSec));
            logicPmma->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 1.0)));

            for (int iHole = 0; iHole <= iRowCounter; iHole++)
            {
                double x = (iHole - 0.5*iRowCounter) * 2.0*holeDiameter;
                if (sqrt(x*x + y*y) + 0.5*holeDiameter < 0.5*Diameter - Margin)
                {
                    double xr = x*cos(a) - y*sin(a);
                    double yr = x*sin(a) + y*cos(a);
                    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {xr, yr, 0}, logicHole, "Phantom_PV", logicPmma, false, 0);
                }
            }

            iRowCounter++;
        }
    }

    return logicPmma;
}

void PhantomDerenzo::doWriteToJson(json11::Json::object &json) const
{
    json["Diameter"]      = Diameter;
    json["Height"]        = Height;

    json11::Json::array ar;
    for (const double d : HoleDiameters) ar.push_back(d);
    json["HoleDiameters"] = ar;

    json["RadialOffset"]  = RadialOffset;
    json["Margin"]        = Margin;
    json["DPhi"]          = DPhi;
}

void PhantomDerenzo::readFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "Diameter",     Diameter);
    jstools::readDouble(json, "Height",       Height);

    json11::Json::array ar;
    jstools::readArray(json, "HoleDiameters", ar);
    HoleDiameters.clear();
    std::string str;
    for (size_t i = 0; i < ar.size(); i++)
    {
        if (i != 0) str += ", ";
        const double dia = ar[i].number_value();
        HoleDiameters.push_back(dia);
        str += std::to_string(dia);
    }
    out("HoleDiameters:", str);

    jstools::readDouble(json, "RadialOffset", RadialOffset);
    jstools::readDouble(json, "Margin",       Margin);
    jstools::readDouble(json, "DPhi",         DPhi);
}

// ---

#include "ParamTest.hh"
#include "G4PVParameterised.hh"
G4LogicalVolume * PhantomParam::definePhantom(G4LogicalVolume * logicWorld)
{
    SessionManager & SM  = SessionManager::getInstance();
    G4NistManager  * man = G4NistManager::Instance();

    G4Material        * matVacuum = man->FindOrBuildMaterial("G4_Galactic");
    G4VSolid          * solidCyl  = new G4Tubs("Phantom_Cyl", 0, 100.0*mm, 100.0*mm, 0, 360.0*deg);
    G4LogicalVolume   * logicCyl  = new G4LogicalVolume(solidCyl, matVacuum, "Phantom");
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, SM.GlobalZ0}, logicCyl, "Phantom_PV", logicWorld, false, 0);
    logicCyl->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 1.0)));

    G4Material      * matBox   = man->FindOrBuildMaterial("G4_TEFLON");
    G4VSolid        * solidBox = new G4Box("chamber", 5.0*mm, 5.0*mm, 5.0*mm);
    G4LogicalVolume * logicBox = new G4LogicalVolume(solidBox, matBox, "Box", 0, 0, 0);

    G4VPVParameterisation * param = new TestParameterisation({10.0*mm, 10.0*mm, 10.0*mm}, 100.0*mm, -50.0*mm);
    new G4PVParameterised("DicomPhant", logicBox, logicCyl, kUndefined, 317*10, param);

    return logicCyl;
}

// ---

G4LogicalVolume * PhantomEnergyCalibration::definePhantom(G4LogicalVolume * logicWorld)
{
    G4NistManager * man = G4NistManager::Instance();

    G4Material * matW = man->FindOrBuildMaterial("G4_WATER");

    const double halfSizeY  = 200.0*mm;  // along the beam!
    const double halfSizeXZ =  75.0*mm;

    G4VSolid          * solid = new G4Box("wPh", halfSizeXZ, halfSizeY, halfSizeXZ);
    G4LogicalVolume   * logic = new G4LogicalVolume(solid, matW, "PhantWaterLog");
    logic->SetVisAttributes(G4Color::Blue());
    new G4PVPlacement(nullptr, {0, -halfSizeY, 0}, logic, "ltp", logicWorld, false, 0);

    return logic;
}

// ---

G4LogicalVolume *PhantomEspana::definePhantom(G4LogicalVolume * logicWorld)
{
    G4NistManager * man = G4NistManager::Instance();
    G4Material * matVacuum = man->FindOrBuildMaterial("G4_Galactic");

    std::vector<double> weightFrac;
    std::vector<G4String> elements;
    elements.push_back("H"); weightFrac.push_back(14.3);
    elements.push_back("C"); weightFrac.push_back(85.7);
    G4Material * matHDPE = man->ConstructNewMaterial("HDPE", elements, weightFrac, 0.95*g/cm3);

    weightFrac.clear();
    elements.clear();
    elements.push_back("H"); weightFrac.push_back(9.6);
    elements.push_back("C"); weightFrac.push_back(14.9);
    elements.push_back("N"); weightFrac.push_back(1.46);
    elements.push_back("O"); weightFrac.push_back(73.8);
    G4Material * matTissue = man->ConstructNewMaterial("GelTissue", elements, weightFrac, 1.13*g/cm3);

    weightFrac.clear();
    elements.clear();
    elements.push_back("H"); weightFrac.push_back(11.03);
    elements.push_back("C"); weightFrac.push_back(1.04);
    elements.push_back("N"); weightFrac.push_back(0.32);
    elements.push_back("O"); weightFrac.push_back(87.6);
    G4Material * matWater = man->ConstructNewMaterial("GelWater", elements, weightFrac, 1.01*g/cm3);

    G4VSolid * solidCont = new G4Box("Phantom_Box",  60.0*mm, 60.0*mm, 60.0*mm);
    G4VSolid * solidHDPE = new G4Box("Phantom_HDPE", 30.0*mm, 60.0*mm, 60.0*mm);
    G4VSolid * solidTis  = new G4Box("Phantom_Tis",  30.0*mm, 60.0*mm, 30.0*mm);
    G4VSolid * solidWat  = new G4Box("Phantom_Wat",  30.0*mm, 60.0*mm, 30.0*mm);

    G4LogicalVolume * logicCont = new G4LogicalVolume(solidCont, matVacuum, "Phantom");
    G4LogicalVolume * logicHDPE = new G4LogicalVolume(solidHDPE, matHDPE,   "HDPE_L");
        logicHDPE->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0, 0)));
    G4LogicalVolume * logicTis  = new G4LogicalVolume(solidTis, matTissue,  "TisG_L");
        logicTis->SetVisAttributes (G4VisAttributes(G4Colour(0, 1.0, 0)));
    G4LogicalVolume * logicWat  = new G4LogicalVolume(solidWat, matWater,   "WatG_L");
        logicWat->SetVisAttributes (G4VisAttributes(G4Colour(0, 0, 1.0)));

    new G4PVPlacement(nullptr, {0, 0, 0}, logicCont, "Phantom_PV", logicWorld, false, 0);
    new G4PVPlacement(nullptr, {30.0, 0, 0}, logicHDPE, "HDPE_PV", logicCont, false, 0);
    new G4PVPlacement(nullptr, {-30.0, 0, -30.0}, logicTis, "TisG_PV", logicCont, false, 0);
    new G4PVPlacement(nullptr, {-30.0, 0, 30.0}, logicWat, "Wat_PV", logicCont, false, 0);

    return logicCont;
}

G4LogicalVolume * PhantomBauerGel::definePhantom(G4LogicalVolume * logicWorld)
{
    G4NistManager * man = G4NistManager::Instance();

    std::vector<G4int> natoms;
    std::vector<G4String> elements;
    elements.push_back("C"); natoms.push_back(5);
    elements.push_back("H"); natoms.push_back(8);
    elements.push_back("O"); natoms.push_back(2);
    G4Material * matPMMA = man->ConstructNewMaterial("PMMA_phantom", elements, natoms, 1.18*g/cm3);

    std::vector<double> weightFrac;
    elements.clear();
    elements.push_back("H"); weightFrac.push_back(11.03);
    elements.push_back("C"); weightFrac.push_back(1.04);
    elements.push_back("N"); weightFrac.push_back(0.32);
    elements.push_back("O"); weightFrac.push_back(87.6);
    G4Material * matGel = man->ConstructNewMaterial("GelWater", elements, weightFrac, 1.01*g/cm3);

    G4VSolid          * solidPMMA = new G4Box("Phantom_Box", 0.5 * 110 * mm, 0.5 * 370 * mm, 0.5 * 110 * mm);
    G4VSolid          * solidGel  = new G4Box("Phantom_Box", 0.5 * 100 * mm, 0.5 * 350 * mm, 0.5 * 100 * mm);

    G4LogicalVolume   * logicPMMA = new G4LogicalVolume(solidPMMA, matPMMA, "PhantomBox");
    logicPMMA->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));
    G4LogicalVolume   * logicGel  = new G4LogicalVolume(solidGel, matGel, "PhantomGel");
    logicGel->SetVisAttributes(G4VisAttributes(G4Colour(0, 1.0, 1.0)));

    new G4PVPlacement(nullptr, {0, 0, 0}, logicPMMA, "PhantomBox_PV", logicWorld, false, 0);
    new G4PVPlacement(nullptr, {0, 0, 0}, logicGel,  "PhantomGel_PV", logicPMMA,  false, 0);

    return logicPMMA;
}

// ---

G4LogicalVolume * PhantomBauerCa::definePhantom(G4LogicalVolume *logicWorld)
{
    G4NistManager * man = G4NistManager::Instance();

    std::vector<G4int> natoms;
    std::vector<G4String> elements;
    elements.push_back("C"); natoms.push_back(5);
    elements.push_back("H"); natoms.push_back(8);
    elements.push_back("O"); natoms.push_back(2);
    G4Material * matPMMA = man->ConstructNewMaterial("PMMA_phantom", elements, natoms, 1.18*g/cm3);

    natoms.clear();
    elements.clear();
    elements.push_back("Ca"); natoms.push_back(1);
    elements.push_back("C");  natoms.push_back(1);
    elements.push_back("O");  natoms.push_back(3);
    G4Material * matChalk = man->ConstructNewMaterial("GelWater", elements, natoms, 1.22*g/cm3);

    G4VSolid          * solidPMMA  = new G4Box("Phantom_Box", 0.5 * 110 * mm, 0.5 * 370 * mm, 0.5 * 110 * mm);
    G4VSolid          * solidChalk = new G4Box("Phantom_Box", 0.5 * 100 * mm, 0.5 * 350 * mm, 0.5 * 100 * mm);

    G4LogicalVolume   * logicPMMA  = new G4LogicalVolume(solidPMMA, matPMMA, "PhantomBox");
    logicPMMA->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));
    G4LogicalVolume   * logicChalk = new G4LogicalVolume(solidChalk, matChalk, "PhantomChalk");
    logicChalk->SetVisAttributes(G4VisAttributes(G4Colour(0, 1.0, 1.0)));

    new G4PVPlacement(nullptr, {0, 0, 0}, logicPMMA, "PhantomBox_PV",    logicWorld, false, 0);
    new G4PVPlacement(nullptr, {0, 0, 0}, logicChalk, "PhantomChalk_PV", logicPMMA,  false, 0);

    return logicPMMA;
}

G4LogicalVolume * PhantomRT::definePhantom(G4LogicalVolume * logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();
    G4NistManager * man = G4NistManager::Instance();

    std::vector<G4int> natoms;
    std::vector<G4String> elements;
    elements.push_back("C"); natoms.push_back(5);
    elements.push_back("H"); natoms.push_back(8);
    elements.push_back("O"); natoms.push_back(2);
    G4Material * matPMMA = man->ConstructNewMaterial("PMMA_phantom", elements, natoms, 1.18*g/cm3);

    G4VSolid          * solidPmma = new G4Tubs("Phantom_Cyl", 0, 100.0*mm, 100.0*mm, 0, 360.0*deg);
    G4LogicalVolume   * logicPmma = new G4LogicalVolume(solidPmma, matPMMA, "Phantom");
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, SM.GlobalZ0}, logicPmma, "Phantom_PV", logicWorld, false, 0);
    logicPmma->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));
    logicPmma->SetVisAttributes(false);

    G4Material * matW = man->FindOrBuildMaterial("G4_WATER");

    G4VSolid          * solidWB = new G4Tubs("WB", 0, 2.0*mm, 30.0*mm, 0, 360.0*deg);
    G4LogicalVolume   * logicWB = new G4LogicalVolume(solidWB, matW, "WBl");
    logicWB->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0, 0)));
    new G4PVPlacement(new CLHEP::HepRotation(0,  10.0*deg, 0), {20.0,  7.2, 0}, logicWB, "WBp", logicPmma, false, 0);
    new G4PVPlacement(new CLHEP::HepRotation(0, -10.0*deg, 0), {20.0, -7.2, 0}, logicWB, "WBp", logicPmma, false, 0);

    G4VSolid          * solidWS = new G4Tubs("WS", 0, 2.0*mm, 4.7*mm, 0, 360.0*deg);
    G4LogicalVolume   * logicWS = new G4LogicalVolume(solidWS, matW, "WSl");
    logicWS->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0, 0)));
    new G4PVPlacement(new CLHEP::HepRotation(0,  90.0*deg, 0), {20.0, 0, 0}, logicWS, "WSp", logicPmma, false, 0);

    /*
    G4VSolid          * st = new G4Box("t", 20.0*mm, 3.0*mm, 33.0*mm);
    G4LogicalVolume   * lt = new G4LogicalVolume(st, matW, "tl");
    lt->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 1.0)));
    new G4PVPlacement(nullptr, {0, 20.0, 0}, lt, "ltp", logicWorld, false, 0);
    */

    return logicPmma;
}

G4LogicalVolume * PhantomMarekWater::definePhantom(G4LogicalVolume * logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();

    G4NistManager * man = G4NistManager::Instance();

    G4Material * matPMMA  = man->FindOrBuildMaterial("G4_PLEXIGLASS");
    G4Material * matWater = man->FindOrBuildMaterial("G4_WATER");

    G4VSolid          * solidPmma = new G4Tubs("Phantom_Outer", 0, 0.5*37.8, 0.5*104.2, 0, 360.0*deg);
    G4LogicalVolume   * logicPmma = new G4LogicalVolume(solidPmma, matPMMA, "Phantom_Outer");
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, SM.GlobalZ0}, logicPmma, "Phantom_Outer_PV", logicWorld, false, 0);
    logicPmma->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * solidWater = new G4Tubs("Phantom_Inner", 0, 0.5*30.6, 0.5*100.0, 0, 360.0*deg);
    G4LogicalVolume   * logicWater = new G4LogicalVolume(solidWater, matWater, "Phantom_inner");
    new G4PVPlacement(nullptr, {0, 0, -0.8}, logicWater, "Phantom_Inner_PV", logicPmma, false, 0);
    logicWater->SetVisAttributes(G4VisAttributes(G4Colour(0, 1.0, 0)));

    return logicPmma;
}

G4LogicalVolume * PhantomMarekCompartments::definePhantom(G4LogicalVolume *logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();

    G4NistManager * man = G4NistManager::Instance();

    G4Material * matPMMA  = man->FindOrBuildMaterial("G4_PLEXIGLASS");
    G4Material * matWater = man->FindOrBuildMaterial("G4_WATER");
    G4Material * matTi = man->FindOrBuildMaterial("G4_Ti");

    double outerD = 37.77;
    double innerD = 30.56;
    double wallOuterD = 33.0;

    double ringT  = 6.4;
    double centerRingT = 0.8;
    double wallT  = 0.00254;
    double frontT = 2.5;
    double backT  = 76.2;

    G4VSolid          * sRing = new G4Tubs("Phantom_sRing", 0, 0.5*outerD, 0.5*ringT, 0, 360.0*deg);
    G4LogicalVolume   * lRing = new G4LogicalVolume(sRing, matPMMA, "Phantom_lRing");
    lRing->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * sInner = new G4Tubs("Phantom_sRing", 0, 0.5*innerD, 0.5*6.4, 0, 360.0*deg);
    G4LogicalVolume   * lInner = new G4LogicalVolume(sInner, matWater, "Phantom_lInner");
    new G4PVPlacement(nullptr, {0, 0, 0}, lInner, "Phantom_Inner_PV", lRing, false, 0);
    lInner->SetVisAttributes(G4VisAttributes(G4Colour(0, 0, 1.0)));

    G4VSolid          * sWall = new G4Tubs("Phantom_sWall", 0, 0.5*wallOuterD, 0.5*wallT, 0, 360.0*deg);
    G4LogicalVolume   * lWall = new G4LogicalVolume(sWall, matTi, "Phantom_lWall");
    lWall->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0, 0)));

    G4VSolid          * sFront = new G4Tubs("Phantom_sFront", 0.5*innerD, 0.5*outerD, 0.5*frontT, 0, 360.0*deg);
    G4LogicalVolume   * lFront = new G4LogicalVolume(sFront, matPMMA, "Phantom_lFront");
    lFront->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * sBack = new G4Tubs("Phantom_sBack", 0.5*innerD, 0.5*outerD, 0.5*backT, 0, 360.0*deg);
    G4LogicalVolume   * lBack = new G4LogicalVolume(sBack, matPMMA, "Phantom_lBack");
    lBack->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * sCenterRing = new G4Tubs("Phantom_sCenterRing", 0, 0.5*outerD, 0.5*centerRingT, 0, 360.0*deg);
    G4LogicalVolume   * lCenterRing = new G4LogicalVolume(sCenterRing, matPMMA, "Phantom_lCenterRing");
    lCenterRing->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * sCenterInner = new G4Tubs("Phantom_sCenterRing", 0, 0.5*innerD, 0.5*centerRingT, 0, 360.0*deg);
    G4LogicalVolume   * lCenterInner = new G4LogicalVolume(sCenterInner, matWater, "Phantom_lCenterInner");
    new G4PVPlacement(nullptr, {0, 0, 0}, lCenterInner, "Phantom_CenterInner_PV", lCenterRing, false, 0);
    lInner->SetVisAttributes(G4VisAttributes(G4Colour(0, 0, 1.0)));

    double pos = SM.GlobalZ0;
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, pos}, lCenterRing, "Phantom_CenterRing_PV", logicWorld, false, 0);
    pos -= 0.5*centerRingT;

    for (int i = 0; i < 4; i++)
    {
        pos -= 0.5*wallT;
        new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, pos}, lWall, "Phantom_Wall_PV", logicWorld, true, 4-i);
        pos -= 0.5*wallT;

        pos -= 0.5*ringT;
        new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, pos}, lRing, "Phantom_Ring_PV", logicWorld, true, 4-i);
        pos -= 0.5*ringT;
    }

    pos -= 0.5*wallT;
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, pos}, lWall, "Phantom_Wall_PV", logicWorld, true, 0);
    pos -= 0.5*wallT;

    pos -= 0.5*frontT;
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, pos}, lFront, "Phantom_Front_PV", logicWorld, false, 0);
    //pos -= 0.5*frontT;

    pos = SM.GlobalZ0 + 0.5*centerRingT;
    for (int i = 0; i < 4; i++)
    {
        pos += 0.5*wallT;
        new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, pos}, lWall, "Phantom_Wall_PV", logicWorld, true, 5+i);
        pos += 0.5*wallT;

        pos += 0.5*ringT;
        new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, pos}, lRing, "Phantom_Ring_PV", logicWorld, true, 5+i);
        pos += 0.5*ringT;
    }

    pos += 0.5*wallT;
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, pos}, lWall, "Phantom_Wall_PV", logicWorld, true, 9);
    pos += 0.5*wallT;

    pos += 0.5*backT;
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, pos}, lBack, "Phantom_Back_PV", logicWorld, false, 0);
    pos += 0.5*backT;

    return lRing;
}

G4LogicalVolume * PhantomBeamDerenzoAndr::definePhantom(G4LogicalVolume * logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();

    G4NistManager * man = G4NistManager::Instance();

    G4Material * matAir  = man->FindOrBuildMaterial("G4_AIR");
    G4Material * matPMMA  = man->FindOrBuildMaterial("G4_PLEXIGLASS");

    double length = 40.0*mm;

    G4VSolid          * sBox = new G4Box("Phantom_sBox", 30.0, 30.0, 30.0);
    G4LogicalVolume   * lBox = new G4LogicalVolume(sBox, matAir, "Phantom_lBox");
    lBox->SetVisAttributes(G4VisAttributes(G4Colour(0, 0, 1.0)));
    new G4PVPlacement(nullptr, {0, 0, SM.GlobalZ0}, lBox, "Phantom_pBox", logicWorld, false, 0);

    G4VSolid          * s2 = new G4Tubs("Phantom_s2", 0, 0.5*2.0*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l2 = new G4LogicalVolume(s2, matPMMA, "Phantom_l2");
    l2->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * s3 = new G4Tubs("Phantom_s3", 0, 0.5*3.0*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l3 = new G4LogicalVolume(s3, matPMMA, "Phantom_l3");
    l3->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * s4 = new G4Tubs("Phantom_s4", 0, 0.5*4.0*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l4 = new G4LogicalVolume(s4, matPMMA, "Phantom_l4");
    l4->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * s5 = new G4Tubs("Phantom_s5", 0, 0.5*5.0*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l5 = new G4LogicalVolume(s5, matPMMA, "Phantom_l5");
    l5->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * s6 = new G4Tubs("Phantom_s6", 0, 0.5*6.0*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l6 = new G4LogicalVolume(s6, matPMMA, "Phantom_l6");
    l6->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    for (int i = 0; i < 8; i++) new G4PVPlacement(nullptr, {(i-3-0.5) * 4.0*mm,   0.0,  0}, l2, "Phantom_p2", lBox, true, i);
    for (int i = 0; i < 5; i++) new G4PVPlacement(nullptr, {(i-2)     * 6.0*mm,   5.5,  0}, l3, "Phantom_p3", lBox, true, i);
    for (int i = 0; i < 4; i++) new G4PVPlacement(nullptr, {(i-1-0.5) * 8.0*mm,  -6.0,  0}, l4, "Phantom_p4", lBox, true, i);
    for (int i = 0; i < 3; i++) new G4PVPlacement(nullptr, {(i-1)     * 10.0*mm, -14.0, 0}, l5, "Phantom_p5", lBox, true, i);
    for (int i = 0; i < 2; i++) new G4PVPlacement(nullptr, {(i-0.5)   * 12.0*mm,  14.0, 0}, l6, "Phantom_p6", lBox, true, i);

    return lBox;
}

G4LogicalVolume *PhantomBeamDerenzoAndr2::definePhantom(G4LogicalVolume *logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();

    G4NistManager * man = G4NistManager::Instance();

    G4Material * matAir  = man->FindOrBuildMaterial("G4_AIR");
    G4Material * matPMMA  = man->FindOrBuildMaterial("G4_PLEXIGLASS");

    double length = 40.0*mm;
    double halfLengthBox = 0.5*length + 1.0*mm;

    G4VSolid          * sBox = new G4Box("Phantom_sBox", halfLengthBox, halfLengthBox, halfLengthBox);
    G4LogicalVolume   * lBox = new G4LogicalVolume(sBox, matAir, "Phantom_lBox");
    lBox->SetVisAttributes(G4VisAttributes(G4Colour(0, 0, 1.0)));
    new G4PVPlacement(nullptr, {0, 0, SM.GlobalZ0}, lBox, "Phantom_pBox", logicWorld, false, 0);

    G4VSolid          * s2 = new G4Tubs("Phantom_s2", 0, 0.5*2.0*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l2 = new G4LogicalVolume(s2, matPMMA, "Phantom_l2");
    l2->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * s2i5 = new G4Tubs("Phantom_s2i5", 0, 0.5*2.5*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l2i5 = new G4LogicalVolume(s2i5, matPMMA, "Phantom_l2i5");
    l2i5->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * s3 = new G4Tubs("Phantom_s3", 0, 0.5*3.0*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l3 = new G4LogicalVolume(s3, matPMMA, "Phantom_l3");
    l3->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * s3i5 = new G4Tubs("Phantom_s3i5", 0, 0.5*3.5*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l3i5 = new G4LogicalVolume(s3i5, matPMMA, "Phantom_l3i5");
    l3i5->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    placeArray(l2,   4.0, -6.0, -5.0, lBox, "Ph_2");
    placeArray(l2i5, 5.0, -6.0,  5.0, lBox, "Ph_2i5");
    placeArray(l3,   6.0,  6.0,  7.5, lBox, "Ph_3");
    placeArray(l3i5, 7.0,  7.0, -6.5, lBox, "Ph_3i5");

    return lBox;
}

void PhantomBeamDerenzoAndr2::placeArray(G4LogicalVolume * element, double step, double x, double y, G4LogicalVolume * container, const std::string & name)
{
    int iCounter = 0;
    for (int ix = 0; ix < 2; ix++)
        for (int iy = 0; iy < 2; iy++)
            new G4PVPlacement(nullptr, {x*mm + (ix - 0.5) * step*mm, y*mm + (iy - 0.5) * step*mm,  0}, element, name, container, true, iCounter++);
}

G4LogicalVolume *PhantomBeamDerenzoAndr2inv::definePhantom(G4LogicalVolume *logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();

    G4NistManager * man = G4NistManager::Instance();

    G4Material * matAir  = man->FindOrBuildMaterial("G4_AIR");
    G4Material * matPMMA  = man->FindOrBuildMaterial("G4_PLEXIGLASS");

    double length = 40.0*mm;
    double cylDiameter = 40.0 * mm;
    double cylLength   = 60.0 * mm;

    G4VSolid          * sCyl = new G4Tubs("sCyl", 0, 0.5*cylDiameter, 0.5*cylLength, 0, 360.0*deg);
    G4LogicalVolume   * lCyl = new G4LogicalVolume(sCyl, matPMMA, "lCyl");
    lCyl->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, SM.GlobalZ0}, lCyl, "lCyl", logicWorld, false, 0);

    G4VSolid          * s2 = new G4Tubs("Phantom_s2", 0, 0.5*2.0*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l2 = new G4LogicalVolume(s2, matAir, "Phantom_l2");
    l2->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * s2i5 = new G4Tubs("Phantom_s2i5", 0, 0.5*2.5*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l2i5 = new G4LogicalVolume(s2i5, matAir, "Phantom_l2i5");
    l2i5->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * s3 = new G4Tubs("Phantom_s3", 0, 0.5*3.0*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l3 = new G4LogicalVolume(s3, matAir, "Phantom_l3");
    l3->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    G4VSolid          * s3i5 = new G4Tubs("Phantom_s3i5", 0, 0.5*3.5*mm, 0.5*length, 0, 360.0*deg);
    G4LogicalVolume   * l3i5 = new G4LogicalVolume(s3i5, matAir, "Phantom_l3i5");
    l3i5->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    placeArray(l2,   4.0, -6.0, -5.0, lCyl, "Ph_2");
    placeArray(l2i5, 5.0, -6.0,  5.0, lCyl, "Ph_2i5");
    placeArray(l3,   6.0,  6.0,  7.5, lCyl, "Ph_3");
    placeArray(l3i5, 7.0,  7.0, -6.5, lCyl, "Ph_3i5");

    return lCyl;
}
