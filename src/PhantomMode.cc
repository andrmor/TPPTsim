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

    if      (Type == "PhantomNone")      ph = new PhantomNone();
    else if (Type == "PhantomPMMA")      ph = new PhantomPMMA();
    else if (Type == "PhantomTinyCube")  ph = new PhantomTinyCube();
    else if (Type == "PhantomCustomBox") ph = new PhantomCustomBox();
    else if (Type == "PhantomDerenzo")   ph = new PhantomDerenzo(100.0, 100.0, {}, 0, 0, 0);
    else if (Type == "PhantomParam")     ph = new PhantomParam();
    else if (Type == "PhantomModeDICOM") ph = new PhantomModeDICOM(100.0, {0,0,50.0}, "DummyFileName.dat", true);
    else if (Type == "PhantomEspana")    ph = new PhantomEspana();
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

G4LogicalVolume * PhantomPMMA::definePhantom(G4LogicalVolume * logicWorld)
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

    return logicPmma;
}

// ---

PhantomCustomBox::PhantomCustomBox(double sizeX, double sizeY, double sizeZ, EMaterial material) :
    SizeX(sizeX), SizeY(sizeY), SizeZ(sizeZ), Material(material) {}

G4LogicalVolume * PhantomCustomBox::definePhantom(G4LogicalVolume * logicWorld)
{
    G4NistManager * man = G4NistManager::Instance();
    G4Material * mat = nullptr;

    switch (Material)
    {
    case PMMA :
        {
            std::vector<G4int> natoms;
            std::vector<G4String> elements;
            elements.push_back("C"); natoms.push_back(5);
            elements.push_back("H"); natoms.push_back(8);
            elements.push_back("O"); natoms.push_back(2);
            mat = man->ConstructNewMaterial("PMMA_phantom", elements, natoms, 1.18*g/cm3);
        }
        break;
    case HDPE :
        {
            std::vector<double> weightFrac;
            std::vector<G4String> elements;
            elements.push_back("H"); weightFrac.push_back(14.3);
            elements.push_back("C"); weightFrac.push_back(85.7);
            mat = man->ConstructNewMaterial("HDPE", elements, weightFrac, 0.95*g/cm3);
        }
        break;
    case PE :
        {
            std::vector<G4int> natoms;
            std::vector<G4String> elements;
            elements.push_back("C"); natoms.push_back(2);
            elements.push_back("H"); natoms.push_back(4);
            mat = man->ConstructNewMaterial("PE", elements, natoms, 0.96*g/cm3);
        }
        break;
    case Graphite :
        {
            std::vector<G4int> natoms;
            std::vector<G4String> elements;
            elements.push_back("C"); natoms.push_back(1);
            mat = man->ConstructNewMaterial("Graphite", elements, natoms, 1.83*g/cm3);
        }
        break;
    case GelTissue :
        {
            std::vector<double> weightFrac;
            std::vector<G4String> elements;
            elements.push_back("H"); weightFrac.push_back(9.6);
            elements.push_back("C"); weightFrac.push_back(14.9);
            elements.push_back("N"); weightFrac.push_back(1.46);
            elements.push_back("O"); weightFrac.push_back(73.8);
            mat = man->ConstructNewMaterial("GelTissue", elements, weightFrac, 1.13*g/cm3);
        }
        break;
    case GelWater :
        {
            std::vector<double> weightFrac;
            std::vector<G4String> elements;
            elements.push_back("H"); weightFrac.push_back(11.03);
            elements.push_back("C"); weightFrac.push_back(1.04);
            elements.push_back("N"); weightFrac.push_back(0.32);
            elements.push_back("O"); weightFrac.push_back(87.6);
            mat = man->ConstructNewMaterial("GelWater", elements, weightFrac, 1.01*g/cm3);
        }
        break;
    default:;
    }
    if (!mat)
    {
        out("Error in material selection of PhantomCustomMatBox");
        exit(10);
    }
    out("Ionization potential for the phantom material:", mat->GetIonisation()->GetMeanExcitationEnergy()/eV, "eV");

    G4VSolid          * solid = new G4Box("Phantom_Box", 0.5 * SizeX * mm, 0.5 * SizeY * mm, 0.5 * SizeZ * mm);
    G4LogicalVolume   * logic = new G4LogicalVolume(solid, mat, "Phantom");
    new G4PVPlacement(nullptr, {0, 0, 0}, logic, "Phantom_PV", logicWorld, false, 0);
    logic->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0)));

    return logic;
}

void PhantomCustomBox::readFromJson(const json11::Json & json)
{
    jstools::readDouble(json, "SizeX", SizeX);
    jstools::readDouble(json, "SizeY", SizeY);
    jstools::readDouble(json, "SizeZ", SizeZ);

    std::string MatStr;
    jstools::readString(json, "Material", MatStr);

    if      (MatStr == "PMMA")      Material = PMMA;
    else if (MatStr == "HDPE")      Material = HDPE;
    else if (MatStr == "PE")        Material = PE;
    else if (MatStr == "Graphite")  Material = Graphite;
    else if (MatStr == "GelTissue") Material = GelTissue;
    else if (MatStr == "GelWater")  Material = GelWater;
    else
    {
        out("Unknown material for phantom PhantomCustomMat");
        exit(1);
    }
}

void PhantomCustomBox::doWriteToJson(json11::Json::object & json) const
{
    json["SizeX"] = SizeX;
    json["SizeY"] = SizeY;
    json["SizeZ"] = SizeZ;

    std::string MatStr;
    switch (Material)
    {
    case PMMA      : MatStr = "PMMA";      break;
    case HDPE      : MatStr = "HDPE";      break;
    case PE        : MatStr = "PE";        break;
    case Graphite  : MatStr = "Graphite";  break;
    case GelTissue : MatStr = "GelTissue"; break;
    case GelWater  : MatStr = "GelWater";  break;
    default:;
    }

    json["Material"] = MatStr;
}

// ---

G4LogicalVolume * PhantomTinyCube::definePhantom(G4LogicalVolume * logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();
    G4NistManager * man = G4NistManager::Instance();

    std::vector<G4int> natoms;
    std::vector<G4String> elements;
    elements.push_back("C"); natoms.push_back(5);
    elements.push_back("H"); natoms.push_back(8);
    elements.push_back("O"); natoms.push_back(2);
    G4Material * matPMMA = man->ConstructNewMaterial("PMMA_phantom", elements, natoms, 1.18*g/cm3);

    G4VSolid          * solidPmma = new G4Box("Box", 3.0*mm, 3.0*mm, 3.0*mm);
    G4LogicalVolume   * logicPmma = new G4LogicalVolume(solidPmma, matPMMA, "Phantom");
    new G4PVPlacement(nullptr, {0, 0, SM.GlobalZ0}, logicPmma, "Phantom_PV", logicWorld, false, 0);
    logicPmma->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 1.0)));

    return logicPmma;
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
    G4Material * matPMMA = man->ConstructNewMaterial("PMMA_phantom", elements, natoms, 1.18*g/cm3);

    G4Material * matWater = man->FindOrBuildMaterial("G4_WATER");

    G4VSolid          * solidPmma = new G4Tubs("Phantom_Cyl", 0, 0.5*Diameter*mm, 0.5*Height*mm, 0, 360.0*deg);
    G4LogicalVolume   * logicPmma = new G4LogicalVolume(solidPmma, matPMMA, "Phantom");
    new G4PVPlacement(new CLHEP::HepRotation(90.0*deg, 0, 0), {0, 0, SM.GlobalZ0}, logicPmma, "Phantom_PV", logicWorld, false, 0);
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
            G4LogicalVolume * logicHole = new G4LogicalVolume(solidHole, matWater, "Hole" + std::to_string(iSec));
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
