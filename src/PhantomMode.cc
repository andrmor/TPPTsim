#include "PhantomMode.hh"
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
#include "out.hh"
#include "jstools.hh"

#include <math.h>

PhantomModeBase * PhantomModeFactory::makePhantomModeInstance(const json11::Json & json)
{
    out("Reading phantom json");
    std::string Type;
    jstools::readString(json, "Type", Type);

    PhantomModeBase * ph = nullptr;

    if      (Type == "PhantomNone")     ph = new PhantomNone();
    else if (Type == "PhantomPMMA")     ph = new PhantomPMMA();
    else if (Type == "PhantomTinyCube") ph = new PhantomTinyCube();
    else if (Type == "PhantomDerenzo")  ph = new PhantomDerenzo(100.0, 100.0, {}, 0, 0, 0);
    else if (Type == "PhantomParam")    ph = new PhantomParam();
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
    logicPmma->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 1.0)));

    return logicPmma;
}

// ---

G4LogicalVolume * PhantomTinyCube::definePhantom(G4LogicalVolume *logicWorld)
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
