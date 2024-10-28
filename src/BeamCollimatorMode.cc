#include "BeamCollimatorMode.hh"
#include "SessionManager.hh"
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

BeamCollimatorBase * BeamCollimatorFactory::makeBeamCollimatorInstance(const json11::Json & json)
{
    out("Reading beam collimator json");
    if (json.is_null()) return nullptr;

    std::string Type;
    jstools::readString(json, "Type", Type);
    if (Type.empty()) return nullptr;

    BeamCollimatorBase * bc = nullptr;

    if      (Type == "BeamCollimatorNone")  bc = new BeamCollimatorNone();
    else if (Type == "BeamCollimatorMarek") bc = new BeamCollimatorMarek(BeamCollimatorMarek::Blind, {0,0,-75*mm});
    else
    {
        out("Unknown beam collimator type!");
        exit(10);
    }

    if (bc) bc->readFromJson(json);
    return bc;
}

// --------

void BeamCollimatorBase::writeToJson(json11::Json::object & json) const
{
    json["Type"] = getTypeName();
    doWriteToJson(json);
}

// --------

BeamCollimatorMarek::BeamCollimatorMarek(EOpeningOptions openingType, const G4ThreeVector & innerFacePosition, double rotationAngle) :
    OpeningType(openingType), FacePosition(innerFacePosition), Angle(rotationAngle) {}


void BeamCollimatorMarek::defineBeamCollimator(G4LogicalVolume * logicWorld)
{
    SessionManager & SM = SessionManager::getInstance();

    double cylLength    = 30.0*mm;
    double cylDiameter  = 30.0*mm;
    double ringLength   = 10.0*mm;
    double ringDiameter = 45.0*mm;

    G4NistManager * man = G4NistManager::Instance();
    G4Material * matBrass = man->FindOrBuildMaterial("G4_BRASS");

    G4Tubs          * sCyl = new G4Tubs("sCyl", 0, 0.5*cylDiameter,0.5*cylLength, 0, 360.0*deg);
    G4LogicalVolume * lCyl = new G4LogicalVolume(sCyl, matBrass, "lCyl");
    lCyl->SetVisAttributes(G4VisAttributes({1.0, 0, 0}));

    G4Tubs          * sRing = new G4Tubs("sRing", 0.5*cylDiameter, 0.5*ringDiameter, 0.5*ringLength, 0, 360.0*deg);
    G4LogicalVolume * lRing = new G4LogicalVolume(sRing, matBrass, "lRing");
    lRing->SetVisAttributes(G4VisAttributes({1.0, 0, 0}));

    switch (OpeningType)
    {
    case Hole3 :
    {
        G4Tubs          * sChan = new G4Tubs("sChan", 0, 0.5*3.0*mm, 0.5*cylLength, 0, 360.0*deg);
        G4LogicalVolume * lChan = new G4LogicalVolume(sChan, SM.WorldMat, "lChan");
        lChan->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));
        new G4PVPlacement(nullptr, {0,0,0}, lChan,  "pChan",  lCyl, false, 0);
        break;
    }
    case Holes369 :
    {
        G4Tubs          * sChan3 = new G4Tubs("sChan3", 0, 0.5*3.0*mm, 0.5*cylLength, 0, 360.0*deg);
        G4LogicalVolume * lChan3 = new G4LogicalVolume(sChan3, SM.WorldMat, "lChan3");
        lChan3->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));
        new G4PVPlacement(nullptr, {10.0*mm,0,0}, lChan3,  "pChan3",  lCyl, false, 0);

        G4Tubs          * sChan6 = new G4Tubs("sChan6", 0, 0.5*6.0*mm, 0.5*cylLength, 0, 360.0*deg);
        G4LogicalVolume * lChan6 = new G4LogicalVolume(sChan6, SM.WorldMat, "lChan6");
        lChan6->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));
        new G4PVPlacement(nullptr, {0,0,0}, lChan6, "pChan6",  lCyl, false, 0);

        G4Tubs          * sChan9 = new G4Tubs("sChan9", 0, 0.5*9.0*mm, 0.5*cylLength, 0, 360.0*deg);
        G4LogicalVolume * lChan9 = new G4LogicalVolume(sChan9, SM.WorldMat, "lChan9");
        lChan9->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));
        new G4PVPlacement(nullptr, {-10.0*mm,0,0}, lChan9,  "pChan9",  lCyl, false, 0);

        break;
    }
    case Holes123 :
    {
        G4Tubs          * sChan1 = new G4Tubs("sChan1", 0, 0.5*1.0*mm, 0.5*cylLength, 0, 360.0*deg);
        G4LogicalVolume * lChan1 = new G4LogicalVolume(sChan1, SM.WorldMat, "lChan1");
        lChan1->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));
        new G4PVPlacement(nullptr, {10.0*mm,0,0}, lChan1,  "pChan1",  lCyl, false, 0);

        G4Tubs          * sChan2 = new G4Tubs("sChan2", 0, 0.5*2.0*mm, 0.5*cylLength, 0, 360.0*deg);
        G4LogicalVolume * lChan2 = new G4LogicalVolume(sChan2, SM.WorldMat, "lChan2");
        lChan2->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));
        new G4PVPlacement(nullptr, {0,0,0}, lChan2, "pChan2",  lCyl, false, 0);

        G4Tubs          * sChan3 = new G4Tubs("sChan3", 0, 0.5*3.0*mm, 0.5*cylLength, 0, 360.0*deg);
        G4LogicalVolume * lChan3 = new G4LogicalVolume(sChan3, SM.WorldMat, "lChan3");
        lChan3->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));
        new G4PVPlacement(nullptr, {-10.0*mm,0,0}, lChan3,  "pChan3",  lCyl, false, 0);

        break;
    }
    case Cross :
    {
        double boxX = 10.0 - 0.5 * 3.0*mm;
        double boxOffset = 0.5 * 3.0*mm + 0.5 * boxX;

        G4Box           * sCenter = new G4Box("sCenter", 0.5*3.0*mm, 0.5*3.0*mm, 0.5*cylLength);
        G4LogicalVolume * lCenter = new G4LogicalVolume(sCenter, SM.WorldMat, "lCenter");
        lCenter->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));

        G4Tubs          * sDom = new G4Tubs("sDom", 0, 0.5*3.0*mm, 0.5*cylLength, -90.0*deg, 180.0*deg);
        G4LogicalVolume * lDom = new G4LogicalVolume(sDom, SM.WorldMat, "lDom");
        lDom->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));

        G4Box           * sBox = new G4Box("sBox", 0.5*boxX, 0.5*3.0*mm, 0.5*cylLength);
        G4LogicalVolume * lBox = new G4LogicalVolume(sBox, SM.WorldMat, "lBox");
        lCenter->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));

        new G4PVPlacement(nullptr, {0,0,0}, lCenter, "pCenter", lCyl, false, 0);

        new G4PVPlacement(nullptr, {boxOffset,0,0}, lBox,  "pBox",  lCyl, true, 0);
        new G4PVPlacement(nullptr, {10.0*mm,0,0},    lDom,  "pDom",  lCyl, true, 0);

        new G4PVPlacement(new CLHEP::HepRotation(90*deg, 0, 0), {0, boxOffset,0}, lBox, "pBox", lCyl, true, 1);
        new G4PVPlacement(new CLHEP::HepRotation(90*deg, 0, 0), {0, 10.0*mm,0}, lDom,  "pDom",  lCyl, true, 1);

        new G4PVPlacement(new CLHEP::HepRotation(180*deg, 0, 0), {-boxOffset,0,0}, lBox,  "pBox",  lCyl, true, 2);
        new G4PVPlacement(new CLHEP::HepRotation(180*deg, 0, 0), {-10.0*mm,0,0}, lDom,  "pDom",  lCyl, true, 2);

        new G4PVPlacement(new CLHEP::HepRotation(-90*deg, 0, 0), {0, -boxOffset,0}, lBox, "pBox", lCyl, true, 3);
        new G4PVPlacement(new CLHEP::HepRotation(-90*deg, 0, 0), {0, -10.0*mm,0}, lDom,  "pDom",  lCyl, true, 3);

        break;
    }
    case Ring :
    {
        double cutAngle = 3.0 / (2*3.1415926*12.5) * 360.0;
        double arcAngle = 90.0 - cutAngle;
        out("Cut angle:", cutAngle, "Arc angle:", arcAngle, "deg");
        G4Tubs          * sChan = new G4Tubs("sChan", 0.5*10.0, 0.5*12.5*mm, 0.5*cylLength, 0, arcAngle*deg);
        G4LogicalVolume * lChan = new G4LogicalVolume(sChan, SM.WorldMat, "lChan");
        lChan->SetVisAttributes(G4VisAttributes({1.0, 1.0, 1.0}));
        new G4PVPlacement(new CLHEP::HepRotation(   0*deg + 0.5*cutAngle*deg, 0, 0), {0,0,0}, lChan,  "pChan0",  lCyl, true, 0);
        new G4PVPlacement(new CLHEP::HepRotation( -90*deg + 0.5*cutAngle*deg, 0, 0), {0,0,0}, lChan,  "pChan1",  lCyl, true, 1);
        new G4PVPlacement(new CLHEP::HepRotation( +90*deg + 0.5*cutAngle*deg, 0, 0), {0,0,0}, lChan,  "pChan2",  lCyl, true, 2);
        new G4PVPlacement(new CLHEP::HepRotation(+180*deg + 0.5*cutAngle*deg, 0, 0), {0,0,0}, lChan,  "pChan3",  lCyl, true, 3);

    }
    default :;
    case Blind :; // no channels
    }

    double rotAngle = 90.0*deg + Angle;
    double zPosCyl = -0.5*cylLength + FacePosition[2] + SM.GlobalZ0;
    new G4PVPlacement(new CLHEP::HepRotation(rotAngle, 0, 0), {FacePosition[0], FacePosition[1], zPosCyl},  lCyl,  "pCyl",  logicWorld, false, 0);
    double zPosRing = -0.5*ringLength + FacePosition[2] + SM.GlobalZ0;
    new G4PVPlacement(new CLHEP::HepRotation(rotAngle, 0, 0), {FacePosition[0], FacePosition[1], zPosRing}, lRing, "pRing", logicWorld, false, 0);
}

void BeamCollimatorMarek::readFromJson(const json11::Json & json)
{
    std::string sType;
    jstools::readString(json, "OpeningType", sType);

    if      (sType == "Blind")    OpeningType = Blind;
    else if (sType == "Hole3")    OpeningType = Hole3;
    else if (sType == "Holes369") OpeningType = Holes369;
    else if (sType == "Holes123") OpeningType = Holes123;
    else if (sType == "Cross")    OpeningType = Cross;
    else if (sType == "Ring")     OpeningType = Ring;
    else
    {
        out("Unknown type of Marek's collimator opening in BeamCollimatorMarek::readFromJson");
        exit(2345);
    }

    json11::Json::array ar;
    jstools::readArray(json, "FacePosition", ar);
    if (ar.size() >= 3)
        for (size_t i = 0; i < 3; i++)
            FacePosition[i] = ar[i].number_value();
}

void BeamCollimatorMarek::doWriteToJson(json11::Json::object & json) const
{
    std::string sType;
    switch (OpeningType)
    {
    case Hole3:    sType = "Hole3"; break;
    case Holes369: sType = "Holes369"; break;
    case Holes123: sType = "Holes123"; break;
    case Cross:    sType = "Cross"; break;
    case Ring:     sType = "Ring"; break;
    case Blind:    sType = "Blind"; break;
    default:
        out("Unknown type of Marek's collimator opening in BeamCollimatorMarek::doWriteToJson");
        exit(2345);
    }

    json["OpeningType"] = sType;

    json11::Json::array ar;
    for (size_t i = 0; i < 3; i++)
        ar.push_back(FacePosition[i]);
    json["FacePosition"] = ar;
}


