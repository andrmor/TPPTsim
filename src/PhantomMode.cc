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

#include <vector>

void PhantomModePMMA::definePhantom(G4LogicalVolume * logicWorld)
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
}
