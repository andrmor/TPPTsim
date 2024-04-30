/** ----------------------
  Copyright (C): OpenGATE Collaboration
  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See LICENSE.md for further details
  ----------------------*/
#ifndef PositroniumDecayChannel_hh
#define PositroniumDecayChannel_hh

#include "globals.hh"
#include "G4GeneralPhaseSpaceDecay.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//original: GatePositroniumDecayChannel.hh

/** Author: Mateusz Bała
 *  Email: bala.mateusz@gmail.com
 *  Theorem author for oPs decay: Daria Kamińska ( Eur. Phys. J. C (2016) 76:445 )
 *  Organization: J-PET (http://koza.if.uj.edu.pl/pet/)
 *  About class: Implements decay of positronium ( pPs and oPs ). Provides support for polarization.
 *  Andrey Morozov, 2023 -> refactor and simplification for implementation in the TPPTsim framework
 **/

class PositroniumDecayChannel : public G4GeneralPhaseSpaceDecay
{
public:
    PositroniumDecayChannel(bool ortho);

    G4DecayProducts * DecayIt(G4double) override; // Return gammas from positronium decay

protected:
    G4DecayProducts * DecayParaPositronium();
    G4DecayProducts * DecayOrthoPositronium();

    // Calculate cross section Mij matrix element
    //Based on "Quantum electrodynamics" V. B. BERESTETSKY.
    //Chapter: 89. Annihilation of positronium
    // Exquantation: 89.14
    G4double GetOrthoPsM(G4double w1, G4double w2, G4double w3) const;

    G4ThreeVector GetPolarization(const G4ThreeVector & momentum) const;
    G4ThreeVector GetPerpendicularVector(const G4ThreeVector & v) const;

protected:
    bool Ortho = true;

    const G4double kPositroniumMass = 2.0 * electron_mass_c2;

    //This is maxiaml number which can be calculated by function GetOrthoPsM() - based on 10^7 iterations
    const G4double kOrthoPsMMax = 7.65928;

    const G4double kElectronMass = electron_mass_c2; //[MeV]
};

#endif //PositroniumDecayChannel_hh
