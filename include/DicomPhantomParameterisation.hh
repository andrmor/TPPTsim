//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file medical/DICOM/include/DicomPhantomParameterisationColour.hh
/// \brief Definition of the DicomPhantomParameterisationColour class
//

#ifndef DicomPhantomParameterisation_HH
#define DicomPhantomParameterisation_HH

#include <map>

#include "G4PhantomParameterisation.hh"
#include "G4ThreeVector.hh"

class G4VisAttributes;
class DicomPhantomZSliceHeader;


class DicomPhantomParameterisation : public G4VPVParameterisation //G4PhantomParameterisation
{
public:
    typedef std::map<G4String,G4VisAttributes*> ColourMap_t;

    //static G4String defaultColorFile;


public:  // with description
    DicomPhantomParameterisation(std::vector<std::pair<double,double>> coord2D, double zStart);

    ~DicomPhantomParameterisation();

    virtual G4Material* ComputeMaterial(const G4int repNo,
                                        G4VPhysicalVolume *currentVol,
                                        const G4VTouchable *parentTouch=0);

    void ComputeTransformation (const int copyNo, G4VPhysicalVolume * physVol) const;

    const ColourMap_t& GetColourMap() const { return fColours; }
    ColourMap_t& GetColourMap() { return fColours; }

    // ANDR
    void SetVoxelDimensions(double DX, double DY, double DZ) {HalfVoxelX = DX; HalfVoxelY = DY; HalfVoxelZ = DZ;}
    void SetAir(G4Material * mat) {air = mat;}
    void SetAir1(G4Material * mat) {air1 = mat;}
    // ----

protected:
    std::vector<std::pair<double,double>> XY;
    double ZStart;

private:
    void ReadColourData(G4String colourFile);

private:
    ColourMap_t fColours;
    std::map<G4int, G4VisAttributes*> mColours;

    DicomPhantomZSliceHeader* ZSliceHeader;
    //G4int control;

    //ANDR
    double HalfVoxelX;
    double HalfVoxelY;
    double HalfVoxelZ;
    G4Material * air;
    G4Material * air1;
    //----

};

#endif
