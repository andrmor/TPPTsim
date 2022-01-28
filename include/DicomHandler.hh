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
/// \file medical/DICOM/include/DicomHandler.hh
/// \brief Definition of the DicomHandler class
//
// The code was written by :
//      *Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268

//*******************************************************//

#ifndef DicomHandler_h
#define DicomHandler_h 1

#include <cstdio>
#include <map>
#include <fstream>
#include <vector>

#include "globals.hh"

//*******************************************************
/// Dicom Handler class
///        - Handling of DICM images
///        - Transforming *.dcm to *.g4 ( pixels->density )
///        - Reading headers and pixels
///        - Transforming pixel to density and creating *.g4
///          files
///        - Functions are in DicomHandler.cc
///
/// Base on previous code by :
///        Dragan Tubic <tdragan@gel.ulaval.ca>
//*******************************************************

class DicomPhantomZSliceHeader;
class DicomPhantomZSliceMerged;

class DicomHandler
{
public:
    static DicomHandler & getInstance();

private:
    DicomHandler();
    ~DicomHandler();

public:
    void configure(const G4String & path, const G4String & convertionFileName, int lateralCompression,
                   const std::vector<std::pair<std::string, float> > & materialUpperDens,
                   const std::vector<std::string> & sliceFiles);

    void setDriver(const G4String & path, const G4String & driverFileName, const G4String & convertionFileName);
    void CheckFileFormat();

private:
    template <class Type> void GetValue(char *, Type &);
    
    const G4int DATABUFFSIZE;
    const G4int LINEBUFFSIZE;
    const G4int FILENAMESIZE;
    
    G4int ReadFile(FILE *, const char *);
    G4int ReadData(FILE *);
    void ReadCalibration();
    void GetInformation(G4int &, char *);
    G4float Pixel2density(G4int pixel);
    void ReadMaterialIndices( std::ifstream& finData);
    unsigned int GetMaterialIndex( G4float density );
    void StoreData(DicomPhantomZSliceHeader* dcmPZSH);
    G4int read_defined_nested(FILE *, G4int);
    void read_undefined_nested(FILE *);
    void read_undefined_item(FILE *);

    short fCompression;
    G4int fNFiles;
    short fRows;
    short fColumns;
    short fBitAllocated;
    G4int fMaxPixelValue, fMinPixelValue;
    
    G4double fPixelSpacingX, fPixelSpacingY;
    G4double fSliceThickness;
    G4double fSliceLocation;
    
    G4int fRescaleIntercept, fRescaleSlope;
    
    G4bool fLittleEndian, fImplicitEndian;
    short fPixelRepresentation;
    
    G4int** fTab;
    std::map<G4float,G4String> fMaterialIndices;
    
    G4int fNbrequali;
    G4double * fValueDensity = nullptr;
    G4double * fValueCT = nullptr;
    G4bool fReadCalibration;
    DicomPhantomZSliceMerged * fMergedSlices = nullptr;

    G4String driverPath;
    G4String fDriverFile;
    G4String fCt2DensityFile;
};

#endif
