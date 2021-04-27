#ifndef SESSIONMANAGER_H
#define SESSIONMANAGER_H

#include <string>
#include <vector>

#include "Modes.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

class G4Material;
class G4String;
class G4UImanager;
class G4UIExecutive;
class G4RunManager;
class G4VisManager;
class SourceModeBase;
class SimModeBase;
class PhantomModeBase;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4ParticleGun;
class G4ParticleDefinition;
class G4Region;
namespace CLHEP { class RanecuEngine; }

class SessionManager
{
    public:
        static SessionManager& getInstance();

    private:
        SessionManager();
        ~SessionManager();

    public:
        SessionManager(SessionManager const&) = delete;
        void operator=(SessionManager const&) = delete;

        void startSession(int argc, char ** argv);
        void endSession();

        void configureGUI(int argc, char ** argv);
        void startGUI();
        void configureOutput();
        void configureRandomGenerator();
        void initializeSource();
        void configureVerbosity();
        void scanMaterials();
        void registerAcollinearGammaModel(G4Region * region);

        int  countScintillators() const;

        bool detectorContains(DetComp component) const;

     // Main settings
        SourceModeBase   * SourceMode    = nullptr;
        SimModeBase      * SimMode       = nullptr;
        PhantomModeBase  * PhantomMode   = nullptr;

        bool bSimAcollinearity = false;

        std::vector<DetComp> DetectorComposition;

        std::string WorkingDirectory = "Sure+Does+Not+Exist";
        std::string FileName = "TpptSim_DefaultSaveName.txt";

        bool bBinOutput   = false;

        long Seed         = 0;
        bool bG4Verbose   = false;
        bool bDebug       = false;

        std::vector<G4ThreeVector> ScintPositions; //Scintillator positions, calculated during DetectorConstruction
        G4ParticleDefinition * GammaPD = nullptr;

     // Geometry
        int    NumScintX  = 8;
        int    NumScintY  = 8;

        double ScintSizeX = 3.005  * mm;
        double ScintSizeY = ScintSizeX;
        double ScintSizeZ = 15.0 * mm;

        double TeflonThick = 0.195 * mm;

        double ScintPitchX = ScintSizeX + TeflonThick;
        double ScintPitchY = ScintSizeY + TeflonThick;

        double EncapsSizeX = 25.8 * mm;
        double EncapsSizeY = EncapsSizeX;
        double EncapsSizeZ = ScintSizeZ;

        int    NumSegments = 12;
        int    NumRows     = 4;

        double RowGap      = 0.6  * mm;
        double Angle0      = 40.5 * deg; // just a guess, calculate please!
        double AngularStep = 9.0  * deg;

        double InnerDiam   = 335.4 * mm;

        double GlobalZ0    = 58.85 * mm;

     // Internal resources
        std::ofstream       * outStream  = nullptr;

     // External resources
        G4Material          * ScintMat    = nullptr;
        G4LogicalVolume     * logicWorld  = nullptr;
        G4VPhysicalVolume   * physWorld   = nullptr;
        G4LogicalVolume     * logicScint  = nullptr;

        G4ParticleGun       * ParticleGun = nullptr;

        CLHEP::RanecuEngine * randGen     = nullptr;

        G4UIExecutive       * ui          = nullptr;
        G4RunManager        * runManager  = nullptr;
        G4VisManager        * visManager  = nullptr;
};

#endif // SESSIONMANAGER_H
