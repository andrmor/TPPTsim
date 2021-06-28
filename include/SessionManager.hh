#ifndef SESSIONMANAGER_H
#define SESSIONMANAGER_H

#include <string>
#include <vector>

#include "Modes.hh"
#include "ScintRecord.hh"
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
class G4VModularPhysicsList;
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
        void createFastSimulationPhysics(G4VModularPhysicsList * physicsList);
        void registerAcollinearGammaModel(G4Region * region);
        void registerParticleKillerModel(G4Region * region);
        void createPhantomRegion(G4LogicalVolume * logVolPhantom);
        void createScintillatorRegion(G4LogicalVolume * logVolScint);

        int  countScintillators() const;

        int  getNumberNatRadEvents(double timeFromInNs, double timeToInNs) const;

        bool detectorContains(DetComp component) const;

        void saveScintillatorTable(const std::string & fileName);

     // Main settings
        SourceModeBase   * SourceMode    = nullptr;
        SimModeBase      * SimMode       = nullptr;
        PhantomModeBase  * PhantomMode   = nullptr;

        bool bSimAcollinearity = false;
        bool bKillNeutrinos    = false;

        std::vector<DetComp> DetectorComposition;

        std::string WorkingDirectory = "Sure+Does+Not+Exist";
        std::string FileName  = "TpptSim_DefaultSaveName.txt";

        bool bBinOutput       = false;

        long Seed             = 0;
        bool bG4Verbose       = false;
        bool bDebug           = false;
        bool bShowEventNumber = false;
        int  EvNumberInterval = 1000;

        std::vector<ScintRecord> ScintRecords;

        G4ParticleDefinition * GammaPD = nullptr;

     // Misc
        double activityLYSO = 281.0; // decays per second per cm3

     // Geometry - Scintillators + Encapsulation
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

        double GlobalZ0    = 55.675 * mm; // (BaseHeight + SystHeight) / 2

     // Geometry - Base
        double RMin        = 16.3 * cm;
        double RMax        = 25.0 * cm;
        double SystHeight  = EncapsSizeX * 4 + RowGap * 3;
        double BaseHeight  = 6.35 * mm;

     // Internal resources
        std::ofstream       * outStream  = nullptr;

     // External resources
        G4Material          * ScintMat    = nullptr;
        G4VPhysicalVolume   * physWorld   = nullptr;

        G4Region            * regPhantom  = nullptr;
        G4Region            * regScint    = nullptr;

        G4ParticleGun       * ParticleGun = nullptr;

        CLHEP::RanecuEngine * randGen     = nullptr;

        G4UIExecutive       * ui          = nullptr;
        G4RunManager        * runManager  = nullptr;
        G4VisManager        * visManager  = nullptr;
};

#endif // SESSIONMANAGER_H
