#ifndef SESSIONMANAGER_H
#define SESSIONMANAGER_H

#include "DetComp.hh"

#include <string>
#include <vector>

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

        void startSession();
        void endSession();

        void configureGUI();
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
        bool detectorContains(const std::string & component);
        void saveScintillatorTable(const std::string & fileName);
        int  isDirExists(const std::string & dirName);
        void saveConfig(const std::string & fileName) const;
        void loadConfig(const std::string & fileName);

     // Main settings
        SourceModeBase   * SourceMode    = nullptr;
        SimModeBase      * SimMode       = nullptr;
        PhantomModeBase  * PhantomMode   = nullptr;

        bool SimAcollinearity = false;
        bool KillNeutrinos    = false;

        DetComp DetectorComposition;

        std::string WorkingDirectory = "Sure+Does+Not+Exist";
        std::string FileName  = "TpptSim_DefaultSaveName.txt";

        bool bBinOutput       = false;

        int  Seed             = 0;      // long->int because of json11
        bool Verbose       = false;
        bool Debug           = false;
        bool ShowEventNumber = false;
        int  EvNumberInterval = 1000;

        std::vector<ScintRecord> ScintRecords;

        G4ParticleDefinition * GammaPD = nullptr;

     // Cuts
        double CutPhantomGamma    = 10.0*mm;
        double CutPhantomElectron = 10.0*mm;
        double CutPhantomPositron = 0.1 *mm;
        double CutScintGamma      = 0.1 *mm;
        double CutScintElectron   = 0.1 *mm;
        double CutScintPositron   = 0.1 *mm;

     // Misc
        double activityLYSO = 281.0; // decays per second per cm3

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
        G4Region            * regPhantom  = nullptr;
        G4Region            * regScint    = nullptr;

     // External resources
        G4Material          * ScintMat    = nullptr;
        G4VPhysicalVolume   * physWorld   = nullptr;

        G4ParticleGun       * ParticleGun = nullptr;

        CLHEP::RanecuEngine * randGen     = nullptr;

        G4UIExecutive       * ui          = nullptr;
        G4RunManager        * runManager  = nullptr;
        G4VisManager        * visManager  = nullptr;
};

#endif // SESSIONMANAGER_H
