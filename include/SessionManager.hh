#ifndef SESSIONMANAGER_H
#define SESSIONMANAGER_H

#include <string>
#include <vector>

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

class G4Material;
class G4String;
class G4UImanager;
class G4UIExecutive;
class G4RunManager;
class G4VisManager;
class SimModeBase;
class G4LogicalVolume;
namespace CLHEP { class RanecuEngine; }

enum class SourceModeEnum   {GammaPair, C11, C10, O15};
enum class DetectorModeEnum {OnlyScint, WithDetector};
enum class PhantomModeEnum  {PMMA};

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

        void startSession(int argc, char **argv);
        void endSession();

        void configureGUI(int argc, char** argv);
        void startGUI();
        void configureOutput();
        void configureRandomGenerator();
        void configureVerbosity();
        void scanMaterials();

        SimModeBase * SimulationMode = nullptr;

        std::string WorkingDirectory;
        std::string FileName         = "TpptSimTest---123.txt";

        long Seed         = 0;

        //Geometry:
        G4Material      * ScintMat  = nullptr;
        G4LogicalVolume * logicWorld = nullptr;
        G4LogicalVolume * logicScint = nullptr;

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
        double Angle0      = 40.5 * deg; // just a guess!
        double AngularStep = 9.0  * deg;

        double InnerDiam   = 335.4 * mm;

        std::vector<G4ThreeVector> ScintPositions;

        //Tests (Stepping Action):
        bool   bVerbose = false;
        double MaxDelta = 0;
        int    Hits     = 0;
        double SumDelta = 0;

        //Scintillators Data:
        std::vector<G4ThreeVector> ScintData;

        double NumParticles = 1;

        CLHEP::RanecuEngine * randGen = nullptr;

        G4UIExecutive * ui         = nullptr;
        G4RunManager  * runManager = nullptr;
        G4VisManager  * visManager = nullptr;

        std::ofstream * outStream = nullptr;
};

#endif // SESSIONMANAGER_H
