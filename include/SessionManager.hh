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
class G4VisManager;
namespace CLHEP { class RanecuEngine; }

class SessionManager
{
    public:
        static SessionManager& getInstance();

        enum RunModeEnum {GUI, ScintPosTest, Main, ShowEvent};

        enum SourceModeEnum {GammaPair, C11, C10, O15};

    private:
        SessionManager();
        ~SessionManager();

    public:
        SessionManager(SessionManager const&) = delete;
        void operator=(SessionManager const&) = delete;

        void startSession();
        void endSession();

        void runSimulation(int NumRuns);

        void testCoordinates();

        //bool bGuiMode     = false;

        void runGUI(G4UImanager* UImanager, G4UIExecutive * ui, G4VisManager * visManager);

        RunModeEnum runMode = GUI;

        SourceModeEnum SourceMode = GammaPair;

        std::string WorkingDirectory = "/data/margarida/Data";
        std::string BaseFileName     = "Data.txt";

        long Seed         = 0;

        //Geometry:

        G4Material * ScintMat  = nullptr;

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
        double Angle0      = 0    * deg;
        double AngularStep = 9.0 * deg;

        double InnerDiam   = 335.4 * mm;

        std::vector<G4ThreeVector> ScintPositions;

        //Tests (Stepping Action):

        bool bScintPositionTestMode = false;
        int Hits = 0;
        int Errors = 0;

        //Scintillators Data:

        std::vector<G4ThreeVector> ScintData;
        std::string FileName_Output;


        double NumParticles = 1;

        CLHEP::RanecuEngine * randGen = nullptr;

private:
        std::ofstream * outStream = nullptr;

};

#endif // SESSIONMANAGER_H
