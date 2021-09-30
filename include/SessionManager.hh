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
        std::string FileName   = "TpptSim_DefaultSaveName.txt";

        bool bBinOutput        = false;

        long Seed              = 0;
        bool bG4Verbose        = false;
        bool bDebug            = false;
        bool bShowEventNumber  = false;
        int  EvNumberInterval  = 1000;

        std::vector<ScintRecord> ScintRecords;

        G4ParticleDefinition * GammaPD = nullptr;

     // Misc
        double activityLYSO    = 281.0; // decays per second per cm3

     // Geometry - Scintillators + Encapsulation
        int    NumScintX       = 8;
        int    NumScintY       = 8;

        double ScintSizeX      = 3.005  * mm;
        double ScintSizeY      = ScintSizeX;
        double ScintSizeZ      = 15.0 * mm;

        double TeflonThick     = 0.195 * mm;

        double ScintPitchX     = ScintSizeX + TeflonThick;
        double ScintPitchY     = ScintSizeY + TeflonThick;

        double EncapsSizeX     = 25.8 * mm;
        double EncapsSizeY     = EncapsSizeX;
        double EncapsSizeZ     = ScintSizeZ;

        int    NumSegments     = 12;
        int    NumRows         = 4;

        double RowGap          = 0.6  * mm;
        double Angle0          = 40.5 * deg;
        double AngularStep     = 9.0  * deg;

        double InnerDiam       = 335.4 * mm;

     // Geometry - Base
        double BaseRMin        = 162.5 * mm;
        double BaseRMax        = 300.0   * mm;
        double SystHeight      = EncapsSizeY * 4.0 + RowGap * 3.0; //105.0 mm
        double BaseHeight      = 6.35  * mm;
        double BaseSegment     = 120.0 * deg;

        double GlobalZ0        = 0.5 * (BaseHeight + SystHeight); //55.675 mm

     // Geometry - Closed Structure
        double InnerWallThick  = 1.0 * mm;
        double OuterWallThick  = 1.0 * mm;
        double WallsSegment    = 120.0 * deg;
        double SideWallSegment = 1.0 * deg;

     // Geometry - SiPMs
     // https://www.hamamatsu.com/resources/pdf/ssd/s14160_s14161_series_kapd1064e.pdf, page 4
        double SIPMSizeX       = EncapsSizeX;
        double SIPMSizeY       = SIPMSizeX;
        double SIPMSizeZ       = 1.35 * mm;

     // Geometry - PCB
     // PCB1 - closer to SiPMs
        double SIPMPCB1Gap     = 0.010 * mm;
        double PCB1SizeX       = 28.2 * mm;
        double PCB1SizeY       = 52.2 * mm;
        double PCB1SizeZ       = 3.18 * mm;

     // PCB2
        double PCB2SizeX       = 25.2 * mm ;
        double PCB2SizeY       = 1.2 * mm;
        double PCB2SizeZ       = 36.0 * mm;
        double MiddlePCBGap    = (40.96 * mm - PCB2SizeZ) / 2.0 * mm;
        double PCB2Z1          = GlobalZ0 + 36.9 * mm;
        double PCB2Z2          = GlobalZ0 - 36.9 * mm;
        double PCB2Z3          = GlobalZ0 + 15.9 * mm;
        double PCB2Z4          = GlobalZ0 - 15.9 * mm;

     // PCB3
        double PCB3SizeX       = PCB1SizeX;
        double PCB3SizeY       = PCB1SizeY;
        double PCB3SizeZ       = 1.95 * mm;

        int    PCB1PCB3NumRows = 2;
        double PCB1PCB3Pitch   = PCB1SizeY;
        double PCB2Pitch       = SystHeight / 4.0;

     // Geometry - Copper Structure
     // Copper columns - from between the scintillators (Column1) to the copper connectors (Column3)
        double Column1SixeX1   = 3.16  * mm;
        double Column1SizeX2   = 2.21  * mm;
        double Column1SizeZ    = 6.0   * mm;
        double Column1SizeY    = SystHeight;
        double Column1RSurf    = 10.8  * mm; //approximate value. Can't be calculated

        double Column2SixeX1   = 1.25  * mm;
        double Column2SizeX2   = 0.75  * mm;
        double Column2SizeY    = SystHeight;
        double Column2SizeZ    = 3.27  * mm;

        double Column3SixeX1   = 4.85  * mm;
        double Column3SizeX2   = 4.25  * mm;
        double Column3SizeY    = SystHeight;
        double Column3SizeZ    = 3.82  * mm;

     // Copper pieces that connect Column3 to the horizontal copper connectors
        double Piece1SizeX1    = 6.29 * mm;
        double Piece1SizeX2    = 4.85 * mm;
        double Piece1SizeY     = 6.6  * mm;
        double Piece1SizeZ     = 9.0  * mm;
        double Piece1Z1        = GlobalZ0 + 37.4 * mm;
        double Piece1Z2        = GlobalZ0 - 37.4 * mm;

        double Piece2SizeX1    = Piece1SizeX1;
        double Piece2SizeX2    = Piece1SizeX2;
        double Piece2SizeY     = 5.9  * mm;
        double Piece2SizeZ     = Piece1SizeZ;
        double Piece2Z1        = GlobalZ0 + 14.15 * mm;
        double Piece2Z2        = GlobalZ0 - 14.15 * mm;

        double Piece3SizeX1    = 6.4  * mm;
        double Piece3SizeX2    = 6.29 * mm;
        double Piece3SizeY     = Piece1SizeY;
        double Piece3SizeZ     = 1.0  * mm;
        double Piece3Z1        = Piece1Z1;
        double Piece3Z2        = Piece1Z2;

        double Piece4SizeX1    = Piece3SizeX1;
        double Piece4SizeX2    = Piece3SizeX2;
        double Piece4SizeY     = Piece2SizeY;
        double Piece4SizeZ     = Piece3SizeZ;
        double Piece4Z1        = Piece2Z1;
        double Piece4Z2        = Piece2Z2;

    // Connectors - horizontal boxes that conect the external columns to the rest of the copper structure and "hide" the PCB2
        double ConnectorSizeX  = 6.4  * mm;
        double OutConnectSizeY = 6.6  * mm;
        double InConnectSizeY  = 5.9  * mm;
        double ConnectorSizeZ  = 21.6 * mm;

        double ConnectorZ1     = GlobalZ0 + 37.4 * mm;
        double ConnectorZ2     = GlobalZ0 - 37.4 * mm;
        double ConnectorZ3     = GlobalZ0 + 14.15 * mm;
        double ConnectorZ4     = GlobalZ0 - 14.15 * mm;

     // External columns (Column4)
        double ExtColumnSizeX  = 6.4 * mm;
        double ExtColumnSizeY  = 5.0 * mm;
        double ExtColumnSizeZ  = SystHeight;

     // Holders (pieces of copper connected to the external columns)
        double HolderSizeX     = 6.4  * mm;
        double HolderSizeY     = HolderSizeX;
        double HolderSizeZ     = 16.7 * mm;
        double HolderPitch     = 49.3 * mm;

     // Geometry - Cooling Assemblies
     // Copper pipe holder (copper box)
        double CopperBoxRmin   = 200.0 * mm;
        double CopperBoxRmax   = 216.0 * mm;
        double CopperBoxHeight = 6.0   * mm;
        double CopperBoxZ1     = ConnectorZ1 - 0.5 * (OutConnectSizeY + CopperBoxHeight) * mm;
        double CopperBoxZ2     = ConnectorZ2 + 0.5 * (OutConnectSizeY + CopperBoxHeight) * mm;
        double CopperBoxZ3     = ConnectorZ3 + 0.5 * (InConnectSizeY  + CopperBoxHeight) * mm;
        double CopperBoxZ4     = ConnectorZ4 - 0.5 * (InConnectSizeY  + CopperBoxHeight) * mm;

     // Water
        double WaterRmin       = 201.5 * mm;
        double WaterRmax       = 211.6 * mm;
        double WaterHeight     = 4.71  * mm;

        double CoolingSegment  = 120.0 * deg - 2.0 * SideWallSegment;

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
