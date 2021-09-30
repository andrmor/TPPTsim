﻿#include "SessionManager.hh"
#include "SourceMode.hh"
#include "FromFileSource.hh"
#include "TimeGenerator.hh"
#include "DefinedParticles.hh"
#include "SimMode.hh"
#include "PhantomMode.hh"
#include "out.hh"

#include <chrono>

int main(int argc, char** argv)
{
    SessionManager & SM = SessionManager::getInstance();

// --- START of user init ---

  // General settings
    SM.Seed              = 0;
    SM.bSimAcollinearity = true;  // only for the phantom region!
    SM.bKillNeutrinos    = true;

    double timeFrom = 0;
    double timeTo   = 1e-5*s;  // currently implemented only for the natural rad from LYSO!

    //SM.WorkingDirectory  = "/home/andr/WORK/TPPT";
    SM.WorkingDirectory = "/data/margarida/Data";

    SM.bG4Verbose        = false;
    SM.bDebug            = false;
    SM.bShowEventNumber  = true; SM.EvNumberInterval = 10000;

  // Phantom
    //SM.PhantomMode      = new PhantomNone;
    SM.PhantomMode      = new PhantomPMMA;
    //SM.PhantomMode      = new PhantomTinyCube;
    //SM.PhantomMode      = new PhantomDerenzo(200.0, 100.0, {1.8, 2.0, 2.2, 2.5, 3.0, 6.0}, 20.0, 10.0, 45.0);
    //SM.PhantomMode      = new PhantomParam;

  // Detector
    //SM.DetectorComposition = {};
    //SM.DetectorComposition = {DetComp::Scintillators};
    //SM.DetectorComposition = {DetComp::Scintillators, DetComp::FirstStageMonitor};
    //SM.DetectorComposition = {DetComp::Scintillators, DetComp::GDML};
    //SM.DetectorComposition = {DetComp::Scintillators, DetComp::SIPM};
    //SM.DetectorComposition = {DetComp::Scintillators, DetComp::SIPM, DetComp::PCB};
    SM.DetectorComposition = {DetComp::Scintillators, DetComp::SIPM, DetComp::PCB, DetComp::CopperStructure};
    //SM.DetectorComposition = {DetComp::Scintillators, DetComp::SIPM, DetComp::PCB, DetComp::CopperStructure, DetComp::CoolingAssemblies};
    //SM.DetectorComposition = {DetComp::Scintillators, DetComp::SIPM, DetComp::PCB, DetComp::CopperStructure, DetComp::CoolingAssemblies, DetComp::Base};
    //SM.DetectorComposition = {DetComp::Scintillators, DetComp::SIPM, DetComp::PCB, DetComp::CopperStructure, DetComp::CoolingAssemblies, DetComp::Base, DetComp::ClosedStructure};

  // Source
    //SM.SourceMode       = new PointSource(new GammaPair, new ExponentialTime(0, 2.034*60*s), {1.2, 2.3, SM.GlobalZ0+2});
    //SM.SourceMode       = new BlurredPointSource(new GammaPair, new ExponentialTime(0, 2.034*60*s), {0, 0, SM.GlobalZ0}, "/data/margarida/Data/AnnihilTest.txt");
    //SM.SourceMode       = new PointSource(new O15, new ConstantTime(0), {0, 0, SM.GlobalZ0});
    //SM.SourceMode       = new PencilBeam(new GammaPair(511.0*keV, true), new ConstantTime(0), {0, 0, SM.GlobalZ0}, {1.0,0,0});
    SM.SourceMode       = new PencilBeam(new Geantino, new ConstantTime(0), {0, 0, SM.GlobalZ0+12.3*mm}, {1.0,0,0});
    //SM.SourceMode       = new MaterialLimitedSource(new O15, new ConstantTime(0), {0, 0, SM.GlobalZ0}, {200.0,200.0,200.0}, "G4_WATER", "/home/andr/WORK/TPPT/der.txt");
    //SM.SourceMode       = new MaterialLimitedSource(new GammaPair, new ExponentialTime(0, 2.034*60.0*s), {0, 0, SM.GlobalZ0}, {200.0,200.0,200.0}, "G4_WATER", "/home/andr/WORK/TPPT/der.txt");
    //SM.SourceMode       = new NaturalLysoSource(timeFrom, timeTo);
    //SM.SourceMode       = new FromFileSource("/home/andr/WORK/TPPT/FirstStage.bin", true);

  // Operation mode
    SM.SimMode          = new SimModeGui();
    //SM.SimMode          = new SimModeShowEvent(119);
    //SM.SimMode          = new SimModeScintPosTest();
    //SM.SimMode          = new SimModeTracing();
    //SM.SimMode          = new SimModeAcollinTest(10000, 2.0, 100, "AcolTest.txt");
    //SM.SimMode          = new SimModeAnnihilTest(1e6, 10, 1000, "AnnihilTest.txt");
    //SM.SimMode          = new SimModeNatRadTest(1000000, 500, "natRadEnergyDistr.txt");
    //SM.SimMode          = new SimModeSingleEvents();
    //SM.SimMode          = new SimModeMultipleEvents(1000, "SimOutput.txt", false);
    //SM.SimMode          = new SimModeMultipleEvents(1e4, "SimOutput.bin", true);
    //SM.SimMode          = new SimModeMultipleEvents(SM.getNumberNatRadEvents(timeFrom, timeTo), "SimOutput.bin", true);
    //SM.SimMode          = new SimModeFirstStage(1e3, "FirstStage.bin", true);
    //SM.SimMode          = new SimModeMultipleEvents(SM.SourceMode->CountEvents(), "SimOutput.txt", false); // if using FromFileSource to use all events in the file

// --- END of user init ---
    SM.startSession(argc, argv);
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    SM.SimMode->run();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    out("Run time", (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count())*1e-6, "s");
    SM.endSession();
}
