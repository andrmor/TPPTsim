#include "SessionManager.hh"
#include "Modes.hh"
#include "SimMode.hh"
#include "SteppingAction.hh"
#include "SensitiveDetectorScint.hh"
#include "out.hh"

#include "G4RunManager.hh"

SimModeGui::SimModeGui()
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void SimModeGui::run()
{
    SessionManager& SM = SessionManager::getInstance();
    SM.startGUI();
}

// ---

SimModeShowEvent::SimModeShowEvent(int EventToShow) : iEvent(EventToShow)
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void SimModeShowEvent::run()
{
    SessionManager& SM = SessionManager::getInstance();
    if (iEvent != 0) SM.runManager->BeamOn(iEvent);
    SimModeGui::run();
}

// ---

SimModeScintPosTest::SimModeScintPosTest()
{
    bNeedGui    = false;
    bNeedOutput = false;
}

void SimModeScintPosTest::run()
{
    SessionManager& SM = SessionManager::getInstance();

    SM.runManager->BeamOn(1);

    outFlush();
    if (Hits > 1) SumDelta /= Hits;
    out("\n---Test results---\nTotal hits of the scintillators:", Hits, "Max delta:", MaxDelta, " Average delta:", SumDelta, "\n\n");
}

G4UserSteppingAction * SimModeScintPosTest::getSteppingAction()
{
    return new SteppingAction_ScintPosTest;
}

// ---

SimModeSingleEvents::SimModeSingleEvents()
{
    bNeedGui    = false;
    bNeedOutput = true;

    NumEvents   = 10000;

    SessionManager& SM = SessionManager::getInstance();
    SM.FileName = "Coincidence-GammaPairs-Test1.txt";
}

void SimModeSingleEvents::run()
{
    SessionManager& SM = SessionManager::getInstance();

    const double EnergyThreshold = 0.500*MeV;

    const int NumScint = SM.countScintillators();
    for(int i=0; i < NumScint; i++) ScintData.push_back({0,0,0});
    std::vector<int> hits;

    for (int iRun = 0; iRun < NumEvents; iRun++)
    {
        SM.runManager->BeamOn(1);

        hits.clear();
        for (int iScint = 0; iScint < NumScint; iScint++)
            if (ScintData[iScint][1] > EnergyThreshold)
                hits.push_back(iScint);

        if (hits.size() == 2)
        {
            for (int i = 0; i < 2; i++)
            {
                const int iScint = hits[i];
                const double Time   = ScintData[iScint][0] / ns;
                const double Energy = ScintData[iScint][1] / MeV;
                const G4ThreeVector & Pos   = SM.ScintPositions.at(iScint);
                const double X = Pos[0] / mm;
                const double Y = Pos[1] / mm;
                const double Z = Pos[2] / mm;

                out("Scint#",iScint, Time,"ns ", Energy, "MeV  xyz: (",X,Y,Z,")  Run# ",iRun);

                if (SM.outStream)
                    *SM.outStream << X << " " << Y << " " << Z << " " << Time << " " << Energy << std::endl;

            }
            out("---");
        }

        for (int i = 0; i < NumScint; i++) ScintData[i] = {0,0,0};
    }

    outFlush();
    if (!SM.outStream) out("\nOutput stream was not created, nothing was saved");
    else out("Data saved to file:", SM.WorkingDirectory + "/" + SM.FileName);
}

G4VSensitiveDetector * SimModeSingleEvents::getScintDetector()
{
    return new SensitiveDetectorScint_SingleEvents("Scint");
}

// ---

SimModeMultipleEvents::SimModeMultipleEvents(int numEvents, const std::string & FileName, bool bBinary)
{
    bNeedGui    = false;
    bNeedOutput = true;

    NumEvents = numEvents;

    SessionManager & SM = SessionManager::getInstance();
    SM.bBinOutput  = bBinary;
    SM.FileName    = FileName;
    InitialReserve = 10000;

    bDoCluster     = true;
    MaxTimeDif     = 0.1 * ns;
}

void SimModeMultipleEvents::run()
{
    SessionManager& SM = SessionManager::getInstance();

    const int NumScint = SM.countScintillators();
    DepositionData.resize(NumScint);
    for(auto & vec : DepositionData) vec.reserve(InitialReserve);

    SM.runManager->BeamOn(NumEvents);

    saveData();

    outFlush();
    if (!SM.outStream) out("\nOutput stream was not created, nothing was saved");
    else
    {
        out("\nData saved to file:", SM.WorkingDirectory + "/" + SM.FileName);
        if (bDoCluster) out("Depositions were clustered using",MaxTimeDif,"ns time threshold");
    }
}

G4VSensitiveDetector *SimModeMultipleEvents::getScintDetector()
{
    return new SensitiveDetectorScint_MultipleEvents("SD");
}

void SimModeMultipleEvents::saveData()
{
    SessionManager& SM = SessionManager::getInstance();
    const int numScint = SM.countScintillators();

    if(SM.bBinOutput)
    {
        for (int iScint = 0; iScint < numScint; iScint++)
        {
            auto & nodes = DepositionData[iScint];

            if (!nodes.empty())
            {
                *SM.outStream << char(0xee);
                SM.outStream->write((char*)&iScint, sizeof(int));
                SM.outStream->write((char*)&SM.ScintPositions[iScint][0], sizeof(double));
                SM.outStream->write((char*)&SM.ScintPositions[iScint][1], sizeof(double));
                SM.outStream->write((char*)&SM.ScintPositions[iScint][2], sizeof(double));

                for (int iNodes = 0; iNodes < nodes.size(); iNodes++)
                {
                    *SM.outStream << char(0xff);
                    SM.outStream->write((char*)&DepositionData[iScint][iNodes].time,   sizeof(double));
                    SM.outStream->write((char*)&DepositionData[iScint][iNodes].energy, sizeof(double));
                }
            }
        }
    }
    else
    {
        for (int iScint = 0; iScint < numScint; iScint++)
        {
            const G4ThreeVector & sp = SM.ScintPositions[iScint];
            auto & nodes = DepositionData[iScint];

            if (!nodes.empty())
            {
                *SM.outStream << "# " << iScint << " " << sp[0] << " " << sp[1] << " " << sp[2] << std::endl;

                for (const DepositionNodeRecord & n : nodes)
                    *SM.outStream << n.time << " " << n.energy << std::endl;
            }
        }
        for (auto & vec : DepositionData) vec.clear();
    }
}

void DepositionNodeRecord::merge(const DepositionNodeRecord & other)
{
    if (other.energy <= 0) return;

    const double newEnergy = energy + other.energy;
    time = (time * energy  +  other.time * other.energy) / newEnergy;
    energy = newEnergy;
}

bool DepositionNodeRecord::isCluster(const DepositionNodeRecord &other, double maxTimeDelta) const
{
    if ( fabs(time - other.time) > maxTimeDelta ) return false;
    return true;
}

// ---

SimModeTracing::SimModeTracing()
{
    bNeedGui    = true;
    bNeedOutput = false;
}

void SimModeTracing::run()
{
    SimModeGui::run();
}

G4UserSteppingAction * SimModeTracing::getSteppingAction()
{
    return new SteppingAction_Tracing;
}

// ---

SimModeAcollinTest::SimModeAcollinTest(int numRuns, const std::string & fileName) :
    NumRuns(numRuns)
{
    bNeedGui    = false;
    bNeedOutput = true;

    SessionManager & SM = SessionManager::getInstance();
    SM.bBinOutput  = false;
    SM.FileName    = fileName;
}

void SimModeAcollinTest::run()
{
    SessionManager& SM = SessionManager::getInstance();

    Histogram.resize(numBins);
    for (int iBin = 0; iBin < numBins; iBin++) Histogram[iBin] = 0;

    for (int iRun = 0; iRun < NumRuns; iRun++)
    {
        ParentTrackId = -1;
        //out("Run #", iRun);

        SM.runManager->BeamOn(1);

        if (Gammas.size() >= 2)
        {
            double angle = Gammas[0].dir.angle(Gammas[1].dir) / deg - 1e-10; // -1e-10 to get 180 deg in the previous bin
            int index = (angle - angleFrom) / deltaAngle;
            if      (index <  0)
            {
                numUnderflows++;

                if (Gammas[0].energy > 0.511 || Gammas[0].energy < 0.51 ||
                    Gammas[1].energy > 0.511 || Gammas[1].energy < 0.51)
                {
                    numNotTherm++;
                    //out("Not thermalized positron");
                    //out("Undeflow angle:", angle);
                    //out("Dir/Energies:");
                    //for (const DirAndEnergy & de : Gammas) std::cout << de.dir << " " << de.energy << std::endl;
                }
                else
                {
                    out("strange event!");
                    out("Undeflow angle:", angle);
                    out("Dir/Energies:");
                    for (const DirAndEnergy & de : Gammas) std::cout << de.dir << " " << de.energy << std::endl;
                }
            }
            else if (index >= numBins) numOverflows++;
            else Histogram[index]++;
        }
        else
        {
            //out("Unexpected: number of gammas is", Gammas.size());
        }
        Gammas.clear();
    }

    outFlush();
    out("\nDistribution of inter-gamma angles (from", angleFrom,"to 180 deg):");
    int sum = 0;
    for (int iBin = 0; iBin < numBins; iBin++)
    {
        std::cout << (iBin == 0 ? '[' : ',');
        std::cout << Histogram[iBin];
        sum += Histogram[iBin];

        if (SM.outStream)
            *SM.outStream << (angleFrom + deltaAngle * (iBin+0.5)) << " " << Histogram[iBin] << std::endl;
    }
    std::cout << ']' << std::endl;
    out("Distribution sum:", sum);
    out("Underflows:", numUnderflows);
    out("Overflows:", numOverflows);
    out("NotThermalized:", numNotTherm);
}

G4UserSteppingAction *SimModeAcollinTest::getSteppingAction()
{
    return new SteppingAction_AcollinearityTester;
}

void SimModeAcollinTest::addDirection(const G4ThreeVector & v, int parentID, double energy)
{
    // the first gamma to track will be annihilation one, some of the following ones can be secondary ones!
    //out("ParentId:", parentID, "Energy:", energy);
    if (ParentTrackId == -1) ParentTrackId = parentID;
    else
    {
        if (parentID != ParentTrackId) return;
    }

    Gammas.push_back( DirAndEnergy(v, energy) );
}

// ---

SimModeAnnihilTest::SimModeAnnihilTest(int numRuns, const std::string & fileName) :
    NumRuns(numRuns)
{
    bNeedGui    = false;
    bNeedOutput = true;

    SessionManager & SM = SessionManager::getInstance();
    SM.bBinOutput  = false;
    SM.FileName    = fileName;
}

G4UserSteppingAction *SimModeAnnihilTest::getSteppingAction()
{
    return new SteppingAction_AnnihilationTester;
}

void SimModeAnnihilTest::run()
{
    SessionManager & SM = SessionManager::getInstance();

    Histogram.resize(numBins);
    for (int iBin = 0; iBin < numBins; iBin++) Histogram[iBin] = 0;

    SM.runManager->BeamOn(NumRuns);

    outFlush();
    out("\nDistribution of annihilation positions (from", positionFrom,"to 4 mm):");
    int sum = 0;
    for (int iBin = 0; iBin < numBins; iBin++)
    {
        std::cout << (iBin == 0 ? '[' : ',');
        std::cout << Histogram[iBin];
        sum += Histogram[iBin];

        if (SM.outStream)
            *SM.outStream << (positionFrom + deltaPosition * iBin) << " " << Histogram[iBin] << std::endl;
    }
    std::cout << ']' << std::endl;
    out("Distribution sum:", sum);
    out("Underflows:", numUnderflows);
    out("Overflows:", numOverflows);
}

void SimModeAnnihilTest::addPosition(double x)
{
    int index = (x - positionFrom) / deltaPosition;
    if      (index <  0)
    {
        numUnderflows++;
        out("Underflow position:", x);
    }
    else if (index >= numBins) numOverflows++;
    else Histogram[index]++;
}
