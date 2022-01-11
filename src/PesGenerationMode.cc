#include "PesGenerationMode.hh"
#include "SessionManager.hh"
#include "StackingAction.hh"
#include "out.hh"
#include "jstools.hh"

#include "G4Material.hh"
#include "G4Track.hh"
#include "G4RandomTools.hh"
#include "Randomize.hh"

void PesGenerationMode::commonConstructor()
{
    SessionManager & SM = SessionManager::getInstance();

    //loadCrossSections("ProductionCrossSections.txt");
    loadCrossSections(SM.WorkingDirectory + "/SecretFile.txt");

    //bNeedGui    = true; // used only for tests!
}

PesGenerationMode::PesGenerationMode(int numEvents, const std::string & outputFileName, bool binaryOutput) :
    SimModeBase(), NumEvents(numEvents), bDirectMode(false)
{
    commonConstructor();

    SessionManager & SM = SessionManager::getInstance();
    bNeedOutput = true;
    SM.FileName   = outputFileName;
    SM.bBinOutput = binaryOutput;
    SaveDirection[0] = 0; SaveDirection[1] = 0; SaveDirection[2] = 1.0;

    // interpolation test
    //    PesGenRecord r(1,1,"");
    //    r.CrossSection = { {0, 1},{1, 2},{2, 3},{3, 4},{5, 5},{10, 6} };
    //    out(r.getCrossSection(0));
    //    out(r.getCrossSection(100));
    //    out(r.getCrossSection(0.5));
    //    out(r.getCrossSection(6));
    //    exit(0);
}

PesGenerationMode::PesGenerationMode(int numEvents, std::array<double, 3> binSize, std::array<int, 3> numBins, std::array<double, 3> origin) :
    SimModeBase(), NumEvents(numEvents),  bDirectMode(true),
    BinSize(binSize), NumBins(numBins), Origin(origin)
{
    commonConstructor();

    initProbArrays();

    // tests for voxel finding algorithms
    //    G4ThreeVector v(0, 2.6, -1.2);
    //    int index[3];
    //    bool ok = getVoxel(v, index);
    //    out(ok, index[0], index[1], index[2]);
    //    exit(1);

    //    std::vector<std::tuple<int, int, int, double>> Path;
    //    addPath({0,1.2,1}, {2,1,1}, Path);
    //    //addPath({1,1,1}, {2.2,1,1}, Path);
    //    //addPath({1.5,1.8,1}, {-2.5,2.1,1}, Path);
    //    for (size_t i = 0; i < Path.size(); i++)
    //        out(std::get<0>(Path[i]), std::get<1>(Path[i]), std::get<2>(Path[i]), std::get<3>(Path[i])); // --> 1 1 1 0.282843
    //    exit(1);

    //    std::vector<std::tuple<int, int, int, double>> Path;
    //    addPathA({0,0,0}, {2.0,-1.0,0}, Path);
    //    //addPathA({-2.0,-1.0,0}, {0,-2.0,1.0}, Path);
    //    //addPathA({0,0,0}, {2,-2,0}, Path);
    //    for (size_t i = 0; i < Path.size(); i++)
    //        out(std::get<0>(Path[i]), std::get<1>(Path[i]), std::get<2>(Path[i]), std::get<3>(Path[i])); // --> 1 1 1 0.282843
    //    exit(1);
}

void PesGenerationMode::initProbArrays()
{
    for (PesGenRecord & r : BaseRecords)
    {
        r.ProbArray = new std::vector<std::vector<std::vector<double>>>();
        // X
        r.ProbArray->resize(NumBins[0]);
        for (std::vector<std::vector<double>> & ary : *r.ProbArray)
        {
            // Y
            ary.resize(NumBins[1]);
            for (std::vector<double> & arz : ary)
            {
                // Z
                arz.resize(NumBins[2]);
            }
        }
    }
}

void PesGenerationMode::loadCrossSections(const std::string & fileName)
{
    std::ifstream inStream(fileName);
    if (!inStream.is_open())
    {
        out("Cannot open file with PES generation cross-sections:\n", fileName);
        exit(1);
    }

    bool fillingRecord = false;
    PesGenRecord currentRecord;
    for (std::string line; std::getline(inStream, line); )
    {
        //out(">>>",line);
        if (line.empty()) continue; //allow empty lines

        if (line[0] == '#')
        {
            //new reaction
            if (fillingRecord) BaseRecords.push_back(currentRecord);

            line.erase(0, 1);
            std::stringstream ss(line);
            ss >> currentRecord.TargetZ
               >> currentRecord.TargetA
               >> currentRecord.PES;
            if (ss.fail())
            {
                out("Unexpected format of a reaction line in the file with the PES cross-sections");
                exit(2);
            }
            currentRecord.CrossSection.clear();
            fillingRecord = true;
            //out("-->Processing reaction:",currentRecord.TargetZ, currentRecord.TargetA, currentRecord.PES);
        }
        else
        {
            std::stringstream ss(line);  // units in the file are MeV and mbarns
            double E, CS;
            ss >> E >> CS;
            if (ss.fail())
            {
                out("Unexpected format of a data line in the file with the PES cross-sections");
                exit(3);
            }
            //out(E, CS);
            currentRecord.CrossSection.push_back({E*MeV, CS});
        }
    }
    if (!currentRecord.CrossSection.empty()) BaseRecords.push_back(currentRecord);

    out("===== PES production cross-section summary:");
    out("Number of reactions:", BaseRecords.size());
    for (const auto & r : BaseRecords)
    {
        out(">", r.TargetZ, r.TargetA, r.PES, "  CS range from:", r.CrossSection.front().first, "to", r.CrossSection.back().first, "MeV");
    }
    out("\n");
}

G4UserStackingAction * PesGenerationMode::getStackingAction()
{
    return new PesGeneratorStackingAction();
}

void PesGenerationMode::preInit()
{
    SessionManager::getInstance().FastPESGeneration = true;
}

#include "G4MTRunManager.hh"
void PesGenerationMode::run()
{
    SessionManager& SM = SessionManager::getInstance();

    exploreMaterials();

    // this sub-mode is just to debug!
    if (bNeedGui)
    {
        SM.startGUI();
        return;
    }

    SM.runManager->BeamOn(NumEvents);

    if (bDirectMode) saveArrays();
}

void PesGenerationMode::exploreMaterials()
{
    MaterialRecords.clear();

    const G4MaterialTable * theMaterialTable = G4Material::GetMaterialTable();
    const size_t numMat = theMaterialTable->size();
    MaterialRecords.resize(numMat);
    out("Defined", numMat, "materials");

    for (size_t iMat = 0; iMat < numMat; iMat++)
    {
        out("---------");
        const G4Material * mat = (*theMaterialTable)[iMat];
        out(mat->GetName(), " Index:", mat->GetIndex());
        const G4ElementVector * ev = mat->GetElementVector();
        for (size_t iEl = 0; iEl < ev->size(); iEl++)
        {
            const G4Element * element = ev->at(iEl);
            double numberDensity = mat->GetVecNbOfAtomsPerVolume()[iEl]; // mm-3
            out(element->GetName(), " Atomic density:", numberDensity, "mm-3");

            const G4IsotopeVector * isoVec = element->GetIsotopeVector();
            double * relAbVec = element->GetRelativeAbundanceVector();
            for (size_t iIs = 0; iIs < isoVec->size(); iIs++)
            {
                const G4Isotope * is = isoVec->at(iIs);
                int z = is->GetZ();
                int a = is->GetN();
                double frac = relAbVec[iIs];
                double PartialNumDens = numberDensity * frac;
                out("  -> Z =", z, "A =", a, "Fraction =", frac, "PartNumDens = ", PartialNumDens);

                updateMatRecords(iMat, z, a, PartialNumDens);
            }
        }
    }
}

void PesGenerationMode::updateMatRecords(int iMat, int Z, int A, double IsotopeNumberDensity)
{
    for (const PesGenRecord & rec : BaseRecords)
    {
        if (rec.TargetZ != Z || rec.TargetA != A) continue;

        out("    ==> Adding PES generation record!");
        PesGenRecord thisRec = rec;
        thisRec.NumberDensity = IsotopeNumberDensity;
        MaterialRecords[iMat].push_back(thisRec);
    }
}

#include <algorithm>
double PesGenRecord::getCrossSection(double energy) const
{
    //Returns an iterator pointing to the first element in the range [first, last) that is greater than value, or last if no such element is found.
    auto it = std::upper_bound(CrossSection.begin(), CrossSection.end(), energy,
                               [](double one, const std::pair<double,double> & two)
                                 {return one < two.first;}
                              );

    //if (it == CrossSection.begin()) return CrossSection.front().second;
    if (it == CrossSection.begin()) return 0; // assuming the first data point is the threshold
    if (it == CrossSection.end())   return CrossSection.back().second;  // first: energy, second: CS

    // interpolation
    // (e1, A) -> (energy, ?) -> (e2, B)  ==>  ? = A + (B-A)*(energy-e1)/(e2-e1)
    const auto lowIt = it - 1;
    const double val = lowIt->second + (it->second - lowIt->second)*(energy - lowIt->first)/(it->first - lowIt->first);
    //out("Energy:", energy, "XS:", val); exit(111);
    return val;
}

void PesGenerationMode::saveRecord(const std::string & Pes, double X, double Y, double Z, double Time) const
{
    SessionManager & SM = SessionManager::getInstance();

    if (SM.bBinOutput)
    {
        *SM.outStream << (char)0xFF;
        *SM.outStream << Pes << (char)0x00;
        SM.outStream->write((char*)&SaveEnergy,     sizeof(double));
        SM.outStream->write((char*)&X,              sizeof(double));
        SM.outStream->write((char*)&Y,              sizeof(double));
        SM.outStream->write((char*)&Z,              sizeof(double));
        SM.outStream->write((char*)SaveDirection, 3*sizeof(double));
        SM.outStream->write((char*)&Time,           sizeof(double));
    }
    else
    {
        std::stringstream ss;
        ss << Pes << ' ';
        ss << SaveEnergy << ' ';
        ss << X << ' ' << Y << ' ' << Z << ' ';
        ss << SaveDirection[0] << ' ' << SaveDirection[1] << ' ' << SaveDirection[2] << ' ';
        ss << Time;

        *SM.outStream << ss.rdbuf() << '\n';
    }

    //out("->",Pes, "(",X,Y,Z,")", Time);
}

void PesGenerationMode::readFromJson(const json11::Json & json)
{
    jstools::readInt(json, "NumEvents",   NumEvents);
    jstools::readBool(json, "DirectMode", bDirectMode);

    if (bDirectMode)
    {
        // !!!*** control array sizes = 3
        // BinSize
        {
            json11::Json::array ar;
            jstools::readArray(json, "BinSize", ar);
            for (int i = 0; i < 3; i++) BinSize[i] = ar[i].number_value();
        }
        // NumBins
        {
            json11::Json::array ar;
            jstools::readArray(json, "NumBins", ar);
            for (int i = 0; i < 3; i++) NumBins[i] = ar[i].int_value();
        }
        // Origin
        {
            json11::Json::array ar;
            jstools::readArray(json, "Origin", ar);
            for (int i = 0; i < 3; i++) Origin[i] = ar[i].number_value();
        }

        bNeedOutput = false;
        initProbArrays();
    }
    else
    {
        SessionManager & SM = SessionManager::getInstance();
        jstools::readString(json, "OutputFileName", SM.FileName);
        jstools::readBool  (json, "BinaryOutput",   SM.bBinOutput);
    }
}

void PesGenerationMode::doWriteToJson(json11::Json::object & json) const
{
    json["NumEvents"] = NumEvents;
    json["DirectMode"] = bDirectMode;

    if (bDirectMode)
    {
        // BinSize
        {
            json11::Json::array ar;
            for (int i = 0; i < 3; i++) ar.push_back(BinSize[i]);
            json["BinSize"] = ar;
        }
        // NumBins
        {
            json11::Json::array ar;
            for (int i = 0; i < 3; i++) ar.push_back(NumBins[i]);
            json["NumBins"] = ar;
        }
        // Origin
        {
            json11::Json::array ar;
            for (int i = 0; i < 3; i++) ar.push_back(Origin[i]);
            json["Origin"] = ar;
        }
    }
    else
    {
        SessionManager & SM = SessionManager::getInstance();
        json["OutputFileName"] = SM.FileName;
        json["BinaryOutput"]   = SM.bBinOutput;
    }
}

void PesGenerationMode::onEventStarted()
{
    if (!bDirectMode)
    {
        SessionManager & SM = SessionManager::getInstance();
        if (SM.bBinOutput)
        {
            *SM.outStream << char(0xEE);
            SM.outStream->write((char*)&CurrentEvent, sizeof(int));
        }
        else
            *SM.outStream << '#' << CurrentEvent << '\n';
    }

    CurrentEvent++;
}

bool PesGenerationMode::modelTrigger(const G4Track * track)
{
    const int StepNumber = track->GetCurrentStepNumber();
    //out("PES call", StepNumber);

    if (StepNumber == 1)
    {
        LastEnergy      = track->GetKineticEnergy();
        LastTrackLength = track->GetTrackLength();
        LastPosition    = track->GetPosition();
        LastMaterial    = track->GetMaterial()->GetIndex();
        return false;
    }

    const double          Energy   = track->GetKineticEnergy();
    const double          Length   = track->GetTrackLength();
    const G4ThreeVector & Position = track->GetPosition();

    if (LastEnergy > Energy)
    {
        if (bDirectMode) doTriggerDirect(track);
        else
        {
            bool kill = doTriggerMC(track);
            if (kill) return true;
        }
    }

    LastEnergy      = Energy;
    LastTrackLength = Length;
    LastPosition    = Position;
    LastMaterial    = track->GetMaterial()->GetIndex();
    return false;
}

bool PesGenerationMode::doTriggerMC(const G4Track * track)
{
    const double stepLength = track->GetTrackLength() - LastTrackLength;
    const double meanEnergy = 0.5 * (track->GetKineticEnergy() + LastEnergy);
    //out("Step", stepLength, "MeanEenergy", meanEnergy, " Material index", LastMaterial);

    const std::vector<PesGenRecord> & Records = MaterialRecords[LastMaterial];
    if (!Records.empty())
    {
        //probability is proportional to CS * NumberDensity
        ProbVec.clear(); ProbVec.reserve(Records.size());
        double sumProb = 0;
        for (const PesGenRecord & r : Records)
        {
            const double cs = r.getCrossSection(meanEnergy);
            const double relProb = cs * r.NumberDensity;
            sumProb += relProb;
            ProbVec.push_back(relProb);
        }

        if (sumProb > 0)
        {
            size_t index = 0;

            if (ProbVec.size() > 1)
            {
                // selecting the reaction
                double val = sumProb * G4UniformRand();
                //out("Probs:",ProbVec.size(),"Sum:",sumProb,"Random:",val);
                for (; index+1 < ProbVec.size(); index++)
                {
                    if (val < ProbVec[index]) break;
                    val -= ProbVec[index];
                }
                //out("-->Selected:",index);
            }

            const double mfp = 1e25 / sumProb; // millibarn = 0.001e-28m2 -> 0.001e-22mm2 -> 1e-25 mm2

            const double trigStep = -mfp * log(G4UniformRand());
            if (trigStep < stepLength)
            {
                G4ThreeVector TriggerPosition = LastPosition + trigStep/stepLength*(track->GetPosition() - LastPosition);
                //out("Triggered! Position:", TriggerPosition, "  Last/FullStep positions:", LastPosition, Position);
                saveRecord(Records[index].PES, TriggerPosition[0], TriggerPosition[1], TriggerPosition[2], track->GetGlobalTime());
                return true; // safe to return, this proton will be killed
            }
        }
    }
    return false;
}

bool PesGenerationMode::getVoxel(const G4ThreeVector & pos, int * index)
{
    for (int i = 0; i < 3; i++)
        index[i] = floor( (pos[i] - Origin[i]) / BinSize[i] );

    for (int i = 0; i < 3; i++)
        if ( index[i] < 0 || index[i] >= NumBins[i]) return false;
    return true;
}

void PesGenerationMode::addPath(const G4ThreeVector & posFrom, const G4ThreeVector & posTo, std::vector<std::tuple<int, int, int, double>> & path)
{
    // check maybe start and finish are in the same voxel
    int From[3];
    bool ok1 = getVoxel(posFrom, From);
    int To[3];
    bool ok2 = getVoxel(posTo, To);
    if (ok1 && ok2 && From[0] == To[0] && From[1] == To[1] && From[2] == To[2])
    {
        path.push_back( {From[0], From[1], From[2], (posTo-posFrom).mag()} );
        return;
    }

    // general
    const G4ThreeVector vector = posTo - posFrom;
    const double distance = vector.mag();
    const int numSteps = 1 + distance * 10.0; // 0.1 mm step
    const double dStep = distance / numSteps;
    const double factor = 1.0/numSteps;
    for (int iStep = 0; iStep < numSteps; iStep++)
    {
        G4ThreeVector curPos = posFrom + vector * factor * iStep;
        bool ok = getVoxel(curPos, From);
        if (ok) path.push_back( {From[0], From[1], From[2], dStep} );
    }
}

void PesGenerationMode::addPathA(const G4ThreeVector & posFrom, const G4ThreeVector & posTo, std::vector<std::tuple<int, int, int, double>> & path)
{
    int From[3];
    bool ok1 = getVoxel(posFrom, From);
    int To[3];
    bool ok2 = getVoxel(posTo, To);
    const double totLength = (posTo-posFrom).mag();
    //out("Total step length:", totLength);
    if (ok1 && ok2 && From[0] == To[0] && From[1] == To[1] && From[2] == To[2])
    {
        path.push_back( {From[0], From[1], From[2], totLength} );
        return;
    }

    std::vector< std::tuple<double, int, int> > store; // progress factor (0->1), axis index (0,1,2) and current index along this axis
    for (int iAxis = 0; iAxis < 3; iAxis++)
    {
        //out(iAxis, "> from->to:", From[iAxis], To[iAxis]);
        int iv = From[iAxis];
        const int step = (posTo[iAxis] > posFrom[iAxis] ? 1 : -1);
        while (iv != To[iAxis])
        {
            //if (iv >= 0 && iv < NumBins[iAxis])
            //{
                const double k = (BinSize[iAxis]*(iv + (step == 1 ? 1 : 0)) + Origin[iAxis] - posFrom[iAxis]) / (posTo[iAxis] - posFrom[iAxis]);  // from 0 to 1.0
                store.push_back( {k, iAxis, iv + step} );
            //}
            iv += step;
        }
    }

    std::sort(store.begin(), store.end(),
              [](const std::tuple<double, int, int> & lhs, const std::tuple<double, int, int> & rhs)
                {return (std::get<0>(lhs) < std::get<0>(rhs));}  );

//    out("Store:");
//    for (size_t i = 0; i < store.size(); i++)
//        out("  ", std::get<0>(store[i]), std::get<1>(store[i]), std::get<2>(store[i]));
//    out("");

    // starting from "From" and updating indexes there
    double length = 0;
    double newK   = 0;
    for (const auto & rec : store)
    {
        newK = std::get<0>(rec);
        if (isValidVoxel(From))
            path.push_back( {From[0], From[1], From[2], (newK - length) * totLength} );
        length = newK;
        From[std::get<1>(rec)] = std::get<2>(rec);
    }
    if (isValidVoxel(From))
        path.push_back( {From[0], From[1], From[2], (1.0 - newK) * totLength} );
}

bool PesGenerationMode::isValidVoxel(int * coords) const
{
    for (int i = 0; i < 3; i++)
        if (coords[i] < 0 || coords[i] >= NumBins[i]) return false;
    return true;
}

void PesGenerationMode::doTriggerDirect(const G4Track * track)
{
    const std::vector<PesGenRecord> & Records = MaterialRecords[LastMaterial];
    if (Records.empty()) return;

    std::vector<std::tuple<int, int, int, double>> Path;
    addPathA(LastPosition, track->GetPosition(), Path);

    const double meanEnergy = 0.5 * (track->GetKineticEnergy() + LastEnergy);
    for (const PesGenRecord & r : Records)
    {
        const double cs = r.getCrossSection(meanEnergy);
        const double DProbByMM = 1e-25 * cs * r.NumberDensity; // millibarn = 0.001e-28m2 -> 0.001e-22mm2 -> 1e-25 mm2
        for (size_t i = 0; i < Path.size(); i++)
            (*r.ProbArray)[std::get<0>(Path[i])][std::get<1>(Path[i])][std::get<2>(Path[i])] += std::get<3>(Path[i]) * DProbByMM;
    }
}

void PesGenerationMode::saveArrays()
{
    SessionManager & SM = SessionManager::getInstance();
    for (PesGenRecord & r : BaseRecords)
    {
        std::string fileName = std::to_string(r.TargetZ) + '_' + std::to_string(r.TargetA) + '_' + r.PES + ".txt";
        std::string fullFileName = SM.WorkingDirectory + '/' + fileName;

        std::ofstream stream;
        stream.open(fullFileName);

        if (!stream.is_open())
        {
            out("Cannot open file to store output data!");
            outFlush();
            exit(1);
        }
        else out("\nSaving array data to file", fullFileName);

        for (const std::vector<std::vector<double>> & ary : *r.ProbArray)
        {
            for (const std::vector<double> & arz : ary)
            {
                for (const double & val : arz) stream << val << ' ';
                stream << '\n';
            }
        }
    }
}
