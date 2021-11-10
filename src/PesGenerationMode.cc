#include "PesGenerationMode.hh"
#include "SessionManager.hh"
#include "StackingAction.hh"
#include "out.hh"

#include "G4Material.hh"

PesGenerationMode::PesGenerationMode(int numEvents, const std::string & outputFileName, bool binaryOutput) :
    SimModeBase(), NumEvents(numEvents)
{
    PesGenRecord C12C11(6, 12, "C11");
        C12C11.CrossSection.push_back({0,     100.0}); // !!!*** todo loader
        C12C11.CrossSection.push_back({500.0, 100.0}); // !!!*** todo loader
    BaseRecords.push_back(C12C11);

    PesGenRecord O16C11(8, 16, "C11");
        O16C11.CrossSection.push_back({0,     100.0}); // !!!*** todo loader
        O16C11.CrossSection.push_back({500.0, 100.0}); // !!!*** todo loader
    BaseRecords.push_back(O16C11);

    //bNeedGui    = true; // temporary! !!!***
    bNeedOutput = true;
    SessionManager & SM = SessionManager::getInstance();
    SM.FileName   = outputFileName;
    SM.bBinOutput = binaryOutput;
    SaveDir[0] = 0; SaveDir[1] = 0; SaveDir[2] = 1.0;
}

G4UserStackingAction * PesGenerationMode::getStackingAction()
{
    return new PesGeneratorStackingAction();
}

void PesGenerationMode::preInit()
{
    SessionManager::getInstance().FastPESGeneration = true;

    /*
    //interpolation test
    PesGenRecord r(1,1,"");
    r.MFP = { {0, 1},{1, 2},{2, 3},{3, 4},{5, 5},{10, 6} };
    out(r.getMFP(0));
    out(r.getMFP(100));
    out(r.getMFP(0.5));
    out(r.getMFP(6));
    exit(0);
    */
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
                out("  -> Z =", z, "A =", a, "Fraction =", frac);

                updateMatRecords(iMat, z, a, numberDensity * frac);
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

    if (it == CrossSection.begin()) return CrossSection.front().second;  // first: energy, second: CS
    if (it == CrossSection.end())   return CrossSection.back(). second;

    // interpolation
    // (e1, A) -> (energy, ?) -> (e2, B)  ==>  ? = A + (B-A)*(energy-e1)/(e2-e1)
    const auto lowIt = it - 1;
    return lowIt->second + (it->second - lowIt->second)*(energy - lowIt->first)/(it->first - lowIt->first);
}

void PesGenerationMode::saveRecord(const std::string & Pes, double X, double Y, double Z, double Time) const
{
    SessionManager & SM = SessionManager::getInstance();

    if (SM.bBinOutput)
    {
        *SM.outStream << (char)0xFF;
        *SM.outStream << Pes << (char)0x00;
        SM.outStream->write((char*)&SaveEnergy,   sizeof(double));
        SM.outStream->write((char*)&X,            sizeof(double));
        SM.outStream->write((char*)&Y,            sizeof(double));
        SM.outStream->write((char*)&Z,            sizeof(double));
        SM.outStream->write((char*)SaveDir,     3*sizeof(double));
        SM.outStream->write((char*)&Time,         sizeof(double));
    }
    else
    {
        std::stringstream ss;
        ss << Pes << ' ';
        ss << SaveEnergy << ' ';
        ss << X          << ' ' << Y          << ' ' << Z          << ' ';
        ss << SaveDir[0] << ' ' << SaveDir[1] << ' ' << SaveDir[2] << ' ';
        ss << Time;

        *SM.outStream << ss.rdbuf() << '\n';
    }

    out("->",Pes, "(",X,Y,Z,")", Time);
}

void PesGenerationMode::onEventStarted()
{
    SessionManager & SM = SessionManager::getInstance();
    if (SM.bBinOutput)
    {
        *SM.outStream << char(0xEE);
        SM.outStream->write((char*)&CurrentEvent, sizeof(int));
    }
    else
        *SM.outStream << '#' << CurrentEvent << '\n';

    CurrentEvent++;
}
