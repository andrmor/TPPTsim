#include "PesGenerationMode.hh"
#include "SessionManager.hh"
#include "out.hh"

#include "G4Material.hh"

PesGenerationMode::PesGenerationMode() :
    SimModeBase()
{
    PesGenRecord C12C11(6, 12, "C11");
        C12C11.CrossSection.push_back({0,     100.0}); // !!!*** todo loader
        C12C11.CrossSection.push_back({500.0, 100.0}); // !!!*** todo loader
    BaseRecords.push_back(C12C11);

    PesGenRecord O16C11(8, 16, "C11");
        O16C11.CrossSection.push_back({0,     100.0}); // !!!*** todo loader
        O16C11.CrossSection.push_back({500.0, 100.0}); // !!!*** todo loader
    BaseRecords.push_back(O16C11);

    bNeedGui = true; // temporary! !!!***
}

void PesGenerationMode::preInit()
{
    SessionManager::getInstance().FastPESGeneration = true;
}

void PesGenerationMode::run()
{
    exploreMaterials();

    SessionManager& SM = SessionManager::getInstance();
    SM.startGUI();
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

        thisRec.MFP = thisRec.CrossSection;
        for (auto & pair : thisRec.MFP) pair.second = 1e25 / pair.second / IsotopeNumberDensity; // millibarh = 0.001e-28m2 -> 0.001e-22mm2 -> 1e-25 mm2

        MaterialRecords[iMat].push_back(thisRec);
    }
}
