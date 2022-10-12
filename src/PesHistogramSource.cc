#include "PesHistogramSource.hh"
#include "SessionManager.hh"
#include "ActivityLoader.hh"
#include "out.hh"
#include "jstools.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4RandomTools.hh"
#include "G4ThreeVector.hh"

#include <ios>
#include <iostream>
#include <sstream>
#include <fstream>

GeantParticleGenerator::GeantParticleGenerator()
{
    std::vector<std::string> allElements = {"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs"};
    for (size_t i = 0; i < allElements.size(); i++)
        ElementToZ.emplace( std::make_pair(allElements[i], i+1) );
}

G4ParticleDefinition * GeantParticleGenerator::makeGeant4Particle(const std::string & particleName)
{
    G4ParticleDefinition * particle = G4ParticleTable::GetParticleTable()->FindParticle(particleName);

    if (!particle)
    {
        // is it an ion?
        int Z, A;
        double E;
        bool ok = extractIonInfo(particleName, Z, A, E);
        if (!ok)
        {
            out("Found an unknown particle: " + particleName);
            exit(112);
        }

        particle = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z, A, E*keV);

        if (!particle)
        {
            out("Failed to generate ion: " + particleName);
            exit(113);
        }

        //std::cout << particleName << "   ->   " << particle->GetParticleName() << std::endl;
    }

    return particle;
}

bool GeantParticleGenerator::extractIonInfo(const std::string & text, int & Z, int & A, double & E)
{
    size_t size = text.length();
    if (size < 2) return false;

    // -- extracting Z --
    const char & c0 = text[0];
    if (c0 < 'A' || c0 > 'Z') return false;
    std::string symbol;
    symbol += c0;

    size_t index = 1;
    const char & c1 = text[1];
    if (c1 >= 'a' && c1 <= 'z')
    {
        symbol += c1;
        index++;
    }
    try
    {
        Z = ElementToZ.at(symbol);
    }
    catch (...)
    {
        return false;
    }

    // -- extracting A --
    A = 0; E = 0;
    char ci;
    while (index < size)
    {
        ci = text[index];
        if (ci < '0' || ci > '9') break;
        A = A*10 + (int)ci - (int)'0';
        index++;
    }
    if (A == 0) return false;

    if (index == size) return true;

    // -- extracting excitation energy --
    if (ci != '[') return false;
    index++;
    std::stringstream excEnergy;
    while (index < size)
    {
        ci = text[index];
        if (ci == ']')
        {
            excEnergy >> E;
            return !excEnergy.fail();
        }
        excEnergy << ci;
        index++;
    }
    return false;
}

// -------------------------------------

PesHistogramSource::PesHistogramSource(const std::string & directory, double multiplier) :
    SourceModeBase(nullptr, nullptr), Directory(directory), Multiplier(multiplier)
{
    init();
}

PesHistogramSource::PesHistogramSource(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);

    init();
}

void PesHistogramSource::init()
{
    IsotopeBase.push_back({ "C11", {"6_12_C11.txt", "8_16_C11.txt"} });
    IsotopeBase.push_back({ "O15", {"8_16_O15.txt"} });
    IsotopeBase.push_back({ "N13", {"8_16_N13.txt"} });

    checkInputData();

    CurrentRecord = 0;
    CurrentFile = 0;

    ParticleGun->SetParticleEnergy(0*keV);
    ParticleGun->SetParticleMomentumDirection({0,0,1.0});
}

void PesHistogramSource::checkInputData()
{
    for (auto it = IsotopeBase.begin(); it < IsotopeBase.end(); ++it)
    {
        for (auto itOther = IsotopeBase.begin(); itOther < IsotopeBase.end(); ++itOther)
        {
            if (it == itOther) continue;
            if (it->Isotope == itOther->Isotope)
            {
                out("Found non-unique record for PES:", it->Isotope);
                exit(3);
            }

            for (const std::string & fn : it->SpatialFiles)
            {
                std::string fullname = Directory + '/' + fn;
                std::ifstream infile(fullname);
                if (!infile.good())
                {
                    out("File not found:", fullname);
                    exit(3);
                }
            }
        }
    }
}

void PesHistogramSource::GeneratePrimaries(G4Event * anEvent)
{
    if (CurrentRecord >= IsotopeBase.size()) return; // finished with all events

    const PesDataRecord & pes = IsotopeBase[CurrentRecord];
    G4ParticleDefinition * pd = ParticleGenerator.makeGeant4Particle(pes.Isotope);
    ParticleGun->SetParticleDefinition(pd);

    const std::string & fileName = Directory + '/' + pes.SpatialFiles[CurrentFile];
    out("==>", fileName);

    BinningParameters binning;
    std::vector<std::vector<std::vector<double>>> data;
    ActivityLoader::load(fileName, data, binning);

    for (int ix = 0; ix < binning.NumBins[0]; ix++)
        for (int iy = 0; iy < binning.NumBins[1]; iy++)
            for (int iz = 0; iz < binning.NumBins[2]; iz++)
            {
                const int number = data[ix][iy][iz] * Multiplier;
                if (number < 1) continue;

                for (int iPart = 0; iPart < number; iPart++)
                {
                    const double x = binning.Origin[0] + (0.5 + ix) * binning.BinSize[0];
                    const double y = binning.Origin[1] + (0.5 + iy) * binning.BinSize[1];
                    const double z = binning.Origin[2] + (0.5 + iz) * binning.BinSize[2];
                    ParticleGun->SetParticlePosition( {x,y,z} );

                    ParticleGun->SetParticleTime(G4UniformRand() * TimeSpan);

                    ParticleGun->GeneratePrimaryVertex(anEvent);
                }
            }

    CurrentFile++;
    if (CurrentFile >= IsotopeBase[CurrentRecord].SpatialFiles.size())
    {
        CurrentRecord++;
        CurrentFile = 0;
    }
}

double PesHistogramSource::CountEvents()
{
    int numEvents = 0;
    for (const PesDataRecord & rec : IsotopeBase) numEvents += rec.SpatialFiles.size();
    out("Number of events is equal to the number of histograms:", numEvents);
    return numEvents;
}

void PesHistogramSource::doWriteToJson(json11::Json::object & json) const
{
    json["Directory"]  = Directory;
    json["Multiplier"] = Multiplier;
}

void PesHistogramSource::doReadFromJson(const json11::Json & json)
{
    jstools::readString(json, "Directory",  Directory);
    jstools::readDouble(json, "Multiplier", Multiplier);
}
