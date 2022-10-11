#include "PesHistogramSource.hh"
#include "SessionManager.hh"
#include "out.hh"
#include "jstools.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

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

PesHistogramSource::PesHistogramSource(const std::string & directory) :
    SourceModeBase(nullptr, nullptr), Directory(directory)
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

PesHistogramSource::~PesHistogramSource()
{
    //delete inStream;
}

void PesHistogramSource::init()
{
    IsotopeBase.push_back({ "C11", {"6_12_C11.txt", "8_16_C11.txt"} });
    IsotopeBase.push_back({ "O15", {"8_16_O15.txt"} });
    IsotopeBase.push_back({ "N13", {"8_16_N13.txt"} });

    checkInputData();

    CurrentRecord = 0;
    CurrentFile = 0;

    /*
    inStream = new std::ifstream(FileName);
    if (!inStream->is_open())
    {
        out("Cannot open input file with particles:" + FileName);
        exit(1000);
        delete inStream; inStream = nullptr;
    }
    */
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
                std::ifstream infile(fn);
                if (!infile.good())
                {
                    out("File not found:", fn);
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
    PesName = pes.Isotope;

    const std::string & fileName = pes.SpatialFiles[CurrentFile];

    /*

        const size_t numChan = iso.SpatialFiles.size();
        out("\n====>Isotope", iso.Isotope, " -> ", numChan, "channel(s)");
        if (numChan == 0) continue;

        for (size_t iChan = 0; iChan < numChan; iChan++)
        {
            const std::string fn = DataDirectory + '/' + iso.SpatialFiles[iChan];
            out(fn);

            BinningParameters thisMapping;
            thisMapping.read(fn);

            std::ifstream in(fn);
            std::string line;
            std::getline(in, line); // skip the first line which is the mapping json

            for (int ix = 0; ix < Mapping.NumBins[0]; ix++)
                for (int iy = 0; iy < Mapping.NumBins[1]; iy++)
                {
                    std::getline(in, line);
                    size_t sz; // std::string::size_type sz;

                    for (int iz = 0; iz < Mapping.NumBins[2]; iz++)
                    {
                        double val = std::stod(line, &sz);
                        line = line.substr(sz);
                        //DATA[ix][iy][iz] += val;
                    }
                }
        }

    CurrentFile++;
    if (CurrentFile >= pes.SpatialFiles.size())
    {
        CurrentRecord++;
        CurrentFile = 0;
    }
    */

    /*
    if (!inStream)
    {
        out("Stream with particles does not exist: check file name!");
        return;
    }

    if (inStream->eof())
    {
        out("End of file reached!");
        return;
    }
    if (!inStream->good())
    {
        out("Error in reading particle data from file!");
        exit(111);
    }

        for (std::string line; std::getline(*inStream, line); )
        {
            //std::cout << "line=>" << line << "<=" << std::endl;
            if (line.size() < 1) continue; //allow empty lines

            if (line[0] == '#') break; //event finished

            if (inStream->eof())
            {
                out("---End of file reached!");
                return;
            }

            std::stringstream ss(line);  // units in file are mm MeV and ns
            ss >> name
               >> energy
               >> pos[0] >> pos[1] >> pos[2]
               >> dir[0] >> dir[1] >> dir[2]
               >> time;
            if (ss.fail())
            {
                out("Unexpected format of a line in the file with the input particles");
                exit(2111);
            }

            addPrimary(anEvent);
        }
    */
}

void PesHistogramSource::addPrimary(G4Event * anEvent)
{
    pd = ParticleGenerator.makeGeant4Particle(PesName);

    //out("->", pd->GetParticleName(), energy, pos, dir, time);

    ParticleGun->SetParticleDefinition(pd);
    ParticleGun->SetParticleEnergy(energy*keV);
    ParticleGun->SetParticleTime(time*ns);
    ParticleGun->SetParticlePosition(pos);
    ParticleGun->SetParticleMomentumDirection(dir);

    ParticleGun->GeneratePrimaryVertex(anEvent);
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
    json["Directory"] = Directory;
}

void PesHistogramSource::doReadFromJson(const json11::Json & json)
{
    jstools::readString(json, "Directory", Directory);
}
