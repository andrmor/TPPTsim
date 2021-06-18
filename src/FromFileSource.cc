#include "FromFileSource.hh"
#include "SessionManager.hh"
#include "out.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include <ios>
#include <iostream>
#include <sstream>
#include <fstream>

FromFileSource::FromFileSource(const std::string & fileName, bool bBinaryFile) :
    SourceModeBase(nullptr, nullptr), bBinary(bBinaryFile)
{
    if (bBinary) inStream = new std::ifstream(fileName, std::ios::in | std::ios::binary);
    else         inStream = new std::ifstream(fileName);

    if (!inStream->is_open())
    {
        out("Cannot open input file with particles:" + fileName);
        exit(1000);
        delete inStream; inStream = nullptr;
    }
    prepareStream();

    std::vector<std::string> allElements = {"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs"};
    for (size_t i = 0; i < allElements.size(); i++)
        ElementToZ.emplace( std::make_pair(allElements[i], i+1) );
}

FromFileSource::~FromFileSource()
{
    delete inStream;
}

void FromFileSource::prepareStream()
{
    inStream->clear();
    inStream->seekg(0);

    bool bError = true;
    if (bBinary)
    {
        char ch = 0;
        int iEv;
        inStream->get(ch);
        if (ch == (char)0xEE)
        {
            inStream->read((char*)&iEv, sizeof(int));
            bError = false;
        }
    }
    else
    {
        std::string line;
        std::getline(*inStream, line);
        if (line.size() > 1 && line[0] == '#') bError = false;
    }

    if (bError)
    {
        out("Input file for FromFileSource does not start with the event marker!");
        exit(1002);
    }
}

void FromFileSource::GeneratePrimaries(G4Event * anEvent)
{
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

    if (bBinary)
    {
        char ch;
        int iEv;
        while (inStream->get(ch))
        {
            if (inStream->eof()) break;

            if (ch == (char)0xEE)
            {
                inStream->read((char*)&iEv, sizeof(int));
                //event finished
                break;
            }
            else if (ch == (char)0xFF)
            {
                name.clear();
                while (inStream->get(ch))
                {
                    if (ch == (char)0x00) break;
                    name += ch;
                }

                inStream->read((char*)&energy, sizeof(double));
                inStream->read((char*)&pos[0], sizeof(double));
                inStream->read((char*)&pos[1], sizeof(double));
                inStream->read((char*)&pos[2], sizeof(double));
                inStream->read((char*)&dir[0], sizeof(double));
                inStream->read((char*)&dir[1], sizeof(double));
                inStream->read((char*)&dir[2], sizeof(double));
                inStream->read((char*)&time,   sizeof(double));

                if (inStream->eof())
                {
                    out("---End of file reached!");
                    return;
                }

                if (!inStream->good())
                {
                    out("Error in reading particle data from file!");
                    exit(111);
                }

                addPrimary(anEvent);
            }
        }
    }
    else
    {
        for (std::string line; std::getline(*inStream, line); )
        {
            std::cout << "line=>" << line << "<=" << std::endl;
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
    }
}

void FromFileSource::addPrimary(G4Event * anEvent)
{
    pd = makeGeant4Particle(name);

    //out("->", pd->GetParticleName(), energy, pos, dir, time);

    ParticleGun->SetParticleDefinition(pd);
    ParticleGun->SetParticleEnergy(energy*keV);
    ParticleGun->SetParticleTime(time*ns);
    ParticleGun->SetParticlePosition(pos);
    ParticleGun->SetParticleMomentumDirection(dir);

    ParticleGun->GeneratePrimaryVertex(anEvent);
}

int FromFileSource::CountEvents()
{
    if (!inStream) return 0;

    inStream->clear();
    inStream->seekg(0);

    int iCounter = 0;
    if (bBinary)
    {
        char ch;
        int iEv;
        double buf[8]; //en + pos + dir + time
        while (inStream->get(ch))
        {
            if (inStream->eof()) break;

            if (ch == (char)0xEE)
            {
                inStream->read((char*)&iEv, sizeof(int));
                iCounter++;
                continue;
            }
            else if (ch == (char)0xFF)
            {
                while (inStream->get(ch))
                {
                    if (ch == (char)0x00) break;
                }

                inStream->read((char*)buf, 8 * sizeof(double));

                if (inStream->eof()) break;

                if (!inStream->good())
                {
                    out("Format error in binary particle data file!");
                    exit(111);
                }
            }
            else
            {
                out("Data header char error in binary particle data file!");
                exit(111);
            }
        }
    }
    else
    {
        for( std::string line; std::getline( *inStream, line ); )
        {
            if (inStream->fail()) break;
            if (line.length() > 0 && line[0] == '#') iCounter++;
        }
    }

    prepareStream();

    out("Found events:", iCounter);
    return iCounter;
}

G4ParticleDefinition * FromFileSource::makeGeant4Particle(const std::string & particleName)
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

bool FromFileSource::extractIonInfo(const std::string & text, int & Z, int & A, double & E)
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
