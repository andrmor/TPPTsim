#include "SourceParticleListFile.hh"
#include "SessionManager.hh"
#include "out.hh"
#include "jstools.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"

#include <ios>
#include <iostream>
#include <sstream>
#include <fstream>

SourceParticleListFile::SourceParticleListFile(const std::string & fileName, bool bBinaryFile) :
    SourceModeBase(nullptr, nullptr), FileName(fileName), bBinary(bBinaryFile)
{
    init();
}

SourceParticleListFile::SourceParticleListFile(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);

    init();
}

void SourceParticleListFile::init()
{
    if (bBinary) inStream = new std::ifstream(FileName, std::ios::in | std::ios::binary);
    else         inStream = new std::ifstream(FileName);

    if (!inStream->is_open())
    {
        out("Cannot open input file with particles:" + FileName);
        exit(1000);
        delete inStream; inStream = nullptr;
    }
    prepareStream();
}

SourceParticleListFile::~SourceParticleListFile()
{
    delete inStream;
}

void SourceParticleListFile::prepareStream()
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

void SourceParticleListFile::GeneratePrimaries(G4Event * anEvent)
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
                onNewEventMarker();
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
            //std::cout << "line=>" << line << "<=" << std::endl;
            if (line.size() < 1) continue; //allow empty lines

            if (line[0] == '#')
            {
                onNewEventMarker();
                break;
            }

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

void SourceParticleListFile::addPrimary(G4Event * anEvent)
{
    pd = ParticleGenerator.makeGeant4Particle(name);

    //out("->", pd->GetParticleName(), energy, pos, dir, time);

    ParticleGun->SetParticleDefinition(pd);
    ParticleGun->SetParticleEnergy(energy*keV);
    ParticleGun->SetParticleTime(time*ns);
    ParticleGun->SetParticlePosition(pos);
    ParticleGun->SetParticleMomentumDirection(dir);

    ParticleGun->GeneratePrimaryVertex(anEvent);
}

double SourceParticleListFile::CountEvents()
{
    if (!inStream) return 0;

    inStream->clear();
    inStream->seekg(0);

    double iCounter = 0;
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

void SourceParticleListFile::doWriteToJson(json11::Json::object & json) const
{
    json["FileName"] = FileName;
    json["bBinary"]  = bBinary;
}

void SourceParticleListFile::doReadFromJson(const json11::Json & json)
{
    jstools::readString(json, "FileName", FileName);
    jstools::readBool  (json, "bBinary",  bBinary);
}

// ---------------------------

SourceParticleListFileMicroPET::SourceParticleListFileMicroPET(const std::string &fileName, bool bBinaryFile, double pushBackDistance) :
    SourceParticleListFile(fileName, bBinaryFile), PushBackDistance(pushBackDistance) {}

SourceParticleListFileMicroPET::SourceParticleListFileMicroPET(const json11::Json & json) :
    SourceParticleListFile(json) {}

void SourceParticleListFileMicroPET::doWriteToJson(json11::Json::object & json) const
{
    SourceParticleListFile::doWriteToJson(json);
    json["PushBackDistance"] = PushBackDistance;
}

void SourceParticleListFileMicroPET::doReadFromJson(const json11::Json & json)
{
    SourceParticleListFile::doReadFromJson(json);
    jstools::readDouble(json, "PushBackDistance", PushBackDistance);
}

void SourceParticleListFileMicroPET::onNewEventMarker()
{
    FirstParticleThisEvent = true;
}

#include "G4RandomTools.hh"
void SourceParticleListFileMicroPET::addPrimary(G4Event * anEvent)
{
    pd = ParticleGenerator.makeGeant4Particle(name);

    ParticleGun->SetParticleDefinition(pd);
    ParticleGun->SetParticleEnergy(energy*keV);
    ParticleGun->SetParticleTime(time*ns);

    if (FirstParticleThisEvent)
    {
        G4ThreeVector posProj{pos[0], pos[1], 0};
        double Xedge = 140.0;  // 0.5*InnerDiam - CubeSize = 0.5*335.4 - 25.4 - 1.0 --> 141.3 mm
        double Yedge = 12.5;
        G4ThreeVector edge{-Xedge, Yedge, 0};
        double phi0 = posProj.angle(edge); // always positive!
        bool upperHalf = (pos[1] > -(Yedge/Xedge*pos[0]));
        if (!upperHalf) phi0 = -phi0;

        double angleOfView = 2.0 * atan(Yedge / Xedge);
        Phi = phi0 + G4UniformRand() * angleOfView;
        if (G4UniformRand() > 0.5) Phi += 180*deg;

        FirstParticleThisEvent = false;
    }
    pos.rotateZ(Phi);
    dir.rotateZ(Phi);

    //out("->", pd->GetParticleName(), energy, pos, dir, time);

    ParticleGun->SetParticleMomentumDirection(dir);
    if (PushBackDistance != 0) pos -= dir * PushBackDistance;
    ParticleGun->SetParticlePosition(pos);

    ParticleGun->GeneratePrimaryVertex(anEvent);
}
