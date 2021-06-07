#ifndef fromfilesource_h
#define fromfilesource_h

#include "SourceMode.hh"

#include "G4ThreeVector.hh"

#include <string>
#include <map>

class G4ParticleDefinition;

class FromFileSource : public SourceModeBase
{
public:
    FromFileSource(const std::string & fileName, bool bBinaryFile);
    ~FromFileSource();

    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    std::ifstream * inStream = nullptr;
    bool bBinary = false;
    std::map<std::string, int> ElementToZ;

    std::string name;
    G4ParticleDefinition * pd;
    double energy, time;
    G4ThreeVector pos, dir;

    G4ParticleDefinition * makeGeant4Particle(const std::string & particleName);
    bool extractIonInfo(const std::string & text, int & Z, int & A, double & E);
};

#endif // fromfilesource_h
