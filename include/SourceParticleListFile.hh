#ifndef fromfilesource_h
#define fromfilesource_h

#include "SourceMode.hh"
#include "GeantParticleGenerator.hh"

#include "G4ThreeVector.hh"

#include <string>
#include <map>

class G4ParticleDefinition;

class SourceParticleListFile : public SourceModeBase
{
public:
    SourceParticleListFile(const std::string & fileName, bool bBinaryFile);
    SourceParticleListFile(const json11::Json & json);
    ~SourceParticleListFile();

    void GeneratePrimaries(G4Event * anEvent) override;

    double CountEvents() override;

    std::string getTypeName() const override {return "SourceParticleListFile";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    void init();

    void prepareStream();
    void addPrimary(G4Event * anEvent);

    GeantParticleGenerator ParticleGenerator;

    std::string FileName;
    bool bBinary = false;
    std::ifstream * inStream = nullptr;
    std::map<std::string, int> ElementToZ;

    std::string name;
    G4ParticleDefinition * pd;
    double energy, time;
    G4ThreeVector pos, dir;
};

#endif // fromfilesource_h
