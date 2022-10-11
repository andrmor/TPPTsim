#ifndef peshistogramsource_h
#define peshistogramsource_h

#include "SourceMode.hh"
#include "G4ThreeVector.hh"

#include <string>
#include <vector>
#include <map>

class G4ParticleDefinition;

// TODO !!!*** to a separate file, and use it in the FromFileSource
class GeantParticleGenerator
{
public:
    GeantParticleGenerator();

    G4ParticleDefinition * makeGeant4Particle(const std::string & particleName);
    bool extractIonInfo(const std::string & text, int & Z, int & A, double & E);

protected:
    std::map<std::string, int> ElementToZ;
};

// -------

struct PesDataRecord
{
    PesDataRecord(std::string isotope, const std::vector<std::string> & spatialFiles) :
        Isotope(isotope), SpatialFiles(spatialFiles) {}

    std::string Isotope;
    std::vector<std::string> SpatialFiles;
};

class PesHistogramSource : public SourceModeBase
{
public:
    PesHistogramSource(const std::string & directory);
    PesHistogramSource(const json11::Json & json);
    ~PesHistogramSource();

    void GeneratePrimaries(G4Event * anEvent) override;

    double CountEvents() override;

    std::string getTypeName() const override {return "PesHistogramSource";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    void init();

    void prepareStream();
    void addPrimary(G4Event * anEvent);
    void checkInputData();

    std::string Directory;

    std::vector<PesDataRecord> IsotopeBase;

    std::ifstream * inStream = nullptr;

    GeantParticleGenerator ParticleGenerator;

    std::string PesName;
    G4ParticleDefinition * pd;
    double energy, time;
    G4ThreeVector pos, dir;
};

#endif // peshistogramsource_h
