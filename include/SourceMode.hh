#ifndef SourceMode_h
#define SourceMode_h

#include "G4ThreeVector.hh"
#include "json11.hh"

#include <vector>

class ParticleBase;
class TimeGeneratorBase;
class G4ParticleGun;
class G4Event;
class G4Material;
class G4Navigator;
class Hist1D;
class Hist1DSampler;

class SourceModeBase
{
public:
    SourceModeBase(ParticleBase * particle, TimeGeneratorBase * timeGenerator); // transfers ownership
    ~SourceModeBase();

    void initialize();

    virtual void GeneratePrimaries(G4Event * anEvent);

    virtual int CountEvents() {return -1;}

    virtual std::string getTypeName() const = 0;
    void writeToJson(json11::Json::object & json) const;

protected:
    virtual void customPostInit() {}
    virtual void doWriteToJson(json11::Json::object & json) const = 0;

    ParticleBase      * Particle      = nullptr;
    TimeGeneratorBase * TimeGenerator = nullptr;
    G4ParticleGun     * ParticleGun   = nullptr;

    G4ThreeVector Direction = {0, 0, 1.0};  // more general approach will be realized later!
    bool bIsotropicDirection = true;

    bool bGeneratePair  = false;    // for second gamma

protected:
    G4ThreeVector generateDirectionIsotropic();
    void generateSecondGamma(G4Event * anEvent);
};

// ---

class PointSource : public SourceModeBase
{
public:
    PointSource(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin);

    std::string getTypeName() const override {return "PointSource";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    G4ThreeVector Origin;
};

// ---

class LineSource : public SourceModeBase
{
public:
    LineSource(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & startPoint, const G4ThreeVector & endPoint);

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "LineSource";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    G4ThreeVector StartPoint;
    G4ThreeVector EndPoint;
};

// ---

class PencilBeam : public SourceModeBase
{
public:
    PencilBeam(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin, const G4ThreeVector & direction);

    std::string getTypeName() const override {return "PencilBeam";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    G4ThreeVector Origin;
};

// ---

class MaterialLimitedSource : public SourceModeBase
{
public:
    MaterialLimitedSource(ParticleBase * particle,
                          TimeGeneratorBase * timeGenerator,
                          const G4ThreeVector & origin, const G4ThreeVector & boundingBoxFullSize,
                          const G4String & material,
                          G4String fileName_EmissionPositions = "");
    ~MaterialLimitedSource();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "MaterialLimitedSource";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void customPostInit() override;

    G4ThreeVector Origin;
    G4ThreeVector BoundingBox;
    G4String      Material;
    G4String      FileName;

    //run-time
    G4Material    * SourceMat = nullptr;
    G4Navigator   * Navigator = nullptr;
    std::ofstream * Stream = nullptr;
};

// ---

class NaturalLysoSource : public SourceModeBase
{
public:
    NaturalLysoSource(double timeFrom, double timeTo);
    ~NaturalLysoSource();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "NaturalLysoSource";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void customPostInit() override;

    double TimeFrom;
    double TimeTo;

    double ScintMaxRadius = 0;
    std::vector<std::pair<double,double>> ElectronSpectrum; //format: energy[keV] relative_probablility

    //run-time
    G4Navigator   * Navigator = nullptr;
    Hist1DSampler * Sampler   = nullptr;
};

// ---

class BlurredPointSource : public PointSource
{
public:
    BlurredPointSource(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin, G4String fileName);
    ~BlurredPointSource();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "BlurredPointSource";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    std::string FileName;
    Hist1DSampler * Sampler   = nullptr;
};

// --

#endif // SourceMode_h
