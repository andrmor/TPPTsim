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
class SourceModeBase;

class SourceModeFactory
{
public:
    static SourceModeBase * makeSourceModeInstance(const json11::Json & json);
};

class SourceModeBase
{
public:
    SourceModeBase(ParticleBase * particle, TimeGeneratorBase * timeGenerator); // transfers ownership
    virtual ~SourceModeBase();

    void initialize();

    virtual void GeneratePrimaries(G4Event * anEvent);

    virtual int CountEvents() {return -1;}

    virtual std::string getTypeName() const = 0;

    void writeToJson(json11::Json::object & json) const;
    void readFromJson(const json11::Json & json);

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
    PointSource(const json11::Json & json);

    std::string getTypeName() const override {return "PointSource";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    G4ThreeVector Origin = {0,0,0};
};

// ---

class LineSource : public SourceModeBase
{
public:
    LineSource(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & startPoint, const G4ThreeVector & endPoint);
    LineSource(const json11::Json & json);

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "LineSource";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    G4ThreeVector StartPoint = {0,0,0};
    G4ThreeVector EndPoint   = {0,0,0};
};

// ---

class ProfileBase
{
public:
    ProfileBase(){}
    virtual ~ProfileBase(){}

    //void addRotation(double Degrees);

    void setDirection(const G4ThreeVector & dir);

    virtual std::string getTypeName() const = 0;
    virtual void generateOffset(G4ThreeVector & pos) const = 0;

    void writeToJson(json11::Json::object & json) const;

protected:
    G4ThreeVector Direction;
    G4ThreeVector UnitPerp;
    double        Angle = 0;

    virtual void doWriteToJson(json11::Json::object & /*json*/) const {}
    //void readFromJson(const json11::Json & json);
};
class UniformProfile : public ProfileBase
{
public:
    UniformProfile(double dx, double dy) : ProfileBase(), DX(dx), DY(dy) {}
    UniformProfile(const json11::Json & json);

    std::string getTypeName() const override {return "Uniform";}
    void generateOffset(G4ThreeVector & pos) const override;

protected:
    double DX = 0;
    double DY = 0;

    void doWriteToJson(json11::Json::object & json) const override;
    void readFromJson(const json11::Json & json);
};
class GaussProfile : public ProfileBase
{
public:
    GaussProfile(double sigmaX, double sigmaY) : ProfileBase(), SigmaX(sigmaX), SigmaY(sigmaY) {}
    GaussProfile(const json11::Json & json);

    std::string getTypeName() const override {return "Gauss";}
    void generateOffset(G4ThreeVector & pos) const override;

protected:
    double SigmaX = 0;
    double SigmaY = 0;

    void doWriteToJson(json11::Json::object & json) const override;
    void readFromJson(const json11::Json & json);
};

class PencilBeam : public SourceModeBase
{
public:
    PencilBeam(ParticleBase * particle, TimeGeneratorBase * timeGenerator,
               const G4ThreeVector & origin, const G4ThreeVector & direction,
               int numParticles = 1, ProfileBase * spread = nullptr);
    PencilBeam(const json11::Json & json);
    ~PencilBeam();

    std::string getTypeName() const override {return "PencilBeam";}
    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    void update();

    G4ThreeVector Origin       = {0,0,0};
    int           NumParticles = 1;
    ProfileBase * Profile      = nullptr;
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
    MaterialLimitedSource(const json11::Json & json);
    ~MaterialLimitedSource();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "MaterialLimitedSource";}

protected:
    void init();
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);
    void customPostInit() override;

    G4ThreeVector Origin;
    G4ThreeVector BoundingBox;
    G4String      Material;
    G4String      FileName;

    //run-time
    G4Material    * SourceMat = nullptr;
    G4Navigator   * Navigator = nullptr;
    std::ofstream * Stream    = nullptr;
};

// ---

class NaturalLysoSource : public SourceModeBase
{
public:
    NaturalLysoSource(double timeFrom, double timeTo);
    NaturalLysoSource(const json11::Json & json);
    ~NaturalLysoSource();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "NaturalLysoSource";}

protected:
    void init();
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);
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
    BlurredPointSource(const json11::Json & json);
    ~BlurredPointSource();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "BlurredPointSource";}

protected:
    void init();
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    std::string FileName;
    Hist1DSampler * Sampler   = nullptr;
};

// --

#endif // SourceMode_h
