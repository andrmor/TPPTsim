#ifndef SourceMode_h
#define SourceMode_h

#include "G4ThreeVector.hh"
#include "json11.hh"

#include <vector>
#include <string>

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
    static SourceModeBase * makeSourceInstance(const json11::Json & json);
};

class SourceModeBase
{
public:
    SourceModeBase(ParticleBase * particle, TimeGeneratorBase * timeGenerator); // transfers ownership
    virtual ~SourceModeBase();

    virtual void initialize();

    virtual void GeneratePrimaries(G4Event * anEvent);

    virtual double CountEvents() {return 1;}

    virtual std::string getTypeName() const = 0;

    void writeToJson(json11::Json::object & json) const;
    void readFromJson(const json11::Json & json);

    virtual void setParticleEnergy(double energy);

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

class SourceMixer : public SourceModeBase
{
public:
    SourceMixer(std::vector<std::pair<SourceModeBase*,double>> sourcesAndStatWeights);
    SourceMixer(const json11::Json & json);

    void initialize() override;
    void setParticleEnergy(double energy) override;

    std::string getTypeName() const override {return "SourceMixer";}

    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    void init();

    std::vector<std::pair<SourceModeBase*,double>> SourcesAndWeights;
    double SumWeight = 0;
};

// ---

class SourcePoint : public SourceModeBase
{
public:
    SourcePoint(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin);
    SourcePoint(const json11::Json & json);

    std::string getTypeName() const override {return "SourcePoint";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    G4ThreeVector Origin = {0,0,0};
};

// ---

class SourceLine : public SourceModeBase
{
public:
    SourceLine(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & startPoint, const G4ThreeVector & endPoint);
    SourceLine(const json11::Json & json);

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "SourceLine";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    G4ThreeVector StartPoint = {0,0,0};
    G4ThreeVector EndPoint   = {0,0,0};
};

// ---

class SourceCylinder : public SourceModeBase
{
public:
    SourceCylinder(ParticleBase * particle, TimeGeneratorBase * timeGenerator,
                      double radius, const G4ThreeVector & startPoint, const G4ThreeVector & endPoint,
                      const std::string & fileName = "");
    SourceCylinder(const json11::Json & json);
    ~SourceCylinder();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "SourceCylinder";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    double Radius = 1.0;
    G4ThreeVector StartPoint = {0,0,0};
    G4ThreeVector EndPoint   = {0,0,0};
    std::string FileName;

    G4ThreeVector Axis;
    G4ThreeVector UnitNormal1;
    G4ThreeVector UnitNormal2;

    std::ofstream * Stream    = nullptr;

private:
    void init();
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
    UniformProfile(double fullSizeX, double fullSizeY) : ProfileBase(), DX(fullSizeX), DY(fullSizeY) {}
    UniformProfile(const json11::Json & json);

    std::string getTypeName() const override {return "Uniform";}
    void generateOffset(G4ThreeVector & pos) const override;

protected:
    double DX = 0;
    double DY = 0;

    void doWriteToJson(json11::Json::object & json) const override;
    void readFromJson(const json11::Json & json);
};
class RoundProfile : public ProfileBase
{
public:
    RoundProfile(double diameter) : ProfileBase(), Diameter(diameter) {}
    RoundProfile(const json11::Json & json);

    std::string getTypeName() const override {return "Round";}
    void generateOffset(G4ThreeVector & pos) const override;

protected:
    double Diameter = 0;

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
class CustomProfile : public ProfileBase
{
public:
    // file should contain pairs [Position_mm, StatWeight] for the entire range (from minus to plus), recommend to start from and stop with the positions with zero weight
    CustomProfile(std::string fileDistributionX, std::string fileDistributionY);
    CustomProfile(const json11::Json & json);
    ~CustomProfile();

    std::string getTypeName() const override {return "Custom";}
    void generateOffset(G4ThreeVector & pos) const override;

protected:
    const bool DoLogPositions = true;
    std::vector<std::pair<double,double>> DistX;
    std::vector<std::pair<double,double>> DistY;

    void doWriteToJson(json11::Json::object & json) const override;
    void readFromJson(const json11::Json & json);

    // runtime
    Hist1DSampler * XSampler = nullptr;
    Hist1DSampler * YSampler = nullptr;
    std::ofstream * logStream = nullptr;

    void init();
};

class SourceBeam : public SourceModeBase
{
public:
    SourceBeam(ParticleBase * particle, TimeGeneratorBase * timeGenerator,
               const G4ThreeVector & origin, const G4ThreeVector & direction,
               int numParticles = 1, ProfileBase * spread = nullptr);
    SourceBeam(const json11::Json & json);
    ~SourceBeam();

    std::string getTypeName() const override {return "SourceBeam";}
    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    void update();

    G4ThreeVector Origin       = {0,0,0};
    int           NumParticles = 1;
    ProfileBase * Profile      = nullptr;
};

struct BeamRecord
{
    double Energy;
    double XIsoCenter;
    double ZIsoCenter;
    double TimeStart;
    double TimeSpan;
    double StatWeight;

    void writeToJson(json11::Json::object & json) const;
    void readFromJson(const json11::Json & json);

    void print() const;
};

// Beam is aligned with Y axis (downwards) and starts from Y = StartBeamFromY
// Divergence is determined using the apex position = Origin
class SourceMultiBeam : public SourceModeBase
{
public:
    SourceMultiBeam(ParticleBase * particle, const std::vector<BeamRecord> & beams, double totalParticles);
    SourceMultiBeam(ParticleBase * particle, const std::string & beamletFileName, double totalParticles); // NomEnergy[MeV] XIso[mm] ZIso[mm] Time0[ns] TimeSpan[ns] StatWeight
    SourceMultiBeam(const json11::Json & json);

    double CountEvents() override {return NumParticles;}

    std::string getTypeName() const override {return "SourceMultiBeam";}
    void GeneratePrimaries(G4Event * anEvent) override;

    std::vector<std::pair<double,double>> getTimeWindows(double delayAfter, double marginBefore) const; // only for beamlets consecutive in time; returns vector of pairs[from, duration]

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    void loadCalibration();
    void loadBeamletData(const std::string & beamletDataFile);
    void calculateParticlesPerStatWeightUnit();

    const G4ThreeVector Origin  = {0, 2520.0, 0};
    const double StartBeamFromY = 150.0;
    std::vector<BeamRecord> Beams;
    double NumParticles = 0;

    std::vector<std::array<double,3>> Calibration; // Nominal_energy[MeV] True_energy[MeV] SpotSigma[mm]
    const std::string CalibrationFileName = "BeamletCalibration.txt";

    //runtime
    double ParticlesPerStatWeightUnit = 0;
    double iRecord = 0;
    double iParticle = 0;
};

// ---

class SourceMaterialLimited : public SourceModeBase
{
public:
    SourceMaterialLimited(ParticleBase * particle,
                          TimeGeneratorBase * timeGenerator,
                          const G4ThreeVector & origin, const G4ThreeVector & boundingBoxFullSize,
                          const G4String & material,
                          G4String fileName_EmissionPositions = "");
    SourceMaterialLimited(const json11::Json & json);
    ~SourceMaterialLimited();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "SourceMaterialLimited";}

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

class SourceLysoNatural : public SourceModeBase
{
public:
    SourceLysoNatural(double timeFrom, double timeTo);
    SourceLysoNatural(const json11::Json & json);
    ~SourceLysoNatural();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "SourceLysoNatural";}

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

class SourceNa22Point : public SourceModeBase
{
public:
    SourceNa22Point(double timeFrom, double timeTo, const G4ThreeVector & origin);
    SourceNa22Point(const json11::Json & json);

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "SourceNa22Point";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    double TimeFrom;
    double TimeTo;
    G4ThreeVector Origin;
};

// ---

class SourcePointBlurred : public SourcePoint
{
public:
    SourcePointBlurred(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin, G4String distributionFileName, int numBins, double range);
    SourcePointBlurred(const json11::Json & json);
    ~SourcePointBlurred();

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "SourcePointBlurred";}

protected:
    void init();
    void doWriteToJson(json11::Json::object & json) const override;
    void doReadFromJson(const json11::Json & json);

    std::string FileName;
    int NumBins = 21;
    double Range = 10.0;
    Hist1DSampler * Sampler = nullptr;
};

// --

#endif // SourceMode_h
