#ifndef SourceMode_h
#define SourceMode_h

#include "G4ThreeVector.hh"

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

protected:
    ParticleBase      * Particle      = nullptr;
    TimeGeneratorBase * TimeGenerator = nullptr;
    G4ParticleGun     * ParticleGun   = nullptr;

    G4ThreeVector Direction = {0, 0, 1.0};  // more general approach will be realized later!
    bool bIsotropicDirection = true;

    bool bGeneratePair  = false;    // for second gamma
    //bool bAcollinearity = false;

protected:
    G4ThreeVector generateDirectionIsotropic();
    void generateSecondGamma(G4Event * anEvent);

    virtual void customPostInit() {}
};

// ---

class PointSource : public SourceModeBase
{
public:
    PointSource(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin);
};

// ---

class PencilBeam : public SourceModeBase
{
public:
    PencilBeam(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin, const G4ThreeVector & direction);
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

protected:
    G4ThreeVector Origin;
    G4ThreeVector BoundingBox;
    G4String      Material;
    G4String      FileName;

    //run-time
    G4Material    * SourceMat = nullptr;
    G4Navigator   * Navigator = nullptr;
    std::ofstream * Stream = nullptr;

    void customPostInit() override;
};

// ---

class NaturalLysoSource : public SourceModeBase
{
public:
    NaturalLysoSource(double timeFrom, double timeTo);
    ~NaturalLysoSource();

    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    double ScintMaxRadius = 0;
    std::vector<std::pair<double,double>> ElectronSpectrum; //format: energy[keV] relative_probablility

    //run-time
    G4Navigator   * Navigator = nullptr;
    Hist1DSampler * Sampler   = nullptr;

    void customPostInit() override;
};

// ---

class BlurredPointSource : public PointSource
{
public:
    BlurredPointSource(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin, string FileName,const Hist1D & dist, long seed);
    ~BlurredPointSource();

    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    Hist1DSampler * Sampler   = nullptr;
    G4ThreeVector Origin;
};

#endif // SourceMode_h
