#ifndef SourceMode_h
#define SourceMode_h

#include "DefinedParticles.hh"

#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;
class G4Material;
class G4Navigator;

class SourceModeBase
{
public:
    SourceModeBase(ParticleBase * particle, int numPerEvent);
    ~SourceModeBase();

    void initialize();

    virtual void GeneratePrimaries(G4Event * anEvent) = 0;

protected:
    ParticleBase  * Particle    = nullptr;  // owns!

    int             NumPerEvent = 1;

    // Run-time resources
    G4ParticleGun * ParticleGun = nullptr;

protected:
    virtual void customPostInit() {}
};

// ---

class PointSource : public SourceModeBase
{
public:
    PointSource(ParticleBase * particle, const G4ThreeVector & origin, int numPerEvent);

    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    bool bSkipDirection = false;
    bool bGeneratePair  = false;
};

class PointSourceUniformTime : public PointSource
{
public:
    PointSourceUniformTime(ParticleBase * particle, const G4ThreeVector & origin, int numPerEvent, double timeWindow);

    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    double TimeWindow = 0;
};

// ---

class PointSourceExponentialTime : public PointSource
{
public:
    PointSourceExponentialTime(ParticleBase * particle, const G4ThreeVector & origin, int numPerEvent, double decayTime);

    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    double DecayTime = 0;
};

// ---

class PencilBeam : public SourceModeBase
{
public:
    PencilBeam(ParticleBase * particle, const G4ThreeVector & origin, const G4ThreeVector & direction, int numPerEvent);

    void GeneratePrimaries(G4Event * anEvent) override;
};

// ---

class MaterialLimitedSource : public SourceModeBase
{
public:
    MaterialLimitedSource(ParticleBase * particle,
                          const G4ThreeVector & origin, const G4ThreeVector & boundingBoxFullSize,
                          const G4String & material,
                          G4String fileName_EmissionPosition = "");
    ~MaterialLimitedSource();

    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    G4ThreeVector Origin;
    G4ThreeVector BoundingBox;
    G4String      Material;
    G4String      FileName;

    bool bSkipDirection = false;
    bool bGeneratePair  = false;

    //run-time
    G4Material    * SourceMat = nullptr;
    G4Navigator   * Navigator = nullptr;
    std::ofstream * Stream = nullptr;

    virtual void customPostInit();
};

#endif // SourceMode_h
