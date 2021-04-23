#ifndef SourceMode_h
#define SourceMode_h

#include "DefinedParticles.hh" // !!!*** to .cc

#include "G4ThreeVector.hh"

class TimeGeneratorBase;
class G4ParticleGun;
class G4Event;
class G4Material;
class G4Navigator;

class SourceModeBase
{
public:
    SourceModeBase(ParticleBase * particle, TimeGeneratorBase * timeGenerator);
    ~SourceModeBase();

    void initialize();

    virtual void GeneratePrimaries(G4Event * anEvent) = 0;

protected:
    ParticleBase      * Particle      = nullptr; // owns!
    TimeGeneratorBase * TimeGenerator = nullptr; // owns!

    // Run-time own resources
    G4ParticleGun     * ParticleGun   = nullptr;

protected:
    G4ThreeVector generateDirectionIsotropic();

protected:
    virtual void customPostInit() {}
};

// ---

class PointSource : public SourceModeBase
{
public:
    PointSource(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin);

    void GeneratePrimaries(G4Event * anEvent) override;

protected:
    bool bSkipDirection = false;
    bool bGeneratePair  = false;
};

// ---

class PencilBeam : public SourceModeBase
{
public:
    PencilBeam(ParticleBase * particle, TimeGeneratorBase * timeGenerator, const G4ThreeVector & origin, const G4ThreeVector & direction);

    void GeneratePrimaries(G4Event * anEvent) override;
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

    bool bSkipDirection = false;
    bool bGeneratePair  = false;

    //run-time
    G4Material    * SourceMat = nullptr;
    G4Navigator   * Navigator = nullptr;
    std::ofstream * Stream = nullptr;

    virtual void customPostInit();
};

#endif // SourceMode_h
