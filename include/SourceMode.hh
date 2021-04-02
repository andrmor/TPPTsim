#ifndef SourceMode_h
#define SourceMode_h

#include "DefinedParticles.hh"

#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;

class SourceModeBase
{
public:
    SourceModeBase(ParticleBase * particle, int numPerEvent);
    ~SourceModeBase();

    void initialize();
    void setOrigin(const G4ThreeVector & origin);

    virtual void GeneratePrimaries(G4Event * anEvent) = 0;

protected:
    ParticleBase  * Particle    = nullptr;  // owns!

    int             NumPerEvent = 1;
    G4ThreeVector   Origin;

    // Run-time resources
    G4ParticleGun * ParticleGun = nullptr;
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



#endif // SourceMode_h
