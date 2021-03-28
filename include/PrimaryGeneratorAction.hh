#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event * anEvent);

private:
    G4ParticleGun * fParticleGun = nullptr;
};

// to be split later: there will be a point source, a pencil beam and I guess something more advanced

#endif // PrimaryGeneratorAction_h
