#ifndef definedparticles_h
#define definedparticles_h

#include <string>

#include "G4SystemOfUnits.hh"

class G4ParticleDefinition;

class ParticleBase
{
public:
    ParticleBase(const std::string & name, bool SkipDirection = false) : Name(name), bSkipDirection(SkipDirection) {}
    virtual ~ParticleBase(){}

    std::string Name = "undefined";
    bool        bSkipDirection = false;
    double      Energy = 1.0;

    virtual G4ParticleDefinition * getParticleDefinition() const = 0;
};

// ---

class ParticleGeantino : public ParticleBase
{
public:
    ParticleGeantino() : ParticleBase("geantino") {}
    G4ParticleDefinition * getParticleDefinition() const override;
};

class ParticleGamma : public ParticleBase
{
public:
    ParticleGamma(double energy = 0.511*MeV) : ParticleBase("gamma") {Energy = energy;}
    G4ParticleDefinition * getParticleDefinition() const override;
};

class ParticleGammaPair : public ParticleBase
{
public:
    ParticleGammaPair(double energy = 0.511*MeV) : ParticleBase("gammaPair") {Energy = energy;}
    G4ParticleDefinition * getParticleDefinition() const override;
};

// ---

class ParticlePES : public ParticleBase
{
public:
    ParticlePES(int z, int a);
    G4ParticleDefinition * getParticleDefinition() const override;

    int Z, A;
};

class ParticleC10 : public ParticlePES
{
public:
    ParticleC10() : ParticlePES(6, 10) {Name = "C10";}
};

class ParticleC11 : public ParticlePES
{
public:
    ParticleC11() : ParticlePES(6, 11) {Name = "C11";}
};

class ParticleO15 : public ParticlePES
{
public:
    ParticleO15() : ParticlePES(8, 15) {Name = "O15";}
};

class ParticleN13 : public ParticlePES
{
public:
    ParticleN13() : ParticlePES(7, 13) {Name = "N13";}
};

// ---

class ParticleProton : public ParticleBase
{
public:
    ParticleProton(double energy = 100.0*MeV) : ParticleBase("proton") {Energy = energy;}
    G4ParticleDefinition * getParticleDefinition() const override;
};

#endif // definedparticles_h
