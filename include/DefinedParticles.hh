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

class Geantino : public ParticleBase
{
public:
    Geantino() : ParticleBase("geantino") {}
    G4ParticleDefinition * getParticleDefinition() const override;
};

class Gamma : public ParticleBase
{
public:
    Gamma(double energy = 0.511*MeV) : ParticleBase("gamma") {Energy = energy;}
    G4ParticleDefinition * getParticleDefinition() const override;
};

class GammaPair : public ParticleBase
{
public:
    GammaPair(double energy = 0.511*MeV) //, bool acollinearity = false) :
        : ParticleBase("gammaPair")      //, bAcollineraity(acollinearity)
          {Energy = energy; }
    G4ParticleDefinition * getParticleDefinition() const override;

    //bool bAcollineraity = false;
};

// ---

class Isotope : public ParticleBase
{
public:
    Isotope(int z, int a);
    Isotope(int z, int a, double excitationEnergy);

    G4ParticleDefinition * getParticleDefinition() const override;

    int Z, A;
    double ExcitationEnergy = 0;
};

class C10 : public Isotope
{
public:
    C10() : Isotope(6, 10) {Name = "C10";}
};

class C11 : public Isotope
{
public:
    C11() : Isotope(6, 11) {Name = "C11";}
};

class O15 : public Isotope
{
public:
    O15() : Isotope(8, 15) {Name = "O15";}
};

class N13 : public Isotope
{
public:
    N13() : Isotope(7, 13) {Name = "N13";}
};

class N12 : public Isotope
{
public:
    N12() : Isotope(7, 12) {Name = "N12";}
};

// ---

class Lu176 : public Isotope
{
public:
    Lu176() : Isotope(72, 176, 596.820*keV){Name = "Hf176exc";}
};

// ---

class Proton : public ParticleBase
{
public:
    Proton(double energy = 100.0*MeV) : ParticleBase("proton") {Energy = energy;}
    G4ParticleDefinition * getParticleDefinition() const override;
};

#endif // definedparticles_h
