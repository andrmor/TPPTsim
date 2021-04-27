#include "TimeGenerator.hh"

#include "G4RandomTools.hh"

double UniformTime::generateTime()
{
    return TimeFrom + (TimeTo - TimeFrom) * G4UniformRand();
}

double ExponentialTime::generateTime()
{
    return TimeFrom + G4RandExponential::shoot(DecayTime);
}
