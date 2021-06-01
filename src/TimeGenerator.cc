#include "TimeGenerator.hh"

#include "G4RandomTools.hh"

double UniformTime::generateTime()
{
    return TimeFrom + (TimeTo - TimeFrom) * G4UniformRand();
}

ExponentialTime::ExponentialTime(double timeFrom, double HalfLife) :
    TimeFrom(timeFrom), DecayTime(HalfLife/log(2)) {}

double ExponentialTime::generateTime()
{
    return TimeFrom + G4RandExponential::shoot(DecayTime);
}
