#include "TimeGenerator.hh"

#include "G4RandomTools.hh"

#include "Hist1D.hh"

double UniformTime::generateTime()
{
    return TimeFrom + (TimeTo - TimeFrom) * G4UniformRand();
}

ExponentialTime::ExponentialTime(double timeFrom, double decayTime) : TimeFrom(timeFrom), DecayTime(decayTime)
{
    Hist = new Hist1D(100,0,2e11);
}

ExponentialTime::~ExponentialTime()
{
    Hist->save("/data/margarida/Data/gammasExp.txt");
}

double ExponentialTime::generateTime()
{
    return TimeFrom + G4RandExponential::shoot(DecayTime);
    Hist->fill(TimeFrom + G4RandExponential::shoot(DecayTime));
}
