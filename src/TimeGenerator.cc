#include "TimeGenerator.hh"

#include "G4RandomTools.hh"

#include "Hist1D.hh"

double UniformTime::generateTime()
{
    return TimeFrom + (TimeTo - TimeFrom) * G4UniformRand();
}

ExponentialTime::ExponentialTime(double timeFrom, double HalfLife) :
    TimeFrom(timeFrom), DecayTime(HalfLife/log(2.0))
{
    //Hist = new Hist1D(100,0,2e11);
}

ExponentialTime::~ExponentialTime()
{
    //Hist->save("/data/margarida/Data/gammasExp.txt");
    //delete Hist;
}

double ExponentialTime::generateTime()
{
    const double time = TimeFrom + G4RandExponential::shoot(DecayTime);
    //Hist->fill(time);
    return time;
}
