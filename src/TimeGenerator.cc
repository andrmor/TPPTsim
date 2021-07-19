#include "TimeGenerator.hh"

#include "G4RandomTools.hh"

#include "Hist1D.hh"

void TimeGeneratorBase::writeToJson(json11::Json::object & json) const
{
    json["Type"] = getTypeName();
    doWriteToJson(json);
}

// ---

void ConstantTime::doWriteToJson(json11::Json::object & json) const
{
    json["Time"] = Time;
}

// ---

double UniformTime::generateTime()
{
    return TimeFrom + (TimeTo - TimeFrom) * G4UniformRand();
}

void UniformTime::doWriteToJson(json11::Json::object &json) const
{
    json["TimeFrom"] = TimeFrom;
    json["TimeTo"]   = TimeTo;
}

// ---

ExponentialTime::ExponentialTime(double timeFrom, double halfLife) :
    TimeFrom(timeFrom), HalfLife(halfLife),
    DecayTime(halfLife/log(2.0))
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

void ExponentialTime::doWriteToJson(json11::Json::object &json) const
{
    json["TimeFrom"] = TimeFrom;
    json["HalfLife"] = HalfLife;
}
