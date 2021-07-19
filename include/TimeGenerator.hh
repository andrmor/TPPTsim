#ifndef TimeGenerator_h
#define TimeGenerator_h

#include "json11.hh"

#include <string>

class Hist1D;
class TimeGeneratorBase
{
public:
    virtual ~TimeGeneratorBase(){}

    virtual double generateTime() = 0;

    virtual std::string getTypeName() const = 0;
    void writeToJson(json11::Json::object & json) const;

protected:
    virtual void doWriteToJson(json11::Json::object & json) const = 0;
};

// ---

class ConstantTime : public TimeGeneratorBase
{
public:
    ConstantTime(double time) : Time(time) {}

    double generateTime() override {return Time;}

    std::string getTypeName() const override {return "ConstantTime";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    double Time = 0;
};

// ---

class UniformTime : public TimeGeneratorBase
{
public:
    UniformTime(double timeFrom, double timeTo) : TimeFrom(timeFrom), TimeTo(timeTo) {}

    double generateTime() override;

    std::string getTypeName() const override {return "UniformTime";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    double TimeFrom = 0;
    double TimeTo   = 0;
};

// ---

class ExponentialTime : public TimeGeneratorBase
{
public:
    ExponentialTime(double timeFrom, double halfLife);
    ~ExponentialTime();

    double generateTime() override;
    std::string getTypeName() const override {return "ExponentialTime";}

protected:
    void doWriteToJson(json11::Json::object & json) const override;

    double TimeFrom  = 0;
    double HalfLife  = 0;
    double DecayTime = 0;

    Hist1D * Hist = nullptr;
};

#endif // TimeGenerator_h
