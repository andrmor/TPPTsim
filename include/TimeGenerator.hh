#ifndef TimeGenerator_h
#define TimeGenerator_h

class Hist1D;
class TimeGeneratorBase
{
public:
    virtual ~TimeGeneratorBase(){}

    virtual double generateTime() = 0;
};

// ---

class ConstantTime : public TimeGeneratorBase
{
public:
    ConstantTime(double time) : Time(time) {}

    double generateTime() override {return Time;}

protected:
    double Time = 0;
};

// ---

class UniformTime : public TimeGeneratorBase
{
public:
    UniformTime(double timeFrom, double timeTo) : TimeFrom(timeFrom), TimeTo(timeTo) {}

    double generateTime() override;

protected:
    double TimeFrom = 0;
    double TimeTo   = 0;
};

// ---

class ExponentialTime : public TimeGeneratorBase
{
public:
    ExponentialTime(double timeFrom, double decayTime);
    ~ExponentialTime();

    double generateTime() override;

protected:
    double TimeFrom  = 0;
    double DecayTime = 0;

    Hist1D    * Hist      = nullptr;
};

#endif // TimeGenerator_h
