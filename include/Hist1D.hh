#ifndef hist1d_h
#define hist1d_h

#include <vector>
#include <string>
#include <random>

class Hist1D
{
public:
    friend class Hist1DSampler;

    Hist1D(int numBins, double from, double to);
    virtual ~Hist1D(){}

    void fill(double value, double weight = 1.0);

    void report();
    void save(const std::string & fileName);

protected:
    int NumBins = 100;
    double From = 0;
    double To   = 100.0;

    double Step = 1.0;

    std::vector<double> Data;
    int NumUnderflows = 0;
    int NumOverflows  = 0;
};

// ---

struct SamplerRec
{
      double x;
      double val;

      SamplerRec(double X, double Val) : x(X), val(Val) {}
      SamplerRec(){}
};

class Hist1DSampler
{
public:
    Hist1DSampler(const Hist1D & hist, long seed);
    virtual ~Hist1DSampler(){}

    double getRandom();

protected:
    std::vector<SamplerRec> Cumulative;

    std::mt19937_64 randEngine;
    std::uniform_real_distribution<double> urd;  // [0,1)

};

#endif // hist1d_h
