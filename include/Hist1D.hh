#ifndef hist1d_h
#define hist1d_h

#include <vector>
#include <string>

class Hist1DRegular
{
public:
    friend class Hist1DSamplerRegular;

    Hist1DRegular(int numBins, double from, double to);
    virtual ~Hist1DRegular(){}

    void fill(double value, double weight = 1.0);

    void report();
    void save(const std::string & fileName);

    double & operator[](size_t index) {return Data[index];}
    const double & operator[](size_t index) const {return Data[index];}

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
      size_t index;
      double val;

      SamplerRec(size_t Index, double Val) : index(Index), val(Val) {}
      SamplerRec(){}

      bool operator<(const SamplerRec & other) const {return val < other.val;}
};

//#include <random>
class Hist1DSamplerRegular
{
public:
    //Hist1DSampler(const Hist1D & hist, long seed);
    Hist1DSamplerRegular(const Hist1DRegular & hist);
    virtual ~Hist1DSamplerRegular(){}

    double getRandom();

protected:
    const Hist1DRegular Hist; // make a copy since original might be already destroyed when random s are generated
    std::vector<SamplerRec> Cumulative;

    //std::mt19937_64 randEngine;
    //std::uniform_real_distribution<double> urd;  // [0,1)

};

#include "G4ThreeVector.hh"
class RandomRadialSampler
{
public:
    RandomRadialSampler(const std::vector<std::pair<double,double>> & distribution);

    double getRandom();
    void   generatePosition(G4ThreeVector & pos);

protected:
    std::vector<std::pair<double,double>> Distribution;
    std::vector<SamplerRec> Cumulative;

};

#endif // hist1d_h
