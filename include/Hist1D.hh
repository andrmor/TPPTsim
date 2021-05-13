#ifndef hist1d_h
#define hist1d_h

#include <vector>
#include <string>

class Hist1D
{
public:
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

#endif // hist1d_h
