#include "Hist1D.hh"

#include <iostream>
#include <fstream>
#include <ios>
#include <algorithm>

Hist1DRegular::Hist1DRegular(int numBins, double from, double to) :
    NumBins(numBins), From(from), To(to)
{
    if (NumBins <= 0) NumBins = 1;

    Step = (To - From) / NumBins;

    Data.resize(NumBins);
    for (int iBin = 0; iBin < NumBins; iBin++) Data[iBin] = 0;
}

void Hist1DRegular::fill(double value, double weight)
{
    int index = (value - From) / Step;

    if      (index <  0)       NumUnderflows++;
    else if (index >= NumBins) NumOverflows++;
    else                       Data[index] += weight;
}

void Hist1DRegular::report()
{
    int sum = 0;
    for (int iBin = 0; iBin < NumBins; iBin++)
    {
        std::cout << (iBin == 0 ? '[' : ',');
        std::cout << Data[iBin];
        sum += Data[iBin];
    }
    std::cout << ']' << std::endl;

    std::cout << "Distribution sum: " << sum << std::endl;
    std::cout << "Underflows: " << NumUnderflows << std::endl;
    std::cout << "Overflows: "  << NumOverflows << std::endl;
}

void Hist1DRegular::save(const std::string & fileName)
{
    std::ofstream outStream;

     outStream.open(fileName);

    if (!outStream.is_open() || outStream.fail() || outStream.bad())
    {
        std::cout  << "Cannot open file to store output data: " << fileName << std::endl;
        return;
    }

    for (int iBin = 0; iBin < NumBins; iBin++)
        outStream << (From + Step * iBin) << " " << Data[iBin] << std::endl;

    outStream.close();
}

// ---

//Hist1DSampler::Hist1DSampler(const Hist1D & hist, long seed) : Hist(hist)
Hist1DSamplerRegular::Hist1DSamplerRegular(const Hist1DRegular & hist) : Hist(hist)
{
    //randEngine.seed(seed);

    double acc = 0;
    for (size_t iBin = 0; iBin < hist.Data.size(); iBin++)
    {
        Cumulative.push_back(SamplerRec(iBin, acc));
        acc += hist.Data[iBin];
    }

    if (acc != 0)
        for (SamplerRec & rec : Cumulative)
            rec.val /= acc;
}

#include "G4RandomTools.hh"
double Hist1DSamplerRegular::getRandom()
{
        //const double rndm = urd(randEngine);  // [0, 1)
    const double rndm = G4UniformRand(); // (0,1)
    auto res = std::upper_bound(Cumulative.begin(), Cumulative.end(), SamplerRec(0, rndm)); // iterator to the element with larger val than rndm

    size_t indexAbove = (res == Cumulative.end() ? Cumulative.size()-1
                                                 : res->index);

    const double dBin = (Hist.To - Hist.From) / Hist.NumBins;
        //return Hist.From + dBin * (indexAbove - urd(randEngine));
    return Hist.From + dBin * (indexAbove - G4UniformRand());
}

// ----------------- VARIABLE BIN SAMPLER FOR RADIAL ----------------------

#include "SessionManager.hh"
#include "out.hh"
RandomRadialSampler::RandomRadialSampler(const std::vector<std::pair<double, double>> & distribution) //: Distribution(distribution)
{
    if (distribution.size() < 2)
    {
        out("Distribution for RandomRadialSampler should have at least 2 points");
        exit(22);
    }
    if (distribution.front().first != 0)
    {
        out("First record for the distribution given to RandomRadialSampler should be for zero radial distance");
        exit(22);
    }

    const size_t sizeInterpolated = 100;
    const double step = distribution.back().first / sizeInterpolated;

    Distribution.reserve(sizeInterpolated);

    const size_t sizeOriginal = distribution.size();
    size_t indexInterpolated = 0;
    for (size_t indexOriginal = 0; indexOriginal < sizeOriginal-1; indexOriginal++)
    {
        double r = step * indexInterpolated;
        while (r <= distribution[indexOriginal+1].first)
        {
            const double interpolationFactor = (r - distribution[indexOriginal].first) / ( distribution[indexOriginal+1].first - distribution[indexOriginal].first );
            const double interpolatedValue   = SessionManager::interpolate(distribution[indexOriginal].second, distribution[indexOriginal+1].second, interpolationFactor);

            Distribution.push_back( {r, interpolatedValue} );

            indexInterpolated++;
            r = step * indexInterpolated;
        }
    }

    double acc = 0;
    const size_t size = Distribution.size();
    for (size_t iBin = 0; iBin < size; iBin++)
    {
        Cumulative.push_back(SamplerRec(iBin, acc));

        double areaFactor = 1.0;
        if (iBin != size-1) areaFactor = Distribution[iBin+1].first*Distribution[iBin+1].first - Distribution[iBin].first*Distribution[iBin].first;

        acc += Distribution[iBin].second * areaFactor;
    }

    if (acc != 0)
        for (SamplerRec & rec : Cumulative)
            rec.val /= acc;
}

double RandomRadialSampler::getRandom()
{
    const double rndm = G4UniformRand(); // (0,1)
    auto res = std::upper_bound(Cumulative.begin(), Cumulative.end(), SamplerRec(0, rndm)); // iterator to the element with larger val than rndm

    size_t indexAbove = (res == Cumulative.end() ? Cumulative.size()-1
                                                 : res->index);

    const double from = Distribution[indexAbove-1].first;
    const double to   = Distribution[indexAbove].first;
    return from + (to - from) * G4UniformRand();
}

void RandomRadialSampler::generatePosition(G4ThreeVector & pos)
{
    const double rndm = G4UniformRand(); // (0,1)
    auto res = std::upper_bound(Cumulative.begin(), Cumulative.end(), SamplerRec(0, rndm)); // iterator to the element with larger val than rndm

    size_t indexAbove = (res == Cumulative.end() ? Cumulative.size()-1
                                                 : res->index);

    pos[2] = 0;
    if (indexAbove != 1)
    {
        const double from = Distribution[indexAbove-1].first;
        const double to   = Distribution[indexAbove].first;

        const double radius = from + (to - from) * G4UniformRand();
        const double phi = 2.0 * 3.1415926535 * G4UniformRand();

        pos[0] = radius;
        pos[1] = 0;
        pos.rotateZ(phi);
    }
    else
    {
        const double r  = Distribution[1].first;
        const double r2 = Distribution[1].first * Distribution[1].first;
        do
        {
            pos[0] = -r + 2.0 * r * G4UniformRand();
            pos[1] = -r + 2.0 * r * G4UniformRand();
        }
        while (pos[0]*pos[0] + pos[1]*pos[1] > r2);
    }
}
