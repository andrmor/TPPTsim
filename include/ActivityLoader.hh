#ifndef activityloader_h
#define activityloader_h

#include "jstools.hh"

#include <string>
#include <vector>

namespace ActivityLoader
{
    static void load(const std::string & fileName, std::vector<std::vector<std::vector<double>>> & data, BinningParameters & binning);
};

#endif // activityloader_h
