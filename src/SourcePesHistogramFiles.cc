#include "SourcePesHistogramFiles.hh"
#include "SessionManager.hh"
#include "ActivityLoader.hh"
#include "out.hh"
#include "jstools.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4RandomTools.hh"
#include "G4ThreeVector.hh"

#include <ios>
#include <iostream>
#include <fstream>

SourcePesHistogramFiles::SourcePesHistogramFiles(const std::string & directory, double multiplier, bool generateUniformOverBin) :
    SourceModeBase(nullptr, nullptr), Directory(directory), Multiplier(multiplier), GenerateUniformOverBin(generateUniformOverBin)
{
    init();
}

SourcePesHistogramFiles::SourcePesHistogramFiles(const json11::Json & json) :
    SourceModeBase(nullptr, nullptr)
{
    doReadFromJson(json);
    readFromJson(json);

    init();
}

void SourcePesHistogramFiles::init()
{
    IsotopeBase.push_back({ "C11",  {"6_12_C11.txt", "8_16_C11.txt"} });
    IsotopeBase.push_back({ "O15",  {"8_16_O15.txt"} });
    IsotopeBase.push_back({ "N13",  {"8_16_N13.txt"} });
    IsotopeBase.push_back({ "P31",  {"15_31_P30.txt"} });
    IsotopeBase.push_back({ "Ca40", {"20_40_K38.txt"} });

    checkInputData();

    CurrentRecord = 0;
    CurrentFile = 0;

    ParticleGun->SetParticleEnergy(0*keV);
    ParticleGun->SetParticleMomentumDirection({0,0,1.0});
}

void SourcePesHistogramFiles::checkInputData()
{
    for (auto it = IsotopeBase.begin(); it < IsotopeBase.end(); ++it)
    {
        for (auto itOther = IsotopeBase.begin(); itOther < IsotopeBase.end(); ++itOther)
        {
            if (it == itOther) continue;
            if (it->Isotope == itOther->Isotope)
            {
                out("Found non-unique record for PES:", it->Isotope);
                exit(3);
            }

            for (const std::string & fn : it->SpatialFiles)
            {
                std::string fullname = Directory + '/' + fn;
                std::ifstream infile(fullname);
                if (!infile.good())
                {
                    out("File not found:", fullname);
                    exit(3);
                }
            }
        }
    }
}

void SourcePesHistogramFiles::GeneratePrimaries(G4Event * anEvent)
{
    if (CurrentRecord >= IsotopeBase.size()) return; // finished with all events

    const PesDataRecord & pes = IsotopeBase[CurrentRecord];
    G4ParticleDefinition * pd = ParticleGenerator.makeGeant4Particle(pes.Isotope);
    ParticleGun->SetParticleDefinition(pd);

    const std::string & fileName = Directory + '/' + pes.SpatialFiles[CurrentFile];
    out("==>", fileName);

    BinningParameters binning;
    std::vector<std::vector<std::vector<double>>> data;
    ActivityLoader::load(fileName, data, binning);

    for (int ix = 0; ix < binning.NumBins[0]; ix++)
        for (int iy = 0; iy < binning.NumBins[1]; iy++)
            for (int iz = 0; iz < binning.NumBins[2]; iz++)
            {
                const int number = data[ix][iy][iz] * Multiplier;
                if (number < 1) continue;

                for (int iPart = 0; iPart < number; iPart++)
                {
                    double x, y, z;
                    if (GenerateUniformOverBin)
                    {
                        x = binning.Origin[0] + (G4UniformRand() + ix) * binning.BinSize[0];
                        y = binning.Origin[1] + (G4UniformRand() + iy) * binning.BinSize[1];
                        z = binning.Origin[2] + (G4UniformRand() + iz) * binning.BinSize[2];
                    }
                    else // bin center
                    {
                        x = binning.Origin[0] + (0.5 + ix) * binning.BinSize[0];
                        y = binning.Origin[1] + (0.5 + iy) * binning.BinSize[1];
                        z = binning.Origin[2] + (0.5 + iz) * binning.BinSize[2];
                    }

                    ParticleGun->SetParticlePosition( {x,y,z} );

                    ParticleGun->SetParticleTime(G4UniformRand() * TimeSpan);

                    ParticleGun->GeneratePrimaryVertex(anEvent);
                }
            }

    CurrentFile++;
    if (CurrentFile >= IsotopeBase[CurrentRecord].SpatialFiles.size())
    {
        CurrentRecord++;
        CurrentFile = 0;
    }
}

double SourcePesHistogramFiles::CountEvents()
{
    int numEvents = 0;
    for (const PesDataRecord & rec : IsotopeBase) numEvents += rec.SpatialFiles.size();
    out("Number of events is equal to the number of histograms:", numEvents);
    return numEvents;
}

void SourcePesHistogramFiles::doWriteToJson(json11::Json::object & json) const
{
    json["Directory"]              = Directory;
    json["Multiplier"]             = Multiplier;
    json["GenerateUniformOverBin"] = GenerateUniformOverBin;
}

void SourcePesHistogramFiles::doReadFromJson(const json11::Json & json)
{
    jstools::readString(json, "Directory",              Directory);
    jstools::readDouble(json, "Multiplier",             Multiplier);
    jstools::readBool  (json, "GenerateUniformOverBin", GenerateUniformOverBin);
}
