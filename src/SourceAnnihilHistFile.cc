#include "SourceAnnihilHistFile.hh"
#include "DefinedParticles.hh"
#include "TimeGenerator.hh"
#include "ActivityLoader.hh"
#include "out.hh"

#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomTools.hh"

SourceAnnihilHistFile::SourceAnnihilHistFile(const std::string & histogramFileName, double activityMultiplier, bool generateUniformOverBin,
                                             std::array<double, 3> offset) :
    SourceModeBase(new GammaPair(0.511*MeV), new UniformTime(0, TimeSpan*ns)),
    HistogramFileName(histogramFileName), ActivityMultiplier(activityMultiplier), GenerateUniformOverBin(generateUniformOverBin),
    Offset(offset)
{
    init();
}

SourceAnnihilHistFile::SourceAnnihilHistFile(const json11::Json & json) :
    SourceModeBase(new GammaPair(0.511*MeV), new UniformTime(0, TimeSpan*ns))
{
    doReadFromJson(json);
    readFromJson(json);

    init();
}

void SourceAnnihilHistFile::init()
{
    ActivityLoader::load(HistogramFileName, ActivityData, Binning);
    if (ActivityData.empty())
    {
        out("Activity data is empty!");
        exit(111);
    }
}

void SourceAnnihilHistFile::GeneratePrimaries(G4Event * anEvent)
{
    if (CurrentIy >= ActivityData.front().size()) return;

    for (int ix = 0; ix < Binning.NumBins[0]; ix++)
        for (int iz = 0; iz < Binning.NumBins[2]; iz++)
        {
            const int number = ActivityData[ix][CurrentIy][iz] * ActivityMultiplier;
            if (number < 1) continue;

            for (int iPart = 0; iPart < number; iPart++)
            {
                double x, y, z;
                if (GenerateUniformOverBin)
                {
                    x = Binning.Origin[0] + (G4UniformRand() + ix)        * Binning.BinSize[0];
                    y = Binning.Origin[1] + (G4UniformRand() + CurrentIy) * Binning.BinSize[1];
                    z = Binning.Origin[2] + (G4UniformRand() + iz)        * Binning.BinSize[2];
                }
                else // bin center
                {
                    x = Binning.Origin[0] + (0.5 + ix)        * Binning.BinSize[0];
                    y = Binning.Origin[1] + (0.5 + CurrentIy) * Binning.BinSize[1];
                    z = Binning.Origin[2] + (0.5 + iz)        * Binning.BinSize[2];
                }

                x += Offset[0];
                y += Offset[1];
                z += Offset[2];

                ParticleGun->SetParticlePosition( {x,y,z} );

                SourceModeBase::GeneratePrimaries(anEvent);
            }
        }

    CurrentIy++;
}

double SourceAnnihilHistFile::CountEvents()
{
    return ActivityData.size();
}

void SourceAnnihilHistFile::doWriteToJson(json11::Json::object & json) const
{
    json["HistogramFileName"]      = HistogramFileName;
    json["ActivityMultiplier"]     = ActivityMultiplier;
    json["GenerateUniformOverBin"] = GenerateUniformOverBin;

    json11::Json::array ar;
    for (size_t i = 0; i < 3; i++) ar.push_back(Offset[i]);
    json["Offset"] = ar;
}

void SourceAnnihilHistFile::doReadFromJson(const json11::Json &json)
{
    jstools::readString(json, "HistogramFileName",      HistogramFileName);
    jstools::readDouble(json, "ActivityMultiplier",     ActivityMultiplier);
    jstools::readBool  (json, "GenerateUniformOverBin", GenerateUniformOverBin);

    json11::Json::array ar;
    jstools::readArray(json, "Offset", ar);
    for (size_t i = 0; i < 3; i++) Offset[i] = ar[i].number_value();
}
