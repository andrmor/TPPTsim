#include "GammaPairFromAnnihilHist.hh"
#include "DefinedParticles.hh"
#include "TimeGenerator.hh"
#include "ActivityLoader.hh"
#include "out.hh"

#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"

GammaPairFromAnnihilHist::GammaPairFromAnnihilHist(const std::string & histogramFileName, double activityMultiplier) :
    SourceModeBase(new GammaPair(0.511*MeV), new UniformTime(0, TimeSpan*ns)),
    HistogramFileName(histogramFileName), ActivityMultiplier(activityMultiplier)
{
    init();
}

GammaPairFromAnnihilHist::GammaPairFromAnnihilHist(const json11::Json & json) :
    SourceModeBase(new GammaPair(0.511*MeV), new UniformTime(0, TimeSpan*ns))
{
    doReadFromJson(json);
    readFromJson(json);

    init();
}

void GammaPairFromAnnihilHist::init()
{
    ActivityLoader::load(HistogramFileName, ActivityData, Binning);
    if (ActivityData.empty())
    {
        out("Activity data is empty!");
        exit(111);
    }
}

void GammaPairFromAnnihilHist::GeneratePrimaries(G4Event * anEvent)
{
    if (CurrentIy >= ActivityData.front().size()) return;

    for (int ix = 0; ix < Binning.NumBins[0]; ix++)
        for (int iz = 0; iz < Binning.NumBins[2]; iz++)
        {
            const int number = ActivityData[ix][CurrentIy][iz] * ActivityMultiplier;
            if (number < 1) continue;

            for (int iPart = 0; iPart < number; iPart++)
            {
                const double x = Binning.Origin[0] + (0.5 + ix)        * Binning.BinSize[0];
                const double y = Binning.Origin[1] + (0.5 + CurrentIy) * Binning.BinSize[1];
                const double z = Binning.Origin[2] + (0.5 + iz)        * Binning.BinSize[2];
                ParticleGun->SetParticlePosition( {x,y,z} );

                SourceModeBase::GeneratePrimaries(anEvent);
            }
        }

    CurrentIy++;
}

double GammaPairFromAnnihilHist::CountEvents()
{
    return ActivityData.size();
}

void GammaPairFromAnnihilHist::doWriteToJson(json11::Json::object &json) const
{
    json["HistogramFileName"]  = HistogramFileName;
    json["ActivityMultiplier"] = ActivityMultiplier;
}

void GammaPairFromAnnihilHist::doReadFromJson(const json11::Json &json)
{
    jstools::readString(json, "HistogramFileName",  HistogramFileName);
    jstools::readDouble(json, "ActivityMultiplier", ActivityMultiplier);
}
