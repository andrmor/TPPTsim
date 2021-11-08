#ifndef fastpesgeneratormodel_h
#define fastpesgeneratormodel_h

#include "G4VFastSimulationModel.hh"

class FastPesGeneratorModel : public G4VFastSimulationModel
{
public:
    FastPesGeneratorModel(const G4String & name, G4Region * region) :
        G4VFastSimulationModel(name, region) {}

    G4bool IsApplicable(const G4ParticleDefinition & particle) override;

    G4bool ModelTrigger(const G4FastTrack & fastTrack) override;

    void DoIt(const G4FastTrack & fastTrack, G4FastStep & step) override;

private:
    double LastEnergy;
    double LastTrackLength;
    G4ThreeVector LastPosition;

    double Lambda = 20.0;

    // tmp
    double sumLength = 0;
    int    num = 0;

};

#endif // fastpesgeneratormodel_h
