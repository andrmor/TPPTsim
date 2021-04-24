#ifndef acollineargammamodel_h
#define acollineargammamodel_h

#include "G4VFastSimulationModel.hh"

class AcollinearGammaModel : public G4VFastSimulationModel
{
public:
    AcollinearGammaModel(const G4String & name, G4Region * region) :
        G4VFastSimulationModel(name, region) {}

    G4bool IsApplicable(const G4ParticleDefinition & particle) override;

    G4bool ModelTrigger(const G4FastTrack & fastTrack) override;

    void DoIt(const G4FastTrack & fastTrack, G4FastStep & step) override;

protected:
    int PrevID       = -1;
    int PrevParentID = -1;
};

#endif // acollineargammamodel_h
