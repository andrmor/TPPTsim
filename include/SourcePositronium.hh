#ifndef sourcepositronium_h
#define sourcepositronium_h

#include "SourceMode.hh"

class SourceThreeGammas : public SourceModeBase
{
public:
    SourceThreeGammas(TimeGeneratorBase * timeGenerator);
    //SourceThreeGammas(const json11::Json & json);

    void GeneratePrimaries(G4Event * anEvent) override;

    std::string getTypeName() const override {return "SourceThreeGammas";}

protected:
    void doWriteToJson(json11::Json::object & json) const override {}
    //void doReadFromJson(const json11::Json & json);

};

#endif // sourcepositronium_h
