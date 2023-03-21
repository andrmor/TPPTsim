#ifndef MaterialBuilder_H
#define MaterialBuilder_H

#include "json11.hh"

enum class EMaterial {Undefined,
                      PMMA, HDPE, Graphite,
                      GelTissue, GelWater,
                      Ni400};

class G4Material;

class MaterialBuilder
{
public:
    MaterialBuilder(EMaterial customMaterial) : CustomMaterial(customMaterial) {}
    MaterialBuilder(std::string standardMaterial) : StandardMaterial(standardMaterial) {}

    MaterialBuilder(const json11::Json & json);

    G4Material * build();

    void writeToJson(json11::Json::object & json);
    void readFromJson(const json11::Json & json);

    std::string StandardMaterial;
    EMaterial   CustomMaterial = EMaterial::Undefined;
};

#endif // MaterialBuilder_H
