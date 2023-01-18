#ifndef MaterialBuilder_H
#define MaterialBuilder_H

#include "json11.hh"

enum class EMaterial {PMMA, HDPE, PE, Graphite, GelTissue, GelWater, Bone, Brain, Blood, Muscle, Tissue};

class G4Material;

namespace MaterialBuilder
{
    G4Material * build(EMaterial material);

    void writeToJson(EMaterial material, json11::Json::object & json);
    void readFromJson(const json11::Json & json, EMaterial & material);
}

#endif // MaterialBuilder_H
