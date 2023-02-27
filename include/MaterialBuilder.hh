#ifndef MaterialBuilder_H
#define MaterialBuilder_H

#include "json11.hh"

enum class EMaterial {PMMA, HDPE, Graphite,
                      GelTissue, GelWater,
                      G4_WATER,
                      G4_BONE_COMPACT_ICRU,
                      G4_BRAIN_ICRP, G4_BLOOD_ICRP,
                      G4_MUSCLE_SKELETAL_ICRP, G4_TISSUE_SOFT_ICRP};

class G4Material;

namespace MaterialBuilder
{
    G4Material * build(EMaterial material);

    void writeToJson(EMaterial material, json11::Json::object & json);
    void readFromJson(const json11::Json & json, EMaterial & material);
}

#endif // MaterialBuilder_H
