#ifndef MaterialBuilder_H
#define MaterialBuilder_H

#include "json11.hh"

enum class EMaterial {PMMA, HDPE, PE, Graphite, GelTissue, GelWater, WATER,
                      BONE_COMPACT_ICRU, BRAIN_ICRP, BLOOD_ICRP, MUSCLE_SKELETAL_ICRP, TISSUE_SOFT_ICRP};

class G4Material;

namespace MaterialBuilder
{
    G4Material * build(EMaterial material);

    void writeToJson(EMaterial material, json11::Json::object & json);
    void readFromJson(const json11::Json & json, EMaterial & material);
}

#endif // MaterialBuilder_H
