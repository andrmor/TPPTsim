#include "MaterialBuilder.hh"
#include "jstools.hh"
#include "out.hh"

#include <string>

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

G4Material * MaterialBuilder::build(EMaterial material)
{
    G4NistManager * man = G4NistManager::Instance();
    G4Material * mat = nullptr;

    switch (material)
    {
    case EMaterial::PMMA :
        {
            std::vector<G4int> natoms;
            std::vector<G4String> elements;
            elements.push_back("C"); natoms.push_back(5);
            elements.push_back("H"); natoms.push_back(8);
            elements.push_back("O"); natoms.push_back(2);
            mat = man->ConstructNewMaterial("PMMA_phantom", elements, natoms, 1.18*g/cm3);
        }
        break;
    case EMaterial::HDPE :
        {
            std::vector<double> weightFrac;
            std::vector<G4String> elements;
            elements.push_back("H"); weightFrac.push_back(14.3);
            elements.push_back("C"); weightFrac.push_back(85.7);
            mat = man->ConstructNewMaterial("HDPE", elements, weightFrac, 0.95*g/cm3);
        }
        break;
    case EMaterial::PE :
        {
            std::vector<G4int> natoms;
            std::vector<G4String> elements;
            elements.push_back("C"); natoms.push_back(2);
            elements.push_back("H"); natoms.push_back(4);
            mat = man->ConstructNewMaterial("PE", elements, natoms, 0.96*g/cm3);
        }
        break;
    case EMaterial::Graphite :
        {
            std::vector<G4int> natoms;
            std::vector<G4String> elements;
            elements.push_back("C"); natoms.push_back(1);
            mat = man->ConstructNewMaterial("Graphite", elements, natoms, 1.83*g/cm3);
        }
        break;
    case EMaterial::GelTissue :
        {
            std::vector<double> weightFrac;
            std::vector<G4String> elements;
            elements.push_back("H"); weightFrac.push_back(9.6);
            elements.push_back("C"); weightFrac.push_back(14.9);
            elements.push_back("N"); weightFrac.push_back(1.46);
            elements.push_back("O"); weightFrac.push_back(73.8);
            mat = man->ConstructNewMaterial("GelTissue", elements, weightFrac, 1.13*g/cm3);
        }
        break;
    case EMaterial::GelWater :
        {
            std::vector<double> weightFrac;
            std::vector<G4String> elements;
            elements.push_back("H"); weightFrac.push_back(11.03);
            elements.push_back("C"); weightFrac.push_back(1.04);
            elements.push_back("N"); weightFrac.push_back(0.32);
            elements.push_back("O"); weightFrac.push_back(87.6);
            mat = man->ConstructNewMaterial("GelWater", elements, weightFrac, 1.01*g/cm3);
        }
        break;
    case EMaterial::WATER :
        mat = man->FindOrBuildMaterial("G4_WATER");
        break;
    case EMaterial::BONE_COMPACT_ICRU :
        mat = man->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
        break;
    case EMaterial::BRAIN_ICRP :
        mat = man->FindOrBuildMaterial("G4_BRAIN_ICRP");
        break;
    case EMaterial::BLOOD_ICRP :
        mat = man->FindOrBuildMaterial("G4_BLOOD_ICRP");
        break;
    case EMaterial::MUSCLE_SKELETAL_ICRP :
        mat = man->FindOrBuildMaterial("G4_MUSCLE_SKELETAL_ICRP");
        break;
    case EMaterial::TISSUE_SOFT_ICRP :
        mat = man->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
        break;
    default:;
    }

    if (!mat)
    {
        out("Error in material selection of PhantomCustomMatBox");
        exit(10);
    }
    //out("-->Ionization potential for the phantom material:", mat->GetIonisation()->GetMeanExcitationEnergy()/eV, "eV");
    return mat;
}

void MaterialBuilder::writeToJson(EMaterial material, json11::Json::object & json)
{
    std::string matStr;
    switch (material)
    {
    case EMaterial::PMMA      : matStr = "PMMA";      break;
    case EMaterial::HDPE      : matStr = "HDPE";      break;
    case EMaterial::PE        : matStr = "PE";        break;
    case EMaterial::Graphite  : matStr = "Graphite";  break;
    case EMaterial::WATER     : matStr = "WATER";     break;

    case EMaterial::GelTissue : matStr = "GelTissue"; break;
    case EMaterial::GelWater  : matStr = "GelWater";  break;

    case EMaterial::BONE_COMPACT_ICRU    : matStr = "BONE_COMPACT_ICRU";    break;
    case EMaterial::BRAIN_ICRP           : matStr = "BRAIN_ICRP";           break;
    case EMaterial::BLOOD_ICRP           : matStr = "BLOOD_ICRP";           break;
    case EMaterial::MUSCLE_SKELETAL_ICRP : matStr = "MUSCLE_SKELETAL_ICRP"; break;
    case EMaterial::TISSUE_SOFT_ICRP     : matStr = "TISSUE_SOFT_ICRP";     break;
    default:
        out("Not implemented material in MaterialBuilder::writeToJson");
        exit(11);
    }

    json["Material"] = matStr;
}

void MaterialBuilder::readFromJson(const json11::Json & json, EMaterial & material)
{
    std::string matStr;
    jstools::readString(json, "Material", matStr);

    if      (matStr == "PMMA")      material = EMaterial::PMMA;
    else if (matStr == "HDPE")      material = EMaterial::HDPE;
    else if (matStr == "PE")        material = EMaterial::PE;
    else if (matStr == "Graphite")  material = EMaterial::Graphite;
    else if (matStr == "WATER")     material = EMaterial::WATER;

    else if (matStr == "GelTissue") material = EMaterial::GelTissue;
    else if (matStr == "GelWater")  material = EMaterial::GelWater;

    else if (matStr == "BONE_COMPACT_ICRU")    material = EMaterial::BONE_COMPACT_ICRU;
    else if (matStr == "BRAIN_ICRP")           material = EMaterial::BRAIN_ICRP;
    else if (matStr == "BLOOD_ICRP")           material = EMaterial::BLOOD_ICRP;
    else if (matStr == "MUSCLE_SKELETAL_ICRP") material = EMaterial::MUSCLE_SKELETAL_ICRP;
    else if (matStr == "TISSUE_SOFT_ICRP")     material = EMaterial::TISSUE_SOFT_ICRP;

    else
    {
        out("Unknown material in MaterialBuilder::readFromJson");
        exit(11);
    }
}
