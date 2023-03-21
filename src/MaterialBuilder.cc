#include "MaterialBuilder.hh"
#include "jstools.hh"
#include "out.hh"

#include <string>

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

MaterialBuilder::MaterialBuilder(const json11::Json & json)
{
    readFromJson(json);
}

G4Material * MaterialBuilder::build()
{
    G4NistManager * man = G4NistManager::Instance();

    if ( StandardMaterial.empty() && (CustomMaterial == EMaterial::Undefined) )  // paranoic
    {
        out("Material builder was not configured!");
        exit(10);
    }

    if (CustomMaterial == EMaterial::Undefined)
    {
        G4Material * mat = man->FindOrBuildMaterial(StandardMaterial);
        if (!mat)
        {
            out("Unknown material name:", StandardMaterial);
            exit(10);
        }
        //out("-->Ionization potential for the phantom material:", mat->GetIonisation()->GetMeanExcitationEnergy()/eV, "eV");
        return mat;
    }

    G4Material * mat = nullptr;
    switch (CustomMaterial)
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
            std::vector<G4int> natoms;
            std::vector<G4String> elements;
            elements.push_back("C"); natoms.push_back(2);
            elements.push_back("H"); natoms.push_back(4);
            mat = man->ConstructNewMaterial("HDPE", elements, natoms, 0.95*g/cm3);
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
    default:;
    }

    if (!mat)
    {
        out("Error in material selection (MaterialBuilder::build)");
        exit(10);
    }
    //out("-->Ionization potential for the phantom material:", mat->GetIonisation()->GetMeanExcitationEnergy()/eV, "eV");
    return mat;
}

void MaterialBuilder::writeToJson(json11::Json::object & json)
{
    if ( StandardMaterial.empty() && (CustomMaterial == EMaterial::Undefined) )  // paranoic
    {
        out("Material builder was not configured!");
        exit(10);
    }

    std::string matStr;

    if (CustomMaterial == EMaterial::Undefined) matStr = StandardMaterial;
    else
    {
        switch (CustomMaterial)
        {
        case EMaterial::PMMA      : matStr = "PMMA";      break;
        case EMaterial::HDPE      : matStr = "HDPE";      break;
        case EMaterial::Graphite  : matStr = "Graphite";  break;

        case EMaterial::GelTissue : matStr = "GelTissue"; break;
        case EMaterial::GelWater  : matStr = "GelWater";  break;

        default:
            out("Not implemented material in MaterialBuilder::writeToJson");
            exit(11);
        }
    }

    json["Material"] = matStr;
}

void MaterialBuilder::readFromJson(const json11::Json & json)
{
    std::string matStr;
    jstools::readString(json, "Material", matStr);

    if (matStr.size() > 2)
    {
        std::string start = matStr;
        start.resize(3);
        if (start == "G4_")
        {
            StandardMaterial = matStr;
            CustomMaterial = EMaterial::Undefined;
            return;
        }
    }

    StandardMaterial.clear();
    if      (matStr == "PMMA")      CustomMaterial = EMaterial::PMMA;
    else if (matStr == "HDPE")      CustomMaterial = EMaterial::HDPE;
    else if (matStr == "Graphite")  CustomMaterial = EMaterial::Graphite;
    else if (matStr == "GelTissue") CustomMaterial = EMaterial::GelTissue;
    else if (matStr == "GelWater")  CustomMaterial = EMaterial::GelWater;

    else
    {
        out("Unknown material in MaterialBuilder::readFromJson");
        exit(11);
    }
}
