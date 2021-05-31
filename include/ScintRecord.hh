#ifndef scintrecord_h
#define scintrecord_h

#include <G4ThreeVector.hh>

struct ScintRecord
{
    G4ThreeVector CenterPos;
    G4ThreeVector FacePos;

    double Angle;

    void write(std::ofstream & stream);
};

#endif // scintrecord_h

