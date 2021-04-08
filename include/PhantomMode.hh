#ifndef PhantomMode_h
#define PhantomMode_h

class G4LogicalVolume;

class PhantomModeBase
{
public:
    virtual ~PhantomModeBase(){}

    virtual void definePhantom(G4LogicalVolume * logicWorld) {}
};

// ---

class PhantomModeNone : public PhantomModeBase
{
    //nothing to add :)
};

// ---

class PhantomModePMMA : public PhantomModeBase
{
    void definePhantom(G4LogicalVolume * logicWorld) override;
};

#endif // PhantomMode_h
