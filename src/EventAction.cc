#include "EventAction.hh"
#include "G4Event.hh"
#include "out.hh"

void EventAction::BeginOfEventAction(const G4Event * event)
{
    int id = event->GetEventID();
    if (id % 10000 == 0) out("--- event #", id, "---");
}
