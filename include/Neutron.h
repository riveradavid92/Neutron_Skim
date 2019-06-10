//////////////////////////////////////////
//Neutron subclass of simb::MCParticle
//D.Rivera
//////////////////////////////////////////
#ifndef NEUTRON_H
#define NEUTRON_H

//some standard C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <memory>

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/fwd.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTracker.h"


//convenient for us! let's not bother with art and std namespaces!
using namespace art;
using namespace std;

namespace MCAna {

  class Neutron: public simb::MCParticle {
    simb::MCParticle particle;
    std::unique_ptr <simb::MCParticle> Np;
    std::unique_ptr <simb::MCTrajectory> Nt;
    Double_t depositedEnergy;
    Double_t trackLength;
    unsigned int  nScatters;
    bool isContained;
    bool isValid;
  public:
    Neutron(simb::MCParticle N); //normal ctor
    Neutron(const Neutron &ob);  //copy ctor
    ~Neutron();                  //dtor
  
    void Show();
    unsigned int  NScatters()  { nScatters = Np->NumberTrajectoryPoints(); return nScatters; }
    Double_t DepositedEnergy(); // { return depositedEnergy; }
    Double_t TrackLength();     //{ return trackLength;     }
    bool Contained()           { return isContained;     }
  private: 
    const int NeutronPDG = 2112;
    bool ValidNeutron();
  }; //Neutron class


bool Neutron::
ValidNeutron()
{ 
  //return ( particle.PdgCode() == NeutronPDG ) ? true : false ;
  return ( Np->PdgCode() == NeutronPDG ) ? true : false ;
} 

Neutron::
~Neutron() 
{
  cout << "Neutron dtor called.\n";
}

Neutron::
Neutron( simb::MCParticle N ) : particle{}, 
                                Np{ std::make_unique<simb::MCParticle>(N) },
                                Nt{ std::make_unique<simb::MCTrajectory>(Np->Trajectory()) },
                                depositedEnergy{0.0},
                                trackLength{0.0},
                                nScatters{0},
                                isContained{false},
                                isValid{false}
{ 
  isValid = ValidNeutron(); //check for validity
  //if ( ValidNeutron() ) { //check for validity
  if (!isValid) { //check for validity
    cerr << "Not a valid Neutron!\n";
  } else {
    cout << "Neutron ctor called.\n";
  }
}

Neutron::
Neutron(const Neutron &cN) : particle{}, 
                             Np{ std::make_unique<simb::MCParticle>(*cN.Np) },
                             Nt{ std::make_unique<simb::MCTrajectory>(Np->Trajectory()) },
                             depositedEnergy{0.0},
                             trackLength{0.0},
                             nScatters{0},
                             isContained{false},
                             isValid{false}
{
  isValid = ValidNeutron();    //check for validity
  if (!isValid){
    cerr << "Not a valid Neutron!\n";
    //Np = NULL; //don't think I need this
  } else {
    cout << "Copy ctor called.\n";
  }
}

void Neutron::
Show()
{
  if(!isValid){
    cerr << "Cannot print info for invalid Neutron!\n";
  } else {
    cout << particle; 
  }
}

Double_t Neutron::
TrackLength()
{
  //for (unsigned int s=0; s<nScatters; s++) {
    //Nt = Np->Trajectory();
    return Nt->TotalLength();
  //}
}

} //namespace MCAna
#endif
