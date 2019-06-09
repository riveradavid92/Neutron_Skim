//////////////////////////////////////////
//Neutron subclass of simb::MCParticle
//D.Rivera
//////////////////////////////////////////
//some standard C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <memory>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystemDirectory.h"
#include "TString.h"
#include "TDirectory.h"

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

//"dunetpc" object includes
#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"

#include "Neutron.h"

//convenient for us! let's not bother with art and std namespaces!
using namespace art;
using namespace std;

int 
main (int argc, char* argv[]) 
{
  simb::MCParticle n0;
  auto n1 = std::make_unique<MCAna::Neutron>(n0);
  n1->Show();

  return 0;
}
