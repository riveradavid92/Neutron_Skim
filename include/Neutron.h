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

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

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
    bool Contained();           //{ return isContained;     }
  private: 
    const int NeutronPDG = 2112;
    bool ValidNeutron();
    double ActiveBounds[6]; // Cryostat boundaries ( neg x, pos x, neg y, pos y, neg z, pos z )
    double fFidVolCut = 0.2;

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
    //cout << particle << endl; 
    cout << *Np << endl; 
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

bool Neutron::
Contained()
{
  const double posX = Np->EndX();
  const double posY = Np->EndY();
  const double posZ = Np->EndZ();  

  art::ServiceHandle<geo::Geometry> geom;                                                           
  double vtx[3] = {posX, posY, posZ};                                                               
  bool inside = false;                                                                              
                                                                                                    
  geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);                                                  
                                                                                                    
  if (geom->HasTPC(idtpc))                                                                          
  {                                                                                                 
    const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);                                            
    double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();                                       
    double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();                                       
    double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();                                       
                                                                                                    
    for (size_t c = 0; c < geom->Ncryostats(); c++)                                                 
    {                                                                                               
      const geo::CryostatGeo& cryostat = geom->Cryostat(c);                                         
      for (size_t t = 0; t < cryostat.NTPC(); t++)                                                  
      {                                                                                             
        const geo::TPCGeo& tpcg = cryostat.TPC(t);                                                  
        if (tpcg.MinX() < minx) minx = tpcg.MinX();                                                 
        if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX();                                                 
        if (tpcg.MinY() < miny) miny = tpcg.MinY();                                                 
        if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();                                                 
        if (tpcg.MinZ() < minz) minz = tpcg.MinZ();                                                 
        if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();                                                 
      }                                                                                             
    }                                                                                               
                                                                                                    
                                                                                                    
    //x                                                                                             
    double dista = fabs(minx - posX);                                                               
    double distb = fabs(posX - maxx);                                                               
    if ((posX > minx) && (posX < maxx) &&                                                           
      (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;                                  
                                                                                                    
    //y                                                                                             
    dista = fabs(maxy - posY);                                                                      
    distb = fabs(posY - miny);                                                                      
    if (inside && (posY > miny) && (posY < maxy) &&                                                 
      (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;                                  
    else inside = false;                                                                            
                                                                                                    
    //z                                                                                             
    dista = fabs(maxz - posZ);                                                                      
    distb = fabs(posZ - minz);                                                                      
    if (inside && (posZ > minz) && (posZ < maxz) &&                                                 
      (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;                                  
    else inside = false;                                                                            
  }                                                                                                 
                                                                                                    
  return inside;                                                                                    
}                                                       
  //auto const* geom = lar::providerFrom<geo::Geometry>();
  //auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
/*
  art::ServiceHandle<geo::Geometry> geom;
  //auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // Build my Cryostat boundaries array...Taken from Tyler Alion in Geometry Core.
  ActiveBounds[0] = ActiveBounds[2] = ActiveBounds[4] = DBL_MAX;
  ActiveBounds[1] = ActiveBounds[3] = ActiveBounds[5] = -DBL_MAX;
  // assume single cryostats
  //auto const* geom = lar::providerFrom<geo::Geometry>();
  for (geo::TPCGeo const& TPC: geom->IterateTPCs()) {
    // get center in world coordinates
    double origin[3] = {0.};
    double center[3] = {0.};
    TPC.LocalToWorld(origin, center);
    double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
 
    if( center[0] - tpcDim[0] < ActiveBounds[0] ) ActiveBounds[0] = center[0] - tpcDim[0];
    if( center[0] + tpcDim[0] > ActiveBounds[1] ) ActiveBounds[1] = center[0] + tpcDim[0];
    if( center[1] - tpcDim[1] < ActiveBounds[2] ) ActiveBounds[2] = center[1] - tpcDim[1];
    if( center[1] + tpcDim[1] > ActiveBounds[3] ) ActiveBounds[3] = center[1] + tpcDim[1];
    if( center[2] - tpcDim[2] < ActiveBounds[4] ) ActiveBounds[4] = center[2] - tpcDim[2];
    if( center[2] + tpcDim[2] > ActiveBounds[5] ) ActiveBounds[5] = center[2] + tpcDim[2];
  } // for all TPC
  std::cout << "Active Boundaries: "
            << "\n\tx: " << ActiveBounds[0] << " to " << ActiveBounds[1]
            << "\n\ty: " << ActiveBounds[2] << " to " << ActiveBounds[3]
            << "\n\tz: " << ActiveBounds[4] << " to " << ActiveBounds[5]
            << std::endl;

  return false;
}
*/
/*
auto const* geom = lar::providerFrom<geo::Geometry>();                                            
auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();                    
                                                                                                  
//compute the drift x range                                                                       
double vDrift = detprop->DriftVelocity()*1e-3; //cm/ns                                            
double xrange[2] = {DBL_MAX, -DBL_MAX };                                                          
for (unsigned int c=0; c<geom->Ncryostats(); ++c) {                                               
  for (unsigned int t=0; t<geom->NTPC(c); ++t) {                                                  
    double Xat0 = detprop->ConvertTicksToX(0,0,t,c);                                              
    double XatT = detprop->ConvertTicksToX(detprop->NumberTimeSamples(),0,t,c);                   
    xrange[0] = std::min(std::min(Xat0, XatT), xrange[0]);                                        
    xrange[1] = std::max(std::max(Xat0, XatT), xrange[1]);                                        
  }                                                                                               
}                                                                                                 
                                                                                                  
double result = 0.;                                                                               
TVector3 disp;                                                                                    
bool first = true;                                                                                
                                                                                                  
for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) {                                    
  // check if the particle is inside a TPC                                                        
  if (p.Vx(i) >= ActiveBounds[0] && p.Vx(i) <= ActiveBounds[1] &&                                 
      p.Vy(i) >= ActiveBounds[2] && p.Vy(i) <= ActiveBounds[3] &&                                 
      p.Vz(i) >= ActiveBounds[4] && p.Vz(i) <= ActiveBounds[5]){                                  
    // Doing some manual shifting to account for                                                  
    // an interaction not occuring with the beam dump                                             
    // we will reconstruct an x distance different from                                           
    // where the particle actually passed to to the time                                          
    // being different from in-spill interactions                                                 
    double newX = p.Vx(i)+(p.T(i)*vDrift);                                                        
    if (newX < xrange[0] || newX > xrange[1]) continue;                                           
    TLorentzVector pos(newX,p.Vy(i),p.Vz(i),p.T());                                               
    if(first){                                                                                    
      start = pos;                                                                                
      starti=i;                                                                                   
      first = false;                                                                              
    }                                                                                             
    else {                                                                                        
      disp -= pos.Vect();                                                                         
      result += disp.Mag();                                                                       
    }                                                                                             
    disp = pos.Vect();                                                                            
    end = pos;                                                                                    
    endi = i;                                                                                     
  }                                                                                               
}                                                                                                 
return result;                                                                             
*/

} //namespace MCAna
#endif
