//some standard C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

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

//convenient for us! let's not bother with art and std namespaces!
using namespace art;
using namespace std;

//I like doing this to not get fooled by underflow/overflow
void ShowUnderOverFlow( TH1* h1){
  h1->SetBinContent(1, h1->GetBinContent(0)+h1->GetBinContent(1));
  h1->SetBinContent(0,0);

  int nbins = h1->GetNbinsX();
  h1->SetBinContent(nbins, h1->GetBinContent(nbins)+h1->GetBinContent(nbins+1));
  h1->SetBinContent(nbins+1,0);
 }

//----modified from simb::MCParticle::Trajectory class
unsigned char ProcessToKey(std::string const& process)
{
  int key = 0;
  
  if     (process.compare("hadElastic")            == 0) key = 1;
  else if(process.compare("pi-Inelastic")          == 0) key = 2;
  else if(process.compare("pi+Inelastic")          == 0) key = 3;
  else if(process.compare("kaon-Inelastic")        == 0) key = 4;
  else if(process.compare("kaon+Inelastic")        == 0) key = 5;
  else if(process.compare("protonInelastic")       == 0) key = 6;
  else if(process.compare("neutronInelastic")      == 0) key = 7;
  else if(process.compare("nCapture")              == 0) key = 8;
  else if(process.compare("nKiller")               == 0) key = 9;  
  else if(process.compare("FastScintillation")     == 0) key =10;
  else if(process.compare("CoupledTransportation") == 0) key =11;
  else if(process.compare("Transportation")        == 0) key =12;
  return key;
}

//----modified from simb::MCParticle::Trajectory class
std::string KeyToProcess(unsigned char const& key) 
{
  std::string process("Unknown");
  
  if     (key == 1) process  = "hadElastic";
  else if(key == 2) process  = "pi-Inelastic";
  else if(key == 3) process  = "pi+Inelastic";
  else if(key == 4) process  = "kaon-Inelastic";
  else if(key == 5) process  = "kaon+Inelastic";
  else if(key == 6) process  = "protonInelastic";
  else if(key == 7) process  = "neutronInelastic";
  else if(key == 8) process  = "nCapture";
  else if(key == 9) process  = "nKiller";
  else if(key == 10) process = "FastScintillation";
  else if(key == 11) process = "CoupledTransportation";
  else if(key == 12) process = "Transportation";
  return process;
}


int counter; 

int main(int argc, char* argv[])
{
  TString dirname1 = "/dune/data2/users/drivera/ProtoDUNE";
  TString output_name = "interaction_summary_PDSP";
  TString contains = "";
  int EvtsPerRun = 100;
  bool track_secondaries = false;
  bool min_info = true;

  if (argc < 3){
    cout << "Need directory and output name" << endl;
    return 1;
  }
  else if (argc == 3){
    dirname1 = argv[1];
    output_name = argv[2];
  }
  else if (argc == 4){
    dirname1 = argv[1];
    output_name = argv[2];
    contains = argv[3];
  }
  else if (argc == 5){
    dirname1 = argv[1];
    output_name = argv[2];
    contains = argv[3];
    EvtsPerRun = atoi(argv[4]);
  }
  else if (argc == 6){
    dirname1 = argv[1];
    output_name = argv[2];
    contains = argv[3];
    EvtsPerRun = atoi(argv[4]);
    track_secondaries = argv[5];
  }
  else if (argc == 7){
    dirname1 = argv[1];
    output_name = argv[2];
    contains = argv[3];
    EvtsPerRun = atoi(argv[4]);
    track_secondaries = argv[5];
    min_info = argv[6];
  }
  else {
    cout << "Too many arguments." << endl;
    return 1;
  }


  
  //std::vector<TString> tags = {"15","18"}; --no assumption on the presence of fasthit objects
  //TODO: add check for products and act on them if they are present
  std::vector<TString> tags = {};
  std::vector<std::vector<unsigned int>> hit_channel(tags.size());
  std::vector<std::vector<unsigned int>> hit_peak_time(tags.size());
  std::vector<std::vector<unsigned int>> hit_tot(tags.size());
  std::vector<std::vector<unsigned int>> hit_sum_adc(tags.size());

  //int truth_id;
  int r_it;
  int run;
  int event;
  //list of particles in each event
  std::vector<int> PDG;
  std::vector<bool> isPrimary;

  ///////////////////////////////////////////////PRIMARY INFO
  
  int processkey;    //process key for primary neutron
  int endprocesskey; //end process key for primary neutron
  int NumTrajPoints; //number of scatters for primary neutron

  double TrajLen; //length of primary trajectory
  //initial pos of primary neutron
  double vtx_x;
  double vtx_y;
  double vtx_z;
  //end position of primary neutron
  double End_X;
  double End_Y;
  double End_Z;
  double E0;  //initial kinetic energy of primary neutron
  double KEf; //initial kinetic energy of primary neutron

  std::vector<int> NScatter; //vector of scatter number to compare multiple particles as a function of scatter number
  std::vector<double> deltaE;
  std::vector<double> Frac_E_Loss; 

  //scatter locations for the primary neutron in each evt. 
  std::vector<double> scatter_X;
  std::vector<double> scatter_Y;
  std::vector<double> scatter_Z;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  std::string processname;         //<---The process which created this particle
  std::string endprocessname;      //<---The final process which this particle underwent

  ////////////////////////////////////////////////SECONDARY INFO

  //Secondary neutron info
  int Nsecondaries;                   //<---Number of secondary particles
  std::vector<int> secprocesskey;     //can have more than one secondary per primary
  std::vector<int> secendprocesskey;  //can have more than one secondary per primary
  std::vector<int> NumSecTrajPoints;  //non-primary neutron trajectory points

  std::vector<double> SecTrajLen; //non-primary neutron trajectory lengths
  std::vector<double> Sec_E0;
  std::vector<double> Sec_KEf;

  std::vector<std::vector<int>>    Sec_NScatter; //vector of scatter number vectors
  std::vector<std::vector<double>> secdeltaE; 
  std::vector<std::vector<double>> Sec_Frac_E_Loss; 

  //there can be more than 1 secondary neutron
  //scatter locations for the secondary neutrons
  std::vector<std::vector<double>> sec_scatter_X;
  std::vector<std::vector<double>> sec_scatter_Y;
  std::vector<std::vector<double>> sec_scatter_Z;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<std::string> secprocessname;      //<---The process which created this particle
  std::vector<std::string> secendprocessname;   //<---The final process which this particle underwent

  //double Vis_E;
  double OLD_vis_E;
  double truth_vis_E;
  double test_energy;
  double Missing_E;         //difference between initial KE of primary and visible energy
  double Missing_E_Frac;    //Missing_E divided by initial_KE
  double Missing_OldE_Frac; //calculated using the old def. of visible energy
  double Missing_testE;     //summing over all IDEs no additional checks involving the trackID
  double Missing_testE_Frac;

  //float E_ROI;
  std::vector<int> ID_vec;
  std::vector<sim::IDE> Sim_IDEs;
  simb::MCParticle neutron;
  simb::MCParticle sec_neutron;
  simb::MCTrajectory primTraj;
  simb::MCTrajectory secTraj;

  std::vector<pair<unsigned short, std::vector<sim::IDE, std::allocator<sim::IDE>>>> Sim_pair_vec;

  unsigned int add;  
  bool good = false;
  int APA_channel;
  
  //Histograms
  TH2F* h_Prim_dE_vs_NScatter = new TH2F("h_Prim_dE_vs_NScatter","Primary Neutron Energy Loss vs. Scatter Number;  NScatter; Fractional KE Loss[MeV]",200, 1-0.5, 400-0.5, 200, 0, 1);
  TH2F* h_Sec_dE_vs_NScatter  = new TH2F("h_Sec_dE_vs_NScatter", "Secondary Neutron Energy Loss vs. Scatter Number;  NScatter; Fractional KE Loss[MeV]",200, 1-0.5, 400-0.5, 200, 0, 1);
  TH2F* h_Prim_dE_visE_vs_KE0 = new TH2F("h_Prim_dE_visE_vs_KE0","Missing Energy Fraction vs. Inital KE; Initial_KE[MeV]; (Initial_KE-Vis_E)/Initial_KE",150,1-0.5,150-0.5, 200, 0, 1); 
  TH2F* h_Prim_dE_oldE_vs_KE0 = new TH2F("h_Prim_dE_oldE_vs_KE0","Missing Energy Fraction vs. Inital KE; Initial_KE[MeV]; (Initial_KE-OLD_E)/Initial_KE",150,1-0.5,150-0.5, 200, 0, 1); 
  TH2F* h_Prim_dE_testE_vs_KE0= new TH2F("h_Prim_dE_testE_vs_KE0","Missing Energy Fraction vs. Inital KE; Initial_KE[MeV]; (Initial_KE-test_E)/Initial_KE",150,1-0.5,150-0.5, 200, 0, 1);

  //We specify our files in a list of file names!
  //Note: multiple files allowed. Just separate by comma.

  TFile* output_file = new TFile(output_name+".root", "RECREATE");

  TString tag;

  //shared 
  TTree tree("tree","Skimmed_Neutrons");

  tree.Branch("Run",       &run);
  tree.Branch("Event_Num", &event);
 
  for (unsigned int i=0; i < tags.size(); ++i){
    tree.Branch("Hit_"+tags.at(i)+"_Channels",&hit_channel.at(i));
    tree.Branch("Hit_"+tags.at(i)+"_Peak_Time",&hit_peak_time.at(i));
    tree.Branch("Hit_"+tags.at(i)+"_TOT",&hit_tot.at(i));
    tree.Branch("Hit_"+tags.at(i)+"_Sum_ADC",&hit_sum_adc.at(i));
  }


  tree.Branch("PDG",               &PDG);
  tree.Branch("isPrimary",         &isPrimary);
  tree.Branch("OLD_Vis_E",         &OLD_vis_E);
  tree.Branch("Vis_E",             &truth_vis_E);
  tree.Branch("Test_E",            &test_energy);
  tree.Branch("Missing_E",         &Missing_E);
  tree.Branch("Missing_testE",     &Missing_testE);
  tree.Branch("Missing_E_Frac",    &Missing_E_Frac);
  tree.Branch("Missing_testE_Frac",&Missing_testE_Frac);
  tree.Branch("Missing_OldE_Frac", &Missing_OldE_Frac);

  tree.Branch("Initial_KE",        &E0); 
  tree.Branch("Final_KE",          &KEf); 
  tree.Branch("Process_Key",       &processkey);
  tree.Branch("Process_Name",      &processname);
  tree.Branch("EndProcess_Key",    &endprocesskey);
  tree.Branch("EndProcess_Name",   &endprocessname);
  tree.Branch("N_Traj_Points",     &NumTrajPoints);
  tree.Branch("Prim_Neut_TrajLen", &TrajLen);

  //primary start and end positions
  tree.Branch("vtx_x", &vtx_x);
  tree.Branch("vtx_y", &vtx_y);
  tree.Branch("vtx_z", &vtx_z);
  tree.Branch("End_X", &End_X);
  tree.Branch("End_Y", &End_Y);
  tree.Branch("End_Z", &End_Z);

  tree.Branch("NScatter",    &NScatter);
  tree.Branch("delta_E",     &deltaE);
  tree.Branch("Frac_E_Loss", &Frac_E_Loss);
  
  tree.Branch("scatter_X", &scatter_X);
  tree.Branch("scatter_Y", &scatter_Y);
  tree.Branch("scatter_Z", &scatter_Z);

  //Secondary
  tree.Branch("Nsecondaries",      &Nsecondaries);
  tree.Branch("secprocesskey",     &secprocesskey);
  tree.Branch("secendprocesskey",  &secendprocesskey);
  tree.Branch("secprocessname",    &secprocessname);
  tree.Branch("secendprocessname", &secendprocessname);
  tree.Branch("Sec_E0",            &Sec_E0); 
  tree.Branch("Sec_KEf",           &Sec_KEf); 
  tree.Branch("NumSecTrajPoints",  &NumSecTrajPoints);
  tree.Branch("SecTrajLen",        &SecTrajLen);
  tree.Branch("Sec_NScatter",      &Sec_NScatter);
  tree.Branch("secdeltaE",         &secdeltaE); //keeping populations separate
  tree.Branch("Sec_Frac_E_Loss",   &Sec_Frac_E_Loss);

  tree.Branch("Sec_Scatter_X", &sec_scatter_X);
  tree.Branch("Sec_Scatter_Y", &sec_scatter_Y);
  tree.Branch("Sec_Scatter_Z", &sec_scatter_Z);

  output_file->cd();

  vector<string> filename;
  vector<int> run_numbers; //holds list of process numbers to be used instead of run number which is always 1 for each job submission
  TSystemDirectory dir(dirname1, dirname1);
  TList *files = dir.GetListOfFiles();
  if (files){
    TSystemFile *file;
    TString fname, tmpname;
    TObjArray *substrings;
    int runnum;

    TIter next(files);
    while ((file=(TSystemFile*)next())){
      fname = file->GetName();
      if (!file->IsDirectory() && fname.Contains(contains) && fname.EndsWith(".root")){
        tmpname = dirname1 + "/" + fname;
        filename.push_back((string)tmpname.Data());
       
        //make list of process number to substitute for the run number in the aux event info
        substrings = fname.Tokenize("_");
        const TString runstring = ( (TObjString *)substrings->At(4) )->GetString();
        runnum = ( (TObjString *)runstring.Tokenize(".")->At(0) )->GetString().Atoi();
        cout << runnum << endl;
        run_numbers.push_back(runnum);
      }
    }
  }


  InputTag Truth_tag  { "generator" };
  InputTag Cosmic_tag { "cosmicgenerator" };
  InputTag Sim_tag    { "largeant"  };
  InputTag fast_hit_tag;
 
  int iterator = -1;
  //ok, now for the event loop! Here's how it works.
  //
  //gallery has these built-in iterator things.
  //
  //You declare an event with a list of file names. Then, you
  //move to the next event by using the "next()" function.
  //Do that until you are "atEnd()".
  //
  //In a for loop, that looks like this:
  r_it = 0;
  for (gallery::Event ev(filename) ; !ev.atEnd(); ev.next()) {
    //initialization on per-event-basis
    //Vis_E = 0;
    int ks = 0; //index for secondary tracks 
    
    //scalars
    OLD_vis_E    = 0;
    truth_vis_E  = 0;
    test_energy  = 0;
    Nsecondaries = 0;
    PDG.clear();
    isPrimary.clear();

    //primary vectors
    NScatter.clear();
    deltaE.clear();
    Frac_E_Loss.clear();

    scatter_X.clear();
    scatter_Y.clear();
    scatter_Z.clear();

    //secondary vectors
    secprocesskey.clear();
    secprocessname.clear();
    secendprocesskey.clear();
    secendprocessname.clear();
    NumSecTrajPoints.clear();
    SecTrajLen.clear();
    Sec_E0.clear();
    Sec_KEf.clear();
 
    Sec_NScatter.clear();
    secdeltaE.clear();
    Sec_Frac_E_Loss.clear();

    sec_scatter_X.clear();
    sec_scatter_Y.clear();
    sec_scatter_Z.clear();
    //

    ID_vec.clear();
    ++iterator;

    
    //Now, we want to get a "valid handle" (which is like a pointer to our collection")
    //We use auto, cause it's annoying to write out the fill type. But it's like
    //vector<recob::Hit>* object.
    
    run = ev.eventAuxiliary().run();
    event = ev.eventAuxiliary().event();
    run = run_numbers.at(int(r_it/EvtsPerRun));

    //to get run and event info, you use this "eventAuxillary()" object.
    cout << "" << endl;
    cout << "Processing "
	    << "Run "   << run   << ", "
	    << "Event " << event << endl;


    auto const& Truth_handle        = ev.getValidHandle<vector<simb::MCTruth>>(Truth_tag);
    auto const& Sim_handle          = ev.getValidHandle<vector<sim::SimChannel>>("elecDrift");
    auto const& Parts_handle        = ev.getValidHandle<vector<simb::MCParticle>>(Sim_tag);

    //From FindMany.h :
    //
    //*ProdB and Data are the only template arguments that must be specified when 
    //constructing a FindMany. Any other items are deducible from arguments.
    //
    //*The FindMany needs a source of objects of type A, a data container (e.g. an event)
    //and an input tag corresponding to the underlying association collection from which
    //to create itself
    //
    //*When constructed, the FindMany will obtain and interrogate the correct Assns and 
    //provide access to the B (and/or D) object(s) associated with each supplied A object
    //in the order in which the A objects were specified.
    //
    //*If the specified A does not have an associated B or D then the vector will be empty
    //
    //*If the required association collection has an extra data object D with each 
    //association then it *must* be specified as a template argument, even if it is not
    //relevant to the current query.
    //
    //Constructors.
    //
    //  // From Handle or ValidHandle to collection of A.
    //  FindMany<ProdB>(Handle<ProdAColl> const&,
    //  		DataContainer const&,
    //  		InputTag const&);
    //
    //  FindMany<ProdB, Data>(Handle<ProdAColl> const&,
    //                        DataContainer const&,
    //                        InputTag const&);
    //
    //the neutron samples have: 'art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo>' 
    //need to include the GeneratedParticleInfo, into the FindMany structure, art complains otherwise
    //
    ///FindMany<simb::MCParticle, sim::GeneratedParticleInfo> part_truth(Truth_handle,ev,Sim_tag);
    FindMany<simb::MCParticle> part_truth(Truth_handle,ev,Sim_tag);
    //associations between MCTruth and simb objects
    for (size_t i_part = 0, size_part = Truth_handle->size(); i_part != size_part; ++i_part){
      std::vector<simb::MCParticle const*> truth_vec;
      part_truth.get(i_part,truth_vec);
      for (unsigned int i=0; i < truth_vec.size(); ++i){
        ID_vec.push_back(truth_vec.at(i)->TrackId());
      }
    }

    auto const& Truth_vec(*Truth_handle);
    auto const& Sim_vec(*Sim_handle);
    auto const& Parts_vec(*Parts_handle);

    //loop over fasthit tags
    for (unsigned int j=0; j < tags.size(); ++j){
      tag = "fasthit"+tags.at(j);
      fast_hit_tag = {(string)tag.Data()};

      auto const& fast_hit_handle = ev.getValidHandle<vector<recob::Hit>>(fast_hit_tag);
      auto const& fast_hit_vec(*fast_hit_handle);

      for (unsigned int i=0; i < fast_hit_vec.size(); ++i){
        hit_channel.at(j).push_back(fast_hit_vec.at(i).Channel());
        hit_peak_time.at(j).push_back(fast_hit_vec.at(i).PeakTime());
        hit_tot.at(j).push_back((unsigned int)(fast_hit_vec.at(i).EndTick()-fast_hit_vec.at(i).StartTick()));
        hit_sum_adc.at(j).push_back(fast_hit_vec.at(i).SummedADC());
      }
    }

 
    //should always be code for a neutron 2112
    //truth_id = Truth_vec.at(0).GetParticle(0).PdgCode();
    double Mneutron = Truth_vec.at(0).GetParticle(0).Mass(); 
    E0 = (Parts_vec.at(0).E(0) - Mneutron)*1000; //starting KE [MeV] of primary
    double KEf = 0.0;

    //end position in the generation stage
    //is particle 0 the beam particle in PD?
    vtx_x = Truth_vec.at(0).GetParticle(0).EndX();
    vtx_y = Truth_vec.at(0).GetParticle(0).EndY();
    vtx_z = Truth_vec.at(0).GetParticle(0).EndZ();

    //end position in the geant propagation stage
    End_X = Parts_vec.at(0).EndX();
    End_Y = Parts_vec.at(0).EndY();
    End_Z = Parts_vec.at(0).EndZ();

    //auto const& Traj = Parts_vec.at(0).Trajectory();
    cout << "Size of MCParticle vector = " << Parts_vec.size() << endl;
    for (unsigned int i=0; i < Parts_vec.size(); ++i){
      int pid = Parts_vec.at(i).PdgCode();
      PDG.push_back(pid);

      if (i==0) { //primary neutron ->this does not work for anything other than single neutrons produced with particle gun

        isPrimary.push_back(true);
        neutron = Parts_vec.at(i); 

        std::string proc = neutron.Process();
        processname = proc;
        processkey  = (int)ProcessToKey(proc);

        NumTrajPoints = Parts_vec.at(i).NumberTrajectoryPoints();
        TrajLen = Parts_vec.at(i).Trajectory().TotalLength(); 
        
        cout << "==================Primary Neutron=================" << endl;
        cout << "Particle " << i << " : " << pid << endl;
        cout << Parts_vec.at(0) << endl;
        cout << "Number of Trajectory points for particle " << i << " = " << NumTrajPoints << endl; 
        cout << "Total lengh of primary neutron trajectory = " << TrajLen << endl;
        
        double prevE = Parts_vec.at(i).E(0);
        cout << "Starting Energy of Primary Neutron  = " << prevE*1000 << "MeV" << endl;
        cout << "Starting Kinetic Energy of Primary  = " << E0 << "MeV" << endl;
        
        //Loop over each trajectory and log changes in energy
        double del = 0.0;
        for (int j = 1; j < NumTrajPoints ; ++j) {
          //scatter position
          //primTraj = neutron.Trajectory();
          scatter_X.push_back(neutron.Vx(j));
          scatter_Y.push_back(neutron.Vy(j));
          scatter_Z.push_back(neutron.Vz(j));

          NScatter.push_back(j); //keep track of E loss as a function of scatter number
          del = prevE - Parts_vec.at(i).E(j); 
          deltaE.push_back(del*1000);
          KEf = (Parts_vec.at(i).E(j) - Mneutron)*1000.; //final kinetic energy for particle

          double frac = del*1000./E0;
          Frac_E_Loss.push_back(frac);
          
          cout << "Lost " << del*1000 << " MeV from " << j << "th scatter" << endl;
          prevE = Parts_vec.at(i).E(j);
          h_Prim_dE_vs_NScatter->Fill(j, frac);
        }

        std::string endproc = neutron.EndProcess();
        endprocessname = endproc;
        endprocesskey  = (int)ProcessToKey(endproc);

        //<--cout << "End Process for neutron: " << endproc << endl;
        cout << "End Process for neutron: " << endprocessname << endl;
        cout << "End Process Key: "         << endprocesskey  << endl;
        cout << "Final Kinetic Energy = "   << KEf            << endl;
        cout << "==================================================" << endl;

      } else if (track_secondaries) {

        isPrimary.push_back(false);
        Nsecondaries++;
        double KEi = 0.0;

        std::vector<int> snscatter;
        std::vector<double> ssecdeltaE, ssec_frac_E_loss; 
        std::vector<double> sscatter_x, sscatter_y, sscatter_z; //position vectors
        
        if (pid == 2112) { //secondary neutrons
          
          sec_neutron = Parts_vec.at(i);
          
          std::string proc = sec_neutron.Process();
          secprocessname.push_back(proc);
          secprocesskey.push_back( (int)ProcessToKey(proc) );

          int NumSecTrajPointsVal = Parts_vec.at(i).NumberTrajectoryPoints(); 
          double SecTrajLenVal = Parts_vec.at(i).Trajectory().TotalLength(); 
          NumSecTrajPoints.push_back(NumSecTrajPointsVal); 
          SecTrajLen.push_back(SecTrajLenVal); 

          cout << "--------------Secondary Neutron----------------" << endl;
          cout << "Particle " << i << " : " << pid << endl;
          cout << sec_neutron << endl;
          cout << "Number of Trajectory points for particle " << i << " = " << NumSecTrajPointsVal << endl;
          cout << "Total lengh of secondary neutron trajectory = " << SecTrajLenVal << endl;

          int process = 0;
          double prevE = Parts_vec.at(i).E(0);
          double sec_E0 = (Parts_vec.at(i).E(0)-Mneutron)*1000.; //staring kinetic energy [MeV] of secondary neutron
          Sec_E0.push_back(sec_E0);
          cout << "Starting Energy of Secondary Neutron = " << prevE*1000 << "MeV" << endl;
          cout << "Starting Kinetic Energy of Secondary = " << sec_E0 << "MeV" << endl;
          double del = 0.0;
          for (int j = 1; j < NumSecTrajPointsVal ; ++j) {
            sscatter_x.push_back(sec_neutron.Vx(j));
            sscatter_y.push_back(sec_neutron.Vy(j));
            sscatter_z.push_back(sec_neutron.Vz(j));

            snscatter.push_back(j);
            del = prevE - Parts_vec.at(i).E(j); 
            ssecdeltaE.push_back(del*1000);
            cout << "Lost " << del*1000 << " MeV from " << j << "th scatter" << endl;
            prevE = Parts_vec.at(i).E(j);
            KEf = (Parts_vec.at(i).E(j) - Mneutron)*1000.; //final kinetic energy for particle
            Sec_KEf.push_back(KEf);
   
            double frac = del*1000./sec_E0;
            ssec_frac_E_loss.push_back(frac);

            h_Sec_dE_vs_NScatter->Fill(j, frac);
          }

          std::string endproc = Parts_vec.at(i).EndProcess();
          secendprocessname.push_back(endproc);
          process = (int)ProcessToKey(endproc); 
          secendprocesskey.push_back( process ); 
          //secprocesskey.push_back(process);
          //secprocesskey.push_back(process);

          cout << "End Process for neutron: " << endproc   << endl;
          cout << "End Process Key: "         << process   << endl;
          cout << "Final Kinetic Energy = "   << KEf       << endl;
          cout << "-----------------------------------------------" << endl;
          ks++; //increment secondary index by 1

        } else if (min_info) { //basic info for all other particles including secondary neutrons when not tracking them specifically 
          Nsecondaries++;
          KEi = (Parts_vec.at(i).E() - Parts_vec.at(i).Mass())*1000; 
          cout << "..............Daughter Particle................" << endl;
          cout << "Summary info for particle type = " << pid << endl;
          cout << Parts_vec.at(i) << endl;
          cout << "Starting Kinetic Energy = " << KEi << "MeV" << endl;
          cout << "..............................................." << endl;
        } 

        //append to the global vector of 
        Sec_NScatter.push_back(snscatter);
        secdeltaE.push_back(ssecdeltaE);
        Sec_Frac_E_Loss.push_back(ssec_frac_E_loss);
        sec_scatter_X.push_back(sscatter_x);
        sec_scatter_Y.push_back(sscatter_y);
        sec_scatter_Z.push_back(sscatter_z);


      } else if (min_info) { //basic info for all other particles including secondary neutrons when not tracking them specifically
        
        isPrimary.push_back(false);
        double KEi = (Parts_vec.at(i).E() - Parts_vec.at(i).Mass())*1000; 
        cout << "..............Daughter Particle................" << endl;
        cout << "Summary info for particle type = " << pid << endl;
        cout << Parts_vec.at(i) << endl;
        cout << "Starting Kinetic Energy = " << KEi << "MeV" << endl;
        cout << "..............................................." << endl;

      } else {
          cout << "Not tracking any daughter particles" << endl;
      }

    }
    
    //loop over energy depositions by simulated tracks on readout channels
    for (unsigned int i=0; i < Sim_vec.size(); ++i){
      APA_channel = (int)Sim_vec.at(i).Channel()%2560; //128ch * 20FEMBs = 2560 channels per APA

      if (APA_channel < 1600) continue; //ch 1600 to 2559 are the collection channels

      Sim_IDEs = Sim_vec.at(i).TrackIDsAndEnergies(0,4500); //map of all trackID, Energy pairs between the start_tick and end_tick
                                                            //This method returns the energy deposited on this channel by each track ID active in the specified TDC time interval.
                                                            //Entries are sorted by track ID number.
      
      Sim_pair_vec = Sim_vec.at(i).TDCIDEMap(); // returns list of pairs of TDC values and sim::IDEs, (std::pair< unsigned short, std::vector<sim::IDE> >)
                                                // sorted by increasing TDC tick 

      //loop over ticks with energy depositions
      for (unsigned int j=0; j < Sim_pair_vec.size(); ++j){
        OLD_vis_E += Sim_vec.at(i).Energy(Sim_pair_vec.at(j).first);
      }

      //loop over all energy depositions in the full drift
      for (unsigned int j=0; j < Sim_IDEs.size(); ++j){
        good = false;
        test_energy += Sim_IDEs.at(j).energy; //add up all energy depositions without further checks

        //only add when the trackID matches associated track IDs
        for (unsigned int k=0; k < ID_vec.size(); ++k){
          if (Sim_IDEs.at(j).trackID == ID_vec.at(k)){
            good = true;
            break;
          }
        }
        if (good){
          truth_vis_E += Sim_IDEs.at(j).energy;
        }
      }
    }

    cout << "Test energy val = " << test_energy << endl;
    cout << "Vis  energy val = " << truth_vis_E << endl;

    
    Missing_E = (E0 - truth_vis_E);
    Missing_E_Frac = Missing_E/E0;
    h_Prim_dE_visE_vs_KE0->Fill(E0, Missing_E_Frac);

    Missing_OldE_Frac = (E0 - OLD_vis_E)/E0;
    h_Prim_dE_oldE_vs_KE0->Fill(E0, Missing_OldE_Frac);

    Missing_testE = (E0 - test_energy);                                                             
    Missing_testE_Frac = Missing_testE/E0;                                                          
    h_Prim_dE_testE_vs_KE0->Fill(E0, Missing_testE_Frac); 


    add = 0;
/*
    //fiducial cut
    //bigger than needs to be but maybe ok for capturing all energy deposits
    //if ( (fabs(vtx_x) > 360) | (fabs(vtx_y) > 600) | (vtx_z < 0 | vtx_z > 1396) ) continue;  
    //actual protodune geometry
    //bounding box: (-391.514,-72.3706,-48.9237) -- (463.526,717.869,806.116)
    if ( (fabs(vtx_x) > 360) | (vtx_y < 0) | (vtx_y > 700) | (vtx_z < 0) | (vtx_z > 1396) ) continue;  
    for (unsigned int i=0; i < Sim_vec.size(); ++i){
      APA_channel = (int)Sim_vec.at(i).Channel()%2560;
      if ( APA_channel < 1600 ) continue;
      Sim_pair_vec = Sim_vec.at(i).TDCIDEMap();
      for (unsigned int j=0; j < Sim_pair_vec.size(); ++j){
	      Vis_E += Sim_vec.at(i).Energy(Sim_pair_vec.at(j).first);
      }
    }
*/
    counter += add;

    ++r_it;
    tree.Fill();
    for (unsigned int i=0; i < tags.size(); ++i) {
      hit_channel.at(i).clear();
      hit_peak_time.at(i).clear();
      hit_tot.at(i).clear();
      hit_sum_adc.at(i).clear();
    } 


  } //end loop over events!
  //cout << "Writing" << endl;
  output_file->Write();
  //cout << "Closing" << endl;
  output_file->Close();

  filename.clear();


  //now, we're in a macro: we can just draw the histogram!
  //Let's make a TCanvas to draw our two histograms side-by-side
  TString outputfilename = output_name + "_neutron_hists_PDSP.root";
  TFile* MyFile = new TFile(outputfilename, "RECREATE");

  TCanvas* canvas1 = new TCanvas("canvas1","Primary",800,800);
  canvas1->cd();
  h_Prim_dE_vs_NScatter->SetDirectory(0);
  h_Prim_dE_vs_NScatter->Draw("COLZ");
  gDirectory->WriteObject(h_Prim_dE_vs_NScatter,"h_Prim_dE_vs_NScatter");

  TCanvas* canvas2 = new TCanvas("canvas2","Secondary",800,800);
  canvas2->cd();
  h_Sec_dE_vs_NScatter->SetDirectory(0);
  h_Sec_dE_vs_NScatter->Draw("COLZ");
  gDirectory->WriteObject(h_Sec_dE_vs_NScatter,"h_Sec_dE_vs_NScatter");

  TCanvas* canvas3 = new TCanvas("canvas3","Primary",800,800);
  canvas3->cd();
  h_Prim_dE_visE_vs_KE0->SetDirectory(0);
  h_Prim_dE_visE_vs_KE0->Draw("COLZ");
  gDirectory->WriteObject(h_Prim_dE_visE_vs_KE0,"h_Prim_dE_visE_vs_KE0");

  TCanvas* canvas4 = new TCanvas("canvas4","Primary Old Vis E",800,800);
  canvas4->cd();
  h_Prim_dE_oldE_vs_KE0->SetDirectory(0);
  h_Prim_dE_oldE_vs_KE0->Draw("COLZ");
  gDirectory->WriteObject(h_Prim_dE_oldE_vs_KE0,"h_Prim_dE_oldE_vs_KE0");

  TCanvas* canvas5 = new TCanvas("canvas5","Primary Old Vis E",800,800);                            
  canvas5->cd();                                                                                    
  h_Prim_dE_testE_vs_KE0->SetDirectory(0);                                                          
  h_Prim_dE_testE_vs_KE0->Draw("COLZ");                                                             
  gDirectory->WriteObject(h_Prim_dE_testE_vs_KE0,"h_Prim_dE_testE_vs_KE0");

  //TDirectory->WriteObject(h_dE_vs_NScatter,"h_dE_vs_NScatter");
  MyFile->Write();
  MyFile->Close();

  //h_E->Draw();
  //canvas2->SaveAs("h_dE_vs_NScatter.root");
  //canvas->cd(2);     //moves us to the second canvas
  //ShowUnderOverFlow(h_hit_peaktime); //use this function to move under/overflow into visible bins. 
  //h_hit_peaktime->Draw();
  //canvas->cd(3);     //moves us to the third canvas
  //ShowUnderOverFlow(h_hit_peakamp); //use this function to move under/overflow into visible bins. 
  //ShowUnderOverFlow(h_hit_integral); //use this function to move under/overflow into visible bins. 
  //h_hit_integral->SetLineColor(kRed);
  //h_hit_peakamp->SetLineColor(kBlue);
  //h_hit_peakamp->Draw();
  //h_hit_integral->Draw("same");
  

  //and ... done!
  return 0;
}
