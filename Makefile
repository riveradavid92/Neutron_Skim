#Makefile for gallery c++ programs.
#Note, being all-incllusive here: you can cut out libraries/includes you do not need
#you can also change the flags if you want too (Werror, pedantic, etc.)

CPPFLAGS=-I $(BOOST_INC) \
         -I $(CANVAS_INC) \
         -I $(CETLIB_INC) \
         -I $(FHICLCPP_INC) \
         -I $(GALLERY_INC) \
				 -I $(HEP_CONCURRENCY_INC) \
         -I $(LARCOREOBJ_INC) \
         -I $(LARDATAOBJ_INC) \
         -I $(LARSIM_INC) \
         -I $(NUSIMDATA_INC) \
	 -I $(CETLIB_EXCEPT_INC) \
	 -I $(CANVAS_ROOT_IO_INC) \
         -I $(ROOT_INC) \
	 -I $(FHICLCPP_INC) \
	 -I $(ART_INC) \
	 -I $(NUTOOLS_INC) \
	 -I $(MESSAGEFACILITY_INC) \
	 -I $(LARCOREALG_INC) \
	 -I $(LARDATAALG_INC) \
	 -I $(CLHEP_INC) \
	 -I $(DUNETPC_INC) 

#CXXFLAGS=-std=c++14 -Wall -Werror -pedantic
CXXFLAGS=-std=c++14 -Wall -pedantic
CXX=g++
LDFLAGS=$$(root-config --libs --cflags) \
        -L $(CETLIB_LIB) -l cetlib \
        -L $(GALLERY_LIB) -l gallery \
        -L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
        -L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
        -L $(LARDATAOBJ_LIB) -l lardataobj_Simulation -l lardataobj_RecoBase -l lardataobj_MCBase -l lardataobj_RawData -l lardataobj_OpticalDetectorData -l lardataobj_AnalysisBase \
	-L $(LARCOREALG_LIB) -l larcorealg_Geometry \
	-L $(LARDATAALG_LIB) -l lardataalg_DetectorInfo \
	-L $(FHICLCPP_LIB) -l fhiclcpp \
	-L $(LARSIM_LIB) -l larsim_MCCheater_BackTracker -l larsim_MCCheater_ParticleInventory \
	-L $(NUTOOLS_LIB) -l nutools_ParticleNavigation

BackTracker: BackTracker.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<	

PDSP_BT_Skim: PDSP_BT_Skim.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<	

Neutron: Neutron.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<	

PDSP_PartSkim: PDSP_PartSkim.cc  
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

PartSkim_v803: PartSkim_v803.cc 
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<
#
RadSkim_v803: RadSkim_v803.cc 
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<
#
#NeutronFinder_Skim: NeutronFinder_Skim.cc 
#	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<
#
NeutronFinder_Skim_v2: NeutronFinder_Skim_v2.cc 
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

Neutron_Skim_PDSP: Neutron_Skim_PDSP.cc 
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

Neutron_Skim_v2: Neutron_Skim_v2.cc 
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

Neutron_Skim_PDSP_data: Neutron_Skim_PDSP_data.cc 
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

Neutron_Skim_PDSP_sim: Neutron_Skim_PDSP_sim.cc 
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

Neutron_Simple_Skim_PDSP: Neutron_Simple_Skim_PDSP.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

PDSP_Simple_Skim_Dirs: PDSP_Simple_Skim_Dirs.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

#all: RadsSkim NeutronFinder_Skim NeutronFinder_Skim_v2 Neutron_Skim_PDSP
all: Neutron_Skim_PDSP Neutron_Skim_v2 Neutron_Skim_PDSP_sim Neutron_Skim_PDSP_data NeutronFinder_Skim_v2 PartSkim_v803 RadSkim_v803 Neutron PDSP_Simple_Skim_Dirs

clean:	
	rm -f Neutron_Skim_PDSP Neutron_Skim_v2 Neutron_Skim_PDSP_data Neutron_Skim_PDSP_sim NeutronFinder_Skim_v2 PartSkim_v803 RadSkim_v803 Neutron Neutron_Simple_Skim_PDSP PDSP_Simple_Skim_Dirs
