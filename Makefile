#Makefile for gallery c++ programs.
#Note, being all-incllusive here: you can cut out libraries/includes you do not need
#you can also change the flags if you want too (Werror, pedantic, etc.)

CPPFLAGS= -I $(BOOST_INC) \
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

####################################################################################################
# Main Skim
MAIN_SRC=./src/Neutron_Skim.cc
OUTDIR=./bin
OUTNAME=Neutron_Skim
OUT=$(OUTDIR)/$(OUTNAME)

####################################################################################################
# Test of Neutron class
TEST_SRC=./test/Neutron.cc
INCLUDEDIR=-I./include
TEST_DIR=./test
TEST_NAME=Neutron
TEST_OUT=$(TEST_DIR)/$(TEST_NAME)

####################################################################################################

#Neutron_Skim: Neutron_Skim.cc 
$(OUT): $(MAIN_SRC)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDEDIR) -o $@ $<

#Neutron: Neutron.cc
$(TEST_OUT): $(TEST_SRC)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDEDIR) -o $@ $<	

#all: Neutron_Skim Neutron 
all: $(OUT) $(TEST_OUT) 

clean:	
	rm -f $(OUT) $(TEST_OUT)
