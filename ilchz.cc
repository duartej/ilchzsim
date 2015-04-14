// main42.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Jordi Duarte-Campderros (TAU), jorge.duarte.campderros@cern.ch.
// 
// Simulation of ILC process: ee -> ZH -> q qbar s sbar
//
// Input and output files are specified on the command line, e.g. like
// ./main42.exe main42.cmnd hepmcout42.dat > out
// The main program contains no analysis; this is intended to happen later.
// It therefore "never" has to be recompiled to handle different tasks.

//#include "Pythia8/Pythia.h"
//#include "Pythia8Plugins/HepMC2.h"

//using namespace Pythia8;
#include "ilchz.h"

// ROOT-related
#include "TFile.h"
#include "TTree.h"

// system
#include<string>

struct kin
{
    int   pdgId;
    int   motherId;
    int   isLeading;
    float p;
    float pt;
    float phi;
    float theta;
};

TTree * treeinit(kin & h0, kin & s, kin & sbar, kin & z0, kin & s_z, kin & sbar_z)
{
    const std::string leaflist("pdgId/I:motherId/I:isLeading/I:p/F:pt/F:phi/F:theta/F");

    TTree * thz = new TTree("mctrue","s sbar decayed from Higgs");
    thz->Branch("h0",&h0,leaflist.c_str());
    thz->Branch("s_h",&s,leaflist.c_str());
    thz->Branch("sbar_h",&sbar,leaflist.c_str());
    thz->Branch("z0",&z0,leaflist.c_str());
    thz->Branch("s_z",&s_z,leaflist.c_str());
    thz->Branch("sbar_z",&sbar_z,leaflist.c_str());

  return thz;
}

void filltreevars(const Pythia & pythia, const int & i, kin & m, kin & p, kin & pbar)
{
    const int s1 = pythia.event[i].daughter1();
    const int s2 = pythia.event[i].daughter2();

    int iS = s1;
    int iSbar = s2;
    if( pythia.event[s1].id() < 0)
    {
        iS = s2;
        iSbar = s1;
    }
    // Fill kinematics 
    // resonance
    m.pdgId = pythia.event[i].id();
    m.motherId = 0;
    m.isLeading = -1;
    m.p  =pythia.event[i].pAbs();
    m.pt =pythia.event[i].pT();
    m.phi=pythia.event[i].phi();
    m.theta=pythia.event[i].theta();
    // q (s)
    p.pdgId=pythia.event[iS].id();
    p.motherId=m.pdgId;
    p.isLeading = pythia.event[iS].pAbs() > pythia.event[iSbar].pAbs();
    p.p  =pythia.event[iS].pAbs();
    p.pt =pythia.event[iS].pT();
    p.phi=pythia.event[iS].phi();
    p.theta=pythia.event[iS].theta();
    // sbar
    pbar.pdgId=pythia.event[iSbar].id();
    pbar.motherId=p.motherId;
    pbar.isLeading = not p.isLeading;
    pbar.p  =pythia.event[iSbar].pAbs();
    pbar.pt =pythia.event[iSbar].pT();
    pbar.phi=pythia.event[iSbar].phi();
    pbar.theta=pythia.event[iSbar].theta();
}

int main(int argc, char* argv[]) 
{
  // Check that correct number of command-line arguments
  if(argc != 3) 
  {
      std::cerr << " Unexpected number of command-line arguments. \n You are"
          << " expected to provide one input and one output file name. \n"
          << " Program stopped! " << std::endl;
      return 1;
  }

  // Check that the provided input name corresponds to an existing file.
  std::ifstream is(argv[1]);
  if(!is) 
  {
      std::cerr << " Command-line file " << argv[1] << " was not found. \n"
          << " Program stopped! " << std::endl;
      return 1;
  }

  // Confirm that external files will be used for input and output.
  std::cout << "\n >>> PYTHIA settings will be read from file " << argv[1]
      << " <<< \n >>> HepMC events will be written to file "
      << argv[2] << " <<< \n" << std::endl;

  // Interface for conversion from Pythia8::Event to HepMC event.
  //HepMC::Pythia8ToHepMC ToHepMC;

  // Specify file where HepMC events will be stored.
  //HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);

  // Generator.
  Pythia pythia;

  // Read in commands from external file.
  pythia.readFile(argv[1]);

  // Extract settings to be used in the main program.
  int    nEvent    = pythia.mode("Main:numberOfEvents");
  int    nAbort    = pythia.mode("Main:timesAllowErrors");

  // Initialization.
  pythia.init();

  // ROOT init tree
  kin h0,s,sbar;
  kin z0,s_z,sbar_z;

  TTree * thz = treeinit(h0,s,sbar,z0,s_z,sbar_z);

  // Begin event loop.
  int iAbort = 0;
  for(int iEvent = 0; iEvent < nEvent; ++iEvent) 
  {
      // Generate event.
      if (!pythia.next()) 
      {
          // If failure because reached end of file then exit event loop.
          if (pythia.info.atEndOfFile()) 
          {
              std::cout << " Aborted since reached end of Les Houches Event File\n";
              break;
          }
          
          // First few failures write off as "acceptable" errors, then quit.
          ++iAbort;
          if(iAbort < nAbort) 
          {
              continue;
          }
          std::cout << " Event generation aborted prematurely, owing to error!\n";
          break;
    }

    // Construct new empty HepMC event and fill it.
    // Units will be as chosen for HepMC build, but can be changed
    // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
    //HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    //ToHepMC.fill_next_event( pythia, hepmcevt );

    // Write the HepMC event to file. Done with it.
    //ascii_io << hepmcevt;
    //delete hepmcevt;

    // ROOT filling: loop over the event (only final state particles
    std::map<int,int> iB;
    iB[25] = 0; // Higgs
    iB[23] = 0; // Z
    int gotit = 0;

    for(int i= 0; i < pythia.event.size(); ++i)
    {
        //if(!pythia.event[i].isFinal())
        //{
        //    continue;
        //}
        // Finding s, sbar coming from h0
        const int pdgid = pythia.event[i].idAbs();
        if(iB.find(pdgid) == iB.end())
        {
            continue;
        }
        else
        {
            // Obtain the last Z0 "carbon" copy
            iB[pdgid] = pythia.event[i].iBotCopyId();
            ++gotit;
        }

        if( gotit > 1 )
        {
            break;
        }
    }
    // Now we have the Higgs and its daughters before decay and so on
    filltreevars(pythia,iB[25],h0,s,sbar); 
    // Now we have the Z and its daughters before decay and so on
    filltreevars(pythia,iB[23],z0,s_z,sbar_z); 

    thz->Fill();
  // End of event loop. Statistics.
  }
  pythia.stat();

  // ROOT File
  TFile *froot = new TFile("hzkin.root","RECREATE");
  thz->Write();
  froot->Close();

  delete thz;
  delete froot;

  // Done.
  return 0;
}
