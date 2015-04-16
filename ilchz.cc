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
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

// system
#include<string>
#include<vector>
#include<algorithm>
//#include<set>


// Helper class to persist data
struct ParticleKinRootAux
{
    std::vector<int> * pdgId;
    std::vector<int> * motherId;
    // FIXME:: remove grandmotherid, instead put in isLeading 
    // this information, change the name isLeading to auxvar
    // which now is a dynamic variable: isLeading when s-quarks,
    // grandmotherId when is hadrons form s
    std::vector<int> * grandmotherId;
    std::vector<int> * isLeading;
    std::vector<float> * p;
    std::vector<float> * pmother;
    std::vector<float> * phi;
    std::vector<float> * theta;

    // Auxiliary
    std::vector<int> * _listofusedpartindex;
    std::vector<std::vector<int> **> _auxI;
    std::vector<std::vector<float> **> _auxF;
    // Strange hadron PDG ID concerned by
    std::vector<int> strangehadrons;

    // Constructor
    ParticleKinRootAux(): 
        pdgId(0),
        motherId(0),
        grandmotherId(0),
        isLeading(0),
        p(0),
        pmother(0),
        phi(0),
        theta(0),
       _listofusedpartindex(0) 
    {
        _auxI.push_back(&pdgId);
        _auxI.push_back(&motherId);
        _auxI.push_back(&grandmotherId);
        _auxI.push_back(&isLeading);
        _auxF.push_back(&p);
        _auxF.push_back(&pmother);
        _auxF.push_back(&phi);
        _auxF.push_back(&theta);

        _auxI.push_back(&_listofusedpartindex);

        // Strange hadrons
        strangehadrons.push_back(130);  // K_Long (should I? too far away from the detector)
        strangehadrons.push_back(310);  // K_short
        strangehadrons.push_back(321);  // K+-
        strangehadrons.push_back(3122); // Lambda
    }

    ~ParticleKinRootAux()
    {
        endloop();
    }

    // used particles in the generated event: to speed up the loop not 
    // looking to particles already looked at and stored
    void usedparticle(const int & i)
    {
        _listofusedpartindex->push_back(i);
    }

    const std::vector<int> * getusedid()
    {
        return _listofusedpartindex;
    }

    // preloop initialization
    void initloop()
    {
        pdgId          = new std::vector<int>;
        motherId       = new std::vector<int>;
        grandmotherId  = new std::vector<int>;
        isLeading      = new std::vector<int>;
        p              = new std::vector<float>;
        pmother        = new std::vector<float>;
        phi            = new std::vector<float>;
        theta          = new std::vector<float>;

        _listofusedpartindex = new std::vector<int>;
    }
    
    // Delete before loop
    void endloop()
    {
        for(std::vector<std::vector<int>** >::const_iterator it=_auxI.begin(); it < _auxI.end(); ++it)
        {
            if( *(*it) != 0 )
            {
                delete *(*it);
                *(*it) = 0;
            }
        }
        for(std::vector<std::vector<float>** >::const_iterator it=_auxF.begin(); it < _auxF.end(); ++it)
        {
            if( *(*it) != 0 )
            {
                delete *(*it);
                *(*it) = 0;
            }
        }
    }
    void filltreevariables(const int & particleindex, const int & id, const int & motherid, const int & grmotherid, const int & leading, 
        const float & _p, const float & _pmother, const float & _phi, const float & _theta)
    {
        this->pdgId->push_back(id);
        this->motherId->push_back(motherid);
        this->grandmotherId->push_back(grmotherid);
        this->isLeading->push_back(leading);
        this->p->push_back(_p);
        this->pmother->push_back(_pmother);
        this->phi->push_back(_phi);
        this->theta->push_back(_theta);

        this->usedparticle(particleindex);
    }


    // Inititalization of the ttree
    void inittree(TTree * t)
    {
        gROOT->ProcessLine("#include <vector>");

        t->Branch("pdgId",&pdgId);
        t->Branch("motherId",&motherId);
        t->Branch("grandmotherId",&grandmotherId);
        t->Branch("isLeading",&isLeading);
        t->Branch("p",&p);
        t->Branch("pmother",&pmother);
        t->Branch("phi",&phi);
        t->Branch("theta",&theta);
    }
};

void fillstrangehadron(const Particle & quarkbeforerad, const Pythia & pythia, const RotBstMatrix & restframe, ParticleKinRootAux & p)
{
    // Note that the quark particle should be obtained throught the method iBotCopyId,
    // The bottom of the chain, just before hadronization (and after radiation)
    const Particle & quark = pythia.event[quarkbeforerad.iBotCopyId()];
    // Getting the original resonance
    int gmId = pythia.event[quarkbeforerad.mother1()].id();
    if(gmId == 0)
    {
        gmId = pythia.event[quarkbeforerad.mother2()].id();
    }
    // Getting the quark in its restframe:
    //Particle  quarkatrest(quark); // quark after radiation (bottom of the chain)
    Particle  quarkatrest(quarkbeforerad);  // quark BEFORE radiation (
    quarkatrest.rotbst(restframe);

    // Get list of quark-daughters which are final and are strange hadrons
    const int ndaughters = quark.daughterList().size(); 
    for(int k = 0; k < ndaughters; ++k)
    {
        const int currI = quark.daughterList()[k];
        const Particle & strhad = pythia.event[currI];
        if( !pythia.event[currI].isFinal() ||
                (std::find(p.strangehadrons.begin(),p.strangehadrons.end(),strhad.idAbs()) == p.strangehadrons.end()) )
        {
            continue;
        }
        Particle restframehadron(strhad);
        restframehadron.rotbst(restframe);
        p.filltreevariables(currI,strhad.id(),quark.id(),gmId,-1,
                restframehadron.pAbs(),quarkatrest.pAbs(),restframehadron.phi(),restframehadron.theta());
    }
}



void fillresonancechain(const int & i, const Pythia & pythia, ParticleKinRootAux & p)
{
    // Fill the resonance relative info: Lab frame
    // Note that we want to store info before showering,
    // directly from Hard-Scattering
    const int iHS = pythia.event[i].iTopCopy();
    const Particle & resonanceHS = pythia.event[iHS];
    const int resId = resonanceHS.id();
    p.filltreevariables(i,resId,0,0,-2,
            resonanceHS.pAbs(),-1,resonanceHS.phi(),resonanceHS.theta());

    // Before dealing with the daughters, recover the lowest copy (already 
    // radiated, therefore, the daughters are coming from a decay process)
    const Particle & resonance = pythia.event[i];
    // Getting the s-quarks
    const int s1 = resonance.daughter1();
    const int s2 = resonance.daughter2();

    // Note that we want to get the info of the first "carbon copy" or
    // same  quark (before the radiation steps)
    int iS    = s1;
    int iSbar = s2;
    if( pythia.event[s1].id() < 0)
    {
        iS    = s2;
        iSbar = s1;
    }
    const Particle & s    = pythia.event[iS];
    const Particle & sbar = pythia.event[iSbar];
    
    // Boost to rest-frame of the s-sbar system (pz is the only non-zero coordinate): 
    RotBstMatrix restframe;
    restframe.toCMframe(s.p(),sbar.p());

    // Save info of the s-quark
    const bool isLeading_s = (s.pAbs() > sbar.pAbs());
    // Note p in the rest-frame of the mother:
    //Particle srestframe(s);
    //srestframe.rotbst(restframe);
    //const float s_prf = srestframe.pAbs();

    p.filltreevariables(pythia.event[iS].iBotCopyId(),s.id(),resId,0,(int)isLeading_s,
            s.pAbs(),resonance.pAbs(),s.phi(),s.theta());
    
    // Save info of the sbar-quark
    // Note p in the rest-frame of the mother:
    //Particle sbarrestframe(s);
    //sbarrestframe.rotbst(restframe);
    //const float sbar_prf = sbarrestframe.pAbs();
    
    p.filltreevariables(iSbar,sbar.id(),resId,0,(int)(not isLeading_s),
            sbar.pAbs(),resonance.pAbs(),sbar.phi(),sbar.theta());
    
    // Track-down the kaons and lambdas from the strange quarks
    fillstrangehadron(s,pythia,restframe,p);
    fillstrangehadron(sbar,pythia,restframe,p);
}


int main(int argc, char* argv[]) 
{
    // Check that correct number of command-line arguments
    if(argc != 2) 
    {
        std::cerr << " Unexpected number of command-line arguments. \n You are"
            << " expected to provide one input file name. \n"
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
    ParticleKinRootAux particles;

    TTree * thz = new TTree("mctrue","ee -> H Z -> s sbar s sbar");
    particles.inittree(thz);

    // Some pseudo-constants initializations
    // resonances PDG IDs
    std::vector<int> idResonance;
    idResonance.push_back(23);  // Z0
    idResonance.push_back(25);  // Higgs (h0)
      
    std::vector<int> idPAbs;
    // resonances (note that the s-quarks are captured by the resonance
    // see "fillresonancechain" method)
    idPAbs.insert(idPAbs.end(),idResonance.begin(),idResonance.end());
    // Strange hadrons 
    idPAbs.insert(idPAbs.end(),particles.strangehadrons.begin(),particles.strangehadrons.end());
      

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
        // HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
        //ToHepMC.fill_next_event( pythia, hepmcevt );

        // Write the HepMC event to file. Done with it.
        //ascii_io << hepmcevt;
        //delete hepmcevt;

        // Allocate variables
        particles.initloop();

        // ROOT filling: loop over the event
        // Particles already checked
        //std::set<int> usedPart;

        for(int i= 0; i < pythia.event.size(); ++i)
        {
            const int pdgid = pythia.event[i].idAbs();
            // "Interesting" particles only!
            //usedPart.insert(currI);
            //if( (std::find(particles.getusedid()->begin(),particles.getusedid()->end(),i) == particles.getusedid()->end())
            if( std::find(idPAbs.begin(),idPAbs.end(),pdgid) == idPAbs.end() ) // is an "interesting" particle ?
            {
                continue;
            }

            // Obtain the last "carbon" copy of the particle to work with it
            const int currI = pythia.event[i].iBotCopy();
            // If was already used, don't duplicate
            if( (std::find(particles.getusedid()->begin(),particles.getusedid()->end(),currI) != particles.getusedid()->end() ) )
            {
                continue;
            }

            std::vector<int> partusedbynexts;
            // Check if it's a resonance, then fill the full chain (up to the s-quarks)
            if( std::find(idResonance.begin(),idResonance.end(),pdgid) != idResonance.end() )
            {
                fillresonancechain(currI,pythia,particles);
            }
            // If not, means that it is a strange hadron, but just want the final state
            else if( pythia.event[currI].isFinal() )
            {
                const Particle & had = pythia.event[currI]; 
                particles.filltreevariables(currI,had.id(),0,0,-1,
                        had.pAbs(),-1,had.phi(),had.theta());
            }
        }
        thz->Fill();
        // deallocate variables after the filling
        particles.endloop();
    }
    // End of event loop. Statistics.
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
