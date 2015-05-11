// Simulation of ILC process: ee -> ZH -> q qbar s sbar
// 
// Author: Jordi Duarte-Campderros (TAU), jorge.duarte.campderros@cern.ch.
// (based on the main42.cc which is a part of the PYTHIA 
// event generator. // Copyright (C) 2015 Torbjorn Sjostrand.)
// PYTHIA is licenced under the GNU GPL version 2.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

//
// Input and output files are specified on the command line, e.g. like
// ./ilchz ilchz.cmnd > out.txt
// The main program contains the generation of a n-tuple (ROOT)
// for further analysis. 

#include "ilchz.h"

// ROOT-related
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

// system
#include<string>
#include<vector>
#include<algorithm>
#include<utility>
//#include<set>


// Helper class to declare hadrons-id
struct FinalStateHadrons
{
    const int K_LONG   =  130;
    const int K_SHORT  =  310;
    const int K_0      =  311;
    const int K_PLUS   =  321;
    const int K_MINUS  = -321;
    const int PI_0     =  111;
    const int PI_PLUS  =  211;
    const int PI_MINUS = -211;

    std::string _selectedtype;
    std::vector<int> _selected;

    // Constructor
    FinalStateHadrons(const std::string & hadrontype)
    {
        if(hadrontype == "kaons")
        {
            _selectedtype = "kaons";
            _selected = getkaons();
        }
        else if(hadrontype == "pions")
        {
            _selectedtype = "pions";
            _selected = getpions();
        }
        else
        {
            _selectedtype = "NONE";
            std::cerr << "FinalStateHadrons() Initialization ERROR:"
                << " Hadron type '" << hadrontype << "' not valid" << std::endl;
        }
    }

    std::vector<int> getpions() const
    {
        std::vector<int> _pions;
        _pions.push_back(PI_0);
        _pions.push_back(PI_PLUS);

        return _pions;
    }
    
    std::vector<int> getkaons() const
    {
        std::vector<int> _kaons;
        _kaons.push_back(K_LONG);
        _kaons.push_back(K_SHORT);
        _kaons.push_back(K_PLUS);

        return _kaons;
    }

    const std::vector<int> getIDs() const
    {
        return _selected;
    }


    const std::string gettype() const
    {
        return _selectedtype;
    }
};


// Helper class to persist data 
struct ParticleKinRootAux
{
    std::vector<int> * pdgId;
    std::vector<int> * motherindex; //n-tuple index: -1 for the strange hadrons case
    // Change meaning depending of the type of the particle:
    //    catchall = 0                  for the resonance (H,Z)
    //    catchall = isLeading (higher p) for the s-squark resonance daughters
    //    catchall = grandmother PDG_ID for the strange hadrons promptly decayed from the resonance
    //    catchall = grandmother PDG_ID for the strange hadrons not prompt
    std::vector<int> * catchall;
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
    ParticleKinRootAux(const std::vector<int> & finalstatehadrons): 
        pdgId(0),
        motherindex(0),
        catchall(0),
        p(0),
        pmother(0),
        phi(0),
        theta(0),
       _listofusedpartindex(0) 
    {
        _auxI.push_back(&pdgId);
        _auxI.push_back(&motherindex);
        _auxI.push_back(&catchall);
        _auxF.push_back(&p);
        _auxF.push_back(&pmother);
        _auxF.push_back(&phi);
        _auxF.push_back(&theta);

        _auxI.push_back(&_listofusedpartindex);

        // User define the hadrons
        strangehadrons.insert(strangehadrons.end(),finalstatehadrons.begin(),
                finalstatehadrons.end());
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
        motherindex    = new std::vector<int>;
        catchall       = new std::vector<int>;
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
    int filltreevariables(const int & pythiaindex, const int & id, const int & _motherindex, 
            const int & _catchall, const float & _p, const float & _pmother, const float & _phi, 
            const float & _theta)
    {
        this->pdgId->push_back(id);
        this->motherindex->push_back(_motherindex);
        this->catchall->push_back(_catchall);
        this->p->push_back(_p);
        this->pmother->push_back(_pmother);
        this->phi->push_back(_phi);
        this->theta->push_back(_theta);

        this->usedparticle(pythiaindex);

        // The current index on the n-tuple for this 
        // particle
        return (this->pdgId->size()-1);
    }

    // Convert from pythia-index to ParticleRootAux index
    int getlocalindex(const int & pythiaindex)
    {
        auto localindex = std::find(_listofusedpartindex->begin(),_listofusedpartindex->end(),
                pythiaindex);
        if(localindex == _listofusedpartindex->end())
        {
            return -1;
        }

        return *localindex;
    }

    // Convert from ParticleRootAux-local index to pythia index
    int getpythiaindex(const int & localindex)
    {
        return this->_listofusedpartindex->at(localindex);
    }

    // Inititalization of the ttree
    void inittree(TTree * t)
    {
        gROOT->ProcessLine("#include <vector>");

        t->Branch("pdgId",&pdgId);
        t->Branch("motherindex",&motherindex);
        t->Branch("catchall",&catchall);
        t->Branch("p",&p);
        t->Branch("pmother",&pmother);
        t->Branch("phi",&phi);
        t->Branch("theta",&theta);
    }
};


/*void fillstrangehadron(const Particle & quarkbeforerad, const int & quarkIndex, 
        const Pythia & pythia, const RotBstMatrix & pre_restframe, ParticleKinRootAux & p)
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
    quarkatrest.rotbst(pre_restframe);
    
    RotBstMatrix restframe(pre_restframe);
    // Checking the quark is defined in the positive axis
    if( quarkatrest.pz() < 0.0 )
    {
        restframe.rot(M_PI);
    }

    // List of pythia indices for all hadrons (in final state) of the types
    // defined at p.strangehadrons
    std::vector<int> finalhadrons; 
std::cout << "[" << (pythia.info.getCounter(3)-1) << "] Chain for " << quark.name() << " quark (" << quark.index() << ") -only final: " ;
    getfinaldaughters(quark.index(),pythia,p.strangehadrons,finalhadrons);
    for(auto & ipythia: finalhadrons)
    {
        // Not duplicate
        if( (std::find(p.getusedid()->begin(),p.getusedid()->end(),ipythia) != p.getusedid()->end() ) )
        {
            continue;
        }
        const Particle & strhad = pythia.event[ipythia];
 std::cout << strhad.name() << " (" << strhad.index() << ") ";
        Particle restframehadron(strhad);
        restframehadron.rotbst(restframe);
        p.filltreevariables(ipythia,strhad.id(),quarkIndex,gmId,
                restframehadron.pAbs(),quarkatrest.pAbs(),restframehadron.phi(),restframehadron.theta());
    }
    // Get list of quark-daughters which are final and are strange hadrons
    *const int ndaughters = quark.daughterList().size();  
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
        p.filltreevariables(currI,strhad.id(),quarkIndex,gmId,
                restframehadron.pAbs(),quarkatrest.pAbs(),restframehadron.phi(),restframehadron.theta());
    }*
}*/

std::pair<int,RotBstMatrix> fillresonancechain(const int & i, const Pythia & pythia, ParticleKinRootAux & p)
{
    // Fill the resonance relative info: Lab frame
    // Note that we want to store info before showering,
    // directly from Hard-Scattering
    const int iHS = pythia.event[i].iTopCopy();
    const Particle & resonanceHS = pythia.event[iHS];
    const int respdgId = resonanceHS.id();
    // Filling n-tuple and getting the n-tuple index of the Resonance
    const int reslocalIndex = p.filltreevariables(i,respdgId,-1,0,
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

    //const int squarkIndex = 
    p.filltreevariables(pythia.event[iS].iBotCopyId(),
            s.id(),reslocalIndex,(int)isLeading_s,
            s.pAbs(),resonance.pAbs(),s.phi(),s.theta());
    
    // Save info of the sbar-quark
    //const int squarkbarIndex = 
    p.filltreevariables(pythia.event[iSbar].iBotCopy(),
            sbar.id(),reslocalIndex,(int)(not isLeading_s),
            sbar.pAbs(),resonance.pAbs(),sbar.phi(),sbar.theta());
    
    // Returning the PDG ID of the resonance and the rest-frame system of 
    // its quarks daughters
    return std::pair<int,RotBstMatrix>(respdgId,restframe);
}

// Get the pythia index of the ancestor particle (which is defined by its PDG ID in consideredmums)
// and the pythia index of the daughter of the ancestor (the daughter quark, before radiation)
int getancestorindex(const int & currIndex, const Pythia & pythia,const std::vector<int> & consideredmums, int & quarkindex)
{
    // Find the top copy, before calling the mother list
    const int & index = pythia.event[currIndex].iTopCopyId();
    const Particle & hadronbeforerad = pythia.event[index];

    // Check the mother list and found the resonance Id:
//#ifdef DEBUG
    std::cout << "--> " << hadronbeforerad.name() << "(" << index << ")"; 
//#endif
    const std::vector<int> & mums = hadronbeforerad.motherList();
    for(unsigned int k = 0; k < mums.size(); ++k)
    {
//#ifdef DEBUG
    std::cout << "--> " << pythia.event[mums[k]].name() << "(" << mums[k] << ")"; 
//#endif
        if( std::find(consideredmums.begin(),consideredmums.end(),pythia.event[mums[k]].id()) != consideredmums.end() )
        {
            return mums[k];
        }
        else if( mums[k] == 0 )
        {
            return 0;
        }
        
        // updating the daugther index before the recursive call
        quarkindex = mums[k];
        return getancestorindex(mums[k],pythia,consideredmums,quarkindex);
    }
    
    quarkindex = -1;
    return 0;
}

void display_usage()
{
    std::cout << "\033[37musage:\033[m ilchz [OPTIONS] ilchz.cmnd [> out.txt]"
        << std::endl;
	std::cout << std::endl;
    std::cout << "Simulate the generation of the e+ e- --> H0 Z0 --> s sbar s sbar" 
        << " process (defined\nin the 'ilchz.cmnd' input file) using the Pythia8.2 "
        << "library. A n-tuple is created \n(called 'mctrue') containing the following info:\n"
        << "\t'pdgId'      : std::vector<int> of the PDG ID code of the stored particle\n"
        << "\t'motherindex': std::vector<int> of the n-tuple vector index of the mother\n"
        << "\t               Note that -1 is used when the particle is the FS hadrons\n"
        << "\t'catchall'   : std::vector<int> a multi-use variable, changing its meaning\n"
        << "\t               depending the type of the particle:\n"
            << "\t\t\t * 0                  for the resonance (H,Z)\n"
            << "\t\t\t * is higher p quark? for the s-squark resonance daughters\n"
            << "\t\t\t * grandmother PDG_ID for the 'final state' strange hadrons\n"
        << "\t'p'          : std::vector<float> momentum of the particle\n"
        << "\t'pmother'    : std::vector<float> momentum of its mother [to be deprecated]\n"
        << "\t'phi'        : std::vector<float> phi of the particle\n"
        << "\t'theta'      : std::vector<float> theta of the particle\n"
        << "Note that the p,phi,theta variables are respect the rest-frame of the q-qbar system\n"
        << "in the case of the final-state hadrons, as well as the pmother\n";
    std::cout << std::endl;
	std::cout << "[OPTIONS]\n\t-o name of the ROOT output file [hzkin.root]\n"
        << "\t-p flag to keep final state PIONS instead of KAONS (default)\n"
        << "\t-h show this help" << std::endl;
}

int main(int argc, char* argv[]) 
{
    // Check that correct number of command-line arguments
    if(argc < 2 && std::string(argv[1]) != "-h") 
    {
        std::cerr << " Unexpected number of command-line arguments. \n You are"
            << " expected to provide one input file name. \n"
            << " Program stopped! " << std::endl;
        return 1;
    }

    // Declare option-related variables
    std::string outputfilename("hzkin.root");
    std::string strangehadrontype("kaons");
	
    std::string cmndfile;
    // get options
    for(int i = 1; i < argc; ++i)
	{
        if( strcmp(argv[i],"-h") == 0 )
		{
            display_usage();
			return 0;
		}
        else if( strcmp(argv[i],"-o") == 0 )
		{
			outputfilename = argv[i+1];
            ++i;
		}
        else if( strcmp(argv[i],"-p") == 0 )
        {
            strangehadrontype = "pions";
        }
        else
        {
            // Check that the provided input name corresponds to an existing file.
            std::ifstream is(argv[i]);
            if(!is && std::string(argv[i]) != "-h") 
            {
                std::cerr << " Command-line file " << argv[i] << " was not found. \n"
                    << " Program stopped! " << std::endl;
                return 2;
            }
            cmndfile = argv[i];
        }
	}

    // Confirm that external files will be used for input and output.
    std::cout << "\n >>> PYTHIA settings will be read from file " << cmndfile << std::endl;
       // << " <<< \n >>> HepMC events will be written to file "
       // << argv[2] << " <<< \n" << std::endl;
    
    // Interface for conversion from Pythia8::Event to HepMC event.
    //HepMC::Pythia8ToHepMC ToHepMC;

    // Specify file where HepMC events will be stored.
    //HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);

    // Generator.
    Pythia pythia;

    // Read in commands from external file.
    pythia.readFile(cmndfile);

    // Extract settings to be used in the main program.
    int    nEvent    = pythia.mode("Main:numberOfEvents");
    int    nAbort    = pythia.mode("Main:timesAllowErrors");

    // Initialization.
    pythia.init();

    // selected hadrons:
    FinalStateHadrons fshadrons(strangehadrontype);
    // ROOT init tree
    ParticleKinRootAux particles(fshadrons.getIDs());

    TTree * thz = new TTree("mctrue","ee -> H Z -> s sbar s sbar");
    particles.inittree(thz);

    // Some pseudo-constants initializations
    // resonances PDG IDs
    std::vector<int> idResonance;
    idResonance.push_back(23);  // Z0
    idResonance.push_back(25);  // Higgs (h0)
     
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

        // Mapping the rest-frame matrix of the quark systems to 
        // its PDG ID resonance number. Declaration
        std::map<int,RotBstMatrix> restframesmap;

        // Note: the extraction algorithm goes backwards: from a
        // final hadron, finding their parents and keeping the info
        // of the original quark and resonance
        for(int i= 0; i < pythia.event.size(); ++i)
        {
            const int abspdgid = pythia.event[i].idAbs();
            
            if( ! pythia.event[i].isFinal() )
            {
                continue;
            }
            // "Interesting" particles only!
            if( std::find(particles.strangehadrons.begin(),particles.strangehadrons.end(),
                        abspdgid) == particles.strangehadrons.end() )
            {
                continue;
            }

            // Obtain the last "carbon" copy of the particle to work with it
            const int currI = pythia.event[i].iBotCopy();
            // If was already used, don't duplicate (note thet currI is the pythia code)
            if( (std::find(particles.getusedid()->begin(),particles.getusedid()->end(),currI) != 
                        particles.getusedid()->end() ) )
            {
                continue;
            }

            Particle & had = pythia.event[currI];

            // Get the resonance and the resonance-daughter quark pythia-index
            int quarkindex = currI;
            const int ancestorindex = getancestorindex(currI,pythia,idResonance,quarkindex);
            // protecting the case of initial state radiation or decays not from the
            // resonance
            if( ancestorindex == 0 && quarkindex == -1)
            {
                continue;
            }
            const int ancestorID = pythia.event[ancestorindex].id();
            // Store the info if isn't, the fillresonance function is also getting
            // the restframe system of the qqbar system 
            if( std::find(particles.getusedid()->begin(),particles.getusedid()->end(),
                        ancestorindex) == particles.getusedid()->end() )
            { 
                // Fill the resonance relative info: plus quarks involved
                restframesmap.insert(fillresonancechain(ancestorindex,pythia,particles));
            }

            // Correct the frame definition (in order to define pz-defined
            // positive for the quark). 
            RotBstMatrix restframe(restframesmap[ancestorID]);
            // Checking the quark is defined in the positive axis
            Particle & quarkatrest = pythia.event[quarkindex];
  std::cout << "[ " << iEvent << " ::: " << quarkindex << " -- " << quarkatrest.name() << "] " << std::endl;
            if( quarkatrest.pz() < 0.0 )
            {
                restframe.rot(M_PI);
            }
            // And convert to qqbar system reference frame
            had.rotbst(restframe);
            
            // storing info in the rest-frame of the quark-bquark ref. system
            particles.filltreevariables(currI,had.id(),particles.getlocalindex(quarkindex),
                    ancestorID,had.pAbs(),quarkatrest.pAbs(),had.phi(),had.theta());
//#ifdef DEBUG
            std::cout << std::endl;
//#endif
            //}
        }
        thz->Fill();
        // deallocate variables after the filling
        particles.endloop();
    }
    // End of event loop. Statistics.
    pythia.stat();

    // ROOT File
    TFile *froot = new TFile(outputfilename.c_str(),"RECREATE");
    thz->Write();
    froot->Close();

    delete thz;
    delete froot;

    // Done.
    return 0;
}
