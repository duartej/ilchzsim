// Simulation of ILC process: ee -> ZH -> q qbar s sbar
// 
// Author: Jordi Duarte-Campderros (TAU), jorge.duarte.campderros@cern.ch.
// (based on the main42.cc which is a part of the PYTHIA 
// event generator. // Copyright (C) 2015 Torbjorn Sjostrand.)
// PYTHIA is licenced under the GNU GPL version 2.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// e+e --> HZ --> qqbar qqbar (q-is defined in the ilchz.cmnd file)
// simulation process. 
//
// Input and output files are specified on the command line, e.g. like
// ./ilchz ilchz.cmnd [OPTIONS] > out.txt
// The main program contains the generation of a n-tuple (ROOT) where
// final state hadrons (kaons or pions, depending the option choosen) 
// relating variables are kept for further analysis. 

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
#include<iterator>


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

    const int B_0      =  511;
    const int B_PLUS   =  521;
    const int B_MINUS  = -521;
    const int B_0_S    =  531;

    const int D_PLUS   =  411;
    const int D_MINUS  = -411;
    const int D_0      =  421;
    const int D_PLUS_S =  431;
    const int D_MINUS_S= -431;

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

    std::vector<int> getcharmed() const
    {
        std::vector<int> _charmed;
        _charmed.push_back(D_0);
        _charmed.push_back(D_PLUS);
        _charmed.push_back(D_MINUS);
        _charmed.push_back(D_PLUS_S);
        _charmed.push_back(D_MINUS_S);

        return _charmed;
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
    
    std::vector<int> getbottoms() const
    {
        std::vector<int> _bottoms;
        _bottoms.push_back(B_0);
        _bottoms.push_back(B_PLUS);
        _bottoms.push_back(B_MINUS);
        _bottoms.push_back(B_0_S);

        return _bottoms;
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
    std::vector<int> * isBCdaughter;
    std::vector<int> * isPrimaryHadron;
    std::vector<float> * p;
    std::vector<float> * pmother;
    std::vector<float> * phi;
    std::vector<float> * theta;
    std::vector<float> * vx;
    std::vector<float> * vy;
    std::vector<float> * vz;

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
        isBCdaughter(0),
        isPrimaryHadron(0),
        p(0),
        pmother(0),
        phi(0),
        theta(0),
        vx(0),
        vy(0),
        vz(0),
       _listofusedpartindex(0) 
    {
        _auxI.push_back(&pdgId);
        _auxI.push_back(&motherindex);
        _auxI.push_back(&catchall);
        _auxI.push_back(&isBCdaughter);
        _auxI.push_back(&isPrimaryHadron);
        _auxF.push_back(&p);
        _auxF.push_back(&pmother);
        _auxF.push_back(&phi);
        _auxF.push_back(&theta);
        _auxF.push_back(&vx);
        _auxF.push_back(&vy);
        _auxF.push_back(&vz);

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
        isBCdaughter   = new std::vector<int>;
        isPrimaryHadron= new std::vector<int>;
        p              = new std::vector<float>;
        pmother        = new std::vector<float>;
        phi            = new std::vector<float>;
        theta          = new std::vector<float>;
        vx             = new std::vector<float>;
        vy             = new std::vector<float>;
        vz             = new std::vector<float>;

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
            const int & _catchall, const int & _isHF, const int & _isPH,
            const float & _p, const float & _pmother, const float & _phi, 
            const float & _theta, const float & _vx, const float & _vy, const float & _vz)
    {
        this->pdgId->push_back(id);
        this->motherindex->push_back(_motherindex);
        this->catchall->push_back(_catchall);
        this->isBCdaughter->push_back(_isHF);
        this->isPrimaryHadron->push_back(_isPH);
        this->p->push_back(_p);
        this->pmother->push_back(_pmother);
        this->phi->push_back(_phi);
        this->theta->push_back(_theta);
        this->vx->push_back(_vx);
        this->vy->push_back(_vy);
        this->vz->push_back(_vz);

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

        return std::distance(_listofusedpartindex->begin(),localindex);
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
        t->Branch("isBCdaughter",&isBCdaughter);
        t->Branch("isPrimaryHadron",&isPrimaryHadron);
        t->Branch("p",&p);
        t->Branch("pmother",&pmother);
        t->Branch("phi",&phi);
        t->Branch("theta",&theta);
        t->Branch("vx",&vx);
        t->Branch("vy",&vy);
        t->Branch("vz",&vz);
    }
};

// Storing variables used related with the resonances (the variable 'i' contains the pythia-index
// of the resonance). The resonance and daughter quarks are persistified. Note that
// the function returns a pair (PDG-ID resonance, RotBstMatrix), where the RotBstMatrix is 
// a Lorentz transformation transforming any 4-vector to the Center of Mass of thee quark-antiquark
// system
std::pair<int,RotBstMatrix> fillresonancechain(const int & i, const Pythia & pythia, ParticleKinRootAux & p)
{
    // Fill the resonance relative info: Lab frame
    // Note that we want to store info before showering,
    // directly from Hard-Scattering
    const int iHS = pythia.event[i].iTopCopy();
    const Particle & resonanceHS = pythia.event[iHS];
    const int respdgId = resonanceHS.id();
    // Filling n-tuple and getting the n-tuple index of the Resonance
    const int reslocalIndex = p.filltreevariables(i,respdgId,-1,0,0,-1,
            resonanceHS.pAbs(),-1,resonanceHS.phi(),resonanceHS.theta(),
            resonanceHS.xProd(),resonanceHS.yProd(),resonanceHS.zProd());

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

    p.filltreevariables(iS,s.id(),reslocalIndex,(int)isLeading_s,0,-1,
            s.pAbs(),resonance.pAbs(),s.phi(),s.theta(),
            s.xProd(),s.yProd(),s.zProd());
    
    // Save info of the sbar-quark
    p.filltreevariables(iSbar,sbar.id(),reslocalIndex,(int)(not isLeading_s),0,-1,
            sbar.pAbs(),resonance.pAbs(),sbar.phi(),sbar.theta(),
            sbar.xProd(),s.yProd(),s.zProd());
    
    // Returning the PDG ID of the resonance and the rest-frame system of 
    // its quarks daughters
    return std::pair<int,RotBstMatrix>(respdgId,restframe);
}

// Get the pythia-index of the ancestor particle (which is defined by its PDG ID in consideredmums)
// and the pythia-index of the daughter of the ancestor (the daughter quark, before radiation). 
// This quark-index is returned by reference (quarkindex variable)
int getancestorindex(const int & currIndex, const Pythia & pythia,const std::vector<int> & consideredmums, int & quarkindex)
{
    // Find the top copy, before calling the mother list
    const int & index = pythia.event[currIndex].iTopCopyId();
    const Particle & hadronbeforerad = pythia.event[index];

    // Check the mother list and found the resonance Id:
#ifdef DEBUG
    std::cout << "--> " << hadronbeforerad.name() << "(" << index << ") :::"; 
#endif
    const std::vector<int> & mums = hadronbeforerad.motherList();
    for(unsigned int k = 0; k < mums.size(); ++k)
    {
#ifdef DEBUG
    std::cout << "--> " << pythia.event[mums[k]].name() << "(" << mums[k] << ")"; 
#endif
        if( std::find(consideredmums.begin(),consideredmums.end(),pythia.event[mums[k]].id()) 
                != consideredmums.end() )
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
        << "\t'pdgId'        : std::vector<int> of the PDG ID code of the stored particle\n"
        << "\t'motherindex'  : std::vector<int> of the n-tuple vector index of the mother\n"
        << "\t                 Note that -1 is used when the particle is the FS hadrons\n"
        << "\t'catchall'     : std::vector<int> a multi-use variable, changing its meaning\n"
        << "\t                 depending the type of the particle:\n"
            << "\t\t\t * 0                  for the resonance (H,Z)\n"
            << "\t\t\t * is higher p quark? for the s-squark resonance daughters\n"
            << "\t\t\t * grandmother PDG_ID for the 'final state' strange hadrons\n"
        << "\t'isBoCdaughter': std::vector<int>   describes if the hadrons is coming from\n"
        << "                   a Bottom or Charm hadron.\n"
        << "\t'p'            : std::vector<float> momentum of the particle\n"
        << "\t'pmother'      : std::vector<float> momentum of its mother [to be deprecated]\n"
        << "\t'phi'          : std::vector<float> phi of the particle\n"
        << "\t'theta'        : std::vector<float> theta of the particle\n"
        << "\t'vx'           : std::vector<float> production vertex, x\n"
        << "\t'vy'           : std::vector<float> production vertex, y\n"
        << "\t'vz'           : std::vector<float> production vertex, z\n"
        << "Note that the p,phi,theta variables are respect the rest-frame of the q-qbar system\n"
        << "in the case of the final-state hadrons, as well as the pmother."
        << "However, the production vertex is respect the Laboratory frame.\n";
    std::cout << std::endl;
	std::cout << "[OPTIONS]\n\t-o name of the ROOT output file [hzkin.root]\n"
        << "\t-b flag to keep track if the final hadrons provenance is from charmed or bottom hadrons\n"
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
    bool accountforheavyflavoured = false;

	
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
        else if( strcmp(argv[i],"-b") == 0 )
        {
            accountforheavyflavoured = true;
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

    // If the user want to keep track of the final hadrons created
    // from heavy flavour hadrons (B or D)
    std::vector<int> * hfhadrons = 0;
    if( accountforheavyflavoured )
    {
        hfhadrons = new std::vector<int>;
        std::vector<int> _provb = fshadrons.getbottoms();
        hfhadrons->insert(hfhadrons->end(),_provb.begin(),_provb.end());
        std::vector<int> _provc= fshadrons.getcharmed();
        hfhadrons->insert(hfhadrons->end(),_provc.begin(),_provc.end());
    }
     
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
        // Allocate variables
        particles.initloop();

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

            Particle had = pythia.event[currI];

            // Get the resonance and the resonance-daughter quark pythia-index
            int pre_quarkindex = currI;
            const int ancestorindex = getancestorindex(currI,pythia,idResonance,pre_quarkindex);
            // protecting the case of initial state radiation or decays not from the
            // resonance
            if( ancestorindex == 0 && pre_quarkindex == -1)
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

            // Get the quark before radiation, so be sure by calling iTopCopyId method
            const int quarkindex = pythia.event[pre_quarkindex].iTopCopyId();
            // Correct the frame definition (in order to define pz-defined
            // positive for the quark).
            RotBstMatrix restframe(restframesmap[ancestorID]);
            Particle quarkatrest = pythia.event[quarkindex];
            quarkatrest.rotbst(restframe);
            // NOTE:: Not use the flipping depending of the quark direction,
            // due to the fact we want to deal with the two highest-pt hadrons
            // in the opposite hemispheres of the CM-qqbar reference system
            //Checking the quark is defined in the positive axis
            //if( quarkatrest.pz() < 0.0 )
            //{
            //    restframe.rot(M_PI);
            //}
            // And convert to qqbar system reference frame
            had.rotbst(restframe);
         
            // Check if the hadron come from a Bottom or charmed hadron (if relevant)
            int isBCdaughter = 0;
            pre_quarkindex   = currI;
            if(hfhadrons != 0)
            {
               if( getancestorindex(currI,pythia,*hfhadrons,pre_quarkindex) != 0 )
               {
                   isBCdaughter = 1;
               }
            }

            // Check if it is a primary hadron: Pythia code 81-89
            int isPrimaryHadron = 0;
            if( had.status() >= 81 && had.status() <= 89)
            {
                isPrimaryHadron = 1;
            }
            // storing info in the rest-frame of the quark-bquark ref. system (except
            // for the vertex)
            const Particle & hadatlab = pythia.event[currI];
            particles.filltreevariables(currI,had.id(),particles.getlocalindex(quarkindex),
                    ancestorID,isBCdaughter,isPrimaryHadron,
                    had.pAbs(),quarkatrest.pAbs(),had.phi(),had.theta(),
                    hadatlab.xProd(),hadatlab.yProd(),hadatlab.zProd());
#ifdef DEBUG
            std::cout << std::endl;
#endif
        }
        thz->Fill();
        // deallocate variables after the filling
        particles.endloop();
    }
    // Deallocate
    if( hfhadrons != 0 )
    {
        delete hfhadrons;
        hfhadrons=0;
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
