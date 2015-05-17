#!/usr/bin/env python
""":script:`processhzroot` -- Process the ROOT files obtained from ilchz executable                                                            
===================================================================================

.. script:: processhzroot     
      :platform: Unix
      :synopsis: Process the root files created by the ilchz executable in order
                 to obtain histogram and efficiency objects

     .. moduleauthor:: Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
"""
NBINS = 50
XMIN  = 0
XMAX  = 60


def gethadronsind(iEvent):
    """.. function:: gethadronsind(iEvent) -> [ i1, i2, ... ]

    return the local indices of all the particles which are
    hadrons (in practical, because of 
         *  motherindex: -1           for H,Z
                         local index  remaining
    
         *  catchall   : 0            for H,Z  
                         0,1          q,qbar
                         PDGID        hadron
    )

    :param iEvent: current tree-entry
    :type  iEvent: TTree
    
    :return: local indices list corresponding to hadrons
    :rtype : list(int)
    """
    hadronsind = {}
    for i,pdgid in enumerate(iEvent.catchall):
        if pdgid > 1:
            try:
                hadronsind[pdgid].append(i)
            except KeyError:
                hadronsind[pdgid] = [i]
    return hadronsind

def evalefficiency(h,cutlist,name):
    """.. function:: evalefficiency(h,cutlist,name) -> ROOT.TEfficiency

    Obtain the efficiency object given a simple momentum (radial) cut
    in the p_{||}^1 vs. p_{||}^2 space

    :param h: histogram with the number of entries for the paralel momentum
              (with respect the quark-antiquark system) of the leading hadrons
    :type  h: ROOT.TH2F
    :param cutlist: list of the momentum cut [GeV]
    :type  cutlist: list(float)
    :param name: Name of the inputfile root used to extract h
    :type  name: str
    """
    import ROOT
    from math import sqrt

    # Count how many events where selected above the
    # radial cut (which are inside the cutlist), 
    # i.e. how many events fulfil:  sqrt(p**2+p**2) > cut
    integral = dict([(x,0) for x in cutlist])
    for ibin in xrange(1,h.GetNbinsX()+1):
        p1 = h.GetXaxis().GetBinCenter(ibin)
        for jbin in xrange(1,h.GetNbinsY()+1):
            p2 = h.GetYaxis().GetBinCenter(jbin)
            content=h.GetBinContent(ibin,jbin)
            r = sqrt(p1**2.0+p2**2.0)
            for cut in cutlist:
                if r > cut:
                    integral[cut]+=content
    total = h.Integral()
    
    # Build the efficiency object (efficiency vs. cut) 
    # and fill it the Npass/Ntotal
    n = len(integral)
    eff = ROOT.TEfficiency(name,"Efficiency simple cut (r > cut);\
                cut (< #sqrt{p_{||1}^{2}+p_{||}^{2}}) [GeV];#varepsilon",50,0,60)
    for (i,(p,npassed)) in enumerate(integral.iteritems()):
        #e = float(npassed)/float(total)
        ibin=eff.FindFixBin(p)
        eff.SetTotalEvents(ibin,int(total))
        eff.SetPassedEvents(ibin,int(npassed))
    return eff


def process(inputfile,outputfile):
    """..function:: process(inputfile) 

    Process the 'inputfile' ROOT (obtained from the ilchz executable)
    in order to obtain several ROOT objects:
      * TH2F: Paralel momentum of the leading hadrons (one per resonance)
      * TEfficiency: efficiency (radial) cut in p1 vs. p2
    These objects are written to the ROOT file: 'processed.root'. If
    any object already exists in the ROOT file will be overwritten.

    Note that the name of these objects are highly dependent of the name
    of the input file. It is assumed that the input file name contains the 
    process (quark-antiquark) and the type of the analysed hadrons (kaons,
    pions,...)

    :param inputfile: name of the input file (should be processed by ilchz)
    :type inputfile : str
                    
    """
    import sys
    import ROOT
    from math import pi,cos

    # ROOT File
    f = ROOT.TFile(inputfile)
    if f.IsZombie():
        mes= "ROOT file %s not found" % inputfile
        raise IOError(mes)
    t = f.Get('mctrue')

    # Histogram definitions
    namehistos = '_th2_'+inputfile.split('/')[-1].replace('.root','')
    hleading ={ 25: ROOT.TH2F('H'+namehistos,'paralel momentum leading hadrons',NBINS,XMIN,XMAX,NBINS,XMIN,XMAX),
            23: ROOT.TH2F('Z'+namehistos,'paralel momentum leading hadrons',NBINS,XMIN,XMAX,NBINS,XMIN,XMAX) }
    # Some cosmethics
    for h in hleading.values():
        h.SetTitle('')
        h.SetXTitle('p_{||}^{1} [GeV]')
        h.SetYTitle('p_{||}^{2} [GeV]')

    # Event loop
    noOpposite = 0
    nOpposite = 0
    msg = "Evaluating %s..." % inputfile
    pointpb = float(t.GetEntries())/100.0
    _i = 0
    for iEvent in t:
        # Progress bar
        sys.stdout.write("\r\033[1;34m+-- \033[1;m"+msg+\
                "[ " +"\b"+str(int(float(_i)/pointpb)+1).rjust(3)+"%]")
        sys.stdout.flush()
        _i+=1
        # Get the local indices of the leading kaons
        # -- store histograms of p (paralel component with respect the quark mother)
        #    for leading kaons in both hemispheres (in the CM-quark-antiquark system
        #    of reference)
        hadronsind    = gethadronsind(iEvent)
        for res,hadlist in hadronsind.iteritems():
            # One line: order hadrons by paralel momentum (or should I do it by momentum?)
            #           and return the local index (just for the hadrons coming from the
            #           resonance 'res'
            sortedhadrind = filter(lambda x: x in hadronsind[res], 
                    map(lambda (x,y): x, sorted(enumerate(iEvent.p),reverse=True,\
                            key=(lambda (k,y):y*abs(cos(iEvent.theta[k]))) ) 
                       )
                    )
            pup   = []
            pdown = []
            # Get the leading kaons in the opposite hemispheres (in the CM q-qbar 
            # reference system
            for k in sortedhadrind:
                if iEvent.theta[k] < pi/2.0:
                    pup.append(iEvent.p[k]*abs(cos(iEvent.theta[k])))
                else:
                    pdown.append(iEvent.p[k]*abs(cos(iEvent.theta[k])))
            try:
                hleading[res].Fill(pup[0],pdown[0])
                nOpposite += 1
            except IndexError:
                noOpposite =+ 1
    print "\nNot found opposite hadrons in %i (of a total of %i) events" % \
            (noOpposite,(noOpposite+nOpposite))
    
    # Persistency
    # Evaluate efficiency of a simple (radial) cut
    efile = ROOT.TFile(outputfile,'UPDATE')
    eff = {}
    resn = { 23: 'Z', 25: 'H'}
    for res in resn.keys():
        name = resn[res]+'_eff_'+inputfile.split('/')[-1].replace('.root','')
        eff[res] = evalefficiency(hleading[res],xrange(0,60,1),name)
        # and storing it
        eff[res].Write("",ROOT.TObject.kOverwrite)
        # also store the TH2F
        hleading[res].Write("",ROOT.TObject.kOverwrite)
    efile.Close()
    del efile
    
if __name__ == '__main__':
    from optparse import OptionParser,OptionGroup
    import os
    
    #Opciones de entrada
    parser = OptionParser()
    parser.set_defaults(inputfile='hzkin.root',outputfile='processed.root')   
    parser.add_option( '-i', '--inputfile', action='store', type='string', dest='inputfile',\
            help="input root filename [hzkin.root]")
    parser.add_option( '-o', '--outputfile', action='store', type='string', dest='outputfile',\
            help="output root filename [processed.root]")
    
    (opt,args) = parser.parse_args()

    process(os.path.abspath(opt.inputfile),os.path.abspath(opt.outputfile))

