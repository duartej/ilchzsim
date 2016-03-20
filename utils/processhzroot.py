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

KAON_ID = 321
PION_ID = 211

#FIXME: NEED TO ADD A WAY TO USE ONLY KAONS, PIONS OR BOTH (a new option)


def gethadronsind(iEvent):
    """.. function:: gethadronsind(iEvent) -> { resPDGID: [ i1, i2, ... ], .. }

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
    
    :return: local indices list corresponding to hadrons (per resonance)
    :rtype : dict(list(int))
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

def settitles(h,totalmomentum):
    """.. function:: settitles(h,ztotalmomentum) 
    
    put the titles in the histograms. The totalmomentum
    flag is indicating to use the total momentum instead
    the paralel to the quark axis
    """
    h.SetTitle('')
    if totalmomentum:
        h.SetXTitle('p^{1} [GeV]')
        h.SetYTitle('p^{2} [GeV]')
    else:
        h.SetXTitle('p_{||}^{1} [GeV]')
        h.SetYTitle('p_{||}^{2} [GeV]')


def process(inputfile,outputfile,d0cut,trackHF):
    """..function:: process(inputfile) 

    Process the 'inputfile' ROOT (obtained from the ilchz executable)
    in order to obtain several ROOT objects:
      * TH2F: Paralel momentum of the leading hadrons (one per resonance)
      * TEfficiency: efficiency (radial) cut in p1 vs. p2
      * TH2F: Number of events vs. d0. pL circular cut
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
    from math import pi,cos,sqrt

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
    # More histograms if the user wants to keep track of the heavy flavour decays
    histodictslist = [hleading] 
    if trackHF:
        hleadingHF = { 25: ROOT.TH2F('H'+namehistos+'_HF','paralel momentum leading hadrons',NBINS,XMIN,XMAX,NBINS,XMIN,XMAX),
            23: ROOT.TH2F('Z'+namehistos+'_HF','paralel momentum leading hadrons',NBINS,XMIN,XMAX,NBINS,XMIN,XMAX) }
        histodictslist.append(hleadingHF)
    # Some cosmethics
    for histodict in histodictslist:
        for res,h in histodict.iteritems():
            dototalmomentum = False
            if res == 23:
                dototalmomentum=True
            settitles(h,dototalmomentum)
    
    # sqrt(d0_1^2+d0^2_2), d0_1,d0^2
    # for 3 different cases: kaons_kaons, kaons_pions and pions_pions
    # title
    TITLE = 'p_{||} circular cut and d_{0}; p_{||} cut [GeV];'\
            'd_{0}^{1} [mm];  d_{0}^{2} [mm]'
    D0MAX = 1.0
    # quick way to initialize the histograms: hd0PL_kk, hd0PL_kp and hd0PL_pp
    _prehists = {}
    for dictname,suffix in  [('hd0PL_kk','_KK_PLd0s'),('hd0PL_kp','_KP_PLd0s'),('hd0PL_pp','_PP_PLd0s')]:      
        _prehists[dictname] = \
                {25:  ROOT.TH3F('H'+namehistos.replace('th2','th3')+suffix,
                    TITLE,NBINS,XMIN,XMAX,
                    NBINS*2,-1.*D0MAX,D0MAX,NBINS*2,-1.*D0MAX,D0MAX),
                 23: ROOT.TH3F('Z'+namehistos.replace('th2','th3')+suffix,
                    TITLE,NBINS,XMIN,XMAX,
                    NBINS*2,-1.*D0MAX,D0MAX,NBINS*2,-1.*D0MAX,D0MAX)
                 }
    hd0PL_kk = _prehists['hd0PL_kk']
    hd0PL_kp = _prehists['hd0PL_kp']
    hd0PL_pp = _prehists['hd0PL_pp']

    histodictslist.append(hd0PL_kk)
    histodictslist.append(hd0PL_kp)
    histodictslist.append(hd0PL_pp)

    # d0-histograms : -->  TO BE DEPRECATED, ioncluded in the TH3
    extrahist = { 25: ROOT.TH2F('H'+namehistos+'_d0','',NBINS*2,-5,5,NBINS*2,-5.,5.),
            23: ROOT.TH2F('Z'+namehistos+'_d0','',NBINS*2,-5,5,NBINS*2,-5.,5.) 
            }
    th1list = extrahist.values()
    if trackHF:
        extrahistHF = { 25: ROOT.TH1F('H'+namehistos+'_d0_HF','',NBINS*2,-5,5,NBINS*2,-5.,5),
                23: ROOT.TH1F('Z'+namehistos+'_d0_HF','',NBINS*2,-5,5,NBINS*2,-5.,5) 
                }
        th1list += extrahistHF.values()
    # cosmethics
    for h in th1list:
        h.SetXTitle('d_{0}^{1} [mm]')
        h.SetYTitle('d_{0}^{2} [mm]')


    # Event loop
    noOpposite = { 23: 0, 25: 0 }; noOppositeHF = {23: 0, 25: 0} ;
    nOpposite = {23: 0, 25: 0 } ; nOppositeHF = {23: 0, 25: 0};
    msg = "Evaluating {0}...".format(inputfile)
    if len(msg) > 100:
        shorten_name = "{0}.../{1}".format(inputfile[:50],
                os.path.basename(inputfile))
        msg = "Evaluating {0}...".format(shorten_name)
    pointpb = float(t.GetEntries())/100.0
    for _i,iEvent in enumerate(t):
        # Progress bar
        sys.stdout.write("\r\033[1;34m+-- \033[1;m"+msg+\
                "[ " +"\b"+str(int(float(_i)/pointpb)+1).rjust(3)+"%]")
        sys.stdout.flush()
        # Get the local indices of the leading kaons
        # -- store histograms of p (paralel component with respect the quark mother)
        #    for leading kaons in both hemispheres (in the CM-quark-antiquark system
        #    of reference)
        hadronsind    = gethadronsind(iEvent)
        for res,hadlist in hadronsind.iteritems():
            # One line: order hadrons by paralel momentum (or should I do it by momentum?)
            #           and return the local index (just for the hadrons coming from the
            #           resonance 'res'
            if res == 23:
                fcos = lambda k: 1.0
            else:
                fcos = lambda k: abs(cos(iEvent.theta[k]))
                
            keycmp = lambda (x,y): y
            sortedhadrind = filter(lambda x: x in hadronsind[res], 
                    map(lambda (x,y): x, sorted(enumerate(iEvent.p),reverse=True,\
                            key=lambda (x,y): y*fcos(x) ) 
                       )
                    )
            # Get the leading kaons in the opposite hemispheres (in the CM q-qbar 
            # reference system):
            # Fill a list and take first elements of the list
            up_List   = []
            down_List = []
            upHF_List   = []
            downHF_List = []
            upPDGID_List   = []
            downPDGID_List = []

            for k in sortedhadrind:
                # Check the impact parameter if the cut is activated
                if d0cut:
                    ## See below [1]
                    d0 = (iEvent.vy[k]-iEvent.vx[k]*iEvent.phi_lab[k])*cos(iEvent.phi_lab[k])
                    #d0 = sqrt(iEvent.vx[k]**2.+iEvent.vy[k]**2.)
                    if d0 > float(d0cut):
                        continue
                momentum = iEvent.p[k]*fcos(k)
                # [1] --> this is not d0, is radius on the transverse plane !!
                # d0       = sqrt(iEvent.vx[k]**2.0+iEvent.vy[k]**2.0)
                # A rough estimation assuming straight line trajectory: use 
                # the phi angle at the lab frame (define the vector director 
                # of a straight line) and the minimum distance of the line 
                # with respect the IP 
                d0 = (iEvent.vy[k]-iEvent.vx[k]*iEvent.phi_lab[k])*cos(iEvent.phi_lab[k])
                if iEvent.theta[k] < pi/2.0:
                    if trackHF and iEvent.isBCdaughter[k]:
                        upHF_List.append( (momentum,d0) )
                    else:
                        up_List.append((momentum,d0))
                        upPDGID_List.append(iEvent.pdgId[k])
                else:
                    if trackHF and iEvent.isBCdaughter[k]:
                        downHF_List.append((momentum,d0))
                    else:
                        down_List.append((momentum,d0))
                        downPDGID_List.append(iEvent.pdgId[k]) 
                # Check if already filled, so break the loop
                if len(up_List) > 0 and len(down_List) > 0:
                    if trackHF:
                        if len(up_List) > 0 and len(down_List) > 0:
                            break
                    else: 
                        break
            # Just getting the highest momentum particles: index 0 of the lists
            try:
                # --- the 
                pSqrt  = sqrt(up_List[0][0]**2.0+down_List[0][0]**2.0)

                hleading[res].Fill(up_List[0][0],down_List[0][0])
                extrahist[res].Fill(up_List[0][1],down_List[0][1])
                nOpposite[res] += 1

                # keeping info regarding the hadron content
                if abs(int(upPDGID_List[0]*downPDGID_List[0])) == KAON_ID*KAON_ID:
                    hd0PL_kk[res].Fill(pSqrt,up_List[0][1],down_List[0][1])
                elif abs(int(upPDGID_List[0]*downPDGID_List[0])) == KAON_ID*PION_ID :
                    hd0PL_kp[res].Fill(pSqrt,up_List[0][1],down_List[0][1])
                elif abs(int(upPDGID_List[0]*downPDGID_List[0])) == PION_ID*PION_ID :
                    hd0PL_pp[res].Fill(pSqrt,up_List[0][1],down_List[0][1])
            except IndexError:
                noOpposite[res] =+ 1
            if trackHF:
                try:
                    hleadingHF[res].Fill(upHF_List[0][0],downHF_List[0][0])
                    extrahistHF[res].Fill(upHF_List[0][1],downHF_List[0][1])
                    nOppositeHF[res] += 1
                except IndexError:
                    noOppositeHF[res] += 1
    for res,opp in nOpposite.iteritems():
        print "\n[%i] Not found opposite hadrons in %i (of a total of %i) events" % \
                (res,noOpposite[res],(noOpposite[res]+opp))
        if trackHF:
            print "[%i] Not found opposite hadrons (decaying from heavy flavour mesons) "\
                    "in %i (of a total of %i) events" % \
                    (res,noOppositeHF[res],(noOppositeHF[res]+nOppositeHF[res]))
    
    # Persistency
    # Evaluate efficiency of a simple (radial) cut
    efile = ROOT.TFile(outputfile,'UPDATE')
    eff = {}
    resn = { 23: 'Z', 25: 'H'}
    for res in resn.keys():
        name = resn[res]+'_eff_'+inputfile.split('/')[-1].replace('.root','')
        eff[res] = evalefficiency(hleading[res],xrange(0,60,1),name)
        if trackHF:
            eff[str(res)+"_HF"] = evalefficiency(hleadingHF[res],xrange(0,60,1),name+'_HF')
            eff[str(res)+"_HF"].Write("",ROOT.TObject.kOverwrite)
        # and storing it
        eff[res].Write("",ROOT.TObject.kOverwrite)
        # also store the TH2F
        for histodict in histodictslist:
            histodict[res].Write("",ROOT.TObject.kOverwrite)
        # and the TH1F
        for _hth1 in th1list:
            _hth1.Write("",ROOT.TObject.kOverwrite)
    efile.Close()
    del efile
    
if __name__ == '__main__':
    from optparse import OptionParser,OptionGroup
    import os
    
    #Opciones de entrada
    parser = OptionParser()
    parser.set_defaults(inputfile='hzkin.root',outputfile='processed.root',d0cut=False,hf=False)  
    parser.add_option( '-i', '--inputfile', action='store', type='string', dest='inputfile',\
            help="input root filename [hzkin.root]")
    parser.add_option( '-o', '--outputfile', action='store', type='string', dest='outputfile',\
            help="output root filename [processed.root]")
    parser.add_option( '-t', '--trackHF', action='store_true', dest='trackhf',\
            help="do the plots split the final hadrons depending whether their provenance"\
            " are heavy flavoured mesons")
    parser.add_option( '-d', '--d0cut', metavar="D0CUT",action='store', dest='d0cut',\
            help="activate the d0 cut (d0<D0CUT) for the hadrons to be considered")
    
    (opt,args) = parser.parse_args()

    process(os.path.abspath(opt.inputfile),os.path.abspath(opt.outputfile),opt.d0cut,opt.trackhf)

