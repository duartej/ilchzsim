#!/usr/bin/env python
""":script:`doplots` -- Plotting stuff created with processhzroot script
======================================================================

.. script:: pion_misid [OPTIONS]    
      :platform: Unix
      :synopsis: Evaluate pion misid...
	  .. moduleauthor:: Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
"""
#_COLOR = [ ROOT.kCyan+2, ROOT.kOrange+5,ROOT.kAzure-7,ROOT.kGreen+2,ROOT.kRed-2, ROOT.kBlue-3,
#        ROOT.kBlack, ROOT.kRed+4]
def getcolor():
    import ROOT
    return [ ROOT.kRed+4, ROOT.kAzure+3, ROOT.kOrange-2, ROOT.kGreen-5, ROOT.kYellow+2, \
        ROOT.kCyan-2, ROOT.kOrange+5,ROOT.kAzure-7,ROOT.kGreen-2,ROOT.kRed-4, ROOT.kGray-3 ]

SUFFIXPLOTS='.pdf'

# Static class to deal with inputs
class higgsinputs:
    # Data obtained from 
    # Branching ratios: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3
    # Cross section   : 
    def __init__(self):
        self.mH         = 125.7 # (GeV)
        self.eeToHat250 = 2.5e2 # (fb) (at 250 GeV center of mass, is almost equivalent to ee->HZ)
        # Generic, only dependent of the higgs mass
        # therefore, following values at self.mH higgs
        # mass
        self.brbbbar    = 5.66e-1
        self.brccbar    = 2.85e-2
        self.brssbar    = 2.41e-4
        self.bruubar    = 0.0  # just below 1e-4
        self.brddbar    = 0.0  # just below 1e-4
        self.Lint       = 0 #(fb)
        self.higgsproduced = 0

    def setLint(self,_Lint):
        self.Lint = float(_Lint)
        self.higgsproduced = self.eeToHat250*self.Lint

    def getEvents(self,decaychannel):
        """
        """
        if self.Lint == 0:
            raise RuntimeError("First set the integral luminosity"\
                    " using the 'setLint(Lint) method")
        try:
            br = getattr(self,'br'+decaychannel)
        except AttributeError:
            raise AttributeError('There is no br "%s" '\
                    'defined in this class' % decaychannel)
        return self.higgsproduced*br

hiInstance = higgsinputs()
hiInstance.setLint(500) # (fb-1)

# Branching ratios with respect ccbar (m_q^2/m_c^2)
# using the pole mass (measured) not the mass at the higgs scale
BRuu_cc = 3.3e-6
BRdd_cc = 1.4e-5
BRss_cc = 0.006
BRbb_cc = 10.75

BRuu_ss = 5.9e-4
BRdd_ss = 2.6e-3
BRcc_ss = 180.1
BRbb_ss = 1936.0

def getbr(name):
    """.. function:: getbr(name) -> br

    given the name of a sample returns the branching
    reation involved (with respect the ccbar)
    """
    if name.find("bbbar") != -1:
        return BRbb_ss
    elif name.find("uubar") != -1:
        return BRuu_ss
    elif name.find("ddbar") != -1:
        return BRdd_ss
    elif name.find("ssbar") != -1:
        return 1.0
    elif name.find("ccbar") != -1:
        return BRcc_ss

def parseprocess(name):
    """
    """
    availpr = [ 'bbbar', 'uubar', 'ddbar', 'ssbar', 'ccbar' ]
    try:
        process = filter(lambda x: name.find(x) != -1,availpr)[0]
    except IndexError:
        raise RuntimeError("Not correctly parsed: '%s'" % name)
      
    return process


def get_filtered_eff(_obj,decay_channel,hadrons,pLcut,d0cut):
    """Extract the efficiency given exist an histogram related with the
    decay channel and hadrons passing the pl and d0 cut

    Parameters
    ----------
    _obj : dict(str,TH3F)
        the histograms with the proper names containing the decay channel and the
        hadrons involved
    decay_channel: str
    hadrons: str, 
        the hadrons to look at, which can be pions, kaons or 
        pions_kaons (kaons_pions)
    pLcut: float 
    d0cut: float

    Returns
    -------
    efficiency: float
    """
    try:
        histo_name = filter(lambda x: 
                x.find('H_th3_hz{0}'.format(decay_chanel)) ==0 and 
                    x.find('{0}_PLd0s'.format(hadrons)) != -1,
                _obj.keys())[0]
    except KeyError:
        raise RuntimeError('Not found the histogram "H_th3_hz{0]_*_{1}_PLd0s"'\
                ' in the root file'.format(decay_channel,hadrons))

        h = _obj['H_th3_hz{0}_{1}_PLd0s'.format(decay_channel, hadrons)]

    # TH1 histos to obtain bin numbers and other manipulations
    pL_h  = h.ProjectionX()
    d01_h = h.ProjectionY()
    d02_h = h.ProjectionZ()

    # First obtaining the bin where the cut are located
    pL_bin = pL_h.FindBin(pLcut)
    d01_bin= d01_h.FindBin(d0cut)
    d02_bin= d02_h.FindBin(d0cut)
    # counting all the events, including those in the overflow bins
    pL_maxBin = h.GetNbinsX()+1

    evts = h.Integral(pL_bin,pL_maxBin,1,d01_bin,1,d02_bin)

    # Entries takes into account overflow bin
    return float(evts)/float(h.GetEntries())

def plots():
    """
    """
    suffixplots = globals()['SUFFIXPLOTS']

    try:
        from PyAnUtils.plotstyles import squaredStyle,setpalette
        lstyle = squaredStyle()
        lstyle.cd()
        ROOT.gROOT.ForceStyle()
        ROOT.gStyle.SetOptStat(0)

        setpalette("darkbody")
    except ImportError:
        pass



def main(rootfile):
    """Main function gathering all the plots to be performed:
     * Plot of all the efficiencies vs. momentum cut in the
     same canvas (including the total background efficiency
     which is calculated in the function gettotalbkgeff)
     * Significance plots (eff_S/sqrt(e_totalbkg))
     * ROC curves (signal efficiency vs. total bkg efficiency)
    The significance plots show as well some working points
    where the signal eff. and total bck. efficiency are shown
    
    """
    import ROOT
    from math import sqrt

    
    # Get the root file with and the TH3 object (only)
    f = ROOT.TFile(rootfile)
    _preobj = dict(map(lambda x:(x.GetName(), f.Get(x.GetName())),
         filter(lambda x: x.GetClassName().find('TH3') == 0,f.GetListOfKeys())))
    
    # Split by resonance
    _obj = { 'H' : dict(filter(lambda (x,y): x.find('H') == 0,_preobj.iteritems())),
            'Z' : dict(filter(lambda (x,y): x.find('Z') == 0,_preobj.iteritems()))
            }

    
    # Plotting
    ROOT.gROOT.SetBatch()
    
    # just take care in the H-ressonance...

    # XXX: can be entered by option
    d0cuts = [0.3, 0.5]
    pLcuts = [20.0]
    eff = {}
    # Prepare efficiencies (changing pL cut) per each d0-cut
    for d0 in d0cuts:
        for pL in pLcuts:
            eff_ssbar_pions = get_filtered_eff(_obj['H'],'ssbar','pions',pL,d0)
            eff_ssbar_kaons = get_filtered_eff(_obj['H'],'ssbar','kaons',pL,d0)
 
            eff_bbbar_pions = get_filtered_eff(_obj['H'],'bbbar','pions',pL,d0)
            eff_bbbar_kaons  = get_filtered_eff(_obj['H'],'bbbar','kaons',pL,d0)

            print '='*80
            print 'd0 = {0:.1f} [mm] , pL = {1:.1f} [GeV] ======='.format(d0,pL)
            print 'Pion efficiency in ssbar:',eff_ssbar_pions
            print 'Pion efficiency in bbbar:',eff_bbbar_pions
            #print 
            #print 'N_ssbar={0:.2f}    N_bbbar={1:.2f}'.format(n_ssbar_pions,n_bbbar_pions)
            #print 'N_bbbar/N_ssbar={0:f}'.format(float(n_bbbar_pions)/float(n_ssbar_pions))
            #print            
            #print 'N_ssbar={0:.2f}    N_bbbar={1:.2f}'.format(n_ssbar_kaons,n_bbbar_pions)
            #print 'N_bbbar/N_ssbar={0:f}'.format(float(n_bbbar_pions)/float(n_ssbar_kaons))
    #    eff = getefficiencies(_obj['H'],d0)



if __name__ == '__main__':
    from optparse import OptionParser,OptionGroup
    import os

    #Opciones de entrada
    parser = OptionParser()
    parser.set_defaults(inputfile='processed.root')    
    parser.add_option( '-i', '--inputfile', action='store', type='string', dest='inputfile',\
            help="input root filename [processed.root]")
    parser.add_option( '-s', '--suffix', action='store', type='string', dest='suffixout',\
            help="output suffix for the plots [.pdf]")
    
    (opt,args) = parser.parse_args()

    if opt.suffixout:
        suff = opt.suffixout.replace('.','')
        globals()['SUFFIXPLOTS'] ='.'+suff
    main(os.path.abspath(opt.inputfile))
