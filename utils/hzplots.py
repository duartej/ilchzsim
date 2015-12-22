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
SPINNING = [ "-","\\","|","/"]

def getleg(**kwd):
    """.. function:: getleg(**kwd) -> ROOT.TLegend()
    
    Return a ROOT.TLegend object with some cosmethics already filled

    kwd accepts the following keys: x0, x1, y0, y1
    """
    import ROOT
    class coord:
        def __init__(self):
            self.x0=0.2
            self.x1=0.45
            self.y0=0.7
            self.y1=0.9
    c = coord()

    for key,val in kwd.iteritems():
        setattr(c,key,val)

    leg = ROOT.TLegend(c.x0,c.y0,c.x1,c.y1)
    leg.SetBorderSize(0)
    leg.SetFillColor(10)
    leg.SetTextSize(0.04)
    leg.SetTextFont(112)

    return leg

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
    reation involved 
    """
    return getattr(hiInstance,'br{0}'.format(name))

def getbr_cc(name):
    """.. function:: getbr(name) -> relative_br

    given the name of a sample returns the branching
    reation involved (with respect the ccbar)
    """
    if name.find("bbbar") != -1:
        return BRbb_cc
    elif name.find("uubar") != -1:
        return BRuu_cc
    elif name.find("ddbar") != -1:
        return BRdd_cc
    elif name.find("ssbar") != -1:
        return BRss_cc
    elif name.find("ccbar") != -1:
        return 1.0

def getbr_ss(name):
    """.. function:: getbr(name) -> relative_br

    given the name of a sample returns the branching
    reation involved (with respect the ssbar)
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

def get_histo_name(th3names,decay_channel,hadron_state):
    """Given a standarized (TH3) histogram list of names (obtained from the 
    processedhzroot.py script), it returns the name of the histogram matching
    the final state

    Parameters
    ----------
    th3names: list(str)
        the names of all the TH3 dicts found in a root processed file by the
        processedhzroot script
    hadron_state: str
        the two opposite-hemisphere final state hadrons: KK, KP, PP (K-kaon,
        P-pion)

    Return
    ------
    str, the name of the input list matching the final state
    """
    try:
        histo_name = filter(lambda x: 
                x.find('H_th3_hz{0}'.format(decay_channel)) ==0 and 
                    x.find('{0}_PLd0s'.format(hadron_state)) != -1,
                th3names)[0]
    except IndexError:
        raise RuntimeError('Not found the histogram "H_th3_hz{0}_*_{1}_PLd0s"'\
                ' in the root file'.format(decay_channel,hadron_state))
    return histo_name

def get_final_state_pr(_obj,decay_channel,hadrons):
    """Extract the probabiliy of having two hadrons of the asked type in each 
    hemisphere given a decay channel

    .. math::f_{AB}^{q\bar{q}} = P( L_{AB} | (ee\rightarrow HZ\rightarrow q\bar{q}) I_0 )

    Parameters
    ----------
    _obj : dict(str,TH3F)
        the histograms with the proper names containing the decay channel and the
        hadrons involved
    decay_channel: str
    hadrons: str, 
        the hadrons to look at, which can be the: KK, KP, PP (K-kaon, P-pion)

    Returns
    -------
    efficiency: float
    """
    
    # Get the final state hadron histo
    histo_name = get_histo_name(_obj.keys(),decay_channel,hadrons)
    n_current_hadrons = _obj[histo_name].GetEntries()

    # get the list of TH3 histos related to this decay_channel
    histos = filter(lambda h: 
            h.GetName().find('H_th3_hz{0}'.format(decay_channel)) == 0, \
                    _obj.values())
    
    n_total_decay = sum(map(lambda x: x.GetEntries(), histos))

    return float(n_current_hadrons)/float(n_total_decay)


def get_cuts_eff(_obj,decay_channel,hadrons,pLcut,d0cut):
    """Extract the efficiency of the cuts (d0 and pL-circular) given a final 
    couple of hedrons in each hemisphere, and a decay channel

    .. math::\varepsilon_{AB}^{q\bar{q}} = P( d0^{c} p_{||}^c |
                   L_{AB} (ee\rightarrow HZ\rightarrow q\bar{q}) I_0 )

    Parameters
    ----------
    _obj : dict(str,TH3F)
        the histograms with the proper names containing the decay channel and the
        hadrons involved
    decay_channel: str
    hadrons: str, 
        the hadrons to look at, which can be the: KK, KP, PP (K-kaon, P-pion)
    pLcut: float 
    d0cut: float

    Returns
    -------
    efficiency: float
    """
    histo_name = get_histo_name(_obj.keys(),decay_channel,hadrons)
    h = _obj[histo_name]

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

    if h.GetEntries() == 0: 
        return 0.0

    # Entries takes into account overflow bin
    return float(evts)/float(h.GetEntries())
    

class eff_cut_hadron(object):
    """
    .. math::\varepsilon_{AB}^{q\bar{q}} = P( d0^{c} p_{||}^c L_{AB} | 
                (ee\rightarrow HZ\rightarrow q\bar{q}) I_0 )
    """
    def __init__(self,decay_channel):
        """
        """
        self.decay_channel = decay_channel
        self.final_state_hadrons = ['KK','KP','PP']
        self.initialized = False
        self.current_d0cut = None
        self.current_pLcut = None

    def __str__(self):
        out = 'eff_cut_hadron instance: \n'
        out+= '  decay channel: {0} \n'.format(self.decay_channel)
        
        eff = '  d0,pL cut eff:            '.format(self.decay_channel)
        p   = '  H1_H2 final states prob:  '.format(self.decay_channel)
        if self.initialized:
            for i in self.final_state_hadrons:
                eff+= '{0}={1:.4f} '.format(i,getattr(self,'eff_cut_{0}'.format(i)))
                p  += '{0}={1:.4f} '.format(i,getattr(self,'p_{0}'.format(i)))
            out += eff+'\n'
            out += p+'\n'
        return out


    def set_total_eff(self,histodict,pLcut,d0cut):
        """
        """
        self.current_d0cut = d0cut
        self.current_pLcut = pLcut
        for i in self.final_state_hadrons:
            # p ( d0 pL | L_AB decay_channel I0 )
            self.__setattr__('eff_cut_{0}'.format(i),
                    get_cuts_eff(histodict,self.decay_channel,i,pLcut,d0cut))
            # p ( L_AB | decay_channel I0 )
            self.__setattr__('p_{0}'.format(i),
                    get_final_state_pr(histodict,self.decay_channel,i))
        self.initialized=True

    def get_total_eff(self,hadrons):
        """
        """
        return getattr(self,'eff_cut_{0}'.format(hadrons))*\
                getattr(self,'p_{0}'.format(hadrons))


    def get_events(self,hadrons):
        """
        """
        if not self.initialized:
            raise AttributeError('eff_total ERROR: You need to call the '\
                    'set_total_eff method before try to calculate the total fraction')
        
        return hiInstance.getEvents(self.decay_channel.split('_')[0])*\
                self.get_total_eff(hadrons)

    def get_total_events(self):
        """
        """
        if not self.initialized:
            raise AttributeError('eff_total ERROR: You need to call the '\
                    'set_total_eff method before try to calculate the total fraction')
        
        N = 0
        for hadrons in self.final_state_hadrons:
            N += hiInstance.getEvents(self.decay_channel.split('_')[0])*\
                getattr(self,'eff_cut_{0}'.format(hadrons))*\
                getattr(self,'p_{0}'.format(hadrons))
        return N


class eff_final_state_hadrons(object):
    """
    .. math::\varepsilon_{AB}^{q\bar{q}} = P( d0^{c} p_{||}^c |
                   L_{AB} (ee\rightarrow HZ\rightarrow q\bar{q}) I_0 )
    """
    def __init__(self,decay_channel):
        """
        """
        self.decay_channel = decay_channel
        self.final_state_hadrons = ['KK','KP','PP']

    def set_total_eff(self,histodict,pLcut,d0cut):
        """
        """
        for i in self.final_state_hadrons:
            self.__setattr__('eff_cut_{0}'.format(i),
                    get_cuts_eff(histodict,self.decay_channel,i,pLcut,d0cut))
            self.__setattr__('p_{0}'.format(i),
                    get_final_state_pr(histodict,self.decay_channel,i))

    def get_total_fraction(self):
        """
        """
        if not hasattr(self,'eff_cut_{0}'.format(self.final_state_hadrons[0])):
            raise AttributeError('eff_total ERROR: You need to call the '\
                    'set_total_eff method before try to calculate the total fraction')
        
        denominator = sum(map(lambda x: 
            self.__getattribute__('eff_cut_{0}'.format(x))*
                         self.__getattribute__('p_{0}'.format(x)), 
            self.final_state_hadrons))

        return ((self.eff_cut_KK)*self.p_KK)/denominator

#def print_events(eff,decay_channel):
#    print 'Total events:',
#    for h in eff[decay_channel].final_state_hadrons:
#        current_n =eff[decay_channel].get_events(h)        
#        print ' {0}: {1:.4f}'.format(h,current_n),
#        print ' Total: {0:.4f}'.format(eff[decay_channel].get_total_events())
#        print float(eff[signal].get_total_events())/sqrt(float(sum(map(lambda (x,y): y.get_total_events(),\
#                filter(lambda (x,y): x != signal, eff.iteritems() )))))


def update_limits(d0_pLlist_dict,pl_attr_inst):
    """
    """
    getfunc = lambda func,i: func(map(lambda x: x[i], plList))
    # get all the values
    min_x0 = 1e10
    max_x1 = 0.0
    min_y0 = 1e10
    max_y1 = 0.0
    for plList in d0_pLlist_dict.values():
        x0 = getfunc(min,0)
        if x0 < min_x0:
            min_x0 = x0*0.98
        x1 = getfunc(max,0)
        if x1 > max_x1:
            max_x1 = x1*1.02
        y0 = getfunc(min,1)
        if y0 < min_y0:
            min_y0 = y0*0.98
        y1 = getfunc(max,1)
        if y1 > max_y1:
            max_y1 = y1*1.02

    if not hasattr(pl_attr_inst,'x0'):
        pl_attr_inst.x0 = min_x0
        
    if not hasattr(pl_attr_inst,'x1'):
        pl_attr_inst.x1 = max_x1

    if not hasattr(pl_attr_inst,'y0'):
        pl_attr_inst.y0 = min_y0
        
    if not hasattr(pl_attr_inst,'y1'):
        pl_attr_inst.y1 = max_y1


class plot_attributes(object):
    """
    """
    def __init__(self,plotname,**kwd):
        """
        """
        self.plotname = plotname
        self.xtitle   = ''
        self.xunit    = ''
        self.ytitle   = ''
        self.yunit    = ''
        self.title    = ''
        self.log      = False
        self.opt      = ''
        for key,val in kwd.iteritems():
            setattr(self,key,val)

def make_plot(points_dict,plot_attr,leg_position="UP"):
    """
    """
    import ROOT
    
    suffixplots = globals()['SUFFIXPLOTS']
    
    # Plotting
    ROOT.gROOT.SetBatch()

    try:
        from PyAnUtils.plotstyles import squaredStyle,setpalette
        lstyle = squaredStyle()
        lstyle.cd()
        ROOT.gROOT.ForceStyle()
        ROOT.gStyle.SetOptStat(0)

        setpalette("darkbody")
    except ImportError:
        pass

    c = ROOT.TCanvas()
    if plot_attr.log:
        c.SetLogy()
    # una grafica para el mismo d0
    # Extract the limits
    dummy = update_limits(points_dict,plot_attr)
    
    # Legend
    if leg_position == "DOWN":
        leg = getleg(y0 = 0.2, y1 = 0.34, x0 = 0.4, x1= 0.65)
    else:
        leg = getleg()

    frame = c.DrawFrame(plot_attr.x0,plot_attr.y0,plot_attr.x1,plot_attr.y1)
    frame.SetXTitle(plot_attr.xtitle+' '+plot_attr.xunit)
    frame.SetYTitle(plot_attr.ytitle+' '+plot_attr.yunit)
    k = 0
    d_graphs = []
    for d0cut,pL_list in points_dict.iteritems():
        g = ROOT.TGraph(len(pL_list))
        d_graphs.append(g)
        g.SetLineWidth(2)
        g.SetLineColor(getcolor()[k])
        g.SetLineStyle(3*k)
        g.SetMarkerColor(getcolor()[k])
        g.SetMarkerStyle(20+k)
        leg.AddEntry(g,'d_{0}={1} [mm]'.format('{0}',d0cut),'PL')
        for i,val in enumerate(pL_list):
            g.SetPoint(i,val[0],val[1])
        g.Draw("LSAME {0}".format(plot_attr.opt))
        k+=1
    leg.Draw()
    c.SaveAs('{0}{1}'.format(plot_attr.plotname.split('.')[0],suffixplots))

def main(rootfile,channels):
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
    import sys
    from math import sqrt

    
    # Get the root file with and the TH3 object (only)
    f = ROOT.TFile(rootfile)
    if f.IsZombie():
        raise IOError("ROOT file '{0}' doesn\'t exist".format(rootfile))
    _preobj = dict(map(lambda x:(x.GetName(), f.Get(x.GetName())),
         filter(lambda x: x.GetClassName().find('TH3') == 0,f.GetListOfKeys())))
    
    # Split by resonance
    _obj = { 'H' : dict(filter(lambda (x,y): x.find('H') == 0,_preobj.iteritems())),
            'Z' : dict(filter(lambda (x,y): x.find('Z') == 0,_preobj.iteritems()))
            }

    #signal=filter(lambda x: x.find('ssbar') != -1,channels)[0]
    signal_PID   = filter(lambda x: x.find('ssbar_PID') != -1,channels)[0]
    signal_noPID = filter(lambda x: x.find('ssbar_noPID') != -1,channels)[0]
    
    # just take care in the H-ressonance...

    # XXX: can be entered by option
    d0cuts = [0.1,0.3, 0.5]
    pLcuts = xrange(0,30)

    # the eff. classes 
    eff = dict(map(lambda x: (x,eff_cut_hadron(x)), channels))
    
    message = "\r\033[1;34mhzplots INFO\033[1;m Obtaining the data"
    # List of plots: { 'd0cut1' [ (pLcut1,valueX,valueY), ... ] }
    effsig_vs_pion_reject = {}
    effsig_vs_purity         = {}
    significance_vs_pLcut    = {}
    for d0 in d0cuts:
        d0str = '{0:.1f}'.format(d0)
        effsig_vs_pion_reject[d0str] = []
        effsig_vs_purity[d0str] = []
        significance_vs_pLcut[d0str] = []
        i = 0
        for pL in pLcuts:
            sys.stdout.write( "{0} {1}".format(message,SPINNING[i % len(SPINNING)]))
            sys.stdout.flush()

            #print '='*80
            #print 'd0 = {0:.1f} [mm] , pL = {1:.1f} [GeV] ======='.format(d0,pL)
            # setting the current cuts to all the efficienciesa
            _dummy = map( lambda e: e.set_total_eff(_obj['H'],pL,d0), eff.values())
            
            # Some needed values
            eff_sig = eff['ssbar_PID'].get_total_eff('KK')
            n_KK    = eff['ssbar_PID'].get_total_events()
            bkg_tot_evts = sum(map(lambda (x,y): y.get_total_events(),\
                    filter(lambda (x,y): x != signal_PID or x != signal_noPID, eff.iteritems() )))
            
            # Signal efficiency calculation: assuming 100% p-K separation with the PID 
            # effsig := eff['ssbar_PID']
            # purity := 2*N_KK/(2*N_KK+2*N_KP+2*N_KK) = N_KK/(N_KK+NKP+N_PP)
            # Pion rejection: N_b/N_signal assuming no p-K separation (no PID)
            # Note: sigma_HZ*L_int*Sum_qq BR(H->qq)/Sum_qq BR(h->qq)*(Sum_qq BR(H->qq) eff_qq*f_qq)
            # the Sum_qq BR(H->qq) terms are cancelled... 
            effsig_vs_pion_reject[d0str].append( (eff_sig,
                    float(bkg_tot_evts)/float(n_KK)) )
            effsig_vs_purity[d0str].append( (eff_sig,
                    float(2.0*eff['ssbar_noPID'].get_events('KK'))/\
                            (float(2.0*eff['ssbar_noPID'].get_events('KK'))+\
                            float(eff['ssbar_noPID'].get_events('KP'))+\
                            float(2.0*eff['ssbar_noPID'].get_events('PP')))
                            ) )
            significance_vs_pLcut[d0str].append( (pL,float(n_KK)/sqrt(float(bkg_tot_evts))) )
            i+=1
    # plotting
    print "\033[1;34mhzplots INFO\033[1;m Plotting..."
    pr_attr = plot_attributes('pion_rejection',
            xtitle='#varepsilon_{S}', 
            ytitle='pion rejection (N_{BKG}/N_{s#bar{s}} with no PID)',
            x0 = 0.0, x1 = 1.0 , y0 = 0.0 )
    make_plot(effsig_vs_pion_reject,pr_attr)

    purity_attr = plot_attributes('purity',
            xtitle='#varepsilon_{S}', 
            ytitle='purity',
            x0 = 0.0, x1 = 1.0 , y0 = 0.0, y1= 0.0 )
    make_plot(effsig_vs_purity,purity_attr)

    significance_attr = plot_attributes('significance',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='N_{S}/#sqrt{N_{B}}',
            x0 = 0.0, y0 = 0.0 )
    make_plot(significance_vs_pLcut,significance_attr,"DOWN")


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
    parser.add_option( '-p', '--pid', action='store_true', dest='pid',\
            help="Assume an efficiency of 100% distinguishing between kaons and pions")
    
    (opt,args) = parser.parse_args()

    if opt.suffixout:
        suff = opt.suffixout.replace('.','')
        globals()['SUFFIXPLOTS'] ='.'+suff
    channels = [ 'ssbar_PID', 'ssbar_noPID', 'bbbar', 'ccbar'] #'ddbar', 'uubar']
    #if opt.pid:
    #    channels.append('ssbar_PID')
    #else:
    #    channels.append('ssbar_noPID')

    main(os.path.abspath(opt.inputfile),channels)
