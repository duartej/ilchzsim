#!/usr/bin/env python
""":script:`hzevaleff_and_plots` -- Plotting stuff created with processhzroot script
====================================================================================

.. script:: hzevaleff_and_plots [OPTIONS]    
      :platform: Unix
      :synopsis: Evaluate pion misid...
.. moduleauthor:: Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
"""
KAON = 321
PION = 211
KSHORT = 310

def isKaon(index):
    """check if the gen particle is kaon

    Parameters
    ----------
    index: int
        the index of the gen-particle in the tree
    """
    return "abs(pdgId[{0}]) == {1}".format(index,KAON)

def isKshort(index):
    """check if the gen particle is K_short

    Parameters
    ----------
    index: int
        the index of the gen-particle in the tree
    """
    return "isKshort[{0}] == 1".format(index)

def isPion(index):
    """check if the gen particle is kaon

    Parameters
    ----------
    index: int
        the index of the gen-particle in the tree
    """
    return "abs(pdgId[{0}]) == {1}".format(index,PION)

# split charged-neutral channels in CN and NC to apply different cuts on neutral and charged part
# while keeping the cutting function simple
FS_CONDITION = { 'KK': '{0} && {1}'.format(isKaon(0),isKaon(1)),
                 'KP': '({0} && {1}) || ({2} && {3})'.format(isKaon(0),isPion(1),isPion(0),isKaon(1)),
                 'PP': '{0} && {1}'.format(isPion(0),isPion(1)),
                 'KsKs': '{0} && {1}'.format(isKshort(0),isKshort(1)),
                 'KsK': '{0} && {1}'.format(isKshort(0),isKaon(1)),
                 'KKs': '{0} && {1}'.format(isKaon(0),isKshort(1)),
                 'KsP': '{0} && {1}'.format(isKshort(0),isPion(1)),
                 'PKs': '{0} && {1}'.format(isPion(0),isKshort(1))
                }
SUFFIXPLOTS='.pdf'
SPINNING = [ "-","\\","|","/"]


def getcolor():
    import ROOT
    return [ ROOT.kRed+4, ROOT.kAzure+3, ROOT.kOrange-2, ROOT.kGreen-5, ROOT.kYellow+2, \
        ROOT.kCyan-2, ROOT.kOrange+5,ROOT.kAzure-7,ROOT.kGreen-2,ROOT.kRed-4, ROOT.kGray-3 ]

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
        self.mH         = 125.09 # (GeV)
        self.eeToHat240 = 230 # (fb) (at 240 GeV center of mass, is almost equivalent to ee->HZ)
        # Generic, only dependent of the higgs mass
        # therefore, following values at self.mH higgs
        # mass
        self.brbbbar    = 5.66e-1
        self.brccbar    = 2.85e-2
        self.brssbar    = 2.41e-4
        self.brgg       = 8.50e-2
        self.bruubar    = 1.3e-7  # needed for number of events but actually negligible
        self.brddbar    = 5.8e-8  # needed for number of events but actually negligible
        self.Lint       = 0 #(fb)
        self.higgsproduced = 0

    def setLint(self,_Lint):
        self.Lint = float(_Lint)
        self.higgsproduced = self.eeToHat240*self.Lint

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
hiInstance.setLint(43478) # (fb-1)

# Branching ratios with respect ccbar (m_q^2/m_c^2)
# using the pole mass (measured) not the mass at the higgs scale
BRuu_cc = 3.3e-6
BRdd_cc = 1.4e-5
BRss_cc = 0.006
BRbb_cc = 10.75
BRgg_cc = 2.98 # ??--> just using values from higgsclass, how can I
               #       if m_g = 0?

BRuu_ss = 5.9e-4
BRdd_ss = 2.6e-3
BRcc_ss = 180.1
BRbb_ss = 1936.0
BTgg_ss = 352.7

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
    elif name.find("gg") != -1:
        return BRgg_cc
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
    availpr = [ 'bbbar', 'uubar', 'ddbar', 'ssbar', 'ccbar', 'gg' ]
    try:
        process = filter(lambda x: name.find(x) != -1,availpr)[0]
    except IndexError:
        raise RuntimeError("Not correctly parsed: '%s'" % name)
      
    return process

def get_tree_name(tree_names,decay_channel):
    """Given a standarized (TTree) tree list of names (obtained from the 
    processedhzroot.py script), it returns the name of the histogram matching
    the final state

    Parameters
    ----------
    three_names: list(str)
        the names of all the trees found in a root processed file by the
        processedhzroot script
    
    Return
    ------
    str, the name of the input list matching the final state

    Raises
    ------
    AttributeError, whenever the tree 'mctrue_{channel}_*' is not found 
    """
    try:
        return filter(lambda x: x.find("mctrue_{0}".format(decay_channel)) != -1,\
                    tree_names)[0]
    except IndexError:
        raise AttributeError('Not found the tree "mctree_{0}_*"'\
                ' in the root file. Note that the input file should '\
                ' follows the standard notation: '\
                ' "mctrue_n-n-n-PID_channel_blahblbah" being nnn: mis-ident.'\
                ' probability'.format(decay_channel))

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
        the two opposite-hemisphere final state hadrons: KK, KP, PP, KsKs, KsK, 
        KKS, KsP, PKs (K-kaon, P-pion, K-short)

    Return
    ------
    str, the name of the input list matching the final state

    TO BE DEPRECATED
    """
    try:
        histo_name = filter(lambda x: 
                x.find('H_h3_pL_d0_d0_{0}'.format(hadron_state)) ==0 and  
                    x.find(decay_channel) != -1,th3names)[0]
    except IndexError:
        raise RuntimeError('Not found the histogram "H_h3_pL_d0_d0_{1}_*_{0}"'\
                ' in the root file'.format(decay_channel,hadron_state))
    return histo_name

def get_final_state_pr(tree,decay_channel,hadrons):
    """Extract the probabiliy of having two hadrons of the asked type in each 
    hemisphere given a decay channel

    .. math::f_{AB}^{q\bar{q}} = P( L_{AB} | (ee\rightarrow HZ\rightarrow q\bar{q}) I_0 )

    Parameters
    ----------
    tree : ROOT.TTree
        the tree
    decay_channel: str
    hadrons: str
        the hadrons to look at, which can be the: KK, KP, PP (K-kaon, P-pion)

    Returns
    -------
    efficiency: float

    TO BE DEPRECATED
    """
    n_current_hadrons = tree.GetEntries(FS_CONDITION[hadrons])
    n_total_decay = tree.GetEntries()

    return float(n_current_hadrons)/float(n_total_decay)

def get_bin(h,order,value):
    """Obtain the bin corresponding to the value given by `value`
    for the coordinate given by `order`

    Parameters
    ----------
    h: ROOT.TH3 (ROOT.TH2)
    order: int
        coordinate, valid values are 0 for X, 1 for Y and 2 for Z
    value: float
        value to check

    Return
    ------
    int, the bin corresponding to the `value`

    NOT NEEDED??
    """
    # TH1 histos to obtain bin numbers and other manipulations
    if order == 0:
        projection  = h.ProjectionX()
    elif order == 1:
        projection  = h.ProjectionY()
    elif order == 2:
        projection  = h.ProjectionZ()
    else:
        raise RuntimeError("get_bin: Invalid coordinate '{0}'".format(order))
    
    return projection.FindBin(value)

def get_cuts_eff(_obj,decay_channel,hadrons,\
        pLcut=None,d0cut=None,pLbin=None,d0binL=None,d0binH=None):
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
        the hadrons to look at, which can be the: KK, KP, PP, KsKs, KsK, KKs, KsP,
        PKs (K-kaon, P-pion, Ks-K_short)
    pLcut: float, optional, incompatible with pLbin  
    d0cut: float, optional, incompatible with d0binL, d0binH
    pLbin: int, optional, incompatible with pLcut
    d0binH: int, optional, incompatible with d0cut
    d0binL: int, optional, incompatible with d0cut

    Returns
    -------
    efficiency: float

    TO BE DEPRECATED
    """
    histo_name = get_histo_name(_obj.keys(),decay_channel,hadrons)
    h = _obj[histo_name]

    # Should it get the bin information?
    if pLbin is None and d0binL is None and d0binH is None:
        # First obtaining the bin where the cut are located
        pLbin  = get_bin(h,0,pLcut)
        d01binH= get_bin(h,1,d0cut)
        d01binL= get_bin(h,1,-1*d0cut)
        d02binH= get_bin(h,2,d0cut)
        d02binL= get_bin(h,2,-1*d0cut)
    elif pLcut is None and d0cut is None:
        d01binH = d0binH
        d01binL = d0binL
        d02binH = d0binH
        d02binL = d0binL
    else:
        raise RuntimeError('Unexpected ERROR! "get_cuts_eff" function'\
                ' called unconsistently.')
    
    # counting all the events, including those in the overflow bins
    pL_maxBin = h.GetNbinsX()+1

    evts = h.Integral(pLbin,pL_maxBin,d01binL,d01binH,d02binL,d02binH)

    if h.GetEntries() == 0: 
        return 0.0

    # Entries takes into account overflow bin
    return float(evts)/float(h.GetEntries())

#def smearing_value(func):
#    """
#    """
#    def resolution(*args,**kwargs):
#        import numpy as np 
#        
#        # d0- sigma calculation 
#        sigma_d0 = "5.0+10.0/(p_lab{[0]}*sin(theta_lab[{0}])**(3./2.))".format(args[0])
#
#        d0_smeared_at1 = np.random.normal(1.0, sigma)
#        # decorate the function with the resolution
#        func.resolution = d0_smeared-at1
#        return resolution
#

def d0_explicit_calculation(vx,vy,phi):
    """Impact parameter extrapolation considering a straight line.
    Explicit calculation.

    Parameters
    ----------
    vx: float
        the x-point of the production vertex
    vy: float
        the y-point of the production vertex
    phi: float
        the angle (in the X-Y plane) defined by the momentum 
        of the particle
    """
    import numpy as np
    return np.sin(np.arctan2(vy,vx)-phi)*np.sqrt(vx**2.+vy**2.)

def d0(index):
    """Impact parameter extrapolation considering a straight line

    Parameters
    ----------
    index: int
        the index of the gen-particle in the tree

    See Also
    --------
    d0_explicit_calculation
    """
    #return  "(vy[{0}]-vx[{0}]*tan(phi_lab[{0}]))*cos(phi_lab[{0}])".format(index)
    #return "sin(atan(vy[{0}]/vx[{0}])-phi_lab[{0}])*sqrt(vx[{0}]**2+vy[{0}]**2)".format(index)
    return "d0[{0}]".format(index)

def R(index):
    """decay length

    Parameters
    ----------
    index: int
        the index of the gen-particle in the tree

    """
    return "R[{0}]".format(index)

def z0_explicit_calculation(vx,vy,phi,vz,theta):
    """Longitudinal impact parameter extrapolation considering a straight line.
    Explicit calculation.

    Parameters
    ----------
    vx: float
        the x-point of the production vertex
    vy: float
        the y-point of the production vertex
    phi: float
        the azimuthal angle (in the X-Y plane) defined by the momentum 
        of the particle
    vz: float
        the z-point of the production vertex
    theta: float
        the angle (in the YZ plane) defined by the momentum of the
        particle
    """
    import numpy as np
    return (d0_explicit_calculation(vx,vy,phi)-np.sqrt(vx**2.+vy**2.))/(np.tan(theta)) + \
            vz

def z0(index):
    """Impact parameter extrapolation considering a straight line

    Parameters
    ----------
    index: int
        the index of the gen-particle in the tree
    See Also
    --------
    d0_explicit_calculation
    """
    return "z0[{0}]".format(index)

def get_common_entries(entrylist_list):
    """Return the common entries found in the list of TEntryList

    Parameters
    ----------
    entrylist_list: list(ROOT.TEntryList)

    Returns
    -------
    common_el: ROOT.TEntryList
    """
    #import ROOT

    # convert in sets, then is easy to obtain the intersection of entries
    entries_sets = []
    for _oel in entrylist_list:
        # is it possible to speed up this?
        entries_sets.append( set(map(lambda i: _oel.GetEntry(i),xrange(_oel.GetN()))) )
    common_evt = set.intersection( *entries_sets )
    # and building a new TEntryList --> need the tree..., so promote this function
    # to the eff class
    #ROOT.TEntryList()
    return len(common_evt)
    
class eff_cut_hadron(object):
    """ Encapsulates the efficiency cut of an hadron-hadron event. 
    Incorporates the cut in impact parameter and in parallel momentum

    .. math::\varepsilon_{AB}^{q\bar{q}} = P( d0^{c} p_{||}^c L_{AB} | 
                (ee\rightarrow HZ\rightarrow q\bar{q}) I_0 )
    """
    def __init__(self,decay_channel,z0cut=False,Ntrackcut=None,d0cut_type="circular",
                 pLcut_type="circular",R0=500,charge_combination='all'):
        """Encapsulates the efficieny cut of a Hadron-hadron event
        Incorporates the cut in impact parameter and in parallel momentum

        Parameters
        ----------
        decay_channel: str
            the Higgs hadronic decay (gg, bbbar,ccbar,ssbar, ddbar, uubar)
        z0cut: bool
            whether or not activated an extra circular cut around z0=0.1 mm
        Ntrackcut: int
            whether or not activated an extra cut in the hadron multiplicity 
            and the maximum allowed tracks surrounding a dR cone of the hadron
        d0cut_type: str, [default: circular]
            defines the function to be used to cut in d0, valid values are
            circular, square
        pLcut_type: str, [default: circular]
            defines the function to be used to cut in parallel momentum, 
            valid values are circular, square and line
        R0: float
            the decay length cut for K_s
        charge_combination: str
            the combination of final state charges to be considered

        FIXME: The class should be associated to a tree!
        """
        self.decay_channel = decay_channel
        self.short_name = decay_channel.split('_')[0]
        self.initialized = False
        self.current_d0cut = None
        self.current_pLcut = None
        
        self.final_state_hadrons = ['KK','KP','PP', 'KsKs', 'KsK', 'KKs', 'KsP', 'PKs']
        self.final_state_NN = ['KsKs']
        self.final_state_CN = ['KKs', 'PKs']
        self.final_state_NC = ['KsK', 'KsP']
        self.final_state_CC = ['KK', 'KP', 'PP']
        self.__entrylist_hadrons = dict(map(lambda i: (i,None),self.final_state_hadrons))
        
        # optimization data-members
        self.__tree = None
        self.__entrylist_d0cut = None
        self.__entrylist_pLcut = None
        # { treename: { 'pLcut': { 'actual cut' : ROOT.TEntryList, ...
        self.__entrylists_reservoir = {}
        self.__finalhadrons_entries = { 'KK': None, 'KP': None, 'PP': None,
                                        'KsKs': None, 'KsK': None, 'KsP': None,
                                        'KKs': None, 'PKs': None }


        if charge_combination not in ['CC', 'NC', 'NN', 'all']:
            raise RuntimeError('{0} is not a valid charge combination. Valid values: CC, NC, NN, all'.format(charge_combination))
        self.cuts = {}
        self.cuts['NN'] = '( ({0}) &&  ({1}))'.format(isKshort(0),isKshort(1))
        self.cuts['NC'] = '( ({0}) && !({1}))'.format(isKshort(0),isKshort(1))
        self.cuts['CN'] = '(!({0}) &&  ({1}))'.format(isKshort(0),isKshort(1))
        self.cuts['CC'] = '(!({0}) && !({1}))'.format(isKshort(0),isKshort(1))
        self.cuts['all'] = ''

        if z0cut:
            self.cuts['CC'] += "&& sqrt( ({0})**2.0 + ({1})**2.0) < {2} ".format(z0(0),z0(1),z0cut)
            self.cuts['NC'] += "&& abs({0}) < {1} ".format(z0(1),z0cut)
            self.cuts['CN'] += "&& abs({0}) < {1} ".format(z0(0),z0cut)
            
        if Ntrackcut:
            self.cuts['CC'] += "&& multiplicity[0] <= {0} && multiplicity[1] <= {0} ".format(Ntrackcut)
            self.cuts['NC'] += "&& multiplicity[0] <= {0} && multiplicity[1] <= {0} ".format(Ntrackcut)
            self.cuts['CN'] += "&& multiplicity[0] <= {0} && multiplicity[1] <= {0} ".format(Ntrackcut)
            self.cuts['NN'] += "&& multiplicity[0] <= {0} && multiplicity[1] <= {0} ".format(Ntrackcut)

        # behaviour members
        if pLcut_type == "circular":
            self.pLcut_function = "sqrt( (p[0]*cos(theta[0]))**2. + (p[1]*cos(theta[1]))**2.) > {0}"
        elif pLcut_type == "square":
            self.pLcut_function = "abs(p[0]*cos(theta[0])) > {0} && abs(p[1]*cos(theta[1])) > {0}"
        elif pLcut_type == "line":
            self.pLcut_function = "p[1]*cos(theta[1]) > (-{0}/30.0)*p[0]*cos(theta[0])+{0}"
        else:
            raise RuntimeError("Not a valid pL-cut functional: valid values:"\
                    " circular, square and line")

        if charge_combination == 'CC':
            self.pLcut_function = "({0}) && {1}".format(self.cuts['CC'],self.pLcut_function)
        elif charge_combination == 'NC':
            self.pLcut_function = "( ({0}) || ({1}) ) && {2}".format(self.cuts['NC'],self.cuts['CN'],self.pLcut_function)
        elif charge_combination == 'NN':
            self.pLcut_function = "({0}) && {1}".format(self.cuts['NN'],self.pLcut_function)
        elif charge_combination == 'all':
            self.pLcut_function = '( ({0}) || ({1}) || ({2}) || ({3}) ) && {4}'.format(self.cuts['CC'],self.cuts['NC'],self.cuts['CN'],self.cuts['NN'],self.pLcut_function)

        if d0cut_type == "circular":
            self.cuts['CC'] += "&& sqrt( ({0})**2.0 + ({1})**2.0) < {2} ".format(d0(0),d0(1),"{0}")
        elif d0cut_type == "square":
            self.cuts['CC'] += "&& ( abs({0}) < {2} && abs({1}) < {2}) ".format(d0(0),d0(1),"{0}")
        else:
            raise RuntimeError("Not a valid d0-cut functional: valid values: circular and square ")
        self.cuts['NC'] += "&& ( abs({0}) < {1} && {2} > {3}) ".format(d0(1),"{0}",R(0),R0)
        self.cuts['CN'] += "&& ( abs({0}) < {1} && {2} > {3}) ".format(d0(0),"{0}",R(1),R0)
        self.cuts['NN'] += "&& ( {0} > {2} && {1} > {2} ) ".format(R(0),R(1),R0)
        if charge_combination == 'CC':
            self.d0cut_function = self.cuts['CC']
        elif charge_combination == 'NC':
            self.d0cut_function = "({0}) || ({1})".format(self.cuts['NC'],self.cuts['CN'])
        elif charge_combination == 'NN':
            self.d0cut_function = "({0})".format(self.cuts['NN'])
        elif charge_combination == 'all':
            self.d0cut_function = "({0}) || ({1}) || ({2}) || ({3})".format(self.cuts['CC'],self.cuts['CN'],self.cuts['NC'],self.cuts['NN'])
        else:
            raise RuntimeError("Not a valid charge_combination: valid values: all|NN|NC|CC")

        print "\033[1;34mhzplots INFO\033[1;m Using the following functional cuts for the charge combination: {0}:".format(charge_combination)
        print "   + pL: {0}".format(self.pLcut_function)
        print "   + d0: {0}".format(self.d0cut_function)


    def get_tree(self):
        """Returns the tree associated with this efficiency
        """
        return self.__tree

    def deactivate_cuts(self):
        """De-activate cuts, the Tree has access to the whole number of 
        entries
        """
        self.__tree.SetEntryList(0)

    def activate_cuts(self,pLcut=None,d0cut=None,hadron_pairs=None,**future_cuts):
        """Activate cuts in the associated TTree, having access only to those
        entries fulfilling the cuts

        Parameters
        ----------
        pLcut: float, optional
            the value of the parallel momentum cut
        d0cut: float, optional
            the value of the impact parameter cut
        hadron_pairs: str
            the leading hadrons final states: KK, KP, PP, KsKs, KsK, KsP, KKs, PKs
        future_cuts: dict()
            not used now, but will allow to incorporate extra cuts

        FIXME: JUST one TEntryList per time, do not intersects it
        """
        print "\033[1;33mWARNING activate_cuts\033[1;m function is"\
                " still not allowing more than one entrylist!".format(pLcut)
        treename = self.__tree.GetName()
        # search the TEntryList to be activatead
        if pLcut:
            el_pLname = "{0}_entrylist_pLcut_{1}".format(treename,pLcut)
            try:
                entrylist = self.__entrylists_reservoir[treename][el_pLname]
                self.__tree.SetEntryList(entrylist)
            except KeyError:
                print "\033[1;33mWARNING\033[1;m the pL={0} [GeV]"\
                        " was not used as cut, the TTree remains as it was".format(pLcut)
        if d0cut:
            el_d0name = "{0}_entrylist_d0cut_{1}".format(treename,d0cut)
            try:
                entrylist = self.__entrylists_reservoir[treename][el_d0name]
            except KeyError:
                print "\033[1;33mWARNING\033[1;m the d0={0} [mm]"\
                        " was not used as cut, the TTree remains as it was".format(d0cut)

        if hadron_pairs:
            try:
                entrylist = self.__entrylist_hadrons[hadron_pair]
            except KeyError:
                print "\033[1;33mWARNING\033[1;m the '{0}"\
                        " was not used as cut, the TTree remains as it was".format(hadron_pair)

        self.__tree.SetEntryList(entrylist)

    
    def __str__(self):
        """string representation

        Returns
        -------
        out: str
        """
        out = 'eff_cut_hadron instance: \n'
        out+= '  decay channel: {0} \n'.format(self.decay_channel)
        
        eff = '  d0,pL cut eff:            '
        p   = '  H1_H2 final states prob:  '
        if self.initialized:
            out+= '  d0cut: {0:.1f} [mm];  pLcut:{1:.0f} [GeV]\n'.format(self.current_d0cut,self.current_pLcut)
            for i in self.final_state_hadrons:
                eff+= '{0}={1:.4f} '.format(i,getattr(self,'eff_cut_{0}'.format(i)))
                p  += '{0}={1:.4f} '.format(i,getattr(self,'p_{0}'.format(i)))
            out += eff+'\n'
            out += p+'\n'
        return out

#    def set_total_eff(self,treedict,pLcut,d0cut,verbose=False):
#        """builds the total efficiency and evaluates it for the different 
#        final state cases KK, KP and PP (if exist)
#
#        Parameters
#        ----------
#        treedict: ROOT.TTree
#            the tree where to extract the info
#        pLcut: float
#            the parallel momentum cut
#        d0cut: float
#            the impact parameter cut
#
#        Implementation Notes
#        --------------------
#        Note that the ROOT.TEntryList should have one unique name, otherwise
#        the ROOT.gDirectory.Get method retrieve the first object in memory, 
#        obtaining wrong results
#        """
#        treename = get_tree_name(treedict.keys(),self.decay_channel)
#        # get the tree and check if it was used before
#        tree = treedict[treename]
#        is_new_tree = False
#        if self.__tree != tree:
#            self.__tree = tree
#            self.__tree_entries = self.__tree.GetEntries()
#            # new EntryList reservoir
#            if not self.__entrylists_reservoir.has_key(treename):
#                self.__entrylists_reservoir[treename] = {}
#            # create the entrylist for final state hadrons final state
#            for i in self.final_state_hadrons:
#                el_HSname = "{0}_entrylist_hadrons_{1}".format(treename,i)
#                print el_HSname,FS_CONDITION[i]
#            raise
#                #self.__tree.Draw(">>{0}".format(el_HSname),FS_CONDITION[i],"entrylist")
#                #self.__entrylist_hadrons[i] = ROOT.gDirectory.Get("{0}".format(el_HSname))
#                #self.__finalhadrons_entries[i] = self.__entrylist_hadrons[i].GetN()
#            is_new_tree = True

    def set_total_eff(self,treedict,pLcut,d0cut,verbose=False):
        """builds the total efficiency and evaluates it for the different 
        final state cases KK, KP, PP, KsKs, KsK, KsP, KKs, and PKs (if exist)

        Parameters
        ----------
        treedict: ROOT.TTree
            the tree where to extract the info
        pLcut: float
            the parallel momentum cut
        d0cut: float
            the impact parameter cut

        Implementation Notes
        --------------------
        Note that the ROOT.TEntryList should have one unique name, otherwise
        the ROOT.gDirectory.Get method retrieve the first object in memory, 
        obtaining wrong results
        """
        import ROOT 
        #ROOT.gROOT.SetBatch()

        treename = get_tree_name(treedict.keys(),self.decay_channel)
        # get the tree and check if it was used before
        tree = treedict[treename]
        is_new_tree = False
        if self.__tree != tree:
            self.__tree = tree
            self.__tree_entries = self.__tree.GetEntries()
            # new EntryList reservoir
            if not self.__entrylists_reservoir.has_key(treename):
                self.__entrylists_reservoir[treename] = {}
            # create the entrylist for final state hadrons final state
            for i in self.final_state_hadrons:
                el_HSname = "{0}_entrylist_hadrons_{1}".format(treename,i)
                self.__tree.Draw(">>{0}".format(el_HSname),FS_CONDITION[i],"entrylist")
                self.__entrylist_hadrons[i] = ROOT.gDirectory.Get("{0}".format(el_HSname))
                self.__finalhadrons_entries[i] = self.__entrylist_hadrons[i].GetN()
            el_BCname1 = '{0}_entrylist_BCdaughter1'.format(treename)
            self.__tree.Draw(">>{0}".format(el_BCname1), "isBCdaughter[0] == 1", "entrylist")
            self.__entrylist_BCdaughters1 = ROOT.gDirectory.Get("{0}".format(el_BCname1))
            el_BCname2 = '{0}_entrylist_BCdaughter2'.format(treename)
            self.__tree.Draw(">>{0}".format(el_BCname2), "isBCdaughter[1] == 1", "entrylist")
            self.__entrylist_BCdaughters2 = ROOT.gDirectory.Get("{0}".format(el_BCname2))
            is_new_tree = True

        is_new_d0cut = False
        ## allow optimize loop
        if not self.current_d0cut or self.current_d0cut != d0cut:
            self.current_d0cut = d0cut
            is_new_d0cut = True
        elif is_new_tree:
            is_new_d0cut = True

        if is_new_d0cut:
            el_d0name = "{0}_entrylist_d0cut_{1}".format(treename,self.current_d0cut)
            try:
                self.__entrylist_d0cut = self.__entrylists_reservoir[treename][el_d0name]
            except KeyError:
                # create the entrylist 
                self.__tree.Draw(">>{0}".format(el_d0name),\
                        self.d0cut_function.format(self.current_d0cut),"entrylist")
                self.__entrylists_reservoir[treename][el_d0name]= ROOT.gDirectory.Get("{0}".format(el_d0name))
                self.__entrylist_d0cut = self.__entrylists_reservoir[treename][el_d0name]
                
        is_new_pLcut = False
        if not self.current_pLcut or self.current_pLcut != pLcut:
            self.current_pLcut = pLcut
            is_new_pLcut = True
        elif is_new_tree:
            is_new_pLcut = True
        
        if is_new_pLcut:
            el_pLname = "{0}_entrylist_pLcut_{1}".format(treename,self.current_pLcut)
            try: 
                self.__entrylist_pLcut = self.__entrylists_reservoir[treename][el_pLname]
            except KeyError:
                # create the entrylist (or extracted from memory if is already there)
                self.__tree.Draw(">>{0}".format(el_pLname),\
                        self.pLcut_function.format(self.current_pLcut),"entrylist")
                self.__entrylists_reservoir[treename][el_pLname]=ROOT.gDirectory.Get("{0}".format(el_pLname))
                self.__entrylist_pLcut = self.__entrylists_reservoir[treename][el_pLname]
        
        # main loop
        for i in self.final_state_hadrons:
            cut_entries = get_common_entries( [self.__entrylist_d0cut,\
                  self.__entrylist_pLcut, self.__entrylist_hadrons[i]] )
            cut_BCdaughter1 = get_common_entries( [self.__entrylist_d0cut,\
                self.__entrylist_pLcut, self.__entrylist_hadrons[i],\
                self.__entrylist_BCdaughters1] )
            cut_BCdaughter2 = get_common_entries( [self.__entrylist_d0cut,\
                self.__entrylist_pLcut, self.__entrylist_hadrons[i],\
                self.__entrylist_BCdaughters2] )
            _eff  = 0.0
            _BCdaughterFraction = 0.0
            if self.__tree_entries != 0:
                _eff = float(cut_entries)/float(self.__tree_entries)
            if float(cut_entries) != 0:
                _BCdaughterFraction = (float(cut_BCdaughter1) + float(cut_BCdaughter2))/(2*float(cut_entries))
            # p ( d0 pL | L_AB decay_channel I0 )
            self.__setattr__('eff_cut_{0}'.format(i),_eff)
            self.__setattr__('BCdaughter_fraction_{0}'.format(i),_BCdaughterFraction)
        self.initialized=True

    def _get_total_eff(self):
        """FIXME: USE THE WHOLE FORMULAE (sum of probabilities ?)
        FIXME
        """
        return sum(map(lambda x: self.get_total_eff(x),self.final_state_hadrons))

    def get_total_eff(self,hadrons=None):
        """
        """
        if not hadrons:
            return self._get_total_eff()
        return getattr(self,'eff_cut_{0}'.format(hadrons))

    def _get_BCdaughter_fraction(self):
        """ Return the fraction of BC daughters among the two selected kaon candidates
        """
        total_evts = self.get_total_events()
        if total_evts == 0:
            return 0
        else:
            return sum(map(lambda x: self.get_events(x)*self.get_BCdaughter_fraction(x), self.final_state_hadrons))/\
                self.get_total_events()

    def get_BCdaughter_fraction(self, hadrons=None):
        """
        """
        if not hadrons:
            return self._get_BCdaughter_fraction()
        return getattr(self,'BCdaughter_fraction_{0}'.format(hadrons))

    def get_BCdaughter_fraction_by_FS(self,fs):
        """
        """
        considered_fs = []
        if fs == 'NN':
            considered_fs = self.final_state_NN
        elif fs == 'NC':
            considered_fs = self.final_state_CN + self.final_state_NC
        elif fs == 'CC':
            considered_fs = self.final_state_CC
        total_evts=self.get_total_events_by_FS(fs)
        if total_evts == 0:
            return 0.0
        else:
            return sum(map(lambda x: self.get_events(x)*self.get_BCdaughter_fraction(x), considered_fs))/\
                total_evts
    
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
                getattr(self,'eff_cut_{0}'.format(hadrons))
        return N

    def get_total_events_by_FS(self,fs):
        """
        """
        if not self.initialized:
            raise AttributeError('eff_total ERROR: You need to call the '\
                    'set_total_eff method before try to calculate the total fraction')
        
        N = 0
        considered_fs = []
        if fs == 'NN':
            considered_fs = self.final_state_NN
        elif fs == 'NC':
            considered_fs = self.final_state_CN + self.final_state_NC
        elif fs == 'CC':
            considered_fs = self.final_state_CC
        for hadrons in considered_fs:
            N += hiInstance.getEvents(self.decay_channel.split('_')[0])*\
                 getattr(self,'eff_cut_{0}'.format(hadrons))
        return N

def update_limits(d0_pLlist_dict,(x_ind,y_ind),pl_attr_inst):
    """
    """
    getfunc = lambda func,i: func(map(lambda x: x[i], plList))
    # get all the values
    min_x0 = 1e10
    max_x1 = 0.0
    min_y0 = 1e10
    max_y1 = 0.0
    for plList in d0_pLlist_dict.values():
        x0 = getfunc(min,x_ind)
        if x0 < min_x0:
            min_x0 = x0*0.98
        x1 = getfunc(max,x_ind)
        if x1 > max_x1:
            max_x1 = x1*1.02
        y0 = getfunc(min,y_ind)
        if y0 < min_y0:
            min_y0 = y0*0.98
        y1 = getfunc(max,y_ind)
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
        self.xlog     = False
        self.opt      = ''
        for key,val in kwd.iteritems():
            setattr(self,key,val)

def get_point_graphs(points_dict,obs_indices,wp_index,working_points_list,leg_entry_format):
    """
    Parameters
    ----------
    points_dict: dict(str: [(val1,val2), ... ])
    """
    import ROOT

    # Indices to plot
    xI = obs_indices[0]
    yI = obs_indices[1]

    # indices for the legend points
    wpX_I = leg_entry_format[1][0]
    wpY_I = leg_entry_format[1][1]

    graphs = {}
    leg    = {}
    # legend coordinates
    y0 = 0.8; y1 = 0.94; x0 = 0.01; x1= 0.25;
    for d0cut,val_list in sorted(points_dict.iteritems(),key=lambda (x,y): float(x)):
        selected_points = []
        TOLERANCE = 1e-10
        # Extract the working points from the val_list list, assuming
        # tuple_index as the value to check and to be equal to the 
        # working_points_list values. Using a very small tolerance
        # and degrading resolution if unable to find
        n_iter = 0
        while len(selected_points) != len(working_points_list):
            try:
                selected_points =  map(lambda y: 
                        filter(lambda x: abs(x[wp_index]-y) <= y*TOLERANCE, val_list)[0], 
                        working_points_list)
            except IndexError:
                # not found anything, so continue
                pass

            if TOLERANCE > 1.0 or n_iter > 100:
                raise RuntimeError('Unable to find the Working Points: "{0}"'\
                        ' Some inconsistency in the call of this function...'.\
                        format(working_points_list))
            TOLERANCE *= 1e1
            n_iter += 1
        # found the working points
        graphs[d0cut] = []
        # legend
        leg[d0cut] = getleg(y0 = y0, y1 = y1, x0 = x0, x1= x1)
        y0 -= 0.2
        y1 -= 0.2
        leg[d0cut].SetTextSize(0.035)
        k = 4
        for ntuple in selected_points:
            graphs[d0cut].append( ROOT.TGraph() )
            graphs[d0cut][-1].SetMarkerStyle(33)
            graphs[d0cut][-1].SetMarkerSize(2)
            graphs[d0cut][-1].SetMarkerColor(getcolor()[k])
            graphs[d0cut][-1].SetPoint(0,ntuple[xI],ntuple[yI])
            text = leg_entry_format[0] % (ntuple[wpX_I],ntuple[wpY_I])
            leg[d0cut].AddEntry(graphs[d0cut][-1],'d0={0:.0f} #mum:  {1}'.format(float(d0cut)*1e3,text),'P')
            k+=1

    return graphs,leg


def make_plot(points_dict,(xindex,yindex),plot_attr,g_points_dict=None,leg_position="UP"):
    """Create a plot superimpossing the expected graphs in the points_dict
    dictionary (one graph per each key). If aux_dict is filled, is interpreted
    as working point dictionary and some points will be painted also
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
    dummy = update_limits(points_dict,(xindex,yindex),plot_attr)
    
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
    for d0cut,pL_list in sorted(points_dict.iteritems(),key=lambda (x,y): float(x)):
        g = ROOT.TGraph(len(pL_list))
        d_graphs.append(g)
        g.SetLineWidth(2)
        g.SetLineColor(getcolor()[k])
        #g.SetLineStyle(3*k)
        g.SetMarkerColor(getcolor()[k])
        g.SetMarkerStyle(20+k)
        leg.AddEntry(g,'d_{0}={1:.0f} #mum'.format('{0}',(float(d0cut)*1e3)),'PL')
        for i,val in enumerate(pL_list):
            g.SetPoint(i,val[xindex],val[yindex])
        g.Draw("LSAME {0}".format(plot_attr.opt))
        if g_points_dict:
            # Format: ( dict(str: [ROOT.TGraph, ...],..), ROOT.TLegend() )
            # see get_points_graph to check the format
            for _graph in g_points_dict[0][d0cut]:
                _graph.Draw("PSAME")
            #g_points_dict[1][d0cut].Draw()
        k+=1
    leg.Draw()
    c.SaveAs('{0}{1}'.format(plot_attr.plotname.split('.')[0],suffixplots))
    #if plot_attr.plotname == "significance":
    #    c.SaveAs('{0}{1}'.format(plot_attr.plotname.split('.')[0],".C"))
    # Legend for the points
    if g_points_dict:
        c.Close()
        c.Draw()
        for d0cut,leg in g_points_dict[1].iteritems():
            leg.Draw()
        c.SaveAs('{0}_legends_WP{1}'.format(plot_attr.plotname.split('.')[0],suffixplots))

def get_latex_table(obsList):
    """.. function:: get_latex_table(obsList) -> str(latex table)
    [ (pLcut1, eff_sig, significance, pion_rejection, purity, N_SIG, N_BKG) ]
    """
    columns = ( 'p_{||}^{c}\;[GeV]','\\varepsilon_{ss}', 'N_{S}/\\sqrt{N_{B}}', 'N_{bkg}/N_{s\\bar{s}}', \
            'Kaon Purity', 'N_{s\\bar{s}}', 'N_{bkg}', '\\varepsilon_{bkg}')
    colformat = map(lambda x: ' ${0:'+str(x)+'}$ &' , \
            ( '.0f','.3f','.2f','.0f','.3f', '.0f', '.0f' , '.3f'))
    ncols = len(obsList[0])
    # Consistency
    if len(columns) != ncols: 
        raise RuntimeError('Unconsistency in "observable" list')

    latex = '\\begin{tabular}{c '
    for i in xrange(ncols):
        latex += ' c '
    latex += '}\\hline\\hline\n'
    for i in columns:
        latex += ' ${0}$ &'.format(i)
    latex = latex[:-1]
    latex += '\\\\ \\hline\n'
    for ntuple in obsList:
        for i,val in enumerate(ntuple):
            latex += colformat[i].format(val)
        latex = latex[:-1]+'\\\\\n'
    latex += '\\hline\n'
    latex += '\\end{tabular}\n'

    return latex

def create_histos(suffix,description,res_int,hc=None):
    """Function gathering all the histograms we are interested
    to plot. So far the histograms defined are:
     * h2_pL_suffix : the parallel momentum of the leading 
                      and subleading hadrons
     * h2_[sub]leading_pL_d0_suffix: the impact parameter of the
                      [sub]leading kaon vs its parallel momentum
     * h2_d0_suffix : the impact parameter of the leading 
                      and subleading hadrons (extrapolated as
                      straight lines)
     * h2_Lxy_suffix: the vertex kaon production in the transverse
                      plane
     * h2_R_suffix  : the vertex kaon production
     * h2_pL_R_suffix: the decay length vs the parallel momentum
                      of the K_s
     * h_d0_suffix  : the impact parameter of the leading 
                      and subleading hadrons (extrapolated as
                      straight lines) in the same histogram (1D)

     * h_z0_suffix  : the longitudinal parameter of the leading 
                      and subleading hadrons (extrapolated as
                      straight lines) in the same histogram (1D)
     * h_pL_pions   : the parallel momentum distribution of truth pions
     * h_pL_kaons   : the parallel momentum distribution of truth kaons
     * h_pL_candidates: the parallel momentum distribution of the kaon 
                      candidates
     * h_plab_candidates: the lab frame momentum distribution of the kaon 
                      candidates
     * h_nM_suffix  : the quark multiplicity (related with the number
                      of constituents of a jet) [FIXME: not well defined]
     * XXX: THE MULTIPLICITY:: MISSING!!!

    Parameters
    ----------
    suffix: str
        used to distinguish the histograms of different samples
    description: str
        the legend to be used when several samples are plotted in 
        the same canvas
    res_int: int
        the considered resonance: 25:=higgs, 23:=Z
    hc: PyAnUtils.histocontainer.HistoContainer instance
        the histogram container gathering all the available histograms
    """
    from PyAnUtils.histocontainer import HistoContainer
    from PyAnUtils.plotstyles import squaredStyle
    import ROOT
    lstyle = squaredStyle()
    lstyle.cd()
    ROOT.gROOT.ForceStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gROOT.SetBatch()

    COLOR = { 'ssbar': 46, 'bbbar': 30, 'gg': 38,
              'ccbar': 14, 'uubar': 20,
              'ddbar': 28}
    RES = { 23: 'Z', 25: 'H' }

    resonance = RES[res_int]

    color = filter(lambda (name,c): suffix.find(name) != -1,COLOR.iteritems())[0][1]
    
    # -- create histo container, if there is no one
    if not hc:
        hc = HistoContainer()
    # -- some cosmethics
    # -- populate the hc with the histograms
    hc.create_and_book_histo("{0}_h2_pL_{1}".format(resonance,suffix),\
            "leading kaons parallel momentum",\
            100,0,65,npoints_y=100,ylow=0,yhigh=65,description=description,
            xtitle="leading-p_{||} [GeV]",ytitle='subleading-p_{||} [GeV]',
            color=color)
    hc.create_and_book_histo("{0}_h2_pL_d0_{1}".format(resonance,suffix),\
            "leading kaons parallel momentum vs impact parameter",\
            100,0,65,npoints_y=100,ylow=0,yhigh=0.1,description=description,
            xtitle="p_{||} [GeV]",ytitle='d_{0} [mm]',
            color=color)
    hc.create_and_book_histo("{0}_h2_pL_R_{1}".format(resonance,suffix),\
            "K_s parallel momentum vs decay length",\
            100,0,65,npoints_y=200,ylow=0,yhigh=1000,description=description,
            xtitle="p_{||} [GeV]",ytitle='R [mm]',
            color=color)
    hc.create_and_book_histo("{0}_h_d0_{1}".format(resonance,suffix),\
            "leading kaons impact parameter",100,-5,5,\
            description=description,
            xtitle="d_{0} [mm]",
            ytitle="A.U.",
            color=color)
    hc.create_and_book_histo("{0}_h_absd0_{1}".format(resonance,suffix),\
            "leading kaons impact parameter",20,10**(-3),10,\
            description=description,
            xtitle="|d_{0}| [mm]",
            ytitle="A.U.",
            color=color,
            xlog=True)
    hc.create_and_book_histo("{0}_h_pLcut15_absd0_{1}".format(resonance,suffix),\
            "leading kaons impact parameter",20,10**(-3),10,\
            description=description,
            xtitle="|d_{0}| [mm]",
            ytitle="A.U.",
            color=color,
            xlog=True)
    hc.create_and_book_histo("{0}_h_Resd0_{1}".format(resonance,suffix),\
            "resolution of leading kaons impact parameter",20,0,30,\
            description=description,
            xtitle="#sigma(d_{0}) [#mum]",
            ytitle="A.U.",
            color=color)
    hc.create_and_book_histo("{0}_h_pLcut15_Resd0_{1}".format(resonance,suffix),\
            "Resolution of leading kaons impact parameter",20,0,30,\
            description=description,
            xtitle="#sigma(d_{0}) [#mum]",
            ytitle="A.U.",
            color=color)
    hc.create_and_book_histo("{0}_h_z0_{1}".format(resonance,suffix),\
            "leading kaons impact parameter",100,-10,10,\
            description=description,
            xtitle="z_{0} [mm]",
            ytitle="A.U.",
            color=color)    
    hc.create_and_book_histo("{0}_h_pL_pions_{1}".format(resonance,suffix),\
            "Parallel momentum distribution of truth pions",40,0,40,\
            description=description,
            xtitle="p_{||} [GeV]",
            ytitle="A.U.",
            color=color)
    hc.create_and_book_histo("{0}_h_pL_kaons_{1}".format(resonance,suffix),\
            "Parallel momentum distribution of truth kaons",40,0,40,\
            description=description,
            xtitle="p_{||} [GeV]",
            ytitle="A.U.",
            color=color)
    hc.create_and_book_histo("{0}_h_pL_candidates_{1}".format(resonance,suffix),\
            "Parallel momentum distribution of kaon candidates",40,0,40,\
            description=description,
            xtitle="p_{||} [GeV]",
            ytitle="A.U.",
            color=color)
    hc.create_and_book_histo("{0}_h_plab_candidates_{1}".format(resonance,suffix),\
            "Lab frame momentum distribution of kaon candidates",40,0,50,\
            description=description,
            xtitle="p_{lab} [GeV]",
            ytitle="A.U.",
            color=color)
    hc.create_and_book_histo("{0}_h_Lxy_{1}".format(resonance,suffix),\
            "leading kaons production vertex (transverse plane)",\
            100,0,5,
            description=description,
            xtitle="L_{xy} [mm]",ytitle='A.U.',
            color=color)
    hc.create_and_book_histo("{0}_h_R_KP_{1}".format(resonance,suffix),\
            "leading kaons production vertex",\
            100,0,5,description=description,
            xtitle="R [mm]", ytitle='A.U.',
            color=color)
    hc.create_and_book_histo("{0}_h_R_Ks_{1}".format(resonance,suffix),\
            "leading K_s decay vertex",\
            100,5,1000,description=description,
            xtitle="R_{K_{s}} [mm]", ytitle='A.U.',
            color=color)
    hc.create_and_book_histo("{0}_h_pLcut10_R_Ks_{1}".format(resonance,suffix),\
            "leading K_s decay vertex",100,10,1000,\
            description=description,
            xtitle="R_{K_{s}} [mm]",
            ytitle="A.U.",
            color=color)
    
    hc.create_and_book_histo("{0}_h_nM_{1}".format(resonance,suffix),\
            "charge particle multiplicity (per quark)",\
            17,0,16,description=description,
            xtitle="N_{part}", ytitle='A.U.',
            color=color)
    
    hc.create_and_book_histo("{0}_h_theta_lab_{1}".format(resonance,suffix),\
            "#theta angle in lab. frame for hadrons with |p| > 20 GeV",
            100,0,91,\
            description=description,
            xtitle="|#theta_{lab}| [^{o}]", ytitle="A.U.",
            color=color)
    
    hc.create_and_book_histo("{0}_h2_pL_theta_lab_{1}".format(resonance,suffix),\
            "leading hadrons #theta vs. p_{L}",\
            100,0,91,npoints_y=100,ylow=0,yhigh=65,description=description,
            ytitle="p_{||} [GeV]",
            xtitle='|#theta_{lab}| [^{o}]',
            color=color)
    
    hc.create_and_book_histo("{0}_h2_cosTheta_{1}".format(resonance,suffix),\
            "leading kaons angle (q#bar{q} system)",\
            100,-1,1,npoints_y=100,ylow=-1,yhigh=1,description=description,
            xtitle="leading-kaon cos(#theta_{q#bar{q}})",
            ytitle='subleading-kaon cos(#theta_{q#bar{q}})',
            color=color)
    
    hc.create_and_book_histo("{0}_h2_pL_cosTheta_{1}".format(resonance,suffix),\
            "leading kaons angle (q#bar{q} system)",\
            100,-1,1,npoints_y=100,ylow=-64,yhigh=64,description=description,
            xtitle="cos(#theta_{q#bar{q}})",
            ytitle='p_{||} [GeV])',
            color=color)
    
    hc.create_and_book_histo("{0}_h2_Resd0_theta_{1}".format(resonance,suffix),\
            "leading kaons: #sigma_{d_{0}} vs. #theta_{lab}",\
            100,0,91,npoints_y=100,ylow=0,yhigh=100.,description=description,
            ytitle="#sigma_{d_{0}} [#mum]",
            xtitle='|#theta_{lab}| [^{o}]',
            color=color)

    hc.create_and_book_histo("{0}_h2_pLcut20_Resd0_theta_{1}".format(resonance,suffix),\
            "leading kaons: #sigma_{d_{0}} vs. #theta_{lab}",\
            100,0,91.,npoints_y=100,ylow=0,yhigh=100.,description=description,
            ytitle="#sigma_{d_{0}} [#mum]",
            xtitle='|#theta_{lab}| [^{o}]',
            color=color)
    
    # NEW: multiplicity
    hc.create_and_book_histo("{0}_h2_pL_multiplicity_0_{1}".format(resonance,suffix),\
            "leading hadrons: p_{||} vs. N_{trk} #in dR < 0.4",
            100,0,65.,npoints_y=20,ylow=-0.5,yhigh=19.5,description=description,
            ytitle="N_{trk} #in dR < 0.4",
            xtitle=' p_{||} [GeV]',
            color=color)
    hc.create_and_book_histo("{0}_h2_pL_multiplicity_1_{1}".format(resonance,suffix),\
            "sub-leading hadrons: p_{||} vs. N_{trk} #in dR < 0.4",
            100,0,65.,npoints_y=20,ylow=-0.5,yhigh=19.5,description=description,
            ytitle="N_{trk} #in dR < 0.4",
            xtitle=' p_{||} [GeV]',
            color=color)
    hc.create_and_book_histo("{0}_h2_pL_multiplicity_Add_{1}".format(resonance,suffix),\
            "leading/subleading hadrons: p_{||} vs. N_{trk} #in dR < 0.4",
            100,0,80.,npoints_y=25,ylow=-0.5,yhigh=24.5,description=description,
            ytitle="N_{trk} #in dR < 0.4",
            xtitle=' p_{||} [GeV]',
            color=color)
    
    # TO BE DEPRECATED -- 
    hc.create_and_book_histo("{0}_h2_d0_{1}".format(resonance,suffix),\
            "leading kaons impact parameter",\
            100,-5,5,npoints_y=100,ylow=-5,yhigh=5,description=description,
            xtitle="leading-kaon d_{0} [mm]",
            ytitle='subleading-kaon d_{0} [mm]',
            color=color)
    
    hc.create_and_book_histo("{0}_h2_z0_{1}".format(resonance,suffix),\
            "leading kaons impact parameter",\
            100,-10,10,npoints_y=100,ylow=-10,yhigh=10,description=description,
            xtitle="leading-kaon z_{0} [mm]",
            ytitle='subleading-kaon z_{0} [mm]',
            color=color)
    
    hc.create_and_book_histo("{0}_h2_Lxy_{1}".format(resonance,suffix),\
            "leading kaons production vertex (transverse plane)",\
            100,0,5,npoints_y=100,ylow=0,yhigh=5,description=description,
            xtitle="leading-kaon L_{xy} [mm]",
            ytitle='subleading-kaon L_{xy} [mm]',
            color=color)
    
    hc.create_and_book_histo("{0}_h2_R_{1}".format(resonance,suffix),\
            "leading kaons production vertex",\
            100,0,5,npoints_y=100,ylow=0,yhigh=5,description=description,
            xtitle="leading-kaon R [mm]",
            ytitle='subleading-kaon R [mm]',
            color=color)
    
    #NBINS = 100
    #D0MAX = 5.0
    ## --- The th3 histograms to be used for efficiency calculations
    #typenames = ['h3_pL_d0_d0_KK','h3_pL_d0_d0_KP','h3_pL_d0_d0_PP']
    #for s in map(lambda x: '{0}_{1}_{2}'.format(resonance,x,suffix),typenames):
    #    hc.create_and_book_histo(s, 'p_{||} circular cut and d_{0}; p_{||} cut [GeV];'\
    #            'd_{0}^{1} [mm];  d_{0}^{2} [mm]',\
    #            NBINS,0,65,
    #            npoints_y=NBINS,ylow=-D0MAX,yhigh=D0MAX,
    #            npoints_z=NBINS,zlow=-D0MAX,zhigh=D0MAX,
    #            xtitle = '#sqrt{p_{||,L}^{2}+p_{||,sL}^{2}}',\
    #            ytitle = 'd_{0}^{lead} [mm]',\
    #            ztitle = 'd_{0}^{sublead} [mm]',
    #            description=description,
    #            color=color)
    return hc

def plot(histo,varname,xtitle='',ytitle='',option=''):
    """
    """
    from PyAnUtils.plotstyles import squaredStyle,setpalette
    import ROOT
    import os
    lstyle = squaredStyle()
    lstyle.cd()
    ROOT.gROOT.ForceStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch()
    #setpalette("forest")
    setpalette("darkbody")
    
    # plotting it
    #histo.GetXaxis().SetTitle(xtitle)
    #histo.GetYaxis().SetTitle(ytitle)
    
    c = ROOT.TCanvas()
    histo.Draw(option)
    #suffix = SUFFIXPLOTS
    c.SaveAs("{0}{1}".format(histo.GetName(),SUFFIXPLOTS))

    c.Close()
    del c

def plot_combined(hc,varname,option='',legposition="RIGHT"):
    """
    """
    from PyAnUtils.plotstyles import squaredStyle
    import ROOT
    import os
    lstyle = squaredStyle()
    lstyle.cd()
    ROOT.gROOT.ForceStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gROOT.SetBatch()
    
    # plotting it
    # --- FIXME:: USE A function
    histonames = filter(lambda _k: _k.find(varname) == 0,\
            hc._histos.keys())
    hc.associated(histonames)
    hc.plot(histonames[0],'comb_{0}{1}'.format(varname,SUFFIXPLOTS),log=True,legposition=legposition)

def plot_profile_combined(hc,varname,axis,
                          ytitle='<#sigma_{d_{0}}> [#mum]',
                          options='',legposition="RIGHT"):
    """Plot in the same canvas the plots with the common name 
    """
    from PyAnUtils.plotstyles import squaredStyle
    import ROOT
    import os
    import random
    lstyle = squaredStyle()
    lstyle.cd()
    ROOT.gROOT.ForceStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gROOT.SetBatch()
    
    
    # plotting it
    # --- FIXME:: USE A function
    histonames = filter(lambda _k: _k.find(varname) == 0 and _k.find('_pf_') == -1,\
            hc._histos.keys())
    # create the profiles and include it in hc
    profile_names = []
    for hname in histonames:
        _h = hc._histos[hname]
        # Note that given the exact name for several cases (with/without) cuts
        # a random number is needed
        _profile_name = "{0}_pf_{1}_{2}".format(hname,axis,hash(random.uniform(0,1e5)))
        if axis == "X":
            fProxy = _h.ProfileX
        elif axis == "Y":
            fProxy = _h.ProfileY
        hc.book_histo(fProxy().Clone(_profile_name),\
                description=hc._description[hname],\
                ytitle=ytitle,
                color=_h.GetLineColor())
        profile_names.append( _profile_name )
        #hc.create_and_book_histo("{0}_prf_{1}".format(hname,axis),\
        #    _h.GetTitle(),_h.GetNbinsX(),_h.GetBinLowEdge(1),_h.GetBinLowEdge(_h.GetNbinsX()),\
        #    description=hc._description[hname],
        #    xtitle=_h.GetXaxis().GetTitle(),
        #    ytitle=_h.GetYaxis().GetTitle())
        #    #color=color)
    hc.associated(profile_names)
    hc.plot(profile_names[0],'comb_{0}_{1}{2}'.format(varname,axis,SUFFIXPLOTS),\
            options=options,log=True,legposition=legposition,\
            normalize=False)

def saveallinfo(d,surname):
    """Store the dictionary "d"

    Parameters
    ----------
    d: any serializable object with pickle
    surname: str, the string to be appended to 'd0cut_dict_' filename
    """
    import cPickle as pickle

    output = open('d0cut_dict_{0}.pkl'.format(surname),'wb')
    pickle.dump(d,output)
    output.close()

def getallinfo(filename):
    """Store the dictionary "d"

    Parameters
    ----------
    d: any serializable object with pickle
    surname: str, the string to be appended to 'd0cut_dict_' filename
    """
    import cPickle as pickle

    input_object = pickle.load( open(filename,'rb') )
    d = input_object
    return d


def main_fixed_pid(rootfile,channels,tables,pLMin,pLMax,pLcut_type,d0cuts,d0cut_type,z0cut,
                   Ntrackcut,wp_activated,R0,charge_combination):
    """Main function steering the efficiency calculation and
    the plots creation.
    
    TO BE FILLED

    * Plot of all the efficiencies vs. momentum cut in the
     same canvas (including the total background efficiency
     which is calculated in the function gettotalbkgeff)
     * Significance plots (eff_S/sqrt(e_totalbkg))
     * ROC curves (signal efficiency vs. total bkg efficiency)
    The significance plots show as well some working points
    where the signal eff. and total bck. efficiency are shown

    Parameters
    ----------
    rootfile: str
    channels: list(str)
        the name of the Higgs' decay channels to consider
    tables: bool
        whether or not print-out the latex tables
    pLMin: float
        the minimal parallel momentum cut considered
    pLMax: float
        the maximal parallel momentum cut considered
    pLcut_type: str
        functional of the pL cut(pL1,pL2): circular|square|line
    d0cuts: list(float)
        the list of d0-cuts to consider, only apply to K+- and Pi+-
    d0cut_type: str
        functional of the d0cut(d01,d02): circular|square
    z0: float|None
        whether activate or not z0cut, if yes, the cut number
    Ntrackcut: int|None
        whether activate or not multiplicity cut, if yes, the cut number
    wp_activated: bool
        whether or not to plot on the working points...
    R0: float
        the cut on the decay length of the K_s.
    charge_combination: str
        the charge combination of the final state to be chosen all|CC|NC|NN
    """
    import ROOT
    import sys
    from math import sqrt,pi

    
    # Get the root file with and the TTree object
    # the root file should contain the tree with the leading hadrons (just two hadrons
    # per hemisphere) ordered by p_parallel
    f = ROOT.TFile(rootfile)
    if f.IsZombie():
        raise IOError("ROOT file '{0}' doesn\'t exist".format(rootfile))
    _preobj = dict(map(lambda x:(x.GetName(), f.Get(x.GetName())),
         filter(lambda x: x.GetClassName().find('TTree') == 0,f.GetListOfKeys())))

    # Split by resonance: FIXME-- lost right now, to be re-incorporated
    #_obj = { 'H' : dict(filter(lambda (x,y): x.find('H') == 0,_preobj.iteritems())),
    #        'Z' : dict(filter(lambda (x,y): x.find('Z') == 0,_preobj.iteritems()))
    #        }
    _obj = { 'H': _preobj }

    #-- get the proper name of the signal channel (with the amount of PID)
    signal_channel   = filter(lambda x: x.find('ssbar') != -1,channels)[0]

    # -- just take care in the H-ressonance...
    pLcuts = xrange(pLMin,pLMax+1)
    
    # --- Ready to extract efficiencies
    # the eff. classes 
    eff = dict(map(lambda x: (x,eff_cut_hadron(x,z0cut,Ntrackcut,d0cut_type=d0cut_type,pLcut_type=pLcut_type,R0=R0,charge_combination=charge_combination)), channels))
    
    message = "\r\033[1;34mhzplots INFO\033[1;m Obtaining efficiencies"
    # { 'docut1': [ (pLcut1, eff_sig, significance, pion_rejection, purity, N_sig, N_bkg), ... ],  }
    observables = {}
    # Indices, see in line [MARK-1]
    I_PL = 0; I_EFF_SIG = 1; I_SIGN=2; I_PION_REJEC = 3; I_PURITY = 4; I_N_SIGNAL = 5; I_N_BKG=6; I_EFF_BKG=7; I_N_BB=8;
    I_N_CC=9; I_N_GG=10; I_N_UU=11;  I_N_DD=12; I_N_SS_CC=13; I_N_SS_NC=14; I_N_SS_NN=15; I_N_BB_CC=16; I_N_BB_NC=17;
    I_N_BB_NN=18; I_N_CC_CC=19; I_N_CC_NC=20; I_N_CC_NN=21; I_N_GG_CC=22; I_N_GG_NC=23; I_N_GG_NN=24;  I_N_UU_CC=25;
    I_N_UU_NC=26; I_N_UU_NN=27; I_N_DD_CC=28; I_N_DD_NC=29; I_N_DD_NN=30; I_EFF_BB=31; I_EFF_CC=32; I_EFF_SS=33; I_EFF_UU=34;
    I_EFF_DD=35; I_EFF_GG=36;
    I_BCDF_BB=37; I_BCDF_CC=38; I_BCDF_SS=39; I_BCDF_UU=40; I_BCDF_DD=41; I_BCDF_GG=42;
    I_BCDF_BB_CC=43; I_BCDF_BB_NC=44; I_BCDF_BB_NN=45; I_BCDF_CC_CC=46; I_BCDF_CC_NC=47; I_BCDF_CC_NN=48;
    I_BCDF_SS_CC=49; I_BCDF_SS_NC=50; I_BCDF_SS_NN=51; I_BCDF_UU_CC=52; I_BCDF_UU_NC=53; I_BCDF_UU_NN=54;
    I_BCDF_DD_CC=55; I_BCDF_DD_NC=56; I_BCDF_DD_NN=57; I_BCDF_GG_CC=58; I_BCDF_GG_NC=59; I_BCDF_GG_NN=60;

    for _d0 in d0cuts:
        d0str = '{0}'.format(_d0)
        observables[d0str] = []
        i = 0
        for pL in pLcuts:
            sys.stdout.write( "{0} {1}".format(message,SPINNING[i % len(SPINNING)]))
            sys.stdout.flush()

            # setting the current cuts to all the efficienciesa
            _dummy = map(lambda e: e.set_total_eff(_obj['H'],pLcut=pL,d0cut=_d0), eff.values())

            # Some needed values
            #eff_sig = eff[signal_channel].get_total_eff('KK') ?
            eff_sig = eff[signal_channel].get_total_eff()
            eff_ss  = eff_sig
            n_KK    = eff[signal_channel].get_total_events()
            n_ss_CC = eff[signal_channel].get_total_events_by_FS('CC')
            n_ss_NC = eff[signal_channel].get_total_events_by_FS('NC')
            n_ss_NN = eff[signal_channel].get_total_events_by_FS('NN')
            bcdf_ss = eff[signal_channel].get_BCdaughter_fraction()
            bcdf_ss_CC = eff[signal_channel].get_BCdaughter_fraction_by_FS('CC')
            bcdf_ss_NC = eff[signal_channel].get_BCdaughter_fraction_by_FS('NC')
            bcdf_ss_NN = eff[signal_channel].get_BCdaughter_fraction_by_FS('NN')
            eff_bb  = eff[filter(lambda x: x.find('bbbar') != -1,channels)[0]].get_total_eff()
            n_bb    = eff[filter(lambda x: x.find('bbbar') != -1,channels)[0]].get_total_events()
            n_bb_CC = eff[filter(lambda x: x.find('bbbar') != -1,channels)[0]].get_total_events_by_FS('CC')
            n_bb_NC = eff[filter(lambda x: x.find('bbbar') != -1,channels)[0]].get_total_events_by_FS('NC')
            n_bb_NN = eff[filter(lambda x: x.find('bbbar') != -1,channels)[0]].get_total_events_by_FS('NN')
            bcdf_bb = eff[filter(lambda x: x.find('bbbar') != -1,channels)[0]].get_BCdaughter_fraction()
            bcdf_bb_CC = eff[filter(lambda x: x.find('bbbar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('CC')
            bcdf_bb_NC = eff[filter(lambda x: x.find('bbbar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('NC')
            bcdf_bb_NN = eff[filter(lambda x: x.find('bbbar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('NN')
            eff_cc  = eff[filter(lambda x: x.find('ccbar') != -1,channels)[0]].get_total_eff()
            n_cc    = eff[filter(lambda x: x.find('ccbar') != -1,channels)[0]].get_total_events()
            n_cc_CC = eff[filter(lambda x: x.find('ccbar') != -1,channels)[0]].get_total_events_by_FS('CC')
            n_cc_NC = eff[filter(lambda x: x.find('ccbar') != -1,channels)[0]].get_total_events_by_FS('NC')
            n_cc_NN = eff[filter(lambda x: x.find('ccbar') != -1,channels)[0]].get_total_events_by_FS('NN')
            bcdf_cc = eff[filter(lambda x: x.find('ccbar') != -1,channels)[0]].get_BCdaughter_fraction()
            bcdf_cc_CC = eff[filter(lambda x: x.find('ccbar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('CC')
            bcdf_cc_NC = eff[filter(lambda x: x.find('ccbar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('NC')
            bcdf_cc_NN = eff[filter(lambda x: x.find('ccbar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('NN')
            eff_gg  = eff[filter(lambda x: x.find('gg') != -1,channels)[0]].get_total_eff()
            n_gg    = eff[filter(lambda x: x.find('gg') != -1,channels)[0]].get_total_events()
            n_gg_CC = eff[filter(lambda x: x.find('gg') != -1,channels)[0]].get_total_events_by_FS('CC')
            n_gg_NC = eff[filter(lambda x: x.find('gg') != -1,channels)[0]].get_total_events_by_FS('NC')
            n_gg_NN = eff[filter(lambda x: x.find('gg') != -1,channels)[0]].get_total_events_by_FS('NN')
            bcdf_gg = eff[filter(lambda x: x.find('gg') != -1,channels)[0]].get_BCdaughter_fraction()
            bcdf_gg_CC = eff[filter(lambda x: x.find('gg') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('CC')
            bcdf_gg_NC = eff[filter(lambda x: x.find('gg') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('NC')
            bcdf_gg_NN = eff[filter(lambda x: x.find('gg') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('NN')
            eff_uu  = eff[filter(lambda x: x.find('uubar') != -1,channels)[0]].get_total_eff()
            n_uu    = eff[filter(lambda x: x.find('uubar') != -1,channels)[0]].get_total_events()
            n_uu_CC = eff[filter(lambda x: x.find('uubar') != -1,channels)[0]].get_total_events_by_FS('CC')
            n_uu_NC = eff[filter(lambda x: x.find('uubar') != -1,channels)[0]].get_total_events_by_FS('NC')
            n_uu_NN = eff[filter(lambda x: x.find('uubar') != -1,channels)[0]].get_total_events_by_FS('NN')
            bcdf_uu = eff[filter(lambda x: x.find('uubar') != -1,channels)[0]].get_BCdaughter_fraction()
            bcdf_uu_CC = eff[filter(lambda x: x.find('uubar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('CC')
            bcdf_uu_NC = eff[filter(lambda x: x.find('uubar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('NC')
            bcdf_uu_NN = eff[filter(lambda x: x.find('uubar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('NN')
            eff_dd  = eff[filter(lambda x: x.find('ddbar') != -1,channels)[0]].get_total_eff()
            n_dd    = eff[filter(lambda x: x.find('ddbar') != -1,channels)[0]].get_total_events()
            n_dd_CC = eff[filter(lambda x: x.find('ddbar') != -1,channels)[0]].get_total_events_by_FS('CC')
            n_dd_NC = eff[filter(lambda x: x.find('ddbar') != -1,channels)[0]].get_total_events_by_FS('NC')
            n_dd_NN = eff[filter(lambda x: x.find('ddbar') != -1,channels)[0]].get_total_events_by_FS('NN')
            bcdf_dd = eff[filter(lambda x: x.find('ddbar') != -1,channels)[0]].get_BCdaughter_fraction()
            bcdf_dd_CC = eff[filter(lambda x: x.find('ddbar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('CC')
            bcdf_dd_NC = eff[filter(lambda x: x.find('ddbar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('NC')
            bcdf_dd_NN = eff[filter(lambda x: x.find('ddbar') != -1,channels)[0]].get_BCdaughter_fraction_by_FS('NN')
            bkg_tot_evts = sum(map(lambda (x,y): y.get_total_events(),\
                    filter(lambda (x,y): x != signal_channel, eff.iteritems() )))
            bkg_tot_eff = sum(map(lambda (x,y): y.get_total_eff(),\
                    filter(lambda (x,y): x != signal_channel, eff.iteritems() )))
            # Signal efficiency calculation 
            # effsig := eff['ssbar_nnnPID'] --< ssbar
            # purity := 2*N_KK/(2*N_KK+2*N_KP+2*N_PP) = N_KK/(N_KK+N_KP+N_PP) --> ssbar 
            # Pion rejection: 1/(# accept. as K / # pions generated) -> equivalent
            #     to ?? --> NOT IMPLEMENTED YET
            # Note: sigma_HZ*L_int*Sum_qq BR(H->qq)/Sum_qq BR(h->qq)*(Sum_qq BR(H->qq) eff_qq*f_qq)
            # the Sum_qq BR(H->qq) terms are cancelled... 
            #pion_rejection = float(bkg_tot_evts)/float(n_KK)
            pion_rejection = -1 # dummy NOT IMPLEMNTED
            try:
                purity = (float(eff[signal_channel].get_events('KK'))+\
                          float(eff[signal_channel].get_events('KsK'))+\
                          float(eff[signal_channel].get_events('KKs'))+\
                          float(eff[signal_channel].get_events('KsKs')))/\
                          (float(eff[signal_channel].get_events('KK'))+\
                           float(eff[signal_channel].get_events('KP'))+\
                           float(eff[signal_channel].get_events('PP'))+\
                           float(eff[signal_channel].get_events('KsKs'))+\
                           float(eff[signal_channel].get_events('KsK'))+\
                           float(eff[signal_channel].get_events('KKs'))+\
                           float(eff[signal_channel].get_events('KsP'))+\
                           float(eff[signal_channel].get_events('PKs')))
            except ZeroDivisionError:
                purity = 0.0
            try:
                significance   = float(n_KK)/sqrt(float(bkg_tot_evts))
            except ZeroDivisionError:
                significance   = 0.0
            # store it, note reference above [MARK-1]
            observables[d0str].append( (pL,eff_sig,significance,pion_rejection,purity,n_KK,bkg_tot_evts,bkg_tot_eff,
                                        n_bb,n_cc,n_gg, n_uu, n_dd, n_ss_CC, n_ss_NC, n_ss_NN, n_bb_CC, n_bb_NC,
                                        n_bb_NN, n_cc_CC, n_cc_NC, n_cc_NN, n_gg_CC, n_gg_NC, n_gg_NN, n_uu_CC,
                                        n_uu_NC, n_uu_NN, n_dd_CC, n_dd_NC, n_dd_NN, eff_bb, eff_cc, eff_ss, eff_uu,
                                        eff_dd, eff_gg, bcdf_bb, bcdf_cc, bcdf_ss, bcdf_uu, bcdf_dd, bcdf_gg,
                                        bcdf_bb_CC, bcdf_bb_NC, bcdf_bb_NN, bcdf_cc_CC, bcdf_cc_NC, bcdf_cc_NN,
                                        bcdf_ss_CC, bcdf_ss_NC, bcdf_ss_NN, bcdf_uu_CC, bcdf_uu_NC, bcdf_uu_NN,
                                        bcdf_dd_CC, bcdf_dd_NC, bcdf_dd_NN, bcdf_gg_CC, bcdf_gg_NC, bcdf_gg_NN ) )
            i+=1
    # plotting
    print
    print "\033[1;34mhzplots INFO\033[1;m Plotting..."
    # --- filling uncutted observables
    plotconstraint = "1 == 1"
    if charge_combination == 'CC':
        plotconstraint = 'isKshort[0] + isKshort[1] == 0'
    elif charge_combination == 'NN':
        plotconstraint = 'isKshort[0] + isKshort[1] == 2'
    elif charge_combination == 'NC':
        plotconstraint = 'isKshort[0] + isKshort[1] == 1'
    hc = None
    for effname,e in eff.iteritems():
        hc = create_histos(e.decay_channel,e.short_name,25,hc)
        #e.get_tree().Project("H_h_d0_{0}".format(e.decay_channel),"(vy-vx*tan(phi_lab))*cos(phi_lab)")
        #e.get_tree().Project("H_h_z0_{0}".format(e.decay_channel),"-(vy-vz*tan(theta_lab))/tan(theta_lab)")
        e.get_tree().Project("H_h_d0_{0}".format(e.decay_channel),"d0","isKshort!=1")
        e.get_tree().Project("H_h_absd0_{0}".format(e.decay_channel),"abs(d0)","isKshort!=1")
        e.get_tree().Project("H_h_Resd0_{0}".format(e.decay_channel),"5.+(10/(p_lab*sin(theta_lab)**(3./2.)))","isKshort!=1")
        e.get_tree().Project("H_h_z0_{0}".format(e.decay_channel),"z0","isKshort!=1")
        e.get_tree().Project("H_h_pL_pions_{0}".format(e.decay_channel),"abs(p*cos(theta))","pdgId==211")
        e.get_tree().Project("H_h_pL_kaons_{0}".format(e.decay_channel),"abs(p*cos(theta))","pdgId==321")
        e.get_tree().Project("H_h_pL_candidates_{0}".format(e.decay_channel),"abs(p*cos(theta))")
        e.get_tree().Project("H_h_plab_candidates_{0}".format(e.decay_channel),"abs(p_lab)")
        e.get_tree().Project("H_h_Lxy_{0}".format(e.decay_channel),"sqrt(vx*vx+vy*vy)","isKshort!=1")
        e.get_tree().Project("H_h_R_KP_{0}".format(e.decay_channel),"R",'isKshort==0')
        e.get_tree().Project("H_h_R_Ks_{0}".format(e.decay_channel),"R",'isKshort==1')
        # two-dim
        e.get_tree().Project("H_h2_pL_{0}".format(e.decay_channel),"abs(p[1]*cos(theta[1])):abs(p[0]*cos(theta[0]))","{0}".format(plotconstraint))
        e.get_tree().Project("H_h2_pL_d0_{0}".format(e.decay_channel),"d0:abs(p*cos(theta))","isKshort!=1 && {0}".format(plotconstraint))
        e.get_tree().Project("H_h2_pL_R_{0}".format(e.decay_channel),"R:abs(p*cos(theta))","isKshort==1 && {0}".format(plotconstraint))
        e.get_tree().Project("H_h2_Resd0_theta_{0}".format(e.decay_channel),\
                "5.+(10/(p_lab*sin(theta_lab)**(3./2.))):acos(abs(cos(theta_lab)))*180./{0}".format(pi),"{0}".format(plotconstraint))
        e.get_tree().Project("H_h2_pL_theta_lab_{0}".format(e.decay_channel),\
                "abs(p*cos(theta)):acos(abs(cos(theta_lab)))*180./{0}".format(pi),"{0}".format(plotconstraint))
        # -- The new-multiplicity
        e.get_tree().Project("H_h2_pL_multiplicity_0_{0}".format(e.decay_channel),"multiplicity[0]:abs(p[0]*cos(theta[0]))","{0}".format(plotconstraint))
        e.get_tree().Project("H_h2_pL_multiplicity_1_{0}".format(e.decay_channel),"multiplicity[1]:abs(p[1]*cos(theta[1]))","{0}".format(plotconstraint)) 
        e.get_tree().Project("H_h2_pL_multiplicity_Add_{0}".format(e.decay_channel),\
                "multiplicity[0]+multiplicity[1]:sqrt( (p[0]*cos(theta[0]))**2.0+ (p[1]*cos(theta[1]))**2.0 )","{0}".format(plotconstraint))
        # cut-dependent
        e.activate_cuts(pLcut=20)
        e.get_tree().Project("H_h_theta_lab_{0}".format(e.decay_channel),\
                "acos(abs(cos(theta_lab)))*180.0/{0}".format(pi),"{0}".format(plotconstraint))
        e.get_tree().Project("H_h2_pLcut20_Resd0_theta_{0}".format(e.decay_channel),\
                "5.+(10/(p_lab*sin(theta_lab)**(3./2.))):acos(abs(cos(theta_lab)))*180./{0}".format(pi),"{0}".format(plotconstraint))
        e.deactivate_cuts()
        e.activate_cuts(pLcut=10)
        e.get_tree().Project("H_h_pLcut10_R_Ks_{0}".format(e.decay_channel),"R","isKshort==1")
        e.deactivate_cuts()
        e.activate_cuts(pLcut=15)        
        e.get_tree().Project("H_h_pLcut15_absd0_{0}".format(e.decay_channel),"abs(d0)","isKshort!=1")
        e.get_tree().Project("H_h_pLcut15_Resd0_{0}".format(e.decay_channel),"5.+(10/(p_lab*sin(theta_lab)**(3./2.)))","isKshort!=1")
        e.deactivate_cuts()
    # -- plotting ...
    for k,h in filter(lambda (_k,_h): _k.find('H_h2_pL')==0 and \
            _k.find("cosTheta") == -1,hc._histos.iteritems()):
        plot(h,k,option='COLZ')
    for k,h in filter(lambda (_k,_h): _k.find('H_h2_pL_theta_lab')==0,hc._histos.iteritems()):
        plot(h,k,option='COLZ')
    #for k,h in filter(lambda (_k,_h): _k.find('cosTheta')!=-1,hc._histos.iteritems()):
    #    plot(h,k,option='COLZ')
   
    # and the combined histograms
    plot_combined(hc,'H_h_d0')
    plot_combined(hc,'H_h_absd0')
    plot_combined(hc,'H_h_pLcut15_absd0')
    plot_combined(hc,'H_h_Resd0')
    plot_combined(hc,'H_h_pLcut15_Resd0')
    plot_combined(hc,'H_h_pLcut10_R_Ks')
    plot_combined(hc,'H_h_z0')
    plot_combined(hc,'H_h_pL_pions')
    plot_combined(hc,'H_h_pL_kaons')
    plot_combined(hc,'H_h_pL_candidates')
    plot_combined(hc,'H_h_plab_candidates')
    plot_combined(hc,'H_h_Lxy')
    plot_combined(hc,'H_h_R_KP')
    plot_combined(hc,'H_h_R_Ks')
    try:
        plot_combined(hc,'H_h_theta_lab',legposition="LEFT")
    except ZeroDivisionError:
        pass
    plot_profile_combined(hc,"H_h2_Resd0_theta","X",options="PE")
    plot_profile_combined(hc,"H_h2_pLcut20_Resd0_theta","X",options="PE")
    plot_profile_combined(hc,"H_h2_pL_theta_lab","X",ytitle="<p_{||}> [GeV]",options="PE")
    plot_profile_combined(hc,"H_h2_pL_theta_lab","Y",ytitle="<#theta> [^{0{}]",options="PE")
    # -- new multiplicity
    plot_profile_combined(hc,"H_h2_pL_multiplicity_Add","X",ytitle="<N_{trk}> #in dR < 0.4",options="PE")
    plot_profile_combined(hc,"H_h2_pL_multiplicity_Add","Y",ytitle="<p_{||}> [GeV]",options="PE")
    #plot_combined(hc,'H_h_nM')

    # --- Some extra points (WP)
    if wp_activated:
        #leg_format_pion_rej = ( 'S/#sqrt{B}=%.1f @ p_{||}^{c}=%.0f GeV',(I_SIGN,I_PL) )
        #graphs_leg_pion_rej = get_point_graphs(observables,(I_EFF_SIG,I_PION_REJEC),I_PL,[10,20],leg_format_pion_rej)
        graphs_leg_pion_rej = None
    
        leg_format_pur = ( 'S/#sqrt{B}=%.1f @ p_{||}^{c}=%.0f GeV',(I_SIGN,I_PL) )
        graphs_leg_pur = get_point_graphs(observables,(I_EFF_SIG,I_PURITY),I_PL,[10,20],leg_format_pur)
        
        leg_format_sig = ( '#varepsilon_{signal}=%.2f, #pi-rej.factor=%.0f',(I_EFF_SIG,I_PION_REJEC) )
        graphs_leg_sig = get_point_graphs(observables,(I_PL,I_SIGN),I_EFF_SIG,[0.4,0.8],leg_format_sig)
        
        leg_format_roc = ( 'S/#sqrt{B}=%.2f, p_{||}^{cut}=%.1f',(I_SIGN,I_PL) )
        graphs_leg_roc = get_point_graphs(observables,(I_EFF_BKG,I_EFF_SIG),I_PL,[10,20],leg_format_roc)
    else:
        graphs_leg_pion_rej = None
        graphs_leg_pur      = None
        graphs_leg_sig      = None
        graphs_leg_roc      = None

    #pr_attr = plot_attributes('pion_rejection',
    #        xtitle='#varepsilon_{S}', 
    #        ytitle='pion rejection (N_{BKG}/N_{s#bar{s}} with no PID)',
    #        x0 = 0.0, x1 = 1.0 , y0 = 0.0 )
    #make_plot(observables,(I_EFF_SIG,I_PION_REJEC),pr_attr,g_points_dict=graphs_leg_pion_rej)

    purity_attr = plot_attributes('purity',
            xtitle='#varepsilon_{S}', 
            ytitle='purity',
            x0 = 0.0, x1 = 1.0 , y0 = 0.0, y1= 0.0 )
    make_plot(observables,(I_EFF_SIG,I_PURITY),purity_attr,g_points_dict=graphs_leg_pur)

    significance_attr = plot_attributes('significance',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='N_{S}/#sqrt{N_{B}}',
            x0 = pLMin, y0 = 0.0 )
    make_plot(observables,(I_PL,I_SIGN),significance_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Eff_BB_attr = plot_attributes('Eff_bb',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='efficiency',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_EFF_BB),Eff_BB_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Eff_CC_attr = plot_attributes('Eff_cc',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='efficiency',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_EFF_CC),Eff_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Eff_SS_attr = plot_attributes('Eff_ss',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='efficiency',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_EFF_SS),Eff_SS_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Eff_UU_attr = plot_attributes('Eff_uu',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='efficiency',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_EFF_UU),Eff_UU_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Eff_DD_attr = plot_attributes('Eff_dd',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='efficiency',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_EFF_DD),Eff_DD_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Eff_GG_attr = plot_attributes('Eff_gg',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='efficiency',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_EFF_GG),Eff_GG_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_SS_attr = plot_attributes('N_ss_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_SIGNAL),Nevents_SS_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_SS_CC_attr = plot_attributes('N_ss_cc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_SS_CC),Nevents_SS_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_SS_NC_attr = plot_attributes('N_ss_nc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_SS_NC),Nevents_SS_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_SS_NN_attr = plot_attributes('N_ss_nn_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_SS_NN),Nevents_SS_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_BB_attr = plot_attributes('N_bb_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_BB),Nevents_BB_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_BB_CC_attr = plot_attributes('N_bb_cc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_BB_CC),Nevents_BB_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_BB_NC_attr = plot_attributes('N_bb_nc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_BB_NC),Nevents_BB_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_BB_NN_attr = plot_attributes('N_bb_nn_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_BB_NN),Nevents_BB_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_CC_attr = plot_attributes('N_cc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_CC),Nevents_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_CC_CC_attr = plot_attributes('N_cc_cc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_CC_CC),Nevents_CC_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_CC_NC_attr = plot_attributes('N_cc_nc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_CC_NC),Nevents_CC_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_CC_NN_attr = plot_attributes('N_cc_nn_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_CC_NN),Nevents_CC_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_GG_attr = plot_attributes('N_gg_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_GG),Nevents_GG_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_GG_CC_attr = plot_attributes('N_gg_cc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_GG_CC),Nevents_GG_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_GG_NC_attr = plot_attributes('N_gg_nc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_GG_NC),Nevents_GG_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_GG_NN_attr = plot_attributes('N_gg_nn_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_GG_NN),Nevents_GG_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_UU_attr = plot_attributes('N_uu_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_UU),Nevents_UU_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_UU_CC_attr = plot_attributes('N_uu_cc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_UU_CC),Nevents_UU_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_UU_NC_attr = plot_attributes('N_uu_nc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_UU_NC),Nevents_UU_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_UU_NN_attr = plot_attributes('N_uu_nn_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_UU_NN),Nevents_UU_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_DD_attr = plot_attributes('N_dd_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_DD),Nevents_DD_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_DD_CC_attr = plot_attributes('N_dd_cc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_DD_CC),Nevents_DD_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_DD_NC_attr = plot_attributes('N_dd_nc_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_DD_NC),Nevents_DD_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    Nevents_DD_NN_attr = plot_attributes('N_dd_nn_events',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# events',
            x0 = pLMin, y0 = 1.,
            log=True)
    make_plot(observables,(I_PL,I_N_DD_NN),Nevents_DD_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    BCdf_BB_attr = plot_attributes('BCDF_BB',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_BB),BCdf_BB_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_BB_CC_attr = plot_attributes('BCDF_BB_CC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_BB_CC),BCdf_BB_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_BB_NC_attr = plot_attributes('BCDF_BB_NC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_BB_NC),BCdf_BB_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_BB_NN_attr = plot_attributes('BCDF_BB_NN',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_BB_NN),BCdf_BB_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_CC_attr = plot_attributes('BCDF_CC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_CC),BCdf_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_CC_CC_attr = plot_attributes('BCDF_CC_CC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_CC_CC),BCdf_CC_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_CC_NC_attr = plot_attributes('BCDF_CC_NC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_CC_NC),BCdf_CC_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_CC_NN_attr = plot_attributes('BCDF_CC_NN',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_CC_NN),BCdf_CC_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_SS_attr = plot_attributes('BCDF_SS',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_SS),BCdf_SS_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_SS_CC_attr = plot_attributes('BCDF_SS_CC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_SS_CC),BCdf_SS_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_SS_NC_attr = plot_attributes('BCDF_SS_NC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_SS_NC),BCdf_SS_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_SS_NN_attr = plot_attributes('BCDF_SS_NN',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_SS_NN),BCdf_SS_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_UU_attr = plot_attributes('BCDF_UU',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_UU),BCdf_UU_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_UU_CC_attr = plot_attributes('BCDF_UU_CC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_UU_CC),BCdf_UU_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_UU_NC_attr = plot_attributes('BCDF_UU_NC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_UU_NC),BCdf_UU_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_UU_NN_attr = plot_attributes('BCDF_UU_NN',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_UU_NN),BCdf_UU_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_DD_attr = plot_attributes('BCDF_DD',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_DD),BCdf_DD_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_DD_CC_attr = plot_attributes('BCDF_DD_CC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_DD_CC),BCdf_DD_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_DD_NC_attr = plot_attributes('BCDF_DD_NC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_DD_NC),BCdf_DD_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_DD_NN_attr = plot_attributes('BCDF_DD_NN',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_DD_NN),BCdf_DD_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_GG_attr = plot_attributes('BCDF_GG',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_GG),BCdf_GG_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_GG_CC_attr = plot_attributes('BCDF_GG_CC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_GG_CC),BCdf_GG_CC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_GG_NC_attr = plot_attributes('BCDF_GG_NC',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_GG_NC),BCdf_GG_NC_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    BCdf_GG_NN_attr = plot_attributes('BCDF_GG_NN',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='# fraction of BC daughters',
            x0 = pLMin, y0 = 0.,
            log=False)
    make_plot(observables,(I_PL,I_BCDF_GG_NN),BCdf_GG_NN_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")
        
    
    # ROC curve
    roc_attr = plot_attributes("roc",ytitle='#varepsilon_{S}',xtitle='#varepsilon_{B}',
            x0=0.0,x1=1.0,y0=0.0,y1=1.0)
    make_plot(observables,(I_EFF_BKG,I_EFF_SIG),roc_attr,g_points_dict=graphs_leg_roc)

    # Tables
    if tables:
        for _d0,obsList in observables.iteritems():
            print "Table for d0: {0} mm".format(_d0)
            print "----------------------"
            print get_latex_table(obsList)
            print 
        print "Max significance: ",max(obsList,key=lambda x: x[2])

    # d0cut: (pLcut1, eff_sig, significance, pion_rejection, purity, N_sig, N_bkg)
    observables['HEADER'] = ('D0-CUT','MOMENTUM-CUT','EFF_SIGNAL', 'SIGNIFICANCE', \
            'NOT-USED1', 'NOT-USED2','N_SIGNAL','N_BKG','EFF_BKG')
    saveallinfo(observables,signal_channel.replace("ssbar_",""))


def indices_max_hadron(indices_list,pV,thetaV):
    """Obtain the higher parallel momentum on the two hemispheres 
    Note that different hemispheres are described by different
    p*cos(theta) sign

    Return
    ------
    max_up,max_down: int,int
        indices of the two opposite hemisphere high parallel momentum
        hadrons
    """
    from math import cos
    # A dictionary
    d = dict(map(lambda k: (k,pV[k]*cos(thetaV[k])),indices_list))
    # Separate dict for up-hemisphere and down
    dup = filter(lambda (k,parp): parp > 1.0, d.iteritems())
    ddown = filter(lambda (k,parp): parp < 1.0, d.iteritems())
    indices = []
    for _d in [ dup, ddown ]:
        try: 
            indices.append(sorted(_d,reverse=True,key=lambda (k,pc): abs(pc))[0][0])
        except IndexError:
            pass
    return indices

# =============================================================================================
# print_decaychain mode functions
def get_decays(t,hadron,pcut=None,d0cut=None):
    """Obtain the decay histogram for the isBCancestor particles
    
    Parameters
    ----------
    t: ROOT.TTree
    hadron: str,
        the family of ancestors to look at. Valid values are
        B or D
    pcut: float, [Default: None]
        the minimum parallel momentum than the hadrons should 
        carry on
    d0cut: float, [Default: None]
        the maximum (absolute value) impact parameter allowed 
        for the hadron (d0 calculated from the vertex to the 
        point of closest approach using a straight line)

    Return
    ------
    h: list(tuple(int,str))
        the list of decays ordered by times produced, so each entry 
        of the list contains the number the decay mode was produced and the decay (str)
    """
    import sys
    import ROOT
    ROOT.gROOT.SetBatch(True)

    if hadron not in [ "B", "D" ]:
        raise AttributeError("Not a valid hadron to seek: '{0}'"\
                " Accepted 'B' 'D' only".format(hadron))
    # The needed branchesa
    pdgIdV = ROOT.std.vector("int")()
    t.SetBranchAddress("pdgId", pdgIdV)
    pV = ROOT.std.vector("float")()
    t.SetBranchAddress("p", pV)
    if d0cut:
        vtxdict = dict(map(lambda x: (x,ROOT.std.vector("float")()), [ 'vx', 'vy']))
        for _vtx, _vect in vtxdict.iteritems():
            t.SetBranchAddress(_vtx,_vect)
        phi_labV= ROOT.std.vector("float")()
        t.SetBranchAddress("phi_lab",phi_labV)
    thetaV = ROOT.std.vector("float")()
    t.SetBranchAddress("theta",thetaV)
    catchallV = ROOT.std.vector("int")()
    t.SetBranchAddress("catchall", catchallV)
    isBCancestorV = ROOT.std.vector("int")()
    t.SetBranchAddress("isBCancestor", isBCancestorV)
    isBCdaughterV = ROOT.std.vector("int")()
    t.SetBranchAddress("isBCdaughter", isBCdaughterV)
    ancestorBCindexV = ROOT.std.vector("int")()
    t.SetBranchAddress("ancestorBCindex", ancestorBCindexV)
    decay_chainV = ROOT.std.vector('std::string')()
    t.SetBranchAddress("decay_chain", decay_chainV)
    
    # order the final state hadrons need for the pre_cut: 
    indices_fs_hadrons = lambda : map(lambda x: x[0], filter(lambda (i,(catch,isBC)): catch == 25 and isBC != 1,\
            enumerate(zip(catchallV,isBCancestorV))))
    # FIXME: just including mesons right now
    if hadron == "B":
        cut = lambda k: abs(pdgIdV[k]) > 500 and abs(pdgIdV[k]) < 600 
    elif hadron == "D":
        cut = lambda k: abs(pdgIdV[k]) > 400 and abs(pdgIdV[k]) < 500 
    
    decays = {}
    total_fs = 0
    # Filling the histo and getting it
    pointpb = float(t.GetEntries())/100.0
    for i in xrange(t.GetEntries()):
        # Progress bar
        sys.stdout.write("\r\033[1;34m+-- \033[1;mExtracting decay modes for "+hadron+
                " ancestors [ " +"\b"+str(int(float(i)/pointpb)+1).rjust(3)+"%]")
        sys.stdout.flush()
        __dum = t.GetEntry(i)
        # Getting the indices of the final state hadrons (only from H)
        idx_fs = indices_fs_hadrons()
        # Only use the two leading hadrons 
        idx_leading_hadrons = indices_max_hadron(idx_fs,pV,thetaV)
        # passing the d0 cuts and p-cut if any
        if pcut:
            if len(filter(lambda _k: pV[_k] > pcut, idx_leading_hadrons)) != 2:
                continue
        if d0cut:
            if len(filter(lambda _k: 
                abs(d0_explicit_calculation(vtxdict['vx'][_k],vtxdict['vy'][_k],
                    phi_labV[_k])) < d0cut , 
                        idx_leading_hadrons)) != 2:
                continue
        # We are here, so we have two final state hadrons
        total_fs += 2
        # checking if those hadrons are BCdaughters
        idx_bcdaughters = filter(lambda _k: isBCdaughterV[_k],idx_leading_hadrons)
        # storing 
        for k in idx_bcdaughters:
            # Get the index of the BC parent
            BC_k = ancestorBCindexV[k]
            # check if is a B or D ancestor
            if not cut(BC_k):
                continue
            try:
                decays[decay_chainV[BC_k]] += 1
            except KeyError:
                decays[decay_chainV[BC_k]] = 1
    print
    t.ResetBranchAddresses()
    # build 
    return sorted(decays.iteritems(),key=lambda (x,y): y,reverse=True),total_fs

def show_decays(pre_d,nfirst,hadron,total_fs):
    """Print a table with the nfirst frequent produced dacay of the 
    first ancestor from the final hadron (pion/kaon)

    Parameters
    ----------
    pre_d: list(tuple(str,int))
        the input list to extract the info (ordered by number of times happened)
    nfirts: int
        the maximum number of decays to show (taken from the nfirst more
        frequent)
    hadron: str,
        the name of the ancestor (valid only B or D)
    total_fs: int
        the total number of final state hadrons present (and pass the cuts if any)
    """
    print "="*80
    print "{0}-meson decays (first backward ancestor from the two leading final "\
            "state hadron)".format(hadron)
    print "-"*80
    # first check if there is any element
    if len(pre_d) == 0:
        print "NOT FOUND"
        print "="*80
        return
    ntotal = sum(map(lambda (dmode,_n): _n, pre_d))
    print "TOTAL final state hadrons: {0} || with {1} ancestors: {2} ({3:.1f}%)".format(
            total_fs,hadron,ntotal,float(ntotal)/float(total_fs)*100.0)
    print "-"*80
    # convert to frequency (and change the order, frequency first)
    d = map(lambda (_dc,_n): (float(_n)/float(ntotal),_dc), pre_d)

    maxline=max(map(lambda (x,y): len(y),d[:nfirst]))
    fmtstr = "{0}0:{1}{2} {3}".format("{",maxline,"s}","{1:.2f}%")
    for i in d[:nfirst]: 
        print fmtstr.format(i[1],i[0]*100.0)
    print "="*80

def table_latex(outfile,hl,hadron):
    """Create a file containing a latex table with the list of 
    decay and their frequency

    Parameters
    ----------
    outfile: str
        name of the output latex file
    hl: list(tuple(float,str))
        the input list to extract the info
    hadron: str,
        the name of the ancestor (valid only B or D)

    Return
    ------
    the name of the created file
    """
    print "table_latex NOT IMPLEMENTED YET... ignoring"
    
def main_decay_chain(rootfiles,want_latex,nfirst=10,**kwargs):
    """Create a list of most frequent decays for the B/D mesons
    ancestors of the final state hadron 

    Parameters
    ----------
    rootfiles: str,
        the root filenames where to extract the info
    want_latex: bool [NOT IMPLEMENTED YET]
        whether or not a latex table will be printed
    nfirst: int, default: 10
        the number of decays mode to print
    pcut: float, optional 
        a cut in the parallel momentum of the hadrons
    d0cut: float, optional
        a cut in the impact parameters of the hadrons
    """
    import ROOT
    # Only recognized: pcut, d0cut
    if kwargs.has_key("pcut"):
        pcut = kwargs["pcut"]
    if kwargs.has_key("d0cut"):
        d0cut = kwargs["d0cut"]

    needed_branches = [ "pdgId", "p", "catchall", "isBCancestor", \
            "isBCdaughter", "ancestorBCindex", "decay_chain", "theta" ]
    if d0cut:
        needed_branches += [ "vx", "vy", "phi_lab" ]
    for i in rootfiles:
        # Get files and tree
        _froot = ROOT.TFile(i)
        if _froot.IsZombie():
            raise IOError("Cannot open '{0}' root file".format(i))
        t = _froot.Get("mctrue")
        # Speeding up the access: only activated the needed branches
        t.SetBranchStatus("*",0)
        _kk = map(lambda x: t.SetBranchStatus(x,1), needed_branches)
        # -- Get histo for decay for B
        hB,total_fs_B_process = get_decays(t,'B',pcut=pcut,d0cut=d0cut)
        hD,total_fs_C_process = get_decays(t,'D',pcut=pcut,d0cut=d0cut)
        # just a cross-check: the total number of final state hadrons
        # should be the same regardless the B/C ancestor check
        assert total_fs_B_process == total_fs_C_process
        total_fs = total_fs_B_process
        # -- and print the list of the top 10
        show_decays(hB,nfirst,"B",total_fs)
        show_decays(hD,nfirst,"D",total_fs)
        if want_latex:
            # Get a kind of file name from the root filename
            tablename = "table_decay_{0}_"+i.replace(".root",".tex")
            table_latex(tablename.format("B"),hB,"B")
            table_latex(tablename.format("D"),hD,"D")

# =============================================================================================
# compare_pid mode functions
COLORS_PLT = [ 'black','firebrick','seagreen', 'goldenrod', 'blue', 'orange', 'skyblue',
               'black','firebrick','seagreen', 'goldenrod', 'blue', 'orange', 'skyblue']
LINESTYLES = [ '-', '-', '-', '-', '-', '-', '-',
               '--', '--', '--', '--', '--', '--', '--']

def extractlegend(pid):
    splitpid = pid.split('-')
    return '($\epsilon_{{track}},\epsilon_{{K^\pm}},\epsilon_{{\pi^\pm}},\epsilon_{{K^0_s}}$)'\
        '=({0}, {1}, {2}, {3})'.format(splitpid[0], splitpid[1], splitpid[2], splitpid[3])

def ordering(pid):
    splitpid = pid.split('-')
    return float(splitpid[1])

def plot_python(_x,ydict,plotname, x_label, y_label):
    """
    """
    from matplotlib import pyplot as plt
    from scipy.interpolate import spline
    import numpy as np

    # Convert to numpy arrays
    x = np.array(_x)
    # and smooth the final lines
    xnew = np.linspace(x.min(),x.max(),5)
    
    # the figure
    #plt.rc('text', usetex=True)
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)
    ymax = 0.0
    ymin = 0.0
    for k,(pid,signlist) in enumerate(sorted(ydict.iteritems(),key=lambda (x,y): ordering(x))):
        pidname = extractlegend(pid)
        ymax = max(ymax,max(signlist[:20]))
        ymin = min(ymin,min(signlist))
        # Just to smooth a little the output lines
        significance_smooth = spline(x,np.array(signlist),xnew)
        plt.plot(xnew,significance_smooth,
                linewidth=3,linestyle=LINESTYLES[k],color=COLORS_PLT[k], 
                label=pidname)
    ax.set_xlim(x.min(),x.max())
    plt.xlabel(r'{0}'.format(x_label))
    plt.ylabel(r'{0}'.format(y_label))
    ax.set_ylim(ymin,ymax*1.3)
    ax.legend(loc=0,frameon=False)
    plt.savefig(plotname)

COLORS_PLT_DET = [ 'black','indianred', 'goldenrod','darksage', 'royalblue', 'cadetblue',
                   'black', 'black', 'black',
                   'indianred', 'indianred', 'indianred',
                   'goldenrod', 'goldenrod', 'goldenrod',
                   'darksage', 'darksage', 'darksage',
                   'royalblue', 'royalblue', 'royalblue',
                   'cadetblue', 'cadetblue', 'cadetblue',
                  ]
                   
LINESTYLES_DET = ['-', '-', '-', '-', '-', '-',
                  '--', '-.', ':',
                  '--', '-.', ':',
                  '--', '-.', ':',
                  '--', '-.', ':',
                  '--', '-.', ':',
                  '--', '-.', ':' ]
LEGEND_DET     = { 'gg': 'gg', 'bb': 'bbbar', 'cc': 'ccbar',
                   'ss': 'ssbar', 'uu': 'uubar', 'dd': 'ddbar',
                   'ss_CC': 'ss_CC', 'ss_NC': 'ss_NC', 'ss_NN': 'ss_NN',
                   'bb_CC': 'bb_CC', 'bb_NC': 'bb_NC', 'bb_NN': 'bb_NN',
                   'cc_CC': 'cc_CC', 'cc_NC': 'cc_NC', 'cc_NN': 'cc_NN',
                   'gg_CC': 'gg_CC', 'gg_NC': 'gg_NC', 'gg_NN': 'gg_NN',
                   'uu_CC': 'uu_CC', 'uu_NC': 'uu_NC', 'uu_NN': 'uu_NN',
                   'dd_CC': 'dd_CC', 'dd_NC': 'dd_NC', 'dd_NN': 'dd_NN'}
ORDER_DET = { 'ss': 0, 'gg':1, 'bb':2, 'cc':3 , 'uu':4, 'dd':5,
              'ss_CC':6, 'ss_NC':7, 'ss_NN':8,
              'bb_CC':9, 'bb_NC':10, 'bb_NN':11,
              'cc_CC':12, 'cc_NC':13, 'cc_NN':14,
              'gg_CC':15, 'gg_NC':16, 'gg_NN':17,
              'uu_CC':18, 'uu_NC':19, 'uu_NN':20,
              'dd_CC':21, 'dd_NC':22, 'dd_NN':23,}


def plot_python_detailed(_x,ydict,plotname,ylabel,ylog=True):
    """
    """
    from matplotlib import pyplot as plt
    from scipy.interpolate import spline
    import numpy as np

    # Convert to numpy arrays
    x = np.array(_x)
    # and smooth the final lines
    xnew = np.linspace(x.min(),x.max(),300)

    # the figure
    #plt.rc('text', usetex=True)
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)
    ymax = 0.0
    ymin = 0.0

    alreadyplotted=[]
    for k,(pid,signlist) in enumerate(sorted(ydict.iteritems(),key=lambda (x,y): ORDER_DET[x])):
        pidname = LEGEND_DET[pid]
        ymax = max(ymax,max(signlist))
        ymin = min(ymin,min(signlist))
        if ylog:
            if 'eff' in plotname:
                ymin = max(ymin, 0.00001)
            elif 'daughter' in plotname:
                ymin = max(ymin, 0.01)
            else:
                ymin = max(ymin, 0.1)
            if ymax < 10*ymin:
                ymax = 10*ymin
        # Just to smooth a little the output lines
        significance_smooth = spline(x,np.array(signlist),xnew)
        plotcolors = COLORS_PLT_DET[k]
        plotlines  = LINESTYLES_DET[k]
        if 'ss' in pidname:
            plotcolors = 'black'
        elif 'bb' in pidname:
            plotcolors = 'indianred'
        elif 'gg' in pidname:
            plotcolors = 'goldenrod'
        elif 'cc' in pidname:
            plotcolors = 'darksage'
        elif 'uu' in pidname:
            plotcolors = 'royalblue'
        elif 'dd' in pidname:
            plotcolors = 'cadetblue'

        if '_CC' in pidname:
            plotlines = '--'
        elif '_NC' in pidname:
            plotlines = '-.'
        elif '_NN' in pidname:
            plotlines = ':'
        # add only to plot when there is something to plot, i.e. not only 0 or 1
        truefalse=map(lambda i: abs(i)<1e-7 or abs(i-1)<1e-7, signlist)
        if (signlist not in alreadyplotted) and (False in truefalse):
            alreadyplotted.append(signlist)
            plt.plot(xnew,significance_smooth,
                     linewidth=3,
                     linestyle=plotlines,
                     color=plotcolors, 
                     label=pidname)
    ax.set_xlim(x[0],x[-1])
    plt.xlabel(r'Parallel momentum cut [GeV]')
    plt.ylabel(r'{0}'.format(ylabel))
    ax.set_ylim(ymin,ymax*1.3)
    if ylog:
        ax.set_yscale('log')
    #ax.legend(loc=1,frameon=False)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc=2, bbox_to_anchor=(1.05,1))
    #ax.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
    if alreadyplotted != []:
        plt.savefig(plotname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close('all')

def save_divide(numerator,denominator,default):
    result=1;
    try:
        result=1.*numerator/denominator
    except ZeroDivisionError:
        result=1
    return result

def main_cmp_pid(listpklfiles, verbose, pLcut):
    """Steering function to plot equ....
    """
    import os
    
    pid_dict = {}
    for pklfile in listpklfiles:
        # get absolute path and basename
        absname  = os.path.abspath(pklfile)
        basename = os.path.basename(absname)
        # Careful, assuming standard names d0cut_dics_PIDRELATED.pkl
        junk1,junk2,pidname = basename.replace(".pkl","").split("_")
        pid_dict[pidname] = getallinfo(absname)
    # Choose a d0cut to be plotted  XXX
    # FIXME use numpy arrays
    x = map(lambda x: x[0],pid_dict.values()[0].values()[0])
    d0values_tmp = filter(lambda x: x != 'HEADER',pid_dict.values()[0].keys())
    d0values = map(lambda x: float(x), d0values_tmp)
    d0values.sort()
    y = {}
    n = {}
    r = {}
    eff = {}
    sigd0 = {}
    total_bcdf = {}
    bcdf = {}
    # make empty lists for significances of each pid 
    for pid,d0dict in pid_dict.iteritems():
        sigd0[pid]=[]
        total_bcdf[pid]=[]

    for d0cut in map(lambda x: str(x), d0values):
        if verbose:
            print '\nd0cut: {0}'.format(d0cut)
        for pid,d0dict in pid_dict.iteritems():
            if verbose:
                print 'pid: {0}'.format(pid)
            d0list = d0dict[d0cut]
            # the figure
            # Create the TGraphs/THistos
            y[pid] = []
            n['bb']=[]
            n['bb_CC']=[]
            n['bb_NC']=[]
            n['bb_NN']=[]
            n['cc']=[]
            n['cc_CC']=[]
            n['cc_NC']=[]
            n['cc_NN']=[]
            n['gg']=[]
            n['gg_CC']=[]
            n['gg_NC']=[]
            n['gg_NN']=[]
            n['uu']=[]
            n['uu_CC']=[]
            n['uu_NC']=[]
            n['uu_NN']=[]
            n['dd']=[]
            n['dd_CC']=[]
            n['dd_NC']=[]
            n['dd_NN']=[]
            n['ss']=[]
            n['ss_CC']=[]
            n['ss_NC']=[]
            n['ss_NN']=[]
            r['bb_CC']=[]
            r['bb_NC']=[]
            r['bb_NN']=[]
            r['cc_CC']=[]
            r['cc_NC']=[]
            r['cc_NN']=[]
            r['gg_CC']=[]
            r['gg_NC']=[]
            r['gg_NN']=[]
            r['uu_CC']=[]
            r['uu_NC']=[]
            r['uu_NN']=[]
            r['dd_CC']=[]
            r['dd_NC']=[]
            r['dd_NN']=[]
            r['ss_CC']=[]
            r['ss_NC']=[]
            r['ss_NN']=[]
            eff['bb']=[]
            eff['cc']=[]
            eff['ss']=[]
            eff['uu']=[]
            eff['dd']=[]
            eff['gg']=[]
            bcdf['bb']=[]
            bcdf['cc']=[]
            bcdf['ss']=[]
            bcdf['uu']=[]
            bcdf['dd']=[]
            bcdf['gg']=[]
            bcdf['bb_CC']=[]
            bcdf['bb_NC']=[]
            bcdf['bb_NN']=[]
            bcdf['cc_CC']=[]
            bcdf['cc_NC']=[]
            bcdf['cc_NN']=[]
            bcdf['ss_CC']=[]
            bcdf['ss_NC']=[]
            bcdf['ss_NN']=[]
            bcdf['uu_CC']=[]
            bcdf['uu_NC']=[]
            bcdf['uu_NN']=[]
            bcdf['dd_CC']=[]
            bcdf['dd_NC']=[]
            bcdf['dd_NN']=[]
            bcdf['gg_CC']=[]
            bcdf['gg_NC']=[]
            bcdf['gg_NN']=[]
            for p,e_signal,significance,_x1,_x2,n_signal,n_bkg,e_bkg,n_bb,n_cc,n_gg,n_uu,n_dd,n_ss_cc,n_ss_nc,n_ss_nn,n_bb_cc,n_bb_nc,n_bb_nn,n_cc_cc,n_cc_nc,n_cc_nn,n_gg_cc,n_gg_nc,n_gg_nn,n_uu_cc,n_uu_nc,n_uu_nn,n_dd_cc,n_dd_nc,n_dd_nn, eff_bb, eff_cc, eff_ss, eff_uu, eff_dd, eff_gg, bcdf_bb, bcdf_cc, bcdf_ss, bcdf_uu, bcdf_dd, bcdf_gg, bcdf_bb_CC, bcdf_bb_NC, bcdf_bb_NN, bcdf_cc_CC, bcdf_cc_NC, bcdf_cc_NN, bcdf_ss_CC, bcdf_ss_NC, bcdf_ss_NN, bcdf_uu_CC, bcdf_uu_NC, bcdf_uu_NN, bcdf_dd_CC, bcdf_dd_NC, bcdf_dd_NN, bcdf_gg_CC, bcdf_gg_NC, bcdf_gg_NN in d0list:
                y[pid].append(significance)
                n['bb'].append(n_bb)
                n['bb_CC'].append(n_bb_cc)
                n['bb_NC'].append(n_bb_nc)
                n['bb_NN'].append(n_bb_nn)
                n['cc'].append(n_cc)
                n['cc_CC'].append(n_cc_cc)
                n['cc_NC'].append(n_cc_nc)
                n['cc_NN'].append(n_cc_nn)
                n['gg'].append(n_gg)
                n['gg_CC'].append(n_gg_cc)
                n['gg_NC'].append(n_gg_nc)
                n['gg_NN'].append(n_gg_nn)
                n['uu'].append(n_uu)
                n['uu_CC'].append(n_uu_cc)
                n['uu_NC'].append(n_uu_nc)
                n['uu_NN'].append(n_uu_nn)
                n['dd'].append(n_dd)
                n['dd_CC'].append(n_dd_cc)
                n['dd_NC'].append(n_dd_nc)
                n['dd_NN'].append(n_dd_nn)
                n['ss'].append(n_signal)
                n['ss_CC'].append(n_ss_cc)
                n['ss_NC'].append(n_ss_nc)
                n['ss_NN'].append(n_ss_nn)
                bcdf['bb'].append(bcdf_bb)
                bcdf['cc'].append(bcdf_cc)
                bcdf['gg'].append(bcdf_gg)
                bcdf['uu'].append(bcdf_uu)
                bcdf['dd'].append(bcdf_dd)
                bcdf['ss'].append(bcdf_ss)
                bcdf['bb_CC'].append(bcdf_bb_CC)
                bcdf['bb_NC'].append(bcdf_bb_NC)
                bcdf['bb_NN'].append(bcdf_bb_NN)
                bcdf['cc_CC'].append(bcdf_cc_CC)
                bcdf['cc_NC'].append(bcdf_cc_NC)
                bcdf['cc_NN'].append(bcdf_cc_NN)
                bcdf['gg_CC'].append(bcdf_gg_CC)
                bcdf['gg_NC'].append(bcdf_gg_NC)
                bcdf['gg_NN'].append(bcdf_gg_NN)
                bcdf['uu_CC'].append(bcdf_uu_CC)
                bcdf['uu_NC'].append(bcdf_uu_NC)
                bcdf['uu_NN'].append(bcdf_uu_NN)
                bcdf['dd_CC'].append(bcdf_dd_CC)
                bcdf['dd_NC'].append(bcdf_dd_NC)
                bcdf['dd_NN'].append(bcdf_dd_NN)
                bcdf['ss_CC'].append(bcdf_ss_CC)
                bcdf['ss_NC'].append(bcdf_ss_NC)
                bcdf['ss_NN'].append(bcdf_ss_NN)
                eff['bb'].append(eff_bb)
                eff['cc'].append(eff_cc)
                eff['ss'].append(eff_ss)
                eff['uu'].append(eff_uu)
                eff['dd'].append(eff_dd)
                eff['gg'].append(eff_gg)
                if p==pLcut:
                    sigd0[pid].append(significance)
                    total_bcdf[pid].append(
                        save_divide(sum(map(lambda ch: n[ch][-1] * bcdf[ch][-1], ['bb','cc','gg','uu','dd'])),
                                    sum(map(lambda ch: n[ch][-1],  ['bb','cc','gg','uu','dd'])),0))

            r['bb_CC']=map(lambda i: save_divide(n['bb_CC'][i],n['bb'][i],1), xrange(0,len(n['bb'])))
            r['bb_NC']=map(lambda i: save_divide(n['bb_NC'][i],n['bb'][i],1), xrange(0,len(n['bb'])))
            r['bb_NN']=map(lambda i: save_divide(n['bb_NN'][i],n['bb'][i],1), xrange(0,len(n['bb'])))
            r['cc_CC']=map(lambda i: save_divide(n['cc_CC'][i],n['cc'][i],1), xrange(0,len(n['cc'])))
            r['cc_NC']=map(lambda i: save_divide(n['cc_NC'][i],n['cc'][i],1), xrange(0,len(n['cc'])))
            r['cc_NN']=map(lambda i: save_divide(n['cc_NN'][i],n['cc'][i],1), xrange(0,len(n['cc'])))
            r['gg_CC']=map(lambda i: save_divide(n['gg_CC'][i],n['gg'][i],1), xrange(0,len(n['gg'])))
            r['gg_NC']=map(lambda i: save_divide(n['gg_NC'][i],n['gg'][i],1), xrange(0,len(n['gg'])))
            r['gg_NN']=map(lambda i: save_divide(n['gg_NN'][i],n['gg'][i],1), xrange(0,len(n['gg'])))
            r['uu_CC']=map(lambda i: save_divide(n['uu_CC'][i],n['uu'][i],1), xrange(0,len(n['uu'])))
            r['uu_NC']=map(lambda i: save_divide(n['uu_NC'][i],n['uu'][i],1), xrange(0,len(n['uu'])))
            r['uu_NN']=map(lambda i: save_divide(n['uu_NN'][i],n['uu'][i],1), xrange(0,len(n['uu'])))
            r['dd_CC']=map(lambda i: save_divide(n['dd_CC'][i],n['dd'][i],1), xrange(0,len(n['dd'])))
            r['dd_NC']=map(lambda i: save_divide(n['dd_NC'][i],n['dd'][i],1), xrange(0,len(n['dd'])))
            r['dd_NN']=map(lambda i: save_divide(n['dd_NN'][i],n['dd'][i],1), xrange(0,len(n['dd'])))
            r['ss_CC']=map(lambda i: save_divide(n['ss_CC'][i],n['ss'][i],1), xrange(0,len(n['ss'])))
            r['ss_NC']=map(lambda i: save_divide(n['ss_NC'][i],n['ss'][i],1), xrange(0,len(n['ss'])))
            r['ss_NN']=map(lambda i: save_divide(n['ss_NN'][i],n['ss'][i],1), xrange(0,len(n['ss'])))
            plot_python_detailed(x,n,'nevents_cmp_{0}_{1}{2}'.format(d0cut,pid,SUFFIXPLOTS), '# events')
            plot_python_detailed(x,eff,'eff_cmp_{0}_{1}{2}'.format(d0cut,pid,SUFFIXPLOTS), 'efficiency')
            plot_python_detailed(x,r,'revents_cmp_{0}_{1}{2}'.format(d0cut,pid,SUFFIXPLOTS), 'fraction of events',False)
            plot_python_detailed(x,bcdf,'BCdaughterfraction_{0}_{1}{2}'.format(d0cut,pid,SUFFIXPLOTS), 'BC daughter fraction',False)
            plot_python_detailed(x,bcdf,'BCdaughterfraction_{0}_{1}_log{2}'.format(d0cut,pid,SUFFIXPLOTS), 'BC daughter fraction')

            #write efficiencies to file
            outfile='efficiencies_{0}_{1}.txt'.format(d0cut,pid)
            f = open(outfile, 'w')
            for i in range(0,len(x)):
                f.write('{0:3}   {1:.4e}   {2:.4e}   {3:.4e}   {4:.4e}   {5:.4e}   {6:.4e}\n'.format(
                    x[i], eff['bb'][i], eff['cc'][i], eff['ss'][i], eff['uu'][i],
                    eff['dd'][i], eff['gg'][i]))
            f.close()
            
            if verbose:
                print 'x={0}'.format(x)
                if False in map(lambda x: x==0., n['bb_CC']):
                    print 'n[bb_CC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['bb_CC'])).replace("'",""))
                if False in map(lambda x: x==0., n['gg_CC']):
                    print 'n[gg_CC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['gg_CC'])).replace("'",""))
                if False in map(lambda x: x==0., n['uu_CC']):
                    print 'n[uu_CC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['uu_CC'])).replace("'",""))
                if False in map(lambda x: x==0., n['dd_CC']):
                    print 'n[dd_CC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['dd_CC'])).replace("'",""))
                if False in map(lambda x: x==0., n['cc_CC']):
                    print 'n[cc_CC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['cc_CC'])).replace("'",""))
                if False in map(lambda x: x==0., n['ss_CC']):
                    print 'n[ss_CC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['ss_CC'])).replace("'",""))
                if False in map(lambda x: x==0., n['bb_NC']):
                    print 'n[bb_NC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['bb_NC'])).replace("'",""))
                if False in map(lambda x: x==0., n['gg_NC']):
                    print 'n[gg_NC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['gg_NC'])).replace("'",""))
                if False in map(lambda x: x==0., n['uu_NC']):
                    print 'n[uu_NC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['uu_NC'])).replace("'",""))
                if False in map(lambda x: x==0., n['dd_NC']):
                    print 'n[dd_NC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['dd_NC'])).replace("'",""))
                if False in map(lambda x: x==0., n['cc_NC']):
                    print 'n[cc_NC]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['cc_NC'])).replace("'",""))
                if False in map(lambda x: x==0., n['ss_NC']):
                    print 'n[ss_NC]={0}'.format(str(map(lambda i: float("{0:.2e}".format(i)), n['ss_NC'])).replace("'",""))
                if False in map(lambda x: x==0., n['bb_NN']):
                    print 'n[bb_NN]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['bb_NN'])).replace("'",""))
                if False in map(lambda x: x==0., n['gg_NN']):
                    print 'n[gg_NN]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['gg_NN'])).replace("'",""))
                if False in map(lambda x: x==0., n['uu_NN']):
                    print 'n[uu_NN]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['uu_NN'])).replace("'",""))
                if False in map(lambda x: x==0., n['dd_NN']):
                    print 'n[dd_NN]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['dd_NN'])).replace("'",""))
                if False in map(lambda x: x==0., n['cc_NN']):
                    print 'n[cc_NN]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), n['cc_NN'])).replace("'",""))
                if False in map(lambda x: x==0., n['ss_NN']):
                    print 'n[ss_NN]={0}'.format(str(map(lambda i: float("{0:.2e}".format(i)), n['ss_NN'])).replace("'",""))
                print 'eff[bb]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), eff['bb'])).replace("'",""))
                print 'eff[cc]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), eff['cc'])).replace("'",""))
                print 'eff[ss]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), eff['ss'])).replace("'",""))
                print 'eff[uu]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), eff['uu'])).replace("'",""))
                print 'eff[dd]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), eff['dd'])).replace("'",""))
                print 'eff[gg]={0}'.format(str(map(lambda i: "{0:.2e}".format(i), eff['gg'])).replace("'",""))
        plot_python(x,y,'significance_cmp_{0}{1}'.format(d0cut,SUFFIXPLOTS),'Parallel momentum cut [GeV]','Significance')
    plot_python(map(lambda x: 1000*x,d0values), sigd0, 'significance_cmp_d0_{0}{1}'.format(pLcut,SUFFIXPLOTS), 'd0 cut in $\mu$m','Significance')
    plot_python(map(lambda x: 1000*x,d0values), total_bcdf, 'total_BCdf_cmp_d0_{0}{1}'.format(pLcut,SUFFIXPLOTS), 'd0 cut in $\mu$m','BC daughter fraction in bkg')

if __name__ == '__main__':
    from argparse import ArgumentParser
    import os

    #Opciones de entrada
    usage = "Plot processing from a input file created with processhzroot script"
    parser = ArgumentParser(prog='hzplots',description=usage)

    # Two use modes: 
    # 1. using as inputs a root file with the same PID definitions for 
    # signal and backgrounds: will create efficiency curves (ROC), significance 
    # curves, and the pickle dictionary
    # 2. using as input a pickle file containing the efficiencies per d0,pL cut, etc... 
    # will create purity curves, pion rejection curves and significance curves comparative
    # using different PID
    subparsers = parser.add_subparsers(title='subcommands',description='valid subcommands',
            help='additional help')


    # just general option to be appearing in both commands
    #parser.add_argument( '-s', '--suffix', action='store', dest='suffixout',\
    #        help="output suffix for the plots [.pdf]")
    
    # 1. same PID
    parser_pid = subparsers.add_parser("fixed_pid",help="Use the same PID from the root input"\
            " file to build a dictionary with efficiencies, number of events, etc per d0 and pL-cut"\
            " and stores this information in the 'd0cuts_dict_PIDUSED.pkl' file. Plus some plots"\
            " will be created: efficiency curves (ROC), significance curves, ...")
    parser_pid.add_argument('channel_mode',nargs=1,help='The PID tree to extract in the ROOTFILE.'\
            ' The tree \'middle\' name of the pre-processed file (with'\
            ' processedhz script) is then used. Note the tree name is created in the'\
            ' \'processhz script\' and should follow the standard: '\
            ' "mctree_n-n-n-PID_channel_blahblah.root"')
    parser_pid.add_argument('rootfile',nargs=1,help='Input root file, '\
            'created with processedhz script')

    parser_pid.add_argument( '-s', '--suffix', action='store', dest='suffixout',\
            help="output suffix for the plots [.pdf]")
    
    parser_pid.add_argument( '-d', '--d0', action='store', nargs='+', dest='d0',\
            metavar='d01[,d02,...]',help="Impact parameters cuts [0.013 0.02 0.05]")
    parser_pid.add_argument( '-z', '--z0', action='store',  dest='z0',\
            help="activate the z0-cut and the value to cut [False]")
    parser_pid.add_argument( '-m', '--multiplicity', action='store',  dest='Ntrack',\
            help="activate the multiplicity cut (the max. number of allowed tracks"\
            " surrounding the final hadron at a dR-cone) and the value to cut [False]")
    parser_pid.add_argument( '--d0cut-type', action='store', dest='d0cut_type',\
            metavar='circular|square',help="functional of the d0 cut [circular]")
    parser_pid.add_argument( '-p', '--pL-cut-min', action='store', dest='pLMin',\
            help="Momentum minimum cut [10 GeV]")
    parser_pid.add_argument( '-P', '--pL-cut-max', action='store', dest='pLMax',\
            help="Momentum maximum cut [20 GeV]")
    parser_pid.add_argument( '--pLcut-type', action='store', dest='pLcut_type',\
            metavar="circular|square|line",help="functional of the pL cut [circular]")
    parser_pid.add_argument( '-R', '--R0', action='store', dest='R0',\
            help="minimum decay length of the K_s [500]")
    parser_pid.add_argument( '-t', '--tables', action='store_true', dest='tables',\
            help="whether or not print the latex tables")
    parser_pid.add_argument( '-w','--working-points', action='store_true', dest='wp_activate',\
            help="whether or not plot the some working points in the plots")
    parser_pid.add_argument( '-l', '--ligth-channels',action='store_true',dest='light_channels',\
            help='whether or not add also the uubar and ddbar channels')
    parser_pid.add_argument( '-c', '--charge_combination', action='store', dest='charge_combination',\
            metavar="all|CC|NC|NN",help="Which final states to choose [all]")
    
    parser_pid.set_defaults(which='fixed_pid',
            pLMin=10,
            pLMax=20,
            d0=[0.013,0.02,0.05],z0=None,
            pLcut_type='circular',d0cut_type='circular',
            channel_mode='0.5-0.08-0.75-PID',
            charge_combination='all',R0=500)    
    
    # 2. decay chain for B/D ancestors
    parser_decaychain = subparsers.add_parser("decay_chain",help='Print the decay of the first'\
            'B or D-hadron ancestor (backward from the final hadron)' )
    parser_decaychain.add_argument('rootfile',nargs='+',action='store',help='The input root file')
    parser_decaychain.add_argument('-p','--pcut',action='store',dest='pcut',type=float,
            help='The parallel momentum cut to be applied for both final hadrons (in GeV)')
    parser_decaychain.add_argument('-d','--d0cut',action='store',dest='d0cut',type=float,
            help='The maximum impact parameter to be allowed for both final hadrons (in mm)')
    parser_decaychain.add_argument('-n','--nfirst',action='store',dest='nfirst',type=int,
            help='Maximum number of decays to show [Default: 10]')
    parser_decaychain.add_argument('--latex',action='store_true',dest='want_latex',help='Print also a latex table')
    parser_decaychain.set_defaults(which='decay_chain',nfirst=10)

    # 3. input pickle file
    parser_cmp_pid = subparsers.add_parser("compare_pid",help='Compare the performance for the'\
            ' different PID assumptions (significance plots and tables for the different PIDs)')
    parser_cmp_pid.add_argument('pickle_files',nargs='+',help='List of pickles files (created'\
            ' previously by the command "fixed_pid" to be compared')
    parser_cmp_pid.add_argument( '-s', '--suffix', action='store', dest='suffixout',\
            help="output suffix for the plots [.pdf]")
    parser_cmp_pid.add_argument( '-v', '--verbose', action='store_true', dest='verbose',
                                 help='Print lists that are plotted')
    parser_cmp_pid.add_argument( '-p', action='store', dest='pLcut', type=float,
            help='the value of the p_|| cut in GeV used for significance plot when looping over'\
            ' d0 cuts [Default:15]')
    parser_cmp_pid.set_defaults(which='compare_pid',
                                suffixout='.pdf',
                                verbose=False,
                                pLcut=15)
    
    args = parser.parse_args()
    
    if args.which in [ 'fixed_pid', 'compare_pid' ] and args.suffixout:
        suff = args.suffixout.replace('.','')
        globals()['SUFFIXPLOTS'] ='.'+suff
    
    if args.which == 'fixed_pid':
        # Which trees should be used?
        pre_channels = [ 'ssbar','bbbar','ccbar', 'gg' ]
        if args.light_channels:
            pre_channels += [ 'uubar', 'ddbar' ]
        channels = [ "{0}_{1}".format(x,args.channel_mode[0]) for x in pre_channels ]
        main_fixed_pid(os.path.abspath(args.rootfile[0]),channels,
                args.tables,
                int(args.pLMin),
                int(args.pLMax),
                args.pLcut_type,
                args.d0,
                args.d0cut_type,
                args.z0,
                args.Ntrack,
                args.wp_activate,
                args.R0,
                args.charge_combination)
    elif args.which == 'decay_chain':
        main_decay_chain(args.rootfile,args.want_latex,args.nfirst,
                pcut=args.pcut,d0cut=args.d0cut)
    elif args.which == 'compare_pid':
        main_cmp_pid(args.pickle_files, args.verbose, args.pLcut)


