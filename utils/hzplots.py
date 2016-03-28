#!/usr/bin/env python
""":script:`doplots` -- Plotting stuff created with processhzroot script
======================================================================

.. script:: pion_misid [OPTIONS]    
      :platform: Unix
      :synopsis: Evaluate pion misid...
	  .. moduleauthor:: Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
"""
KAON = 321
PION = 211
def isKaon(index):
    """check if the gen particle is kaon

    Parameters
    ----------
    index: int
        the index of the gen-particle in the tree
    """
    return "abs(pdgId[{0}]) == {1}".format(index,KAON)

def isPion(index):
    """check if the gen particle is kaon

    Parameters
    ----------
    index: int
        the index of the gen-particle in the tree
    """
    return "abs(pdgId[{0}]) == {1}".format(index,PION)

FS_CONDITION = { 'KK': '{0} && {1}'.format(isKaon(0),isKaon(1)), \
        'KP': '({0} && {1} || {2} && {3})'.format(isKaon(0),isPion(1),isPion(0),isKaon(1)),
                'PP': '{0} && {1}'.format(isPion(0),isPion(1)) 
                }

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
        self.mH         = 125.09 # (GeV)
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

def get_tree_name(tree_names,decay_channel,hadron_state):
    """Given a standarized (TTree) tree list of names (obtained from the 
    processedhzroot.py script), it returns the name of the histogram matching
    the final state

    Parameters
    ----------
    three_names: list(str)
        the names of all the trees found in a root processed file by the
        processedhzroot script
    hadron_state: str
        the two opposite-hemisphere final state hadrons: KK, KP, PP (K-kaon,
        P-pion)

    Return
    ------
    str, the name of the input list matching the final state
    """
    try:
        tree_name = filter(lambda x: 
                x.find('mctree_{0}'.format(hadron_state)) ==0 and  
                    x.find(decay_channel) != -1,tree_names)[0]
    except IndexError:
        raise RuntimeError('Not found the tree "mctree_{0}_*_{1}"'\
                ' in the root file'.format(decay_channel,hadron_state))
    return tree_name

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
        the hadrons to look at, which can be the: KK, KP, PP (K-kaon, P-pion)
    pLcut: float, optional, incompatible with pLbin  
    d0cut: float, optional, incompatible with d0binL, d0binH
    pLbin: int, optional, incompatible with pLcut
    d0binH: int, optional, incompatible with d0cut
    d0binL: int, optional, incompatible with d0cut

    Returns
    -------
    efficiency: float
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

def d0(index):
    """Impact parameter extrapolation considering a straight line

    Parameters
    ----------
    index: int
        the index of the gen-particle in the tree
    """
    return  "(vy[{0}]-vx[{0}]*phi_lab[{0}])*cos(phi_lab[{0}])".format(index)

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
    def __init__(self,decay_channel):
        """Encapsulates the efficieny cut of a Hadron-hadron event
        Incorporates the cut in impact parameter and in parallel momentum

        Parameters
        ----------
        decay_channel: str
            the Higgs hadronic decay (bbbar,ccbar,ssbar, ddbar, uubar)
        """
        self.decay_channel = decay_channel
        self.initialized = False
        self.current_d0cut = None
        self.current_pLcut = None
        
        self.final_state_hadrons = ['KK','KP','PP']
        self.__entrylist_hadrons = dict(map(lambda i: (i,None),self.final_state_hadrons))
        
        self.pLcut_function="sqrt( (p[0]*cos(theta[0]))**2. + (p[1]*cos(theta[1]))**2.)"
        self.d0cut_function="sqrt( {0}**2.0 + {1}**2.0)".format(d0(0),d0(1))

        # optimization data-members
        self.__tree = None
        self.__entrylist_d0cut = None
        self.__entrylist_pLcut = None
        self.__finalhadrons_entries = { 'KK': None, 'KP': None, 'PP': None }

    
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


    def set_total_eff(self,treedict,pLcut,d0cut,verbose=False):
        """builds the total efficiency and evaluates it for the different 
        final state cases KK, KP and PP (if exist)

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
        the ROOt.gDirectory.Get method retrieve the first object in memory, 
        obtaining wrong results
        """
        import ROOT 
        #ROOT.gROOT.SetBatch()

        try:
            treename = filter(lambda x: x.find("mctrue_{0}".format(self.decay_channel)) != -1,\
                    treedict.keys())[0]
        except IndexError:
            raise AttributeError("\033[1;31mERROR:\033[1;m not found "\
                    "decay channel '{0}' tree".format(self.decay_channel))
        # get the tree and check if it was used before
        tree = treedict[treename]
        is_new_tree = False
        if self.__tree != tree:
            self.__tree = tree
            self.__tree_entries = self.__tree.GetEntries()
            # create the entrylist for final state hadrons final state
            for i in self.final_state_hadrons:
                el_HSname = "{0}_entrylist_hadrons_{1}".format(treename,i)
                self.__tree.Draw(">>{0}".format(el_HSname),FS_CONDITION[i],"entrylist")
                self.__entrylist_hadrons[i] = ROOT.gDirectory.Get("{0}".format(el_HSname))
                self.__finalhadrons_entries[i] = self.__entrylist_hadrons[i].GetN()
            is_new_tree = True

        is_new_d0cut = False
        ## allow optimize loop
        if not self.current_d0cut or self.current_d0cut != d0cut:
            self.current_d0cut = d0cut
            is_new_d0cut = True
        elif is_new_tree:
            is_new_d0cut = True

        if is_new_d0cut:
            el_d0name = "{0}_entrylist_d0cut_{1:.1f}".format(treename,self.current_d0cut)
            # create the entrylist 
            self.__tree.Draw(">>{0}".format(el_d0name),"{0} < {1}".format(\
                    self.d0cut_function,self.current_d0cut),"entrylist")
            self.__entrylist_d0cut = ROOT.gDirectory.Get("{0}".format(el_d0name))

        is_new_pLcut = False
        if not self.current_pLcut or self.current_pLcut != pLcut:
            self.current_pLcut = pLcut
            is_new_pLcut = True
        elif is_new_tree:
            is_new_pLcut = True
        
        if is_new_pLcut:
            el_pLname = "{0}_entrylist_pLcut_{1:.1f}".format(treename,self.current_pLcut)
            # create the entrylist (or extracted from memory if is already there)
            #self.__entrylist_pLcut = ROOT.gDirectory.Get("{0}".format(el_pLname))
            #if isinstance(k,ROOT.TObject):
                # there wasn't there, so build it
            self.__tree.Draw(">>{0}".format(el_pLname),"{0} > {1}".format(\
                    self.pLcut_function,self.current_pLcut),"entrylist")
            self.__entrylist_pLcut = ROOT.gDirectory.Get("{0}".format(el_pLname))

        # return if everything is done
        ### NEW if not is_new_tree and not is_new_d0cut and not is_new_pLcut:
        ### NEW   return
        
        # main loop
        for i in self.final_state_hadrons:
            #cut_entries_OLD = tree.GetEntries(\
            #        "{0} > {1} && {2} < {3} && {4}".format(\
            #        self.pLcut_function,pLcut,\
            #        self.d0cut_function,d0cut,\
            #        FS_CONDITION[i]))
            #n_hadrons_entries = tree.GetEntries(FS_CONDITION[i])
            cut_entries = get_common_entries( [self.__entrylist_d0cut,\
                  self.__entrylist_pLcut, self.__entrylist_hadrons[i]] )
            #if cut_entries != cut_entries_OLD:
            #    print 
            #    print "\033[1;33mWARNING\033[1;m"
            #    print self.__tree.GetName(),self.current_d0cut,self.current_pLcut
            #    print cut_entries,cut_entries_OLD,"+"*20
            #    print self
            _eff  = 0.0
            _prob = 0.0
            if self.__tree_entries != 0:
                _eff = float(cut_entries)/float(self.__tree_entries)
                _prob= float(self.__finalhadrons_entries[i])/float(self.__tree_entries)
            # p ( d0 pL | L_AB decay_channel I0 )
            self.__setattr__('eff_cut_{0}'.format(i),_eff)
            # p ( L_AB | decay_channel I0 )
            self.__setattr__('p_{0}'.format(i),_prob)
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
    for d0cut,val_list in points_dict.iteritems():
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
            leg[d0cut].AddEntry(graphs[d0cut][-1],'d0={0} mm:  {1}'.format(d0cut,text),'P')
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
    for d0cut,pL_list in points_dict.iteritems():
        g = ROOT.TGraph(len(pL_list))
        d_graphs.append(g)
        g.SetLineWidth(2)
        g.SetLineColor(getcolor()[k])
        g.SetLineStyle(3*k)
        g.SetMarkerColor(getcolor()[k])
        g.SetMarkerStyle(20+k)
        leg.AddEntry(g,'d_{0}={1} mm'.format('{0}',d0cut),'PL')
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
    columns = ( 'p_{||}^{c}\;[GeV]','\\varepsilon_{ss}', 'S/\\sqrt{B}', 'N_{bkg}/N_{s\\bar{s}}', \
            'Kaon Purity', 'N_{s\\bar{s}}', 'N_{bkg}')
    colformat = map(lambda x: ' ${0:'+str(x)+'}$ &' , \
            ( '.0f','.3f','.2f','.0f','.3f', '.0f', '.0f' ))
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

def main(rootfile,channels,tables,pLMax,d0cuts,wp_activated):
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
    pLMax: float
        the parallel momentum cut considered
    d0cuts: list(float)
        the list of d0-cuts to consider
    wp_activated: bool
        whether or not to plot on the working points...
    """
    import ROOT
    import sys
    from math import sqrt

    
    # Get the root file with and the TTree object
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

    signal_PID   = filter(lambda x: x.find('ssbar_PID') != -1,channels)[0]
    signal_noPID = filter(lambda x: x.find('ssbar_noPID') != -1,channels)[0]

    # check if the noPID histograms are present
    try:
        _test = filter(lambda name: name.find('ssbar_noPID') != -1,_obj['H'].keys())[0]
        purity_evaluation = True
    except IndexError:
        # not found the related histograms, we cannot calculate purity
        # just warn, but continue with the plots
        print "\033[1;33mhzplots WARNING\033[1;m No 'ssbar_noPID' available." \
                " The purity terms are not evaluated"
        purity_evaluation = False
        # and remove the signal_noPID from the channels list
        channels.remove(signal_noPID)
    
    # -- just take care in the H-ressonance...
    pLcuts = xrange(0,pLMax+1)

    # the eff. classes 
    eff = dict(map(lambda x: (x,eff_cut_hadron(x)), channels))
    
    message = "\r\033[1;34mhzplots INFO\033[1;m Obtaining the data"
    # { 'docut1': [ (pLcut1, eff_sig, significance, pion_rejection, purity, N_sig, N_bkg), ... ],  }
    observables = {}
    # Indices
    I_PL = 0; I_EFF_SIG = 1; I_SIGN=2; I_PION_REJEC = 3; I_PURITY = 4; I_N_SIGNAL = 5; I_N_BKG=6;

    for d0 in d0cuts:
        d0str = '{0:.1f}'.format(d0)
        observables[d0str] = []
        i = 0
        for pL in pLcuts:
            sys.stdout.write( "{0} {1}".format(message,SPINNING[i % len(SPINNING)]))
            sys.stdout.flush()

            # setting the current cuts to all the efficienciesa
            _dummy = map(lambda e: e.set_total_eff(_obj['H'],pLcut=pL,d0cut=d0), eff.values())
            
            # Some needed values
            eff_sig = eff['ssbar_PID'].get_total_eff('KK')
            n_KK    = eff['ssbar_PID'].get_total_events()
            bkg_tot_evts = sum(map(lambda (x,y): y.get_total_events(),\
                    filter(lambda (x,y): x != signal_PID or x != signal_noPID, eff.iteritems() )))
            # --- > probably interesting also...
            #bkg_eff = sum(map(lambda (x,y): y.get_total_eff(),\
            #        filter(lambda (x,y): x != signal_PID or x != signal_noPID, eff.iteritems() )))
            
            # Signal efficiency calculation: assuming 100% p-K separation with the PID 
            # effsig := eff['ssbar_PID']
            # purity := 2*N_KK/(2*N_KK+2*N_KP+2*N_KK) = N_KK/(N_KK+NKP+N_PP)
            # Pion rejection: N_b/N_signal assuming no p-K separation (no PID)
            # Note: sigma_HZ*L_int*Sum_qq BR(H->qq)/Sum_qq BR(h->qq)*(Sum_qq BR(H->qq) eff_qq*f_qq)
            # the Sum_qq BR(H->qq) terms are cancelled... 
            pion_rejection = float(bkg_tot_evts)/float(n_KK)
            purity = -1
            if purity_evaluation:
                purity         = float(eff['ssbar_noPID'].get_events('KK'))/\
                        (float(eff['ssbar_noPID'].get_events('KK'))+\
                        float(eff['ssbar_noPID'].get_events('KP'))+\
                        float(eff['ssbar_noPID'].get_events('PP')))
                purity_evaluation = True
            significance   = float(n_KK)/sqrt(float(bkg_tot_evts))

            # store it
            observables[d0str].append( (pL,eff_sig,significance,pion_rejection,purity,n_KK,bkg_tot_evts) )
            i+=1
    # plotting
    print
    import ROOT
    print ROOT.gDirectory.Print()
    print "\033[1;34mhzplots INFO\033[1;m Plotting..."
    # --- Some extra points (WP)
    if wp_activated:
        leg_format_pion_rej = ( 'S/#sqrt{B}=%.1f @ p_{||}^{c}=%.0f GeV',(I_SIGN,I_PL) )
        graphs_leg_pion_rej = get_point_graphs(observables,(I_EFF_SIG,I_PION_REJEC),I_PL,[10,20],leg_format_pion_rej)
    
        leg_format_pur = None
        graphs_leg_pur = None
        if purity_evaluation:
            leg_format_pur = ( 'S/#sqrt{B}=%.1f @ p_{||}^{c}=%.0f GeV',(I_SIGN,I_PL) )
            graphs_leg_pur = get_point_graphs(observables,(I_EFF_SIG,I_PURITY),I_PL,[10,20],leg_format_pur)
        
        leg_format_sig = ( '#varepsilon_{signal}=%.2f, #pi-rej.factor=%.0f',(I_EFF_SIG,I_PION_REJEC) )
        graphs_leg_sig = get_point_graphs(observables,(I_PL,I_SIGN),I_EFF_SIG,[0.4,0.8],leg_format_sig)
    else:
        graphs_leg_pion_rej = None
        graphs_leg_pur      = None
        graphs_leg_sig      = None

    pr_attr = plot_attributes('pion_rejection',
            xtitle='#varepsilon_{S}', 
            ytitle='pion rejection (N_{BKG}/N_{s#bar{s}} with no PID)',
            x0 = 0.0, x1 = 1.0 , y0 = 0.0 )
    make_plot(observables,(I_EFF_SIG,I_PION_REJEC),pr_attr,g_points_dict=graphs_leg_pion_rej)

    if purity_evaluation:
        purity_attr = plot_attributes('purity',
                xtitle='#varepsilon_{S}', 
                ytitle='purity',
                x0 = 0.0, x1 = 1.0 , y0 = 0.0, y1= 0.0 )
        make_plot(observables,(I_EFF_SIG,I_PURITY),purity_attr,g_points_dict=graphs_leg_pur)

    significance_attr = plot_attributes('significance',
            xtitle='p_{||}^{c}', xunit = '[GeV]', 
            ytitle='N_{S}/#sqrt{N_{B}}',
            x0 = 0.0, y0 = 0.0 )
    make_plot(observables,(I_PL,I_SIGN),significance_attr,g_points_dict=graphs_leg_sig,leg_position="DOWN")

    # Tables
    if tables:
        for d0,obsList in observables.iteritems():
            print "Table for d0: {0} mm".format(d0)
            print "----------------------"
            print get_latex_table(obsList)


if __name__ == '__main__':
    from optparse import OptionParser,OptionGroup
    import os

    #Opciones de entrada
    usage = "usage: hzplots INPUTFILE [options]"
    parser = OptionParser(usage=usage)
    parser.set_defaults(pLMax=30,d0="0.1,0.3,0.5")    
    parser.add_option( '-s', '--suffix', action='store', type='string', dest='suffixout',\
            help="output suffix for the plots [.pdf]")
    parser.add_option( '-d', '--d0', action='store', type='string', dest='d0',\
            metavar='d01[,d02,...]',help="Impact parameters cuts [0.1,0.3,0.5]")
    parser.add_option( '-p', '--pL-cut', action='store', type='string', dest='pLMax',\
            help="Circular momentum maximum cut [30 GeV]")
    parser.add_option( '-t', '--tables', action='store_true', dest='tables',\
            help="whether or not print the latex tables")
    parser.add_option( '-w','--working-points', action='store_true', dest='wp_activate',\
            help="whether or not plot the some working points in the plots")
    parser.add_option( '-l','--leading-hadrons', action='store_true', dest='doTH2Plots',\
            help="whether or not activate the paralel momentum 2-dim plots of the"\
            " leading hadrons")
    
    (opt,args) = parser.parse_args()
    
    if len(args) < 1:
        message = "\033[31;1mhzroot ERROR\033[1;m Missing input file(s)"
        raise RuntimeError(message)

    if opt.suffixout:
        suff = opt.suffixout.replace('.','')
        globals()['SUFFIXPLOTS'] ='.'+suff
    channels = [ 'ssbar_PID', 'ssbar_noPID', 'bbbar', 'ccbar'] #'ddbar', 'uubar']
    
    d0cuts = []
    for d0cut in sorted(opt.d0.split(','),key=lambda x: float(x)):
        d0cuts.append(float(d0cut))

    main(os.path.abspath(args[0]),channels,
            opt.tables,
            int(opt.pLMax),
            d0cuts,
            opt.wp_activate)
