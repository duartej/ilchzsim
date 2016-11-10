#!/usr/bin/env python
""":script:`processhzroot` -- Process the ROOT files obtained from ilchz executable                                                            
===================================================================================

.. script:: processhzroot     
      :platform: Unix
      :synopsis: Process the root files created by the ilchz executable in order
                 to obtain histogram and efficiency objects

.. moduleauthor:: Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
FIXME: MISSING DOCUMENTATION
"""
SUFFIXPLOTS='.pdf'
KAON_ID = 321
PION_ID = 211

import functools

@functools.total_ordering
class hadron(object):
    """Encapsulating the hadron
    """
    def __init__(self,**kwd):
        """
        """
        for key,val in kwd.iteritems():
            setattr(self,key,val)
        if self.p < 0:
            self.hemisphere = -1
            self.p = abs(self.p)
        else:
            self.hemisphere = 1

    def __eq__(self,other):
        return self.p == other.p
    
    def __lt__(self,other):
        return self.p < other.p
    
    def __gt__(self,other):
        return self.p > other.p

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

def get_leading_kaons(tree,applycharge):
    """
    """
    from math import cos,sqrt
    import os
    import sys

    # auxiliar function to obtain the signed (+1 top hemisphere, 
    # -1 bottom hemispher) parallel momentum
    signed_pm = lambda _k: tree.p[_k]*cos(tree.theta[_k])
    d0_f      = lambda _k: (tree.vy[_k]-tree.vx[_k]*tree.phi_lab[_k])*cos(tree.phi_lab[_k])
    z0_f      = lambda _k: (tree.vy[_k]-tree.vz[_k]*tree.theta_lab[_k])*cos(tree.theta_lab[_k])
    L_f       = lambda _k: sqrt(tree.vx[_k]**2.0+tree.vy[_k]**2.0)
    R_f       = lambda _k: sqrt(tree.vx[_k]**2.0+tree.vy[_k]**2.0+tree.vz[_k]**2.0)
    nentries = tree.getentries()
    leading_kaons = {}
    ## count the particle multiplicity per quark
    nM = []
    msg = "Evaluating {0}...".format(tree._rootfiles[0])
    if len(msg) > 100:
        shorten_name = "{0}.../{1}".format(tree._rootfiles[0][:50],
                os.path.basename(tree._rootfiles[0]))
        msg = "Evaluating {0}...".format(shorten_name)
    pointpb = float(nentries)/100.0
    for i in xrange(nentries):
        # Progress bar
        sys.stdout.write("\r\033[1;34m+-- \033[1;m"+msg+\
                "[ " +"\b"+str(int(float(i)/pointpb)+1).rjust(3)+"%]")
        sys.stdout.flush()
        _dummy=tree.getentry(i)
        ## obtain just higgs daughters
        kaons_pm = []
        ## multiplicity
        _nM_evt = {}
        for k in xrange(tree.catchall.size()):
            if tree.catchall[k] != 25:
                continue
            # parallel momentum with sign, and charge
            ###kaons_pm.append((signed_pm(k),abs(tree.pdgId[k])/tree.pdgId[k],k))
            kaons_pm.append(  hadron(p=signed_pm(k),
                                charge=abs(tree.pdgId[k])/tree.pdgId[k],
                                d0 = d0_f(k),
                                z0 = z0_f(k),
                                L  = L_f(k),
                                R  = R_f(k),
                                cosTheta = cos(tree.theta[k]),
                                pdgId = tree.pdgId[k],
                                index = k,
                                evt   = i,
                                ) )
            # FIXME: WRONG! this is considering only the hadrons we decide to 
            #        keep at the ilchz exec. (kaons or kaons-pions), this 
            #        variable should be included in the c++ code. 
            #        Just by now be careful in the meaning: it is counting
            #        the hadron multiplicity of a quark (but just counting
            #        the kaons or kaons and pions, depending on the launched
            #        mode)
            # count how many hadrons proceed from the same quark
            # 
            try:
                _nM_evt[tree.motherindex] += 1
            except KeyError:
                _nM_evt[tree.motherindex] = 1

        # sort in decrease order, the first is the highest p in
        # the top hemisphere, and last one, the highest p in the
        # bottom hemisphere
        sorted_kaons= sorted(kaons_pm,reverse=True)
        # separate between both hemispheres (using the 
        # sign of the parallel momentum)
        ###u_h = filter(lambda x: x[0] > 0,sorted_kaons)
        # Re-order again because we're recovering the positive
        # sign for the parallel momentum
        ###d_h = sorted(map(lambda y: (abs(y[0]),y[1],y[2]),\
        ###        filter(lambda x: x[0] < 0,sorted_kaons)),reverse=True)
        u_h = filter(lambda x: x.hemisphere == 1, sorted_kaons)
        d_h = filter(lambda x: x.hemisphere == -1, sorted_kaons)
        # -- continue if none was found
        if len(u_h) < 1 or len(d_h) < 1:
            #leading_kaons[i] = ()
            continue
        # -- check charge if needed
        if applycharge:
            for ku in u_h:
                try:
                    #opposite_down = filter(lambda kd: ku[1]*kd[1] < 1,d_h)[0]
                    opposite_down = filter(lambda kd: ku.charge*kd.charge < 1,d_h)[0]
                    leading_kaons[i] = (ku, opposite_down)
                    break
                except IndexError:
                    pass
        else:
            leading_kaons[i] = (u_h[0],d_h[0])

        # multiplicity: added to the general counter
        if len(nM) < 2:
            # just taking into account cases with no charge particles
            nM.append(0)
            if nM < 1:
                nM.append(0)
        _dummy = map(lambda x: nM.append(x),_nM_evt.values())
    print
    return leading_kaons,nM

def init_tree(filename):
    """
    """
    from PyAnUtils.retrievetrees import plaintree
    # dummy class inheriting from stored tree

    return plaintree(filename,'mctrue')

def momentum_1d(p1,p2):
    """ momentum_1d: R x R --> R
    """
    from math import sqrt
    return sqrt(p1*p1+p2*p2)

def pythonize(x):
    import ROOT

    # Fix the string case    
    return eval('ROOT.std.'+x.replace('<','(').replace('>',')').replace('string','str')+'()')

def store_hadrons(outname,hadronlist,old_tree,treename):
    """
    """
    import ROOT
    import sys

    f = ROOT.TFile.Open(outname,'UPDATE')
    if f.IsZombie():
        raise IOError('Problems opening ROOT file {0}'.format(outname))

    # Tree creation
    tree = ROOT.TTree(treename,'leading and subleading hadrons')
    #tree = old_tree.CloneTree(0)
    # -- initialize the containers to be filled
    #pythonize = lambda x:  eval('ROOT.std.'+x.replace('<','(').replace('>',')')+'()')
    containers = dict(map(lambda x: (x.GetName(),pythonize(x.GetClassName())),old_tree.GetListOfBranches()))
    # -- setting branch addresses
    _dumm = map(lambda (bname,vobject): tree.Branch(bname,vobject),sorted(containers.iteritems()))

    # -- old tree, set containers
    oldcont = dict(map(lambda x: (x.GetName(),pythonize(x.GetClassName())),old_tree.GetListOfBranches()))
    _dumm = map(lambda (bname,vobject): old_tree.SetBranchAddress(bname,vobject),oldcont.iteritems())
    #_dumm = map(lambda (bname,vobject): tree.SetBranchAddress(bname,vobject),containers.iteritems())


    msg = "Copying leading hadrons ..."
    pointpb = float(len(hadronlist))/100.0

    for (_ip,(h_l,h_sl)) in enumerate(hadronlist):
        # Progress bar
        sys.stdout.write("\r\033[1;34m+-- \033[1;m"+msg+\
                "[ " +"\b"+str(int(float(_ip)/pointpb)+1).rjust(3)+"%]")
        sys.stdout.flush()
        # get the indices of the elements
        k_l,k_sl = h_l.index,h_sl.index
        # ordered
        ind_keep = sorted([k_l,k_sl])
        # fill the containers with the event entries
        _dumm = old_tree.GetEntry(h_l.evt)
        # fill the branches, removing all the others items 
        # but the selected hadrons
        for (branch_name, vobject) in containers.iteritems():
            vobject.clear()
            vobject.reserve(2)
            # using the two leading hadrons (index) and filling it
            _dumm = map(lambda k: vobject.push_back(oldcont[branch_name][k]),[k_l,k_sl])
        tree.Fill()
    print
    f.Write()
    f.Close()
        

def main(args,suffixout,is_charge_considered,outfilename):
    """
    """
    import os

    if suffixout:
        suff = suffixout.replace('.','')
        globals()['SUFFIXPLOTS'] ='.'+suff
    # evaluating a list of files
    absfilenames = map(lambda fname: os.path.abspath(fname),args)
    hc = None
    for fname in absfilenames:
        sname = os.path.basename(fname).replace(".root","").\
                replace("hz","")#.replace("_PID_","").replace('kaons_only','')
        # -- initialize file
        t = init_tree(fname)
        # get the leading kaons dict { event#: ((up_pm,k),(down,k)), ... } 
        # and multiplicity of hadrons per quark
        leading_kaons,nM = get_leading_kaons(t,is_charge_considered)
        # put always higher pL in position 0
        ordered_lk = map(lambda x: sorted(x, reverse=True),leading_kaons.values())
        # persistency, copy of the original tree but keeping 
        #  only the leading hadrons
        store_hadrons(outfilename,ordered_lk,t._tree,"mctrue_"+sname)

if __name__ == '__main__':
    from optparse import OptionParser,OptionGroup

    #Opciones de entrada
    usage = "usage: processhzroot INPUTFILE1 [INPUTFILE2 ...] [options]"
    parser = OptionParser(usage=usage)
    parser.set_defaults(hadrons='kaons',notcharge=False,outfname='processed.root')    
    parser.add_option( '-o', '--outfile', action='store', type='string', dest='outfname',\
            help="output filename to persistify the created histograms [processed.root]")
    parser.add_option( '-s', '--suffix', action='store', type='string', dest='suffixout',\
            help="output suffix for the plots [.pdf]")
    parser.add_option( '--no-charge', action='store_true',  dest='notcharge',\
            help="not applying the opposite charge requirement"\
            "between leading kaons")
    
    (opt,args) = parser.parse_args()

    if len(args) < 1:
        message = "\033[31mprocesshzroot ERROR\033[m Missing input file(s)"
        raise RuntimeError(message)
    
    main(args,opt.suffixout,(not opt.notcharge),opt.outfname)

