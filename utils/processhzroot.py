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
KSHORT_ID = 310
USEKSHORTS=0

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

    def __add__(self,other):
        return self.p + other.p
    
    def __radd__(self,other):
        return self.p + other
    
def get_opposite_charge(ref_hadron, the_other_list):
    for kref in ref_hadron:
        try:
            opposite_other = filter(lambda kd: abs(kref.charge+kd.charge) < 1.5,the_other_list)[0]
            return (kref,opposite_other)
        except IndexError:
            pass
    # didn't found it
    return (-1,-1)

def smear(value, sigma):
    """ Smear value with gaussian
    """
    import numpy as np

    return np.random.normal(value, sigma)

def get_leading_kaons(tree,applycharge,useKshorts):
    """Obtain the leading kaons
    """
    from math import cos,sqrt,tan,sin,atan2
    import os
    import sys
    
    # resoltuion of the IP reconstruction
    IP_resolution = 0.000 # mm

    # FIXME: As the d0 and z0 are calculated here, it could be useful to
    #        included them in the trees later

    # auxiliar function to obtain the signed (+1 top hemisphere, 
    # -1 bottom hemispher) parallel momentum
    signed_pm = lambda _k: tree.p[_k]*cos(tree.theta[_k])
    #d0_f      = lambda _k: (tree.vy[_k]-tree.vx[_k]*tan(tree.phi_lab[_k]))*cos(tree.phi_lab[_k])
    #z0_f      = lambda _k: -(tree.vy[_k]-tree.vz[_k]*tan(tree.theta_lab[_k]))/tan(tree.theta_lab[_k])

    # distance in transverse plane of the decay vertex to the IP
    L_f       = lambda _k: sqrt(tree.vx[_k]**2.0+tree.vy[_k]**2.0)
    # distance in 3D of the decay vertex to the IP
    R_f       = lambda _k: sqrt(tree.vx[_k]**2.0+tree.vy[_k]**2.0+tree.vz[_k]**2.0)
    # impact parameter in transverse plane
    d0_f      = lambda _k: sin(atan2(tree.vy[_k],tree.vx[_k])-tree.phi_lab[_k])*L_f(_k)
    # longitudinal impact parameter
    z0_f      = lambda _k: (d0_f(_k)-L_f(_k))/tan(tree.theta_lab[_k])+tree.vz[_k]
    # resolution of impact parameter d0 in mm
    # ILC
    d0_resolution = lambda _k: 0.001*sqrt(5.0**2+(15.0/(tree.p_lab[_k]*sin(tree.theta_lab[_k])**(3./2.)))**2)
    # CLIC 
    # d0_resolution = lambda _k: 0.001*sqrt(5.0**2+(15.0/(tree.p_lab[_k]*sin(tree.theta_lab[_k])**(3./2.)))**2)
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
            # just final hadrons
            if tree.catchall[k] != 25 or tree.isBCancestor[k] == 1:
                continue
            # only use the K-shorts if we want to
            if tree.isKshort[k] == 1 and useKshorts == -1:
                continue
            # only use the K+- if we want to
            if tree.isKshort[k] != 1 and useKshorts == 1:
                continue
            # get charge of the particle
            if tree.isKshort[k] == 1:
                particlecharge = 0;
            else:
                particlecharge = abs(tree.pdgId[k])/tree.pdgId[k];
            # parallel momentum with sign, and charge
            kaons_pm.append(  hadron(p=signed_pm(k),
                                charge=particlecharge,
                                d0 = smear(d0_f(k), sqrt(d0_resolution(k)**2 + IP_resolution**2)),
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
        u_h = filter(lambda x: x.hemisphere == 1, sorted_kaons)
        d_h = filter(lambda x: x.hemisphere == -1, sorted_kaons)
        # -- continue if none was found
        if len(u_h) < 1 or len(d_h) < 1:
            #leading_kaons[i] = ()
            continue
        # -- check charge if needed
        if applycharge:
            # double loop and get the higher of the two
            up_down_pair = get_opposite_charge(u_h,d_h)
            down_up_pair = get_opposite_charge(d_h,u_h)

            # check that get_opposite_charge does not return (-1, -1) for both
            # combinations (happens when no match is found). Then discard event,
            # otherwise pick the one with the largest summed momentum
            if  sum(up_down_pair) == sum(down_up_pair) < 0:
                continue
            elif  sum(up_down_pair) > sum(down_up_pair):
                leading_kaons[i] = up_down_pair
            else:
                leading_kaons[i] = down_up_pair
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
    containers = dict(map(lambda x: (x.GetName(),pythonize(x.GetClassName())),old_tree.GetListOfBranches()))
    # -- setting branch addresses
    _dumm = map(lambda (bname,vobject): tree.Branch(bname,vobject),sorted(containers.iteritems()))
    # -- and the new d0, z0 for hadrons
    pseudoimpactpar = { 'd0': ROOT.std.vector("float")(), 'z0' : ROOT.std.vector("float")(), 'R' : ROOT.std.vector("float")() }
    _dumm = map(lambda (bname,vobject): tree.Branch(bname,vobject),sorted(pseudoimpactpar.iteritems()))

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
        # and the new d0,z0 
        for (bname, ipvector) in pseudoimpactpar.iteritems():
            ipvector.clear()
            ipvector.reserve(2)
            _dumm = map(lambda _hadron: ipvector.push_back(getattr(_hadron,bname)),[h_l,h_sl])
        tree.Fill()
    # clear the vectors for filling the tree with empty events (to match number of generated events)
    _dummy = map(lambda _v: _v.clear(), containers.values())
    _dummy = map(lambda _v: _v.clear(), pseudoimpactpar.values())
    # Filling empty events
    _dummy = map(lambda k: tree.Fill(), xrange(old_tree.GetEntries()-len(hadronlist)))
    print
    f.Write()
    f.Close()
        

def main(args,suffixout,kshorts_considered,is_charge_considered,outfilename):
    """
    """
    import os

    if suffixout:
        suff = suffixout.replace('.','')
        globals()['SUFFIXPLOTS'] ='.'+suff
    if kshorts_considered == None:
        kshorts_considered = USEKSHORTS
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
        leading_kaons,nM = get_leading_kaons(t,is_charge_considered,kshorts_considered)
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
                       " between leading kaons. When applying, the requirement is that"\
                       " the absolute value of the summed charges is less than 2. In"\
                       " this way the possible presence of K_shorts is taken into account.")
    parser.add_option( '-k', '--kshorts', action='store', type='int', dest='usekshorts',\
                       help="-1: do not use K_shorts; 0: use K_shorts and K+-; +1 use only K_shorts"\
                       " and no K+-. [0]")
    
    (opt,args) = parser.parse_args()

    if len(args) < 1:
        message = "\033[31mprocesshzroot ERROR\033[m Missing input file(s)"
        raise RuntimeError(message)
    
    main(args,opt.suffixout,opt.usekshorts,(not opt.notcharge),opt.outfname)

