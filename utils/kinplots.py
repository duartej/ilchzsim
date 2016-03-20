#!/usr/bin/env python
""":script:`kinplots` -- Plotting stuff created with processhzroot script
=========================================================================

.. script:: doplots [OPTIONS]    
      :platform: Unix
      :synopsis: Do some plots of the objects previously obtained by the 
                 processhzroot script
.. moduleauthor:: Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>

FIXME:: MISSING DOCUMENTATION
"""
import functools

@functools.total_ordering
class kaon(object):
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

SUFFIXPLOTS='.pdf'
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
        for k in xrange(tree.catchall.size()):
            if tree.catchall[k] != 25:
                continue
            # parallel momentum with sign, and charge
            ###kaons_pm.append((signed_pm(k),abs(tree.pdgId[k])/tree.pdgId[k],k))
            kaons_pm.append(  kaon(p=signed_pm(k),
                                charge=abs(tree.pdgId[k])/tree.pdgId[k],
                                d0 = d0_f(k),
                                z0 = z0_f(k),
                                L  = L_f(k),
                                R  = R_f(k),
                                index = k,
                                ) )
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
    print
    return leading_kaons

def init_tree(filename):
    """
    """
    from PyAnUtils.retrievetrees import plaintree
    # dummy class inheriting from stored tree

    return plaintree(filename,'mctrue')

def create_histos(suffix,description,hc=None):
    """Function gathering all the histograms we are interested
    to plot. So far the histograms defined are:
     * h2_pL_suffix : the parallel momentum of the leading 
                      and subleading hadrons
     * h2_d0_suffix : the impact parameter of the leading 
                      and subleading hadrons (extrapolated as
                      straight lines)
     * h2_Lxy_suffix: the vertex kaon production in the transverse
                      plane
     * h2_R_suffix  : the vertex kaon production

     * h_d0_suffix  : the impact parameter of the leading 
                      and subleading hadrons (extrapolated as
                      straight lines) in the same histogram (1D)


    Parameters
    ----------
    suffix: str
        used to distinguish the histograms of different samples
    description: str
        the legend to be used when several samples are plotted in 
        the same canvas
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

    COLOR = { 'ssbar': 46, 'bbbar': 12,
            'ccbar': 14, 'uubar': 16,
            'ddbar': 18}
    
    # -- create histo container, if there is no one
    if not hc:
        hc = HistoContainer()
    # -- some cosmethics
    # -- populate the hc with the histograms
    hc.create_and_book_histo("h2_pL_{0}".format(suffix),\
            "leading kaons parallel momentum",\
            100,0,65,npoints_y=100,ylow=0,yhigh=65,description=description,
            xtitle="leading-p_{||} [GeV]",ytitle='subleading-p_{||} [GeV]',
            color=COLOR[suffix])
    hc.create_and_book_histo("h_d0_{0}".format(suffix),\
            "leading kaons impact parameter",100,-5,5,\
            description=description,
            xtitle="d_{0} [mm]",
            ytitle="A.U.",
            color=COLOR[suffix])
    # TO BE DEPRECATED
    hc.create_and_book_histo("h2_d0_{0}".format(suffix),\
            "leading kaons impact parameter",\
            100,-5,5,npoints_y=100,ylow=-5,yhigh=5,description=description,
            xtitle="leading-kaon d_{0} [mm]",
            ytitle='subleading-kaon d_{0} [mm]',
            color=COLOR[suffix])
    
    hc.create_and_book_histo("h_z0_{0}".format(suffix),\
            "leading kaons impact parameter",100,-10,10,\
            description=description,
            xtitle="z_{0} [mm]",
            ytitle="A.U.",
            color=COLOR[suffix])
    # TO BE DEPRECATED
    hc.create_and_book_histo("h2_z0_{0}".format(suffix),\
            "leading kaons impact parameter",\
            100,-10,10,npoints_y=100,ylow=-10,yhigh=10,description=description,
            xtitle="leading-kaon z_{0} [mm]",
            ytitle='subleading-kaon z_{0} [mm]',
            color=COLOR[suffix])
    
    hc.create_and_book_histo("h_Lxy_{0}".format(suffix),\
            "leading kaons production vertex (transverse plane)",\
            100,0,5,
            description=description,
            xtitle="L_{xy} [mm]",ytitle='A.U.',
            color=COLOR[suffix])
    # TO BE DEPRECATED
    hc.create_and_book_histo("h2_Lxy_{0}".format(suffix),\
            "leading kaons production vertex (transverse plane)",\
            100,0,5,npoints_y=100,ylow=0,yhigh=5,description=description,
            xtitle="leading-kaon L_{xy} [mm]",
            ytitle='subleading-kaon L_{xy} [mm]',
            color=COLOR[suffix])
    
    hc.create_and_book_histo("h_R_{0}".format(suffix),\
            "leading kaons production vertex",\
            100,0,5,description=description,
            xtitle="R [mm]", ytitle='A.U.',
            color=COLOR[suffix])
    # TO BE DEPRECATED
    hc.create_and_book_histo("h2_R_{0}".format(suffix),\
            "leading kaons production vertex",\
            100,0,5,npoints_y=100,ylow=0,yhigh=5,description=description,
            xtitle="leading-kaon R [mm]",
            ytitle='subleading-kaon R [mm]',
            color=COLOR[suffix])

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
    c.SaveAs("{0}.{1}".format(histo.GetName(),'png'))

    c.Close()
    del c

def plot_combined(hc,varname,option=''):
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
    hc.plot(histonames[0],'comb_{0}.png'.format(varname),log=True)

def main(args,suffixout,hadrons,is_charge_considered):
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
                replace("hz","").replace("_PID_","").replace('kaons_only','')
        # create the histos
        hc = create_histos(sname,sname,hc)
        # -- initialize file
        t = init_tree(fname)
        # get the leading kaons dict { event#: ((up_pm,k),(down,k)), ... } 
        leading_kaons = get_leading_kaons(t,is_charge_considered)
        # put always higher pL in position 0
        ordered_lk = map(lambda x: sorted(x, reverse=True),leading_kaons.values())
        # and filling the histos
        for (x_h,x_l) in ordered_lk:
            _dummy = hc.fill("h2_pL_{0}".format(sname),x_h.p,x_l.p)
            _dummy = hc.fill("h2_d0_{0}".format(sname),x_h.d0,x_l.d0)
            _dummy = hc.fill("h_d0_{0}".format(sname),x_h.d0)
            _dummy = hc.fill("h_d0_{0}".format(sname),x_l.d0)
            _dummy = hc.fill("h2_z0_{0}".format(sname),x_h.z0,x_l.z0)
            _dummy = hc.fill("h_z0_{0}".format(sname),x_h.z0)
            _dummy = hc.fill("h_z0_{0}".format(sname),x_l.z0)
            _dummy = hc.fill("h2_Lxy_{0}".format(sname),x_h.L,x_l.L)
            _dummy = hc.fill("h_Lxy_{0}".format(sname),x_h.L)
            _dummy = hc.fill("h_Lxy_{0}".format(sname),x_l.L)
            _dummy = hc.fill("h2_R_{0}".format(sname),x_h.R,x_l.R)
            _dummy = hc.fill("h_R_{0}".format(sname),x_h.R)
            _dummy = hc.fill("h_R_{0}".format(sname),x_l.R)
    # plotting 
    # FIXME-- a function: plot those histos with a reg_expr 
    for k,h in filter(lambda (_k,_h): _k.find('h2_pL')==0,hc._histos.iteritems()):
        plot(h,k,option='COLZ')
    for k,h in filter(lambda (_k,_h): _k.find('h2_d0')==0,hc._histos.iteritems()):
        plot(h,k,option='COLZ')
    for k,h in filter(lambda (_k,_h): _k.find('h2_z0')==0,hc._histos.iteritems()):
        plot(h,k,option='COLZ')
    for k,h in filter(lambda (_k,_h): _k.find('h2_Lxy')==0,hc._histos.iteritems()):
        plot(h,k,option='COLZ')
    for k,h in filter(lambda (_k,_h): _k.find('h2_R')==0,hc._histos.iteritems()):
        plot(h,k,option='COLZ')
    # plot the combined histograms
    plot_combined(hc,'h_d0')
    plot_combined(hc,'h_z0')
    plot_combined(hc,'h_Lxy')
    plot_combined(hc,'h_R')

if __name__ == '__main__':
    from optparse import OptionParser,OptionGroup

    #Opciones de entrada
    usage = "usage: kinplots INPUTFILE1 [INPUTFILE2 ...] [options]"
    parser = OptionParser(usage=usage)
    parser.set_defaults(hadrons='kaons',notcharge=False)    
    parser.add_option( '-s', '--suffix', action='store', type='string', dest='suffixout',\
            help="output suffix for the plots [.pdf]")
    parser.add_option( '--no-charge', action='store_true',  dest='notcharge',\
            help="not applying the opposite charge requirement"\
            "between leading kaons")
    parser.add_option( '--hadron', action='store', type='string', dest='hadrons',\
            help="final state hadron type [kaons]")
    
    (opt,args) = parser.parse_args()

    if len(args) < 1:
        message = "\033[31mkinplots ERROR\033[m Missing input file(s)"
        raise RuntimeError(message)
    
    main(args,opt.suffixout,opt.hadrons,(not opt.notcharge))

