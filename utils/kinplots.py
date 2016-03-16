#!/usr/bin/env python
""":script:`kinplots` -- Plotting stuff created with processhzroot script
=========================================================================

.. script:: doplots [OPTIONS]    
      :platform: Unix
      :synopsis: Do some plots of the objects previously obtained by the 
                 processhzroot script
	  .. moduleauthor:: Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
"""

SUFFIXPLOTS='.pdf'
def get_leading_kaons(tree,applycharge):
    """
    """
    from math import cos
    import sys

    # auxiliar function to obtain the signed (+1 top hemisphere, 
    # -1 bottom hemispher) parallel momentum
    signed_pm = lambda _k: tree.p[_k]*cos(tree.theta[_k])
    nentries = tree.getentries()
    leading_kaons = {}    
    msg = "Evaluating %s..." % tree._rootfiles[0]
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
            kaons_pm.append((signed_pm(k),abs(tree.pdgId[k])/tree.pdgId[k],k))
        # sort in decrease order, the first is the highest p in
        # the top hemisphere, and last one, the highest p in the
        # bottom hemisphere
        sorted_kaons= sorted(kaons_pm,reverse=True)
        # separate between both hemispheres (using the 
        # sign of the parallel momentum)
        u_h = filter(lambda x: x[0] > 0,sorted_kaons)
        # Re-order again because we're recovering the positive
        # sign for the parallel momentum
        d_h = sorted(map(lambda y: (abs(y[0]),y[1],y[2]),\
                filter(lambda x: x[0] < 0,sorted_kaons)),reverse=True)
        # -- continue if none was found
        if len(u_h) < 1 or len(d_h) < 1:
            #leading_kaons[i] = ()
            continue
        # -- check charge if needed
        if applycharge:
            for ku in u_h:
                try:
                    opposite_down = filter(lambda kd: ku[1]*kd[1] < 1,d_h)[0]
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

def plots(filename,hadrontype,applycharge):
    """
    """
    from PyAnUtils.histocontainer import HistoContainer
    from PyAnUtils.plotstyles import squaredStyle,setpalette
    import ROOT
    import os

    lstyle = squaredStyle()
    lstyle.cd()
    ROOT.gROOT.ForceStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch()
    setpalette("inv_darkbody")

    # plotting machinery
    h = HistoContainer()
    h.create_and_book_histo("h2_pL","leading kaons parallel momentum",\
            50,0,65,npoints_y=50,ylow=0,yhigh=65)

    t = init_tree(filename)
    # get the leading kaons dict { event#: ((up_pm,k),(down,k)), ... } 
    leading_kaons = get_leading_kaons(t,applycharge)
    # put always higher pL in position 0
    ordered_lk = map(lambda x: sorted(x,key=lambda x: abs(x[0]),reverse=True),\
            leading_kaons.values())
    # And filling the histos
    _dummy = map(lambda (x_h,x_l): h.h2_pL.Fill(x_h[0],x_l[0]),ordered_lk)

    # plotting it
    h.h2_pL.GetXaxis().SetTitle("high-p_{||} [GeV]")
    h.h2_pL.GetYaxis().SetTitle("low-p_{||} [GeV]")

    c = ROOT.TCanvas()
    h.h2_pL.Draw("COLZ")
    c.SaveAs("{0}_pL.{1}".format(
        os.path.basename(filename).replace(".root",""),'png'))



if __name__ == '__main__':
    from optparse import OptionParser,OptionGroup
    import os

    #Opciones de entrada
    usage = "usage: kinplots INPUTFILE1 [INPUTFILE2 ...] [options]"
    parser = OptionParser(usage=usage)
    parser.set_defaults(hadrons='kaons',notcharge=False)    
    parser.add_option( '-s', '--suffix', action='store', type='string', dest='suffixout',\
            help="output suffix for the plots [.pdf]")
    parser.add_option( '--not-charge', action='store_true',  dest='notcharge',\
            help="not applying the opposite charge requirement"\
            "between leading kaons")
    parser.add_option( '--hadron', action='store', type='string', dest='hadrons',\
            help="final state hadron type [kaons]")
    
    (opt,args) = parser.parse_args()

    if len(args) < 1:
        message = "\033[31mkinplots ERROR\033[m Missing input file(s)"
        raise RuntimeError(message)

    if opt.suffixout:
        suff = opt.suffixout.replace('.','')
        globals()['SUFFIXPLOTS'] ='.'+suff
    absfilenames = map(lambda fname: os.path.abspath(fname),args)
    for fname in absfilenames:
        plots(fname,opt.hadrons,(not opt.notcharge))
