#!/usr/bin/env python
""":script:`doplots` -- Plotting stuff created with processhzroot script
======================================================================

.. script:: doplots [OPTIONS]    
      :platform: Unix
      :synopsis: Do some plots of the objects previously obtained by the 
                 processhzroot script
	  .. moduleauthor:: Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
"""
def setupefficiencies(rootobjdict,hadrons):
    """ ..function ::setupefficiencies(rootobjdict) -> {str: ROOT.TCanvas(), .. }

    extract the available efficiencies in the rootobjdict and set up one canvas
    for each resonance (Z,H) where the frame is used to set the axis titles 
    and so on.

    :param rootobjdict: dictionary with the key name,class name and object availables
                        in the open TFile
    :type  rootobjdict: { (str,str): ROOT.TObject, ... }

    :return: the canvas where to draw the efficiencies
    :rtype : {str: ROOT.TCanvas(), ... }
    """
    import ROOT
    # Get the TEfficiency objects available
    effpredict = dict(map(lambda ((x,y),z): (x,z),filter(lambda ((x,classname),y): \
            classname.find('TEfficiency') == 0, rootobjdict.iteritems())))

    # Just the selected hadrons
    effdict = dict(filter(lambda (x,y): x.find(hadrons) != -1,effpredict.iteritems()))
    
    # cosmethics
    canvasdict = { 'Z': ROOT.TCanvas(), 'H': ROOT.TCanvas() }
    for name,canvas in canvasdict.iteritems():
        canvas.cd()
        frame = canvas.DrawFrame(0,0,60,1)
        frame.GetXaxis().SetTitle("cut ( < #sqrt{p_{1}^{2}+p_{2}^{2}} ) [GeV]")
        frame.GetYaxis().SetTitle("#varepsilon")
        frame.Draw()

    return canvasdict,effdict

def ploteffs(canvas,eff,color):
    """ .. function:: ploteffs()
    """
    # Plotting eff
    canvas.cd()
    eff.SetLineColor(color)
    eff.SetLineWidth(2)
    eff.SetMarkerColor(color)
    eff.SetMarkerStyle(20)
    eff.SetMarkerSize(0.5)
    eff.SetFillColor(color)
    eff.Draw("CZSAME")

def getleg(**kwd):
    """.. function:: getleg(**kwd) -> ROOT.TLegend()
    
    Return a ROOT.TLegend object with some cosmethics already filled

    kwd accepts the following keys: x0, x1, y0, y1
    """
    import ROOT
    class coord:
        def __init__(self):
            self.x0=0.4
            self.x1=0.65
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

def geteffsignalandbkgs(effdict):
    """.. function:: geteffsignalandbkgs(effdict) -> { 'res': [eff1,eff2,..], .. },
                                                          { str: {str: [bkg1,..], ...}, ...}

    Extracts from the dictionary input a list of the efficiency
    per cut, for the signal and backgrounds involved. Returns
    signal and backgrounds
    """
    sigdict = {}
    bkgdict = {}
    for name,eff in effdict.iteritems():
        if name.find('Z') == 0:
            res = 'Z'
        elif name.find('H') == 0:
            res = 'H'
        # Storing the values for 
        effvallist =  [eff.GetEfficiency(i) \
                for i in xrange(1,eff.GetTotalHistogram().GetNbinsX()+1)]
        if name.find('ssbar') != -1:
            sigdict[res] = effvallist
        else:
            try:
                bkgdict[res][name] = effvallist
            except KeyError:
                bkgdict[res] = { name: effvallist}

    return sigdict,bkgdict

# Branching ratios with respect ccbar
BRuu_cc = 0.002
BRbb_cc = 3.32
BRss_cc = 0.075
BRdd_cc = 0.004

def getbr(name):
    """.. function:: getbr(name) -> br

    given the name of a sample returns the branching
    reation involved (with respect the ccbar)
    """
    global BRuu_cc
    global BRbb_cc
    global BRss_cc
    global BRdd_cc 

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

def gettotalbkgeff(bkg):
    """.. function:: gettotalbkgeff(bkg) -> [eff1,eff2,...]

    Given the available backgrounds (bkg), the total background efficiency
    is calculated as (using ccbar as reference br):
    
                        sum_{qqbar}(br_qqbar/br^_ccbar)*e^i_{qqbar} 
         e^i_totalB =-------------------------------------------------
                                sum_qqbar br_qqbar/br_ccbar

    where i-stands for the i-bin of momentum (radial cut in the p1,p2
    space); br-stands for branching ratio and e-stands for efficiency    

    :param bkg: dictionary with the efficiency of each background
                (the keys are the names of the background)
    :type  bkg: { str: list(float), .. }

    :return: list of the total background efficiency evaluated by bin
             of momentum (radial cut in the p1,p2 space)
    :rtype: list(float)
    """
    # Get the backgrounds
    totalbkgdenominator = sum([ getbr(bkgname) for bkgname in bkg.keys()])
    ilist = xrange(len(bkg.values()[0]))
    totalbkgeff = []
    for i in ilist:
        # Adding up all the background efficiency
        # for the momentum cut i-essim
        ibkg = 0
        for bkgname,bkgefflist in bkg.iteritems():
            ibkg += getbr(bkgname)*bkgefflist[i]
        totalbkgeff.append( float(ibkg)/float(totalbkgdenominator) )
    return totalbkgeff 

def draweffpointsin(usefuldict,outname):
    """.. draweffpointsin(usefuldict,outname) -> graph,leg

    Draws a graph made of relavant working points in the significance
    which relates with the signal and backgrounds efficiency

    :param usefuldict. 

    :return: The graph and the legend
    :rtype:  (ROOT.TGraph,ROOT.TLegend)
    """
    import ROOT
    # Get the resonance
    x0 = 0.2
    x1 = 0.3
    if outname.find('H') == 0:
        y0 = 0.45
        y1 = 0.7
    else:
        y0 = 0.2
        y1 = 0.55
    # Prepare a new graph to be include in
    # the canvas c, (containing the significance curve)
    _COLOR = [ ROOT.kCyan+2, ROOT.kOrange+5,ROOT.kAzure-7,ROOT.kGreen+2,ROOT.kRed-2, ROOT.kBlue-3 ]
    textpos = {}
    _g = {}
    leg = getleg(x0=x0,y0=y0,x1=x1,y1=y1)
    leg.SetTextSize(0.03)
    for (k,(seff,(p,significance,sigeff,bkgeff))) in \
            enumerate(sorted(usefuldict.iteritems(),key=lambda (k,(p,sig,sigeff,bkgeff)): p)):
        _g[k] = ROOT.TGraph()
        _g[k].SetMarkerStyle(33)
        _g[k].SetMarkerSize(2)
        _g[k].SetMarkerColor(_COLOR[k])
        _g[k].SetPoint(0,p,significance)
        effstr="%.1f" % (sigeff*100.)
        bkgstr="%.1f" % ((bkgeff**2.0)*100.)
        text = " #varepsilon^{i}_{sig}=%s%s, #varepsilon^{i}_{bkg}=%s%s" % (effstr,"%",bkgstr,"%")
        leg.AddEntry(_g[k],text,'P')
    return _g,leg

def drawgraph(g,**kwd):
    """.. function:: drawgraph(g[,**kwd) 

    Draw a graph using some pre-defined styles.
    The x-y titles, name of the outputfile, etc..
    are got them from the kwd
    """
    import ROOT
    class graphattr:
        def __init__(self):
            self.xtitle = ''
            self.yttile = ''
            self.outname= 'outputplot'
            self.opt    = ''
            self.addtext= None
            self.log    = False
    a = graphattr()
    
    for key,value in kwd.iteritems():
        setattr(a,key,value)

    c = ROOT.TCanvas()
    if a.log:
        c.SetLogy()
    g.SetLineWidth(2)
    g.Draw(a.opt)
    
    _h = g.GetHistogram()
    _h.GetXaxis().SetTitle(a.xtitle)
    _h.GetYaxis().SetTitle(a.ytitle)
    if a.addtext:
        ng,leg = draweffpointsin(a.addtext,a.outname)
        for g in ng.values():
            g.Draw("PSAME")
        leg.Draw()
    c.SaveAs(a.outname+'.pdf')
    c.SaveAs(a.outname+'.png')



def plots(rootfile,hadrons='kaons'):
    """.. function:: plots(rootfile[,hadrons='kaons'])

    Main function gathering all the plots to be performed:
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
    global BRuu_cc
    global BRbb_cc

    try:
        from PyAnUtils.plotstyles import squaredStyle,setpalette
        lstyle = squaredStyle()
        lstyle.cd()
        ROOT.gROOT.ForceStyle()
        ROOT.gStyle.SetOptStat(0)

        setpalette("darkbody")
    except ImportError:
        pass
    # Pseudo-global
    COLOR = [ ROOT.kRed+4, ROOT.kAzure+3, ROOT.kOrange-2, ROOT.kGreen-5, ROOT.kBlue+5, \
            ROOT.kCyan-2, ROOT.kOrange+5,ROOT.kAzure-7,ROOT.kGreen-2,ROOT.kRed-4, ROOT.kGray-3 ]
    
    # Get the root file with the TH2 and TEfficiency objects
    f = ROOT.TFile(rootfile)
    _obj = dict(((x.GetName(),x.GetClassName()), \
            f.Get(x.GetName())) for x in f.GetListOfKeys())
    # Plotting
    ROOT.gROOT.SetBatch()

    # Plotting TH2F (p1 vs. p2 of the leading hadrons)
    for ((name,k),h) in filter(lambda ((x,classname),y): \
            classname.find('TH2') == 0, _obj.iteritems()):
        c = ROOT.TCanvas()
        h.Draw("COLZ")    
        c.SaveAs(name.replace('_th2f_','')+'.pdf')
        c.SaveAs(name.replace('_th2f_','')+'.png')

    # Prepare the efficiencies
    c,effdict = setupefficiencies(_obj,hadrons)

    # Plot efficiencies of different processes
    k = { 'H':0, 'Z':0 }
    leg = {}
    for _k in k.keys():
        leg[_k] = getleg()
    
    for name,eff in sorted(effdict.iteritems()):
        if name.find('Z') == 0:
            res = 'Z'
        elif name.find('H') == 0:
            res = 'H'
        # Plotting efficiency of the given process
        entry = name.replace(res+'_eff_','').replace('hz','').replace('_',' ')
        leg[res].AddEntry(eff,entry,"PL")
        ploteffs(c[res],eff,COLOR[k[res]])
        k[res] += 1
    # Get the total background efficiency to include it
    # in the efficiency plots
    # For this I need to extract signal and backgrounds dictionaries
    # which are obtained when do the significance calculation 
    sigdict,bkgdict = geteffsignalandbkgs(effdict)
        
    totalbkgeff= dict(map(lambda (res,bkgs): (res,gettotalbkgeff(bkgs)),bkgdict.iteritems()))
    # Getting the momentum involved
    p = [effdict.values()[0].GetTotalHistogram().GetBinCenter(i) \
            for i in xrange(1,effdict.values()[0].GetTotalHistogram().GetNbinsX()+1)]
    
    # Include the total background efficiency
    gnew = {}
    for res,bkgefflist in totalbkgeff.iteritems():
        gnew[res] = ROOT.TGraph()
        for i,val in enumerate(bkgefflist):
            gnew[res].SetPoint(i,p[i],val)
        leg[res].AddEntry(gnew[res],"total backgrounds","PL")
        ploteffs(c[res],gnew[res],ROOT.kRed)
        k[res]+=1
    
    for name,canvas in c.iteritems():
        canvas.cd()
        leg[name].Draw()
        canvas.SaveAs('effcmp_'+name+'.pdf')
        canvas.SaveAs('effcmp_'+name+'.png')
    #-- End efficiency plots
    
    # Filling the graphs for the significance and ROC 
    WP = [1,5,20,40,60,80]
    TOLERANCEPERCENT = 0.2
    rocsgf = {}
    for res,sig in sigdict.iteritems():
        ## -- Get total background efficiency
        bkgeff = totalbkgeff[res]
        # -- Get significance
        # Graph for the resonance res
        sgf = ROOT.TGraph()
        roc = ROOT.TGraph()
        rocsgf[res] = {}
        for i,efsignal in enumerate(sig):
            ibkgeff = bkgeff[i]
            try:
                significance = float(efsignal)/sqrt(ibkgeff)
            except ZeroDivisionError:
                significance = 0.0
            sgf.SetPoint(i,p[i],significance)
            roc.SetPoint(i,efsignal,ibkgeff)
            # Storing some working points  (respecting the signal efficiency)
            try: 
                wp = filter(lambda x: abs(efsignal*100.0-x) < x*TOLERANCEPERCENT,WP)[0]
                # Check if we already got it
                if rocsgf[res].has_key(wp):
                    # Only update it if is the lowest distance
                    if abs(wp-efsignal) < abs(wp-rocsgf[res][wp][2]):
                        rocsgf[res][wp]= (p[i],significance,efsignal,ibkgeff)
                else:
                    rocsgf[res][wp]= (p[i],significance,efsignal,ibkgeff)
            except IndexError:
                pass

        # Plot the graphs:
        # -- significance
        xtitle = 'cut ( < #sqrt{p^{2}_{1||}+p^{2}_{2||}} ) [GeV]'
        ytitle = '#varepsilon_{S}/#sqrt{#varepsilon_{B}}'
        outname= res+'_significance'
        drawgraph(sgf,xtitle=xtitle,ytitle=ytitle,outname=outname,\
                addtext=rocsgf[res],log=True,opt='ALC')
        
        # -- roc
        drawgraph(roc,xtitle='#varepsilon_{S}',ytitle='#varepsilon_{B}',outname=res+'_ROC',opt='ALC')
        

if __name__ == '__main__':
    from optparse import OptionParser,OptionGroup
    import os
    
    #Opciones de entrada
    parser = OptionParser()
    parser.set_defaults(inputfile='processed.root')    
    parser.add_option( '-i', '--inputfile', action='store', type='string', dest='inputfile',\
            help="input root filename [processed.root]")
    
    (opt,args) = parser.parse_args()

    plots(os.path.abspath(opt.inputfile))
