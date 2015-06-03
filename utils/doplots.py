#!/usr/bin/env python
""":script:`doplots` -- Plotting stuff created with processhzroot script
======================================================================

.. script:: doplots [OPTIONS]    
      :platform: Unix
      :synopsis: Do some plots of the objects previously obtained by the 
                 processhzroot script
	  .. moduleauthor:: Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
"""
#_COLOR = [ ROOT.kCyan+2, ROOT.kOrange+5,ROOT.kAzure-7,ROOT.kGreen+2,ROOT.kRed-2, ROOT.kBlue-3,
#        ROOT.kBlack, ROOT.kRed+4]
def getcolor():
    import ROOT
    return [ ROOT.kRed+4, ROOT.kAzure+3, ROOT.kOrange-2, ROOT.kGreen-5, ROOT.kYellow+2, \
        ROOT.kCyan-2, ROOT.kOrange+5,ROOT.kAzure-7,ROOT.kGreen-2,ROOT.kRed-4, ROOT.kGray-3 ]

SUFFIXPLOTS='.pdf'


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

# Static class to deal with inputs
class higgsinputs:
    # Data obtained from 
    # Branching ratios: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3
    # Cross section   : 
    def __init__(self):
        self.mH         = 125.7 # (GeV)
        self.eeToHat250 = 2.5e2 # (fb)
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
    global BRuu_ss
    global BRdd_ss 
    global BRcc_ss
    global BRbb_ss

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
    global hiInstance
    # Get the backgrounds
    #totalbkgdenominator = sum([ getbr(bkgname) for bkgname in bkg.keys()])
    bkgparsed = map(lambda x: parseprocess(x), bkg.keys())
    totalbkgdenominator = sum([ getattr(hiInstance,'br'+bkgname) for bkgname in bkgparsed])
    ilist = xrange(len(bkg.values()[0]))
    totalbkgeff = []
    for i in ilist:
        # Adding up all the background efficiency
        # for the momentum cut i-essim
        ibkg = 0
        for bkgname,bkgefflist in bkg.iteritems():
            #ibkg += getbr(bkgname)*bkgefflist[i]
            ibkg += getattr(hiInstance,'br'+parseprocess(bkgname))*bkgefflist[i]
        totalbkgeff.append( float(ibkg)/float(totalbkgdenominator) )
    return totalbkgeff 

def draweffpointsinsignificance(usefuldict,outname):
    """.. draweffpointsinsignificance(usefuldict,outname) -> graph,leg

    Draws a graph made of relevant working points in the significance
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
        y0 = 0.65
        y1 = 0.8
    else:
        y0 = 0.2
        y1 = 0.45
    # Prepare a new graph to be include in
    # the canvas c, (containing the significance curve)
    #_COLOR = [ ROOT.kCyan+2, ROOT.kOrange+5,ROOT.kAzure-7,ROOT.kGreen+2,ROOT.kRed-2, ROOT.kBlue-3,
    #        ROOT.kBlack, ROOT.kRed+4]
    COLOR=getcolor()
    textpos = {}
    _g = {}
    leg = getleg(x0=x0,y0=y0,x1=x1,y1=y1)
    leg.SetTextSize(0.035)
    for (k,(seff,(p,significance,sigeff,bkgeff))) in \
            enumerate(sorted(usefuldict.iteritems(),key=lambda (k,(p,sig,sigeff,bkgeff)): p)):
        _g[k] = ROOT.TGraph()
        _g[k].SetMarkerStyle(33)
        _g[k].SetMarkerSize(2)
        _g[k].SetMarkerColor(COLOR[k])
        _g[k].SetPoint(0,p,significance)
        effstr="%.1f" % (sigeff*100.)
        bkgstr="%.2f" % (bkgeff*100.)
        text = " #varepsilon^{i}_{sig}=%s%s, #varepsilon^{i}_{bkg}=%s%s" % (effstr,"%",bkgstr,"%")
        leg.AddEntry(_g[k],text,'P')
    return _g,leg

def draweffpointsinroc(usefuldict,outname):
    """.. draweffpointsinroc(usefuldict,outname) -> graph,leg

    Draws a graph made of relavant working points in the roc curve
    which relates with the used cut

    :param usefuldict. 

    :return: The graph and the legend
    :rtype:  (ROOT.TGraph,ROOT.TLegend)
    """
    import ROOT
    # Get the resonance
    x0 = 0.3
    x1 = 0.4
    y0 = 0.2
    y1 = 0.45
    # Prepare a new graph to be include in
    # the canvas c, (containing the significance curve)
    #_COLOR = [ ROOT.kCyan+2, ROOT.kOrange+5,ROOT.kAzure-7,ROOT.kGreen+2,ROOT.kRed-2, ROOT.kBlue-3,
    #        ROOT.kBlack, ROOT.kRed+4]
    COLOR=getcolor()

    textpos = {}
    _g = {}
    leg = getleg(x0=x0,y0=y0,x1=x1,y1=y1)
    leg.SetTextSize(0.035)
    for (k,(seff,(p,significance,sigeff,bkgeff))) in \
            enumerate(sorted(usefuldict.iteritems(),key=lambda (k,(p,sig,sigeff,bkgeff)): p)):
        _g[k] = ROOT.TGraph()
        _g[k].SetMarkerStyle(33)
        _g[k].SetMarkerSize(2)
        _g[k].SetMarkerColor(COLOR[k])
        _g[k].SetPoint(0,bkgeff,sigeff)
        text = " S/#sqrt{B}=%.1f (@ p=%.1f GeV)" % (significance,p)
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
        if a.outname.find('significance') != -1:
            ng,leg = draweffpointsinsignificance(a.addtext,a.outname)
        elif a.outname.find('ROC') != -1:
            ng,leg = draweffpointsinroc(a.addtext,a.outname)
        for g in ng.values():
            g.Draw("PSAME")
        leg.Draw()
    c.SaveAs(a.outname+'.pdf')
    c.SaveAs(a.outname+'.root')

def getlatextable(cutdict):
    """.. function:: getlatextable(cutdict) -> str(latex table)
    """
    # Set an order 
    _setorder = { 's': 1, 'u': 2, 'd': 3, 'c': 4, 'b': 5 }
    
    ncols = len(cutdict.values()[0])
    latex = '\\begin{tabular}{c '
    for i in xrange(ncols):
        latex += ' c '
    latex += '}\n'
    latex += ' & ' 
    for i in sorted(cutdict.values()[0].keys()):
        latex += ' $d_{0}$ < %.1f [mm] &' % i
    latex = latex[:-1]
    latex += '\\\\ \\hline\\hline\n'
    for process,cuts in sorted(cutdict.iteritems(),key=lambda (x,y): _setorder[x[0]]):
        pre_latexify = process.replace('bar','')
        latexify = '$'+pre_latexify[0]+'\\bar{'+pre_latexify[1]+'}$'+\
                pre_latexify[2:]+' &'
        latex += latexify
        for cut,val in sorted(cuts.iteritems()):
            latex += ' %.3f &' % val
        latex = latex[:-1]+'\\\\\n'
    latex += '\\hline\\hline\n'
    latex += '\\end{tabular}\n'

    return latex

def saveallinfo(d):
    """
    """
    import pickle

    output = open('allinfo.pkl','wb')
    pickle.dump(d,output)
    output.close()

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
    global hiInstance

    try:
        from PyAnUtils.plotstyles import squaredStyle,setpalette
        lstyle = squaredStyle()
        lstyle.cd()
        ROOT.gROOT.ForceStyle()
        ROOT.gStyle.SetOptStat(0)

        setpalette("darkbody")
    except ImportError:
        pass
    
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
        c.SaveAs(name.replace('_th2_','_')+'.pdf')
    
    # Prepare the efficiencies
    c,effdict = setupefficiencies(_obj,hadrons)

    # Plot efficiencies of different processes
    k = { 'H':0, 'Z':0 }
    leg = {}
    for _k in k.keys():
        leg[_k] = getleg()
    
    COLOR = getcolor()
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
    #-- End efficiency plots
    
    # Filling the graphs for the significance and ROC 
    WP = [20,40,60]
    SIGNIFICANCE = [1,3,5]
    TOLERANCEPERCENT = 0.02
    rocsgf = {}
    sigsgf = {}
    allinfodict = {}
    for res,sig in sigdict.iteritems():
        ## -- Get total background efficiency
        bkgeff = totalbkgeff[res]
        # -- Get significance
        # Graph for the resonance res
        sgf = ROOT.TGraph()
        roc = ROOT.TGraph()
        rocsgf[res] = {}
        sigsgf[res] = {}
        allinfodict[res] = []
        for i,efsignal in enumerate(sig):
            ibkgeff = bkgeff[i]
            try:
                #significance = float(efsignal)/sqrt(ibkgeff)
                ###correct = getbr('uubar')+getbr('ddbar')+getbr('ccbar')+getbr('bbbar')
                #m_c=1.275
                #m_s=0.095
                #correct *= (m_c/m_s)**2.0
                ###significance = float(efsignal)/(ibkgeff*correct)
                Nbkgevents = sum(map(lambda decay: hiInstance.getEvents(decay),
                    ['bbbar','ccbar','uubar','ddbar']))
                _B = Nbkgevents*ibkgeff
                Nsignalevents = hiInstance.getEvents('ssbar')
                _S = Nsignalevents*efsignal
                significance = float(_S)/sqrt(float(_B))

            except ZeroDivisionError:
                significance = 0.0
            sgf.SetPoint(i,p[i],significance)
            roc.SetPoint(i,ibkgeff,efsignal)
            # extra info
            allinfodict[res].append( (p[i],_S,_B,significance) )
            # Storing some working points  (respect the signal efficiency)
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
            #Doing the same for the significance
            try:
                wps = filter(lambda x: abs(significance-x) < x*1.,SIGNIFICANCE)[0]
                # Check if we already got it
                if sigsgf[res].has_key(wps):
                    # Only update it if is the lowest distance
                    if abs(wps-significance) < abs(wps-sigsgf[res][wps][1]):
                        sigsgf[res][wps]= (p[i],significance,efsignal,ibkgeff)
                else:
                    sigsgf[res][wps]= (p[i],significance,efsignal,ibkgeff)
            except IndexError:
                pass
        # Adding the significance working points to the roc dict
        rocsgf[res].update(sigsgf[res])

        # Plot the graphs:
        # -- significance
        if res == 'Z':
            xtitle = 'cut ( < #sqrt{p^{2}_{1}+p^{2}_{2}} ) [GeV]'
        else:
            xtitle = 'cut ( < #sqrt{p^{2}_{1||}+p^{2}_{2||}} ) [GeV]'
        #ytitle = '#varepsilon_{S}/#sqrt{#varepsilon_{B}}'
        #ytitle = '#varepsilon_{S}/#sum_{q}(m_{q}/m_{s}}^{2}#varepsilon_{B}'
        ytitle = 'S/#sqrt{B}'
        outname= res+'_significance'
        drawgraph(sgf,xtitle=xtitle,ytitle=ytitle,outname=outname,\
                addtext=rocsgf[res],opt='ALC')
        
        # -- roc
        outnameroc=res+'_ROC'
        drawgraph(roc,ytitle='#varepsilon_{S}',xtitle='#varepsilon_{B}',outname=outnameroc,\
                addtext=rocsgf[res],opt='ALC')
        
    # Finally, plotting TH1F (d0 of the leading hadrons)
    _pre_th1hists = filter(lambda ((x,classname),y): classname.find('TH1') == 0,\
                _obj.iteritems())
    th1hists = { 'H': map(lambda ((z,x),h): h, \
                  filter(lambda ((name,classname),h): name.find('H') == 0,\
                    _pre_th1hists)),
                  'Z': map(lambda ((z,x),h): h, \
                     filter(lambda ((name,x),y): name.find('Z') == 0, \
                     _pre_th1hists))
               }
    for res,histlist in th1hists.iteritems():
        c = ROOT.TCanvas()
        c.SetLogy()
        #c.SetLogx()
        leg = getleg(x0=0.6,x1=0.75)
        cutdict = {}
        for i,h in enumerate(sorted(histlist)):
            h.SetLineColor(COLOR[i])
            h.SetLineWidth(2)
            h.SetMarkerStyle(20)
            h.SetMarkerSize(0.5)
            h.SetMarkerColor(COLOR[i])
            h.SetNormFactor(1.0/h.Integral())
            if i == 0:
                h.GetYaxis().SetTitle('A.U./'+str(h.GetXaxis().GetBinWidth(1)))
                h.Draw("PE")
            else:
                h.Draw("PESAME") 
            name = h.GetName().replace(res+'_th1_hz','').replace('_kaons_d0',' ').replace('_','')
            leg.AddEntry(h,name,'PL')
            # Getting the values of some cuts
            _int  = h.Integral()
            cutdict[name] = {}
            for cut in [ 0.1,0.5,0.7,1.0]: 
                cutdict[name][cut] = h.Integral(1,h.FindBin(cut))/_int
        leg.Draw()
        c.SaveAs(res+'_d0.pdf')
        
        print "===  d0 cuts |%s| ==========================" % (res)
        print getlatextable(cutdict)
    
    allinfodict['HEADER'] = ('MOMENTUM-CUT','EFF_SIGNAL','EFF_BKG','SIGNIFICANCE')
    saveallinfo(allinfodict)




if __name__ == '__main__':
    from optparse import OptionParser,OptionGroup
    import os

    #Opciones de entrada
    parser = OptionParser()
    parser.set_defaults(inputfile='processed.root')    
    parser.add_option( '-i', '--inputfile', action='store', type='string', dest='inputfile',\
            help="input root filename [processed.root]")
    parser.add_option( '-s', '--suffix', action='store', type='string', dest='suffixout',\
            help="output suffix for the plots .pdf]")
    
    (opt,args) = parser.parse_args()

    if opt.suffixout:
        globals()['SUFFIXPLOTS'] ='.'+opt.suffixout

    plots(os.path.abspath(opt.inputfile))