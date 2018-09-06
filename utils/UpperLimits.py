#!/usr/bin/env python3
"""script to obtain the realistic significance and upper limits for different collider scenarios
   =========================

   usage: 
"""
# Final states in the Higgs channel
HiggsChannels=['CC', 'NC', 'NN']

# Higgs branching ratios
Hbbbar = 5.66e-1
Hccbar = 2.85e-2
Hssbar = 2.41e-4
Hgg    = 8.50e-2
Huubar = 1.3e-7
Hddbar = 5.8e-8

HBR={}
HBR['bb']=Hbbbar
HBR['cc']=Hccbar
HBR['ss']=Hssbar
HBR['gg']=Hgg
HBR['uu']=Huubar
HBR['dd']=Hddbar

#Z branching ratios into hadrons
Zbbbar = 0.1512
Zccbar = 0.1203
Zssbar = (3*0.156 - Zbbbar)/2
Zuubar = 2 * 0.116 - Zccbar
Zddbar = (3*0.156 - Zbbbar)/2

total = Zbbbar + Zccbar + Zssbar + Zuubar + Zddbar
Zbbbar = Zbbbar/total
Zccbar = Zccbar/total
Zssbar = Zssbar/total
Zuubar = Zuubar/total
Zddbar = Zddbar/total

    
class progressbar(object):
    def __init__(self, nsteps):
        import sys
        self.nsteps=nsteps
        self.counter=0
        self.progressstring='[{}:{}{}] {}%'.format('{',self.nsteps,'}','{}')
        sys.stdout.write('\r')
        sys.stdout.write(self.progressstring.format('='*0, 0))
        sys.stdout.flush()

    def next(self):
        import sys
        self.counter+=1
        sys.stdout.write('\r')
        sys.stdout.write(self.progressstring.format(
            '='*self.counter, int(round(100./self.nsteps*self.counter))))
        sys.stdout.flush()

    def stop(self):
        import sys
        sys.stdout.write('\n')

def difflist(l,k, changedindex):
    """ check how many different items a list has. Used for checking how different a parameter point is
    """
    v=0
    for i in range(len(l)):
        if i==3: #index 2 and 3 are correlated
            continue
        if l[i]!=k[i]:
            v+=1
            if i==changedindex: # prefer points that changed in the same index as before
                v-=0.1
    return v

class parameterspace(object):
    """ administer a list of all points of the rectangular parameterspace.
    """
    def __init__(self, listofpoints):
        try:
            dummy=listofpoints[0][0]
        except:
            print("no parameter points provided")

        self.parameterpoints=listofpoints
        self.dimension = len(self.parameterpoints[0])
        # get the possible values for all parameters"
        self.values=[sorted(list(set([v[i] for v in self.parameterpoints]))) for i in range(self.dimension)]

    def getneighbors(self, point, depth, oldWP, maxdiff=2):
        """get the neighboring parameter points of the input points in the multidimensional space. See
        below for a sketch in 2D. 'p' is the input point and 'r' the neighbors that are
        returned. oldWP is used to first return neighbours which have a change in the same index as
        in the step from oldWP to point
        -------
        --rrr--
        --rpr--
        --rrr--
        -------

        """
        import itertools

        # get the (first) index that changed from going from oldWP to point
        changedindex=-1
        for i in range(len(point)):
            if point[i]!=oldWP[i]:
                changedindex=i
                break
        
        # return all parameter points if required
        if depth == -1:
            neighbours= [p for p in self.parameterpoints]
            neighbours.remove(point)
            neighbours= list(filter(lambda p: difflist(point, p, changedindex) <= maxdiff, neighbours))
            return sorted(neighbours, key=lambda p: difflist(point, p, changedindex))
        
        # get indexes of the variables of point
        indexlist = list(map(lambda i: self.values[i].index(point[i]), range(self.dimension)))

        # get the corresponding values (+-1 if possible)
        valuelist = []
        for i in range(self.dimension):
            valuelist.append([])
            valuelist[-1].append(self.values[i][indexlist[i]])
            for d in range(1, depth+1):
                if indexlist[i]-d >= 0:
                    valuelist[-1].append(self.values[i][indexlist[i]-d])
                if indexlist[i]+d < len(self.values[i]):
                    valuelist[-1].append(self.values[i][indexlist[i]+d])

        neighbours = []
        for element in itertools.product(*valuelist):
            if list(element) in self.parameterpoints:
                neighbours.append(list(element))
        neighbours.remove(point)
        neighbours= list(filter(lambda p: difflist(point, p, changedindex) <= maxdiff, neighbours))
        return sorted(neighbours, key=lambda p: difflist(point, p, changedindex))
        
def Possible_SubAnalyses(scenario):
    """ Return the possible analysis channels for a given collider scenario
    """

    if scenario == 'CLIC350':
        return ['Ztoinv', 'Ztohad']
    elif scenario in ['ILC250', 'ILC350',]:
        return ['inv', 'had', 'electron', 'muon']
    elif scenario in ['FCCee']:
        return ['inv', 'had', 'electron', 'muon']
    else:
        return [None]


class colliderscenarios(object):
    def __init__(self, collider, analysischannel=None):
        if collider=='CLIC350':
            # Number of events from 1608.07538
            if analysischannel == 'Ztoinv':
                self.analysischannel = 'Ztoinv'
                self.NHiggs = (8000./Hbbbar + 372./Hccbar)/2
                self.NnHiggs= 2100+2090+104+30+1230
            elif analysischannel == 'Ztohad':
                self.analysischannel = 'Ztohad'
                self.NHiggs = (11100./Hbbbar+ 434./Hccbar)/2
                self.NnHiggs= 60+89+9990+11400
            else:
                raise ValueError('For CLIC350 the analysischannel needs to be specified')

        elif collider=='CLIC1400':
            # Number of events from 1608.07538
            self.analysischannel = 'Hnunu'
            self.NHiggs = (65400./Hbbbar + 3790./Hccbar)/2
            self.NnHiggs= 18500+23600+18500+170000+22200

        elif collider=='CLIC3000':
            # Number of events from 1608.07538
            self.analysischannel = 'Hnunu'
            self.NHiggs = (120000./Hbbbar + 6380./Hccbar)/2
            self.NnHiggs= 47400+52200+118000+394000+207000

        elif collider == 'ILC250':
            # numbers from 1207.0300
            if analysischannel == 'inv':
                self.analysischannel = 'inv'
                self.NHiggs  = 6293
                self.NnHiggs = 10940
            elif analysischannel == 'had':
                self.analysischannel = 'had'
                self.NHiggs  = 13726
                self.NnHiggs = 166807
            elif analysischannel == 'electron':
                self.analysischannel = 'electron'
                self.NHiggs  = 1184
                self.NnHiggs = 1607
            elif analysischannel == 'muon':
                self.analysischannel = 'muon'
                self.NHiggs  = 1365
                self.NnHiggs = 983

        elif collider == 'ILC350':
            if analysischannel == 'inv':
                self.analysischannel = 'inv'
                self.NHiggs  = 9543
                self.NnHiggs = 11092
            elif analysischannel == 'had':
                self.analysischannel = 'had'
                self.NHiggs  = 8686
                self.NnHiggs = 25393
            elif analysischannel == 'electron':
                self.analysischannel = 'electron'
                self.NHiggs  = 567
                self.NnHiggs = 590
            elif analysischannel == 'muon':
                self.analysischannel = 'muon'
                self.NHiggs  = 638
                self.NnHiggs = 465
        elif collider == 'CEPC':
            if analysischannel == 'inv':
                self.analysischannel = 'inv'
                self.NHiggs =(69820./Hbbbar+3029./Hccbar + 9522./Hgg)/3
                self.NnHiggs=16031.
        elif collider == 'FCCee':
            # numbers from 1207.0300
            lumiratio=0.95*10**7/77921. # ratio of #events FCCee/ILC
            ILC250=colliderscenarios('ILC250', analysischannel)
            self.NHiggs = ILC250.NHiggs*lumiratio
            self.NnHiggs= ILC250.NnHiggs*lumiratio
            #
            # Nleptons = 0  # required number of leptons
            # MrecoilPMin = 80
            # MrecoilPMax = 100
            # MinvMin = 85  # minimal invariant missing mass
            # HiggsPmin = 55 # minimal pT of the Higgs candidate
            # HiggsPmax = 65 # maximal pT of the Higgs candidate
            # MdijetMin = 95 # minimal dijet mass
            # MdijetMinMax = 120 # maximal value of the minimal dijet mass
            # MdijetDelta = 5 # step size for the variation of the minimal Mdijet value
            # MdijetMax = 127 # maximal dijet mass
            self.roclist=[[1156501,1115275,1057146,969313,835254,601540],
                          [2643139,2366095,2099555,1755290,1326952,817905]]
            # Nleptons = 0  # required number of leptons
            # MrecoilPMin = 80
            # MrecoilPMax = 100
            # MinvMin = 87  # minimal invariant missing mass
            # MinvMax = 110 # maximal invariant missing mass
            # HiggsPLmax = 45 # maximal pL of the Higgs candidate
            # Y12min = 0.3   # maximal value of the 1->2 jet splitting
            # MdijetMin = 111 # minimal dijet mass
            # MdijetMinMax = 123 # maximal value of the minimal dijet mass
            # MdijetDelta = 3 # step size for the variation of the minimal Mdijet value
            # MdijetMax = 126 # maximal dijet mass
            self.roclist[0] += [751280, 693966, 598959, 468650, 284505]
            self.roclist[1] += [1044710, 900730, 679656, 451379, 219288]
            Nleptons = 0  # required number of leptons

            # MrecoilPMin = 85
            # MrecoilPMax = 95
            # MinvMin = 87  # minimal invariant missing mass
            # MinvMax = 110 # maximal invariant missing mass
            # HiggsPLmax = 40 # maximal pL of the Higgs candidate
            # Y12min = 0.3   # maximal value of the 1->2 jet splittin
            # Y12max = 0.8   # maximal value of the 1->2 jet splitting
            # MdijetMin = 111 # minimal dijet mass
            # MdijetMinMax = 123 # maximal value of the minimal dijet mass
            # MdijetDelta = 3 # step size for the variation of the minimal Mdijet value
            # MdijetMax = 126 # maximal dijet mass
            self.roclist[0]+= [538111, 501478, 449817, 374512, 238360]
            self.roclist[1] += [476490, 405028, 312735, 213530, 102018]

            # Nleptons = 0  # required number of leptons
            # MinvMin = 87  # minimal invariant missing mass
            # Y12min = 0.2   # maximal value of the 1->2 jet splittin
            # Y12max = 0.7   # maximal value of the 1->2 jet splitting
            # chisquaredMax = 7
            # chisquaredDelta = 1
            self.chi2list=[[17120, 40492, 63836, 87343, 113813, 140037, 165800],
                           [4098, 15454, 31106, 48505, 73069, 95177, 122374]]

            # additionally promote from chi2 to explicit cut
            # HiggsPLmax = 60
            #self.chi2list[0] += [137673, 222272, 288310, 344129, 391687, 433294, 471557]
            #self.chi2list[1] += [62940, 147056, 233462, 322121, 419652, 515762, 605948]

            # Nleptons = 0  # required number of leptons
            # MinvMin = 87  # minimal invariant missing mass
            # Y12min = 0.4   # maximal value of the 1->2 jet splittin
            # Y12max = 0.8   # maximal value of the 1->2 jet splitting
            # HiggsPLmax = 60
            # chisquaredMax = 7
            # chisquaredDelta = 1
            self.chi2list[0] += [185584, 282521, 355326, 414569, 463568, 504549, 538791]
            self.chi2list[1] += [87539, 191665, 296176, 403129, 517304, 622066, 719059]

        elif collider=='FCCeePerfect':
            self.analysischannel = 'HZ'
            self.NHiggs = 10**7
            self.NnHiggs= 0

        
                
def Expected_UpperLimit(SB, CL=0.95):
    """ Return the expected upper limit on the signal strength
    
    Parameters
    ----------
    SB=[signal, background]
    signal: float
            number of expected signal events
    background: float
            number of expected background events
    CL:     float
            required confidence level. Default: 95%
    """
    from scipy.stats import chi2

    [signal, background] = SB
    df = 2* (signal+background+1)
    absUpperLimit=0.5*chi2.ppf(1-(1-CL)*(1-chi2.cdf(2*background,df)),df)-background
    return absUpperLimit/signal

def Expected_Significance(SB):
    """ Return the expected significance of the signal
        Eq. 97 of 1007.1727

    Parameters
    ----------
    SB=[signal, background]
    signal: float
            number of expected signal events
    background: float
            number of expected background events
    """

    import numpy as np
    [signal, background] = SB
    return np.sqrt(2*((signal+background) * np.log(1+signal/background) - signal) )

# for color printing
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def Print_Fail(message):
    print (bcolors.FAIL + bcolors.BOLD + 'WARNING: ' + message + bcolors.ENDC)

def Print_Warning(message):
    print (bcolors.WARNING + bcolors.BOLD + 'WARNING: ' + message + bcolors.ENDC)

COLORS=['crimson', 'blue', 'forestgreen', 'darkorange', 'skyblue', 'm', 'darkgrey']

def LogPlot(x, ys, xlabel, ylabel, plotname):
    """ make a logplot 

    parameters:
    -----------
    x: [floats]
           the x values of all points
    ys [lists of floats]
           a list containing lists of the y values. the order is bb, cc, ss, uu, dd, gg, WW, (non-Higgs)
    xlabel: str
           the label of the x axis
    ylabel: str
           the label of the y axis
    plotname: str
           the file name under which the plot is saved.
    """
    
    import matplotlib.pyplot as plt

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=13)
    
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)

    if len(ys)==7:
        fs=[r'$b\bar{b}$',r'$c\bar{c}$',r'$s\bar{s}$',r'$u\bar{u}$',r'$d\bar{d}$',r'$gg$', r'$WW^*$']
    elif len(ys)==8:
        fs=[r'$b\bar{b}$',r'$c\bar{c}$',r'$s\bar{s}$',r'$u\bar{u}$',r'$d\bar{d}$',r'$gg$',
            r'$WW^*$', r'non-Higgs']
    else:
        fs=map(lambda x: 'n.a.', ys)
        Print_Warning('in plot {0}, unknown number ({1}) of entries'.format(ylabel), len(ys))

    for (plotcolor, plotlabel, yvalue) in zip(COLORS[0:len(ys)], fs, ys):
        plt.plot(x, yvalue, linewidth=2, linestyle='-', color=plotcolor, label=plotlabel)

    ax.set_xlim(min(x), max(x))

    # remove 0 out of the list before determining the y range and remove empty lists
    dummy=map(lambda sublist: [x for x in sublist if x != 0], ys)
    dummy=filter(None, dummy)
    ymin=0.5*min(map(lambda x: min(x), dummy))
    ymax=1.3*max(map(lambda x: max(x), ys))
    ax.set_ylim(ymin, ymax)
    ax.set_yscale('log')
    ax.tick_params(direction='in', top=True, right=True)

    plt.xlabel(r'{0}'.format(xlabel))
    plt.ylabel(r'{0}'.format(ylabel))
    ax.legend(loc=0, frameon=False)

    plt.savefig(plotname.replace('_None',''))
    plt.close()

    
def LinPlot(x, ys, xlabel, ylabel, plotlabels, plotname):
    """ make a plot 

    parameters:
    -----------
    x: [floats]
           the x values of all points
    ys [lists of floats]
           a list containing lists of the y values. the order is bb, cc, ss, uu, dd, gg, WW, (non-Higgs)
    xlabel: str
           the label of the x axis
    ylabel: str
           the label of the y axis
    plotlabels: [str]
           the labels for the different ys[i] (in the same order)
    plotname: str
           the file name under which the plot is saved.
    """
    import matplotlib.pyplot as plt

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=13)
    
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)

    for (plotcolor, plotlabel, yvalue) in zip(COLORS[0:len(ys)], plotlabels, ys):
        plt.plot(x, yvalue, linewidth=2, linestyle='-', color=plotcolor, label=plotlabel)

    ax.set_xlim(min(x), max(x))

    ymin=0
    ymax=1.3*max(map(lambda x: max(x), ys))
    ax.set_ylim(ymin, ymax)
    ax.tick_params(direction='in', top=True, right=True)

    plt.xlabel(r'{0}'.format(xlabel))
    plt.ylabel(r'{0}'.format(ylabel))
    ax.legend(loc=0, frameon=False)

    plt.savefig(plotname.replace('_None',''))
    plt.close()

def Plot2D(X, Y, Z, xlabel, ylabel, zlabel, plotname, loglog=False):
    """ contour plot

    parameters:
    -----------
    X, Y, Z: mesh grids
           contain the x, y, and z values
    xlabel: str
           the label of the x axis
    ylabel: str
           the label of the y axis
    zlabel: str
           the label of the z axis
    plotname: str
           the file name under which the plot is saved.
    loglog: Boolean
           if it is a loglog plot or not
    """

    import matplotlib.pyplot as plt
    from matplotlib import colors
    import matplotlib as mpl

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=15)
    
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.tick_params(top=True, right=True)
    ax.set_ylim(bottom=100, top=10000000)

    if 'BestMu' in plotname:
        levels=[1,2,5,10,20,50,100,200,500]
        contour=plt.contourf(X, Y, Z, cmap=plt.cm.viridis, levels=levels, norm=colors.LogNorm(), rasterized=True)
        #contour.set_edgecolor('face')
        #contour=plt.pcolormesh(X, Y, Z, cmap=plt.cm.viridis, norm=colors.LogNorm())
        for c in contour.collections:
            c.set_edgecolor("face")
        cbar=plt.colorbar(ticks=levels)
        cbar.set_ticklabels(levels)
    elif 'BestpL' in plotname:
        levels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
        cmap = plt.cm.get_cmap('viridis', len(levels)-1)
        contour = plt.pcolormesh(X, Y, Z,
                                 cmap=cmap,
                                 norm=mpl.colors.BoundaryNorm(levels, ncolors=len(levels)-1, clip=False),
                                 rasterized=True)
        shiftedlevels= list(map(lambda x: x+0.5, levels))
        cbar=plt.colorbar(ticks=shiftedlevels)
        cbar.set_ticklabels(levels)
    elif 'Bestd0' in plotname:
        levels = [0.014, 0.015, 0.016, 0.017, 0.018, 0.019,
                  0.02, 0.021, 0.022, 0.023, 0.024, 0.025]
        levelsmu = list(map(lambda x: int(x*1000), levels))
        cmap = plt.cm.get_cmap('viridis', len(levels)-1)
        contour = plt.pcolormesh(X, Y, Z,
                                 cmap=cmap,
                                 norm=mpl.colors.BoundaryNorm(levels, ncolors=len(levels)-1, clip=False),
                                 rasterized=True)
        shiftedlevels= list(map(lambda x: x+0.0005, levels))
        cbar=plt.colorbar(ticks=shiftedlevels)
        cbar.set_ticklabels(levelsmu)
    elif 'BestPID' in plotname:
        levels = [0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01]
        levelsmu = list(map(lambda x: int(x*1000), levels))
        cmap = plt.cm.get_cmap('viridis', len(levels)-1)
        contour = plt.pcolormesh(X, Y, Z,
                                 cmap=cmap,
                                 norm=mpl.colors.BoundaryNorm(levels, ncolors=len(levels)-1, clip=False),
                                 rasterized=True)
        shiftedlevels= list(map(lambda x: x+0.005, levels))
        cbar=plt.colorbar(ticks=shiftedlevels)
        cbar.set_ticklabels(levels)
    elif 'BestEff' in plotname:
        contour=plt.contourf(X, Y, Z, cmap=plt.cm.viridis)
        cbar=plt.colorbar(contour)
    elif 'Best' in plotname:
        levels = list(set([item for sublist in Z for item in sublist]))
        levels.sort()
        contour=plt.contourf(X, Y, Z, cmap=plt.cm.viridis, levels=levels)
        cbar=plt.colorbar(ticks=levels)
        cbar.set_ticklabels(levels)
    else:
        contour=plt.contourf(X, Y, Z, cmap=plt.cm.viridis)    
        cbar = plt.colorbar(contour)
    if 'Best' in plotname:
        FCCee=colliderscenarios('FCCee', 'inv')
        ILC=colliderscenarios('ILC250', 'inv')
        CEPC=colliderscenarios('CEPC', 'inv')
        plt.plot([100,10**7],[100*ILC.NnHiggs/ILC.NHiggs,10**7*ILC.NnHiggs/ILC.NHiggs],'--k', label='Cut\&Count')
        plt.plot([100,10**7],[100*CEPC.NnHiggs/CEPC.NHiggs,10**7*CEPC.NnHiggs/CEPC.NHiggs],':k', label='BDT')

        plt.plot(ILC.NHiggs*200, ILC.NnHiggs*200, 'r+', markersize=13, label='$\mathcal{L}=50\,$ab$^{-1}$')
        plt.plot(ILC.NHiggs*20, ILC.NnHiggs*20, 'mP', markersize=13, label='$\mathcal{L}=5\,$ab$^{-1}$')
        plt.plot(ILC.NHiggs, ILC.NnHiggs, 'b*', markersize=13, label='$\mathcal{L}=250\,$fb$^{-1}$')

        plt.plot(CEPC.NHiggs*10, CEPC.NnHiggs*10, 'r+', markersize=13)
        plt.plot(CEPC.NHiggs, CEPC.NnHiggs, 'mP', markersize=13)
        plt.plot(CEPC.NHiggs/20., CEPC.NnHiggs/20., 'b*', markersize=13)

        #plt.plot(8000./Hbbbar, 2110+2090+104+30+1230, 'ro', markersize=13, label='CLIC')
        plt.legend(loc=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    cbar.ax.set_ylabel(zlabel)
    if loglog:
        ax.set_yscale('log')
        ax.set_xscale('log')
    plt.savefig(plotname.replace('_None',''))
    plt.close()
    
    
def nonHiggsEff(process, effs):
    """ calculate the efficiency for the non-higgs Background
    
    Parameters
    ----------
    process: string
              the background process in consideration
    effs:     [float*6]
              the efficiencies of [bb, cc, ss, uu, dd, gg, ww] final states
    """

    [bb, cc, ss, uu, dd, gg, ww] = effs
    if process == 'ZZstarInv':
        relativeBRs=np.array([Zbbbar, Zccbar, Zssbar, Zuubar, Zddbar, 0])
        return sum(effs*relativeBRs)
    if process == 'CEPCInv':
        return 0.653*ww + 0.061*uu + 0.06*dd + 0.064*cc + 0.06*ss + 0.098*bb + 0.00*gg
    elif process == 'WWstarInv':
        return -1
    elif process == 'WW1stGen':
        return np.sqrt(uu*dd)
    elif process == 'WW2ndGen':
        return np.sqrt(cc*ss)
    elif process == 'GG':
        return gg
    elif process == 'BB':
        return bb
    else:
        print('unknown process {0}'.format(process))
        relativeBRs=np.array([Zbbbar, Zccbar, Zssbar, Zuubar, Zddbar, 0])
        return sum(effs*relativeBRs)

def evaluateWP(HiggsEvents, NonHiggsEvents, eff, WP, chargechannel, analysischannel):
    """ calculate the upper limit and other values for a given working point
    """

    if chargechannel == '1C':
        [d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC] = WP
        signal=HiggsEvents*Hssbar*eff[chargechannel, 'ss', d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC]
        background=HiggsEvents*sum(list(map(lambda x:
                                            HBR[x]*eff[chargechannel, x, d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC],
                                            ['bb', 'cc', 'uu', 'dd', 'gg'])))
        background+=(NonHiggsEvents *
                     nonHiggsEff(analysischannel, list(map(lambda c:
                                                eff[chargechannel, c, d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC],
                                                    ['bb', 'cc', 'ss', 'uu', 'dd', 'gg', 'ww']))))
    else:
        [d0, etrack, eK, ePi, eK0, pLcut] = WP
        signal=HiggsEvents*Hssbar*eff[chargechannel, 'ss',etrack, eK, ePi, eK0, d0, pLcut]
        background=HiggsEvents*sum(list(map(lambda x:
                                            HBR[x]*eff[chargechannel, x, etrack, eK, ePi, eK0, d0, pLcut],
                                            ['bb', 'cc', 'uu', 'dd', 'gg'])))
        background+=(NonHiggsEvents *
                     nonHiggsEff(analysischannel, list(map(lambda c:
                                                eff[chargechannel, c, etrack, eK, ePi, eK0, d0, pLcut],
                                                    ['bb', 'cc', 'ss', 'uu', 'dd', 'gg', 'ww']))))
                        
    mu = Expected_UpperLimit([signal, background])

    return [signal, background, mu]

def GetBestWP(HiggsEvents, NonHiggsEvents, eff, WP, chargechannel, analysischannel, allparameters):
    """start from the workingpoint 'WP' and check if its neighbours improve the performance, iterate
    until a minimum is found

    """
    ps=parameterspace(allparameters)

    mutmp = evaluateWP(HiggsEvents, NonHiggsEvents, eff, WP, chargechannel, analysischannel)[2]
    if HiggsEvents == 10000000.0 and NonHiggsEvents == 4941713.361323838:
        print(WP)
        print(mutmp)
    change = True
    checked = []
    checked.append(WP)
    oldWP = [x for x in WP]
    while change:
        change = False
        print('checking {0} points'.format(len([p for p in ps.getneighbors(WP,1,oldWP) if p not in checked])))
        for point in [p for p in ps.getneighbors(WP,1,oldWP) if p not in checked]:
            checked.append(point)
            mutest = evaluateWP(HiggsEvents, NonHiggsEvents, eff, point, chargechannel, analysischannel)[2]
            if mutest < mutmp:
                oldWP = [x for x in WP]
                WP = point
                mutmp = mutest
                change = True
                break
        if change:
            continue
        print('checked first neighbours')

        # if nearest neighbours did not succeed, try second neighbours
        print('checking {0} points'.format(len([p for p in ps.getneighbors(WP,2,oldWP,1) if p not in checked])))
        for point in [p for p in ps.getneighbors(WP,2,oldWP,1) if p not in checked]:
            checked.append(point)
            mutest = evaluateWP(HiggsEvents, NonHiggsEvents, eff, point, chargechannel, analysischannel)[2]
            if mutest < mutmp:
                oldWP = [x for x in WP]
                WP = point
                mutmp = mutest
                change = True
                break
        if change:
            continue
        print('checked secon neighbours')
        
        # if nearest neighbours did not succeed, try third neighbours
        print('checking {0} points'.format(len([p for p in ps.getneighbors(WP,3,oldWP,1) if p not in checked])))
        for point in [p for p in ps.getneighbors(WP,3,oldWP,1) if p not in checked]:
            checked.append(point)
            mutest = evaluateWP(HiggsEvents, NonHiggsEvents, eff, point, chargechannel, analysischannel)[2]
            if mutest < mutmp:
                oldWP = [x for x in WP]
                WP = point
                mutmp = mutest
                change = True
                break
        if change:
            continue
        print('checked third neighbours')

        # if nearest neighbours did not succeed, try fourth neighbours
        print('checking {0} points'.format(len([p for p in ps.getneighbors(WP,4,oldWP,1) if p not in checked])))
        for point in [p for p in ps.getneighbors(WP,4,oldWP,1) if p not in checked]:
            checked.append(point)
            mutest = evaluateWP(HiggsEvents, NonHiggsEvents, eff, point, chargechannel, analysischannel)[2]
            if mutest < mutmp:
                oldWP = [x for x in WP]
                WP = point
                mutmp = mutest
                change = True
                break
        if change:
            continue
        print('checked fourth neighbours')

        """
        # if second neighbours did not succeed, try all neighbours
        for point in [p for p in ps.getneighbors(WP,-1,oldWP,5) if p not in checked]:
            checked.append(point)
            mutest = evaluateWP(HiggsEvents, NonHiggsEvents, eff, point, chargechannel, analysischannel)[2]
            if mutest < mutmp:
                Print_Warning('found improvement beyond 4th neighour')
                oldWP = [x for x in WP]
                WP = point
                print('oldWP = {0}  => {1}'.format(oldWP, mutmp))
                print('WP    = {0}  => {1}'.format(WP, mutest))
                mutmp = mutest
                change = True
                break
        """

    return WP
            

def transpose(listoflist):
    return list(map(list, zip(*listoflist)))

def PIDlabel(pid):
    return r'$\epsilon_{0}={1},~\epsilon_{2}={3}$'.format('K^\pm',pid[0],'\pi^\pm',pid[1])

# create static class for collider scenarios
includedscenarios=['CLIC350', 'CLIC1400', 'CLIC3000', 'ILC250', 'ILC350', 'FCCeePerfect', 'FCCee']
    

#############################################################################
#############################################################################
if __name__ == '__main__':
    from argparse import ArgumentParser
    import os
    import sys
    import numpy as np
    import itertools
    
    usage = "Significance and upper limit for collider scenarios"
    parser = ArgumentParser(prog='UpperLimits', description=usage)
    parser.add_argument('-b', '--basedir', action='store', dest='basedir',
                        help='The base directory from where the input files are found [pwd]')
    parser.add_argument( '-s', '--suffix', action='store', dest='suffix',\
                         help="output suffix for the plots [.png]")
    parser.add_argument( '-c', '--colliders', action='store', nargs='*', dest='scenarios',\
                         help=r'collider scenarios which should be analyzed. Possible choices '\
                         'are: \n{}\n default: all'.format(includedscenarios))
    parser.add_argument( '-n', '--nraster', action='store', dest='nraster', type=int,\
                         help='number of raster points per axis [100]')
    parser.add_argument( '--noPID', action='store_true', help='assume no PID is possible')

    pwd=os.getcwd()
    parser.set_defaults(suffix='.png', basedir=pwd, scenarios=includedscenarios, nraster=100)
    args = parser.parse_args()

    if args.scenarios == []:
        Print_Warning('Flag for collider scenarios raised but not used. Run with all scenarios')
        scenarios=includedscenarios
    elif not all([(x in includedscenarios) for x in args.scenarios]):
        Print_Fail('At least one chosen collider scenario is not included in the code.'\
                   ' Possible choices are {}'.format(includedscenarios))
        exit()
    else:
        scenarios = args.scenarios
    basedir   = args.basedir
    suffix    = args.suffix
    if suffix[0] != '.':
        suffix = ".{0}".format(suffix)
    noPID     = args.noPID
    if noPID:
        suffix = '_noPID{0}'.format(suffix)
    nraster   = args.nraster

    print('\nrun in directory {0}'.format(basedir))
    print('run on scenarios {}'.format(scenarios))
    print('save plots as {0}'.format(suffix))
    print('Number of raster points {0}^2'.format(nraster))

    # test if basedir exists
    try:
        basedirentries = os.listdir(basedir)
    except OSError:
        Print_Fail('Could not open {0}'.format(basedir))
        exit()

    # Read all the efficiencies
    parameters = {} # (d0, etrack, eK, ePi, eK0)
    pcutlist   = []
    eff        = {}
    firstFile  = True
    processedChargeChannels = []
    for chargechannel in HiggsChannels:
        if chargechannel in basedirentries and os.path.isdir(basedir+'/'+chargechannel):
            parameters[chargechannel]=[]

            processedChargeChannels.append(chargechannel)
            effFiles = filter(lambda x: 'combinedeffs' in x and 'txt' in x,
                              os.listdir(basedir+'/'+chargechannel))

            for effFile in effFiles:
                # extract information out of the file name
                # print('file: {0}/{1}'.format(chargechannel, effFile))
                try:
                    [dummy2, d0cut, etrack, eK, ePi, eK0, dummy2] = effFile.replace('_','-').split('-')
                    parameter = list(map(lambda x: float(x), [d0cut, etrack, eK, ePi, eK0]))
                    [d0cut, etrack, eK, ePi, eK0] = parameter
                except:
                    Print_Warning('failed with file {0}'.format(effFile))
                    continue

                if noPID and eK != 1.0:
                    continue
                
                # read the efficiencies from the file
                f = open(basedir + '/' + chargechannel + '/' +effFile, 'r')
                for line in f:
                    [pcut, effB, effC, effS, effU, effD, effG, effW200, effW250] = map(lambda x: float(x), line.split())
                    eff[chargechannel, 'bb',etrack, eK, ePi, eK0, d0cut, pcut] = effB
                    eff[chargechannel, 'cc',etrack, eK, ePi, eK0, d0cut, pcut] = effC
                    eff[chargechannel, 'ss',etrack, eK, ePi, eK0, d0cut, pcut] = effS
                    eff[chargechannel, 'uu',etrack, eK, ePi, eK0, d0cut, pcut] = effU
                    eff[chargechannel, 'dd',etrack, eK, ePi, eK0, d0cut, pcut] = effD
                    eff[chargechannel, 'gg',etrack, eK, ePi, eK0, d0cut, pcut] = effG
                    eff[chargechannel, 'ww',etrack, eK, ePi, eK0, d0cut, pcut] = effW250

                    if firstFile:
                        pcutlist.append(pcut)
                    else:
                        if not pcut in pcutlist:
                            Print_Fail('pcut={0} not the same everywhere'.format(pcut))
                f.close()
                firstFile=False

                parameters[chargechannel].append(parameter)

            parameters[chargechannel].sort()

    firstScenario=True
    for scenario in scenarios:
        break # Fixme: Remove

        for subAnalysis in Possible_SubAnalyses(scenario):
            print("")
            print("##########################################")
            print("##########################################")
            print("      {0}           {1}       ".format(scenario, subAnalysis))
            collider=colliderscenarios(scenario, subAnalysis)

            for chargechannel in processedChargeChannels:
                print("==========================================")
                print('            {0}'.format(chargechannel))
                SignalBackground={}
                SignalBackgroundOnlyHiggs={}
                significance={}
                significanceOnlyHiggs={}
                UpperLimit={}
                UpperLimitOnlyHiggs={}
                d0cutlist=[]
                pidlist=[]
                print('getting the efficiency and Nevent plots, calculate upper limits')

                progress1=progressbar(len(parameters))
                
                for (d0, etrack, eK, ePi, eK0) in parameters:
                    # print((d0, etrack, eK, ePi, eK0))

                    d0cutlist.append(d0)
                    pidlist.append([eK, ePi])
                    # make list of efficiencies
                    # [[effB], [effC], [effS], [effU], [effD], [effG], [effW], [effNon-Higgs]]
                    # where each entry is a list as funcion of pcut
                    efflist=[]
                    for pcut in pcutlist:
                        dummy = list(map(lambda x: eff[chargechannel, x, etrack, eK, ePi, eK0, d0, pcut],
                                         ['bb', 'cc', 'ss', 'uu', 'dd', 'gg', 'ww']))
                        dummy.append(nonHiggsEff('generic', dummy))
                        efflist.append(dummy)
                    efflist=np.array(transpose(efflist))
                    if firstScenario:
                        LogPlot(pcutlist, efflist, '$p_{||}^{\mathrm{cut}}$ [GeV]',
                                '$\epsilon_{s\mathrm{-tag}}$',
                                'eff_{0}_{1}_{2}_{3}_{4}_{5}{6}'.format(
                                    chargechannel, d0, etrack, eK, ePi, eK0, suffix))

                    # get number of events
                    NHiggs =collider.NHiggs
                    NnHiggs=collider.NnHiggs
                    effT=transpose(efflist)
                    pureNumbers=np.append(NHiggs*np.array([Hbbbar, Hccbar, Hssbar, Huubar, Hddbar, Hgg]),
                                          NnHiggs)
                    Nevents=transpose(effT*pureNumbers)

                    LogPlot(pcutlist, Nevents,  '$p_{||}^{\mathrm{cut}}$ [GeV]',
                                '\# events',
                                'Nevents_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}{8}'.format(
                                    scenario, subAnalysis, chargechannel, d0, etrack, eK, ePi, eK0, suffix))

                    # get signal and background numbers
                    dummy = transpose(Nevents)
                    SignalBackgroundOnlyHiggs[d0,eK,ePi] = list(map(lambda x: [x[2], sum(x)-x[2]-x[6]], dummy))
                    if NnHiggs != 0:
                        SignalBackground[d0, eK, ePi] = list(map(lambda x: [x[2], sum(x)-x[2]], dummy))

                    significanceOnlyHiggs[d0, eK, ePi] = list(map(lambda x: Expected_Significance(x),
                                                                  SignalBackgroundOnlyHiggs[d0, eK, ePi]))
                    UpperLimitOnlyHiggs[d0, eK, ePi] = list(map(lambda x: Expected_UpperLimit(x),
                                                                SignalBackgroundOnlyHiggs[d0, eK, ePi]))
                    if NnHiggs !=0:
                        UpperLimit[d0, eK, ePi] = list(map(lambda x: Expected_UpperLimit(x),
                                                           SignalBackground[d0, eK, ePi]))
                        significance[d0, eK, ePi] = list(map(lambda x: Expected_Significance(x),
                                                             SignalBackground[d0, eK, ePi]))

                    progress1.next()

                progress1.stop()
                d0cutlist=list(set(d0cutlist))
                d0cutlist.sort()
                pidlist  =sorted(map(list, set(map(tuple, pidlist))))

                print('getting the significance, upper limit and 2D plots')
                progress2=progressbar(len(d0cutlist)+ len(pidlist))
                
                # Plot significance and upper limits
                for d0 in d0cutlist:
                    significancelist=[]
                    significanceOnlyHiggslist=[]
                    upperlimitlist=[]
                    upperlimitOnlyHiggslist=[]
                    plotlabels=[]
                    for pid in pidlist:
                        significanceOnlyHiggslist.append(significanceOnlyHiggs[d0, pid[0], pid[1]])
                        upperlimitOnlyHiggslist.append(UpperLimitOnlyHiggs[d0, pid[0], pid[1]])
                        if NnHiggs != 0:
                            significancelist.append(significance[d0, pid[0], pid[1]])
                            upperlimitlist.append(UpperLimit[d0, pid[0], pid[1]])
                        plotlabels.append(PIDlabel(pid))

                    LinPlot(pcutlist, significanceOnlyHiggslist, '$p_{||}^{\mathrm{cut}}$ [GeV]',
                            'significance', plotlabels, 'SignificanceOnlyHiggs_{0}_{1}_{2}_{3}{4}'.format(
                                scenario, subAnalysis, chargechannel, d0, suffix))
                    LinPlot(pcutlist, upperlimitOnlyHiggslist, '$p_{||}^{\mathrm{cut}}$ [GeV]',
                            '95\% CL on $\mu$', plotlabels,
                            'UpperLimitOnlyHiggs_{0}_{1}_{2}_{3}{4}'.format(
                                scenario, subAnalysis, chargechannel, d0, suffix))
                    if NnHiggs != 0:
                        LinPlot(pcutlist, significancelist, '$p_{||}^{\mathrm{cut}}$ [GeV]', 'significance',
                                plotlabels, 'Significance_{0}_{1}_{2}_{3}{4}'.format(
                                    scenario, subAnalysis, chargechannel, d0, suffix))
                        LinPlot(pcutlist, upperlimitlist, '$p_{||}^{\mathrm{cut}}$ [GeV]', '95\% CL on $\mu$',
                                plotlabels, 'UpperLimit_{0}_{1}_{2}_{3}{4}'.format(
                                    scenario, subAnalysis, chargechannel, d0, suffix))

                    progress2.next()

                for pid in pidlist:
                    contourX=np.array(d0cutlist)
                    contourY=np.array(pcutlist)
                    X,Y = np.meshgrid(contourX, contourY)
                    Z=[]
                    for y in range(0, len(contourY)):
                        dummy=[]
                        for x in range(0, len(contourX)):
                            dummy.append(UpperLimitOnlyHiggs[X[0][x], pid[0], pid[1]][y])
                        Z.append(dummy)
                    
                    Plot2D(X, Y, np.array(Z), '$d_{0}$ [mm]', '$p_{||}^{\mathrm{cut}}$ [GeV]',
                           '95\% CL on $\mu$', '2DupperlimitOnlyHiggs_{0}_{1}_{2}_{3}_{4}.{5}'.format(
                               scenario, subAnalysis, chargechannel, pid[0], pid[1], suffix))

                    if NnHiggs != 0:
                        Z=[]
                        for y in range(0, len(contourY)):
                            dummy=[]
                            for x in range(0, len(contourX)):
                                dummy.append(UpperLimit[X[0][x], pid[0], pid[1]][y])
                            Z.append(dummy)
                    
                        Plot2D(X, Y, np.array(Z), '$d_{0}$ [mm]', '$p_{||}^{\mathrm{cut}}$ [GeV]',
                               '95\% CL on $\mu$', '2Dupperlimit_{0}_{1}_{2}_{3}_{4}.{5}'.format(
                                   scenario, subAnalysis, chargechannel, pid[0], pid[1], suffix))

                    progress2.next()
                progress2.stop()
                
            firstScenario=False


    # Find best cuts for varying S/B numbers
    NHiggs=np.logspace(np.log10(100), np.log10(10**7), num=nraster)
    NnHiggs=np.logspace(np.log10(100), np.log10(10**7), num=nraster)
    NH,NnH = np.meshgrid(NHiggs, NnHiggs)
    
    for analysischannel in ['CEPCInv' ]: #'WWstarInv', 'WW1stGen', 'WW2ndGen', 'GG', 'BB']: #['ZZstarInv', 'WWstarInv']:
        for chargechannel in ['1C']:
            print('{0}   {1}'.format(analysischannel, chargechannel))

            allparameters = []

            print('prepare parameter space')
            if chargechannel == '1C':
                for CCparameter in parameters['CC']: # d0cut, etrack, eK, ePi, eK0
                    [d0cutCC, etrack, eK, ePi, eK0] = CCparameter
                    for NCparameter in parameters['NC']: # d0cut, etrack, eK, ePi, eK0
                        [d0cutNC, etrackNC, eKNC, ePiNC, eK0NC] = NCparameter
                        if not CCparameter[1::]==NCparameter[1::]:
                            continue
                        for pcutCC in pcutlist:
                            for pcutNC in pcutlist:
                                allparameters.append(CCparameter+[d0cutNC, pcutCC, pcutNC])
                for parameter in allparameters:
                    [d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC] = parameter
                    eff['1C', 'bb', d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC] = \
                        eff['CC', 'bb', etrack, eK, ePi, eK0, d0cutCC, pcutCC] + \
                        eff['NC', 'bb', etrack, eK, ePi, eK0, d0cutNC, pcutNC]
                    eff['1C', 'cc', d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC] = \
                        eff['CC', 'cc', etrack, eK, ePi, eK0, d0cutCC, pcutCC] + \
                        eff['NC', 'cc', etrack, eK, ePi, eK0, d0cutNC, pcutNC]
                    eff['1C', 'ss', d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC] = \
                        eff['CC', 'ss', etrack, eK, ePi, eK0, d0cutCC, pcutCC] + \
                        eff['NC', 'ss', etrack, eK, ePi, eK0, d0cutNC, pcutNC]
                    eff['1C', 'uu', d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC] = \
                        eff['CC', 'uu', etrack, eK, ePi, eK0, d0cutCC, pcutCC] + \
                        eff['NC', 'uu', etrack, eK, ePi, eK0, d0cutNC, pcutNC]
                    eff['1C', 'dd', d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC] = \
                        eff['CC', 'dd', etrack, eK, ePi, eK0, d0cutCC, pcutCC] + \
                        eff['NC', 'dd', etrack, eK, ePi, eK0, d0cutNC, pcutNC]
                    eff['1C', 'gg', d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC] = \
                        eff['CC', 'gg', etrack, eK, ePi, eK0, d0cutCC, pcutCC] + \
                        eff['NC', 'gg', etrack, eK, ePi, eK0, d0cutNC, pcutNC]
                    eff['1C', 'ww', d0cutCC, etrack, eK, ePi, eK0, d0cutNC, pcutCC, pcutNC] = \
                        eff['CC', 'ww', etrack, eK, ePi, eK0, d0cutCC, pcutCC] + \
                        eff['NC', 'ww', etrack, eK, ePi, eK0, d0cutNC, pcutNC]
                                
            else:
                for element in itertools.product(*[parameters[chargechannel], pcutlist]):
                    allparameters.append(element[0]+[element[1]])
            print('finished preparing parameter space')
            BestWP      = {}
            Bestd0cut   = [[0 for i in range(nraster)] for j in range(nraster)]
            BestpLcut   = [[0 for i in range(nraster)] for j in range(nraster)]
            Bestpid     = [[0 for i in range(nraster)] for j in range(nraster)]
            Bestmu      = [[0 for i in range(nraster)] for j in range(nraster)]
            Bestd0cutNC = [[0 for i in range(nraster)] for j in range(nraster)]
            BestpLcutNC = [[0 for i in range(nraster)] for j in range(nraster)]
            BestSigEff  = [[0 for i in range(nraster)] for j in range(nraster)]
            BestBkgEff  = [[0 for i in range(nraster)] for j in range(nraster)]

            progress2 =progressbar(len(NnHiggs))

            f = open(basedir + '/' + chargechannel + '_' + analysischannel + '.dat', 'w')
            for [NonHiggsEvents, HiggsEvents] in itertools.product(*[NnHiggs,NHiggs]):
                NonHiggsIndex = np.where(NnHiggs==NonHiggsEvents)[0][0]
                HiggsIndex    = np.where(NHiggs ==HiggsEvents)[0][0]
                print('GetBestWP[{0}, {1}]'.format(NonHiggsIndex, HiggsIndex))
                if NonHiggsIndex == HiggsIndex == 0:
                    BestWP[(NonHiggsEvents, HiggsEvents)] = GetBestWP(HiggsEvents, NonHiggsEvents, eff,
                                                                allparameters[0], chargechannel,
                                                                analysischannel, allparameters)
                elif NonHiggsIndex != 0:
                    BestWP[(NonHiggsEvents, HiggsEvents)] = GetBestWP(HiggsEvents, NonHiggsEvents, eff,
                                                                BestWP[(NnHiggs[NonHiggsIndex-1], HiggsEvents)],
                                                                chargechannel, analysischannel, allparameters)
                elif NonHiggsIndex == 0:
                    BestWP[(NonHiggsEvents, HiggsEvents)] = GetBestWP(HiggsEvents, NonHiggsEvents, eff,
                                                                BestWP[(NonHiggsEvents, NHiggs[HiggsIndex-1])],
                                                                chargechannel, analysischannel, allparameters)
                if chargechannel == '1C':
                    [d0, etrack, eK, ePi, eK0, d0cutNC, pLcut, pcutNC] = BestWP[(NonHiggsEvents, HiggsEvents)]
                    Bestd0cutNC[NonHiggsIndex][HiggsIndex] = d0cutNC
                    BestpLcutNC[NonHiggsIndex][HiggsIndex] = pcutNC
                else:
                    [d0, etrack, eK, ePi, eK0, pLcut] = BestWP[(NonHiggsEvents, HiggsEvents)]
                    d0cutNC = 0
                    pcutNC  = 0
                    Bestd0cutNC[NonHiggsIndex][HiggsIndex] = d0cutNC
                    BestpLcutNC[NonHiggsIndex][HiggsIndex] = pcutNC
                Bestd0cut[NonHiggsIndex][HiggsIndex]  = d0
                BestpLcut[NonHiggsIndex][HiggsIndex]  = pLcut
                Bestpid[NonHiggsIndex][HiggsIndex]    = eK
                [signal, background, mu]              = evaluateWP(HiggsEvents, NonHiggsEvents, eff,
                                                                   BestWP[(NonHiggsEvents, HiggsEvents)],
                                                                   chargechannel, analysischannel)
                sigeff                                = signal/(HiggsEvents*Hssbar)
                bkgeff                                = background/(HiggsEvents*sum(list(map(lambda x: HBR[x],
                                                            ['bb', 'cc', 'uu', 'dd', 'gg'])))+NonHiggsEvents)
                Bestmu[NonHiggsIndex][HiggsIndex]     = mu
                BestSigEff[NonHiggsIndex][HiggsIndex] = sigeff
                BestBkgEff[NonHiggsIndex][HiggsIndex] = bkgeff

                f.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(
                    HiggsEvents, NonHiggsEvents, d0 , pLcut, eK, d0cutNC, pcutNC, mu, sigeff, bkgeff))
                if HiggsIndex == nraster-1:
                    progress2.next()

            progress2.stop()

            Plot2D(NH, NnH, np.array(Bestmu), r'\# $h \to jj$ events', '\# non-$h2j$ events',
                               '95\% CL on $\mu$', 'BestMu_{0}_{1}.{2}'.format(
                                   chargechannel, analysischannel, suffix), True)
            Plot2D(NH, NnH, np.array(Bestd0cut), r'\# $h \to jj$ events', '\# non-$h2j$ events',
                               '$d_0^\mathrm{cut}|_\mathrm{best}$ [$\mu$m]', 'Bestd0_{0}_{1}.{2}'.format(
                                   chargechannel,analysischannel, suffix), True)
            Plot2D(NH, NnH, np.array(BestpLcut), r'\# $h \to jj$ events', '\# non-$h2j$ events',
                               '$p_{||}^\mathrm{cut}|_\mathrm{best}$ [GeV]', 'BestpL_{0}_{1}.{2}'.format(
                                   chargechannel,analysischannel, suffix), True)
            Plot2D(NH, NnH, np.array(Bestd0cutNC), r'\# $h \to jj$ events', '\# non-$h2j$ events',
                               '$d_0^\mathrm{cut}|_\mathrm{best}$ [$\mu$m]', 'Bestd0NC_{0}_{1}.{2}'.format(
                                   chargechannel,analysischannel, suffix), True)
            Plot2D(NH, NnH, np.array(BestpLcutNC), r'\# $h \to jj$ events', '\# non-$h2j$ events',
                               '$p_{||}^\mathrm{cut}|_\mathrm{best}$ [GeV]', 'BestpLNC_{0}_{1}.{2}'.format(
                                   chargechannel,analysischannel, suffix), True)
            Plot2D(NH, NnH, np.array(Bestpid), r'\# $h \to jj$ events', '\# non-$h2j$ events',
                               '$\epsilon_{K^\pm}|_\mathrm{best}$', 'BestPID_{0}_{1}.{2}'.format(
                                   chargechannel,analysischannel, suffix), True)
            Plot2D(NH, NnH, np.array(BestSigEff)/np.array(BestBkgEff), r'\# $h \to jj$ events',
                   '\# non-$h2j$ events', '$\epsilon_s/\epsilon_b$',
                   'BestEffRatio_{0}_{1}.{2}'.format(chargechannel,analysischannel, suffix), True)
            Plot2D(NH, NnH, np.array(BestSigEff)/np.sqrt(np.array(BestBkgEff)), r'\# $h \to jj$ events',
                   '\# non-$h2j$ events', '$\epsilon_s/\sqrt{\epsilon_b}$',
                   'BestEffsqrtRatio_{0}_{1}.{2}'.format(chargechannel,analysischannel, suffix), True)
