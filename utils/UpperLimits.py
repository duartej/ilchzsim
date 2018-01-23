#!/usr/bin/env python3
"""script to obtain the realistic significance and upper limits for different collider scenarios
   =========================

   usage: 
"""
# Final states in the Higgs channel
HiggsChannels=['CC', 'NC', 'NN', '1C']

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
                self.NnHiggs=14106.
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

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=13)
    
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.tick_params(top=True, right=True)
    if 'BestMu' in plotname:
        levels=[1,2,5,10,20,50,100,200,500]
        contour=plt.contourf(X, Y, Z, cmap=plt.cm.viridis, levels=levels, norm=colors.LogNorm())
        cbar=plt.colorbar(ticks=levels)
        cbar.set_ticklabels(levels)
    elif 'BestpL' in plotname:
        levels = [0,2,4,6,8,10,12,14,16,18,20]
        contour=plt.contourf(X, Y, Z, cmap=plt.cm.viridis, levels=levels)
        cbar=plt.colorbar(ticks=levels)
        cbar.set_ticklabels(levels)
    elif 'Bestd0' in plotname:
        levels = [0.016, 0.018, 0.02, 0.022, 0.024]
        levelsmu = list(map(lambda x: x*1000, levels))
        contour=plt.contourf(X, Y, Z, cmap=plt.cm.viridis, levels=levels)
        cbar=plt.colorbar(ticks=levels)
        cbar.set_ticklabels(levelsmu)
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
        plt.plot(FCCee.NHiggs, FCCee.NnHiggs, 'r+', markersize=13, label='FCCee')
        plt.plot(FCCee.roclist[0], FCCee.roclist[1], 'r+', markersize=5)
        plt.plot(FCCee.chi2list[0], FCCee.chi2list[1], 'mP', markersize=5, label='FCCee $\chi^2$')
        CEPC=colliderscenarios('CEPC', 'inv')
        plt.plot(CEPC.NHiggs, CEPC.NnHiggs, 'mP', markersize=13, label='CEPC')
        ILC=colliderscenarios('ILC250', 'inv')
        plt.plot(ILC.NHiggs, ILC.NnHiggs, 'b*', markersize=13, label='ILC 250')
        plt.plot([100,10**7],[100,10**7],'--k')
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
    w = 0.5*np.sqrt(uu*dd) + 0.5*np.sqrt(cc*ss)
    if process == 'ZZstarInv':
        relativeBRs=np.array([Zbbbar, Zccbar, Zssbar, Zuubar, Zddbar, 0])
        return sum(effs*relativeBRs)
    if process == 'CEPCInv':
        return 0.16*ww + 0.06*uu + 0.06*dd + 0.06*cc + 0.06*ss + 0.10*bb + 0.49*w + 0.00*gg
    elif process == 'WWstarInv':
        return w
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

    usage = "Significance and upper limit for collider scenarios"
    parser = ArgumentParser(prog='UpperLimits', description=usage)
    parser.add_argument('-b', '--basedir', action='store', dest='basedir',
                        help='The base directory from where the input files are found [pwd]')
    parser.add_argument( '-s', '--suffix', action='store', dest='suffix',\
                         help="output suffix for the plots [.png]")
    parser.add_argument( '-c', '--colliders', action='store', nargs='*', dest='scenarios',\
                         help=r'collider scenarios which should be analyzed. Possible choices '\
                         'are: \n{}\n default: all'.format(includedscenarios))
                         
    pwd=os.getcwd()
    parser.set_defaults(suffix='.png', basedir=pwd, scenarios=includedscenarios)
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

    print('\nrun in directory {0}'.format(basedir))
    print('run on scenarios {}'.format(scenarios))
    print('save plots as {0}'.format(suffix))

    # test if basedir exists
    try:
        basedirentries = os.listdir(basedir)
    except OSError:
        Print_Fail('Could not open {0}'.format(basedir))
        exit()

    # Read all the efficiencies
    parameters = [] # (d0, etrack, eK, ePi, eK0)
    pcutlist   = []
    eff        = {}
    firstFile  = True
    firstChannel      = True
    processedChargeChannels = []
    for chargechannel in HiggsChannels:
        if chargechannel in basedirentries and os.path.isdir(basedir+'/'+chargechannel):
            processedChargeChannels.append(chargechannel)
            effFiles = filter(lambda x: 'efficiencies' in x and 'txt' in x,
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

                # read the efficiencies from the file
                f = open(basedir + '/' + chargechannel + '/' +effFile, 'r')
                for line in f:
                    [pcut, effB, effC, effS, effU, effD, effG, effW] = map(lambda x: float(x), line.split())
                    eff[chargechannel, 'bb',etrack, eK, ePi, eK0, d0cut, pcut] = effB
                    eff[chargechannel, 'cc',etrack, eK, ePi, eK0, d0cut, pcut] = effC
                    eff[chargechannel, 'ss',etrack, eK, ePi, eK0, d0cut, pcut] = effS
                    eff[chargechannel, 'uu',etrack, eK, ePi, eK0, d0cut, pcut] = effU
                    eff[chargechannel, 'dd',etrack, eK, ePi, eK0, d0cut, pcut] = effD
                    eff[chargechannel, 'gg',etrack, eK, ePi, eK0, d0cut, pcut] = effG
                    eff[chargechannel, 'ww',etrack, eK, ePi, eK0, d0cut, pcut] = effW

                    if firstFile:
                        pcutlist.append(pcut)
                    else:
                        if not pcut in pcutlist:
                            Print_Fail('pcut={0} not the same everywhere'.format(pcut))
                f.close()
                firstFile=False

                
                if firstChannel:
                    parameters.append(parameter)
                else:
                    if not parameter in parameters:
                        Print_Fail('parameters {0} not the same in each chargechannel'.format(parameter))

            if firstChannel:
                firstChannel=False
                parameters.sort()
            
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
                           '95\% CL on $\mu$', '2DupperlimitOnlyHiggs_{0}_{1}_{2}_{3}_{4}{5}'.format(
                               scenario, subAnalysis, chargechannel, pid[0], pid[1], suffix))

                    if NnHiggs != 0:
                        Z=[]
                        for y in range(0, len(contourY)):
                            dummy=[]
                            for x in range(0, len(contourX)):
                                dummy.append(UpperLimit[X[0][x], pid[0], pid[1]][y])
                            Z.append(dummy)
                    
                        Plot2D(X, Y, np.array(Z), '$d_{0}$ [mm]', '$p_{||}^{\mathrm{cut}}$ [GeV]',
                               '95\% CL on $\mu$', '2Dupperlimit_{0}_{1}_{2}_{3}_{4}{5}'.format(
                                   scenario, subAnalysis, chargechannel, pid[0], pid[1], suffix))

                    progress2.next()
                progress2.stop()
                
            firstScenario=False


    # Find best cuts for varying S/B numbers
    nraster=5
    NHiggs=np.logspace(np.log10(100), np.log10(10**7), num=nraster)
    NnHiggs=np.logspace(np.log10(100), np.log10(10**7), num=nraster)
    NH,NnH = np.meshgrid(NHiggs, NnHiggs)

    for analysischannel in ['CEPCInv' ]: #'WWstarInv', 'WW1stGen', 'WW2ndGen', 'GG', 'BB']: #['ZZstarInv', 'WWstarInv']:
        for chargechannel in ['CC', '1C']:
            print('{0}   {1}'.format(analysischannel, chargechannel))
            Bestd0cut=[]
            BestpLcut=[]
            Bestpid  =[]
            Bestmu   =[]
            BestSigEff=[]
            BestBkgEff=[]

            progress2 =progressbar(len(NnHiggs))

            f = open(basedir + '/' + chargechannel + '_' + analysischannel + '.dat', 'w')
            for y in range(len(NnHiggs)):
                NonHiggsEvents=NnHiggs[y]

                d0cutdummy =[]
                pLcutdummy =[]
                piddummy   =[]
                mudummy    =[]
                sigeffdummy=[]
                bkgeffdummy=[]

                for x in range(len(NHiggs)):
                    HiggsEvents=NHiggs[x]

                    currentbest=[10**10, -1, -1, -1, 0, 0]
                    for (d0, etrack, eK, ePi, eK0) in parameters:
                        for pLcut in pcutlist:
                            signal=HiggsEvents*Hssbar*eff[chargechannel, 'ss',etrack, eK, ePi, eK0, d0, pLcut]
                            background=HiggsEvents*sum(list(map(lambda x:
                                                HBR[x]*eff[chargechannel, x, etrack, eK, ePi, eK0, d0, pLcut],
                                                ['bb', 'cc', 'uu', 'dd', 'gg'])))
                            background+=(NonHiggsEvents *
                                         nonHiggsEff(analysischannel, list(map(lambda c:
                                        eff[chargechannel, c, etrack, eK, ePi, eK0, d0, pLcut],
                                                                      ['bb', 'cc', 'ss', 'uu', 'dd', 'gg', 'ww']))))
                        
                            mutmp=Expected_UpperLimit([signal, background])
                            if mutmp < currentbest[0]:
                                currentbest=[mutmp, d0, pLcut, eK, signal/(HiggsEvents*Hssbar),
                                             background/(HiggsEvents*sum(list(map(lambda x: HBR[x],
                                            ['bb', 'cc', 'uu', 'dd', 'gg'])))+NonHiggsEvents)]
                    f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(
                        HiggsEvents, NonHiggsEvents, currentbest[1], currentbest[2],
                        currentbest[3], currentbest[0], currentbest[4], currentbest[5]))
                    d0cutdummy.append(currentbest[1])
                    pLcutdummy.append(currentbest[2])
                    piddummy.append(currentbest[3])
                    mudummy.append(currentbest[0])
                    sigeffdummy.append(currentbest[4])
                    bkgeffdummy.append(currentbest[5])

                Bestd0cut.append(d0cutdummy)
                BestpLcut.append(pLcutdummy)
                Bestpid.append(piddummy)
                Bestmu.append(mudummy)
                BestBkgEff.append(bkgeffdummy)
                BestSigEff.append(sigeffdummy)
                progress2.next()

            f.close()
            progress2.stop()

            Plot2D(NH, NnH, np.array(Bestmu), '\# Higgs events', '\# non-Higgs events',
                               '95\% CL on $\mu$', 'BestMu_{0}_{1}{2}'.format(
                                   chargechannel, analysischannel, suffix), True)
            Plot2D(NH, NnH, np.array(Bestd0cut), '\# Higgs events', '\# non-Higgs events',
                               '$d_0^\mathrm{cut}|_\mathrm{best}$ [$\mu$m]', 'Bestd0_{0}_{1}{2}'.format(
                                   chargechannel,analysischannel, suffix), True)
            Plot2D(NH, NnH, np.array(BestpLcut), '\# Higgs events', '\# non-Higgs events',
                               '$p_{||}^\mathrm{cut}|_\mathrm{best}$ [GeV]', 'BestpL_{0}_{1}{2}'.format(
                                   chargechannel,analysischannel, suffix), True)
            Plot2D(NH, NnH, np.array(Bestpid), '\# Higgs events', '\# non-Higgs events',
                               '$\epsilon_{K^\pm}|_\mathrm{best}$', 'BestPID_{0}_{1}{2}'.format(
                                   chargechannel,analysischannel, suffix), True)
