#!/usr/bin/env python
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
Hbr=[Hbbbar, Hccbar, Hssbar, Huubar, Hddbar, Hgg]

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

    
# create static class for collider scenarios
scenarios=['CLIC350', 'CLIC1400', 'CLIC3000']
scenarios=['FCCeePerfect']
    
def Possible_SubAnalyses(scenario):
    """ Return the possible analysis channels for a given collider scenario
    """

    if scenario == 'CLIC350':
        return ['Ztoinv', 'Ztohad']
    else:
        return [None]

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
        sys.stdout.write(self.progressstring.format('='*self.counter, int(round(100./self.nsteps*self.counter))))
        sys.stdout.flush()

    def stop(self):
        import sys
        sys.stdout.write('\n')

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
           a list containing lists of the y values. the order is bb, cc, ss, uu, dd, gg, (non-Higgs)
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

    if len(ys)==6:
        fs=[r'$b\bar{b}$',r'$c\bar{c}$',r'$s\bar{s}$',r'$u\bar{u}$',r'$d\bar{d}$',r'$gg$']
    elif len(ys)==7:
        fs=[r'$b\bar{b}$',r'$c\bar{c}$',r'$s\bar{s}$',r'$u\bar{u}$',r'$d\bar{d}$',r'$gg$', r'non-Higgs']
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

    plt.savefig(plotname)
    plt.close()

    
def LinPlot(x, ys, xlabel, ylabel, plotlabels, plotname):
    """ make a plot 

    parameters:
    -----------
    x: [floats]
           the x values of all points
    ys [lists of floats]
           a list containing lists of the y values. the order is bb, cc, ss, uu, dd, gg, (non-Higgs)
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

    plt.savefig(plotname)
    plt.close()

def Plot2D(X, Y, Z, xlabel, ylabel, zlabel, plotname):
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
    """

    import matplotlib.pyplot as plt
 
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.tick_params(top=True, right=True)
    contour=plt.contourf(X, Y, Z, cmap=plt.cm.viridis)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    cbar = plt.colorbar(contour)
    cbar.ax.set_ylabel(zlabel)
    plt.savefig(plotname)
    plt.close()
    
    
def nonHiggsEff(scenario, subAnalysis, channel, effs):
    """ calculate the efficiency for the non-higgs Background
    
    Parameters
    ----------
    scenario: string
              the collider scenario in consideration
    subAnalysis: string
              the sub analysis in consideration, e.g. HZ, with Z to invisible or Z to hadronic
    channel:  string
              the channel in consideration, CC, NC, or NN
    effs:     [float*6]
              the efficiencies of [bb, cc, ss, uu, dd, gg] final states
    """

    relativeBRs=np.array([Zbbbar, Zccbar, Zssbar, Zuubar, Zddbar, 0])
    return sum(effs*relativeBRs)


def transpose(listoflist):
    return list(map(list, zip(*listoflist)))

def PIDlabel(pid):
    return r'$\epsilon_{0}={1},~\epsilon_{2}={3}$'.format('K^\pm',pid[0],'\pi^\pm',pid[1])


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

    pwd=os.getcwd()
    parser.set_defaults(suffix='.png', basedir=pwd)
    args = parser.parse_args()

    basedir = args.basedir
    suffix  = args.suffix
    print('\nrun in directory {0}'.format(basedir))
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
    processedChannels = []
    for channel in HiggsChannels:
        if channel in basedirentries and os.path.isdir(basedir+'/'+channel):
            processedChannels.append(channel)
            effFiles = filter(lambda x: 'efficiencies' in x and 'txt' in x, os.listdir(basedir+'/'+channel))

            for effFile in effFiles:
                # extract information out of the file name
                # print('file: {0}/{1}'.format(channel, effFile))
                try:
                    [dummy2, d0cut, etrack, eK, ePi, eK0, dummy2] = effFile.replace('_','-').split('-')
                    parameter = list(map(lambda x: float(x), [d0cut, etrack, eK, ePi, eK0]))
                    [d0cut, etrack, eK, ePi, eK0] = parameter
                except:
                    Print_Warning('failed with file {0}'.format(effFile))
                    continue

                # read the efficiencies from the file
                f = open(basedir + '/' + channel + '/' +effFile, 'r')
                for line in f:
                    [pcut, effB, effC, effS, effU, effD, effG] = map(lambda x: float(x), line.split())
                    eff[channel, 'bb',etrack, eK, ePi, eK0, d0cut, pcut] = effB
                    eff[channel, 'cc',etrack, eK, ePi, eK0, d0cut, pcut] = effC
                    eff[channel, 'ss',etrack, eK, ePi, eK0, d0cut, pcut] = effS
                    eff[channel, 'uu',etrack, eK, ePi, eK0, d0cut, pcut] = effU
                    eff[channel, 'dd',etrack, eK, ePi, eK0, d0cut, pcut] = effD
                    eff[channel, 'gg',etrack, eK, ePi, eK0, d0cut, pcut] = effG

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
                        Print_Fail('parameters {0} not the same in each channel'.format(parameter))

            if firstChannel:
                firstChannel=False
                parameters.sort()
            

    firstScenario=True
    for scenario in scenarios:
        for subAnalysis in Possible_SubAnalyses(scenario):
            print("")
            print("##########################################")
            print("##########################################")
            print("      {0}           {1}       ".format(scenario, subAnalysis))
            collider=colliderscenarios(scenario, subAnalysis)

            for channel in processedChannels:
                print("==========================================")
                print('            {0}'.format(channel))
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
                    # [[effB], [effC], [effS], [effU], [effD], [effG], [effNon-Higgs]]
                    # where each entry is a list as funcion of pcut
                    efflist=[]
                    for pcut in pcutlist:
                        dummy = list(map(lambda x: eff[channel, x, etrack, eK, ePi, eK0, d0, pcut], ['bb', 'cc', 'ss', 'uu', 'dd', 'gg']))
                        dummy.append(nonHiggsEff(scenario, subAnalysis, channel, dummy))
                        efflist.append(dummy)
                    efflist=np.array(transpose(efflist))
                    if firstScenario:
                        LogPlot(pcutlist, efflist, '$p_{||}^{\mathrm{cut}}$ [GeV]',
                                '$\epsilon_{s\mathrm{-tag}}$',
                                'eff_{0}_{1}_{2}_{3}_{4}_{5}{6}'.format(channel, d0, etrack, eK, ePi, eK0, suffix))

                    # get number of events
                    NHiggs =collider.NHiggs
                    NnHiggs=collider.NnHiggs
                    effT=transpose(efflist)
                    pureNumbers=np.append(NHiggs*np.array([Hbbbar, Hccbar, Hssbar, Huubar, Hddbar, Hgg]),
                                          NnHiggs)
                    Nevents=transpose(effT*pureNumbers)

                    LogPlot(pcutlist, Nevents,  '$p_{||}^{\mathrm{cut}}$ [GeV]',
                                '\# events',
                                'Nevents_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}{8}'.format(scenario, subAnalysis, channel, d0, etrack, eK, ePi, eK0, suffix))

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

                    LinPlot(pcutlist, significanceOnlyHiggslist, '$p_{||}^{\mathrm{cut}}$ [GeV]', 'significance',
                            plotlabels, 'SignificanceOnlyHiggs_{0}_{1}_{2}_{3}{4}'.format(
                                scenario, subAnalysis, channel, d0, suffix))
                    LinPlot(pcutlist, upperlimitOnlyHiggslist, '$p_{||}^{\mathrm{cut}}$ [GeV]',
                            '95\% CL on $\mu$', plotlabels,
                            'UpperLimitOnlyHiggs_{0}_{1}_{2}_{3}{4}'.format(
                                scenario, subAnalysis, channel, d0, suffix))
                    if NnHiggs != 0:
                        LinPlot(pcutlist, significancelist, '$p_{||}^{\mathrm{cut}}$ [GeV]', 'significance',
                                plotlabels, 'Significance_{0}_{1}_{2}_{3}{4}'.format(
                                    scenario, subAnalysis, channel, d0, suffix))
                        LinPlot(pcutlist, upperlimitlist, '$p_{||}^{\mathrm{cut}}$ [GeV]', '95\% CL on $\mu$',
                                plotlabels, 'UpperLimit_{0}_{1}_{2}_{3}{4}'.format(
                                    scenario, subAnalysis, channel, d0, suffix))

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
                               scenario, subAnalysis, channel, pid[0], pid[1], suffix))

                    if NnHiggs != 0:
                        Z=[]
                        for y in range(0, len(contourY)):
                            dummy=[]
                            for x in range(0, len(contourX)):
                                dummy.append(UpperLimit[X[0][x], pid[0], pid[1]][y])
                            Z.append(dummy)
                    
                        Plot2D(X, Y, np.array(Z), '$d_{0}$ [mm]', '$p_{||}^{\mathrm{cut}}$ [GeV]',
                               '95\% CL on $\mu$', '2Dupperlimit_{0}_{1}_{2}_{3}_{4}{5}'.format(
                                   scenario, subAnalysis, channel, pid[0], pid[1], suffix))

                    progress2.next()
                progress2.stop()
                
            firstScenario=False
