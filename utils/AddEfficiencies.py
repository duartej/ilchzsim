#!/usr/bin/env python
"""script to add the efficiencies of the different charge channels. Use only if events can appear in
   only one of the channels.

"""

#############################################################################
#############################################################################
if __name__ == '__main__':
    from argparse import ArgumentParser
    import os
    import sys
    import shutil
    import numpy as np

    usage = "Add efficiencies of different charge channels"
    parser = ArgumentParser(prog='AddEfficiencies', description=usage)
    parser.add_argument('-d', '--basedir', action='store', dest='basedir',
                        help='The base directory from where the input files are found [pwd]')
    parser.add_argument('-o', action='store', dest='outchannel',
                        help='name of the resulting combined channel [1C]')
    parser.add_argument('-c', '--channels', action='store', nargs='*', dest='channels',
                        help='channels that will be added [CC, NC]')


    pwd=os.getcwd()
    parser.set_defaults(basedir=pwd, outchannel='1C', channels=['CC','NC'])
    args = parser.parse_args()
    basedir=args.basedir
    outchannel=args.outchannel
    channels=args.channels

    # make sure that channels contains more than 1 channel
    try:
        if len(channels) < 2:
            raise ValueError
    except:
        print('Need to specify more than 1 channel')
        exit()
    
    # test if basedir exists
    try:
        basedirentries = os.listdir(basedir)
    except OSError:
        Print_Fail('Could not open {0}'.format(basedir))
        exit()

    # check that all required channels are present
    try:
        if not set(channels) <= set(basedirentries):
            raise ValueError
    except:
        print('Not all channels are present in this directory')
        exit()

    # check that the combined channel directory is not yet present, otherwise delete.
    # create the outchannel directory
    if outchannel in basedirentries and os.path.isdir(basedir+'/'+outchannel):
        print('The directory {0} already exists. Deleting it'.format(outchannel))
        shutil.rmtree(basedir+'/'+outchannel)
    os.makedirs(basedir+'/'+outchannel)


    effFiles = filter(lambda x: 'efficiencies' in x and 'txt' in x and '~' not in x,
                      os.listdir(basedir+'/'+channels[0]))

    print('\nEfficiency files that are considered:')
    print('\n'.join(effFiles))
    print('\n')
    
    for effFile in effFiles:
        combinedEff={}
        pcuts=[]
        first=True

        for channel in channels:
            try:
                f = open(basedir+'/'+channel+'/'+effFile, 'r')
            except FileNotFoundError:
                print('Did not find {0}/{1}/{2}'.format(basedir,channel,effFile))

            for line in f:
                dummy=line.split()
                pcut, effs = int(dummy[0]), np.array(list(map(lambda x: float(x), dummy[1:])))
                if first:
                    combinedEff[pcut]=effs
                    pcuts.append(pcut)
                else:
                    try:
                        combinedEff[pcut]+=effs
                    except KeyError:
                        print('KeyError for pcut={0} in channel {1}'.format(pcut,channel))
            first=False
            f.close()
            
        result=[]
        
        o = open(basedir+'/'+ outchannel +'/'+effFile, 'w')
        for pcut in pcuts:
            
            o.write('{0:3}   {1}\n'.format(pcut, '   '.join(map('{0:.4e}'.format, combinedEff[pcut]))))
        #np.savetxt(basedir+'/'+outchannel+'/'+effFile, result, fmt='%.4e')
        o.close()

    print('done')
