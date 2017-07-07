#!/usr/bin/env python

# use the separation of two gaussians with unit width and a kaon reconstruction efficiency as input
# and return the corresponding pion reconstruction efficiency. The separation is given in units of
# the combined widths, i.e. in units of sqrt(2)
#
# usage: kaonpionseparation.py kaonefficiency separation
#
# matthias.schlaffer@weizmann.ac.il
#

def effpi(effk, sep):
    """ The efficiency for mis-tagging a charged pion as kaon

    Parameters
    -----------------------
    effk: float
          The efficiency for tagging a charged kaon
    sep:  float
          The separation of the two gaussians corresponding to the kaon and pion distributions of
          some observable that discriminates between them
    """
    from scipy.special import erfc
    from scipy.special import erfcinv

    return 1-0.5*erfc(-sep + erfcinv(2-2*effk))


if __name__ == '__main__':
    from argparse import ArgumentParser
    import os
    import numpy as np

    usage = "Get the efficiency for the reconstruction of charged pions"
    parser = ArgumentParser(prog='kaonpionseparation', description=usage)
    parser.add_argument('eff_and_sep', nargs=2, action='store', type=float, default=[0.9, 1])

    args=parser.parse_args()
    [eff, sep] = args.eff_and_sep

    print '{0:.{1}f}'.format(effpi(eff,sep),4)
    
