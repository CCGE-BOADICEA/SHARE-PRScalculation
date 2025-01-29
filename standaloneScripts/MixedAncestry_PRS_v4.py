import argparse
import csv
import numpy as np

###############
##USER INPUTS##
###############

description = """
\nCommand line tool for calculating PGSs for mixed ancestry individuals.
To calculate an individual's PGS, these inputs must be provided:
\t$ -p <pgs name> - full path to the PGS info file
\t$ -a <ancestry path> - full path to the ancestry proportions file
"""
    
# Set up a parser with the desired options
parser = argparse.ArgumentParser(description=description,
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument(
    '-p',
    '--pgs',
    help='PGS info file contains SNPs frequencies and logOR, and PGS model parameters',
    default='./BCAC_307_PRS-UKB_all.prs',
    metavar='<PRS name>',
    dest='pars_full_path',
    type=str,
    action='store'
    )

parser.add_argument(
    '-a',
    '--ancestries',
    help='Ancestry proportions file contains the raw PGS and the ancestry proportions',
    default='./example_values.csv',
    metavar='<Input full path>',
    dest='props_full_path',
    type=str,
    action='store'
    )

# Complete paths from user inputs
File_LociPars = parser.parse_args().pars_full_path
###File_AncestryPars = ????
File_Ancestries = parser.parse_args().props_full_path


#################
##INPUT READING##
#################

# sample-specific parameters; raw PRS and population proportion, in one line
tmpTab = np.genfromtxt(File_Ancestries, delimiter=',',skip_header=1)
rawPRS = tmpTab[0]
sampleProp = tmpTab[1:5]
if np.sum(sampleProp) != 1:
    raise ValueError("Error in the input file!")

# SNPs log(OR) and EAFs
tmpTab = np.genfromtxt(File_LociPars, delimiter=',',skip_header=9,usecols = (4,6,7,8,9),comments='#')
logORs = tmpTab[:,0]
EAFtab = tmpTab[:,1:]

# ethnic parameters; europeans, african, eastAsian, southAsian
tmpTab = np.genfromtxt(File_LociPars, delimiter=',',skip_header=0,max_rows=4,comments='#')
POPmeans = tmpTab[1,1: ]
POPsigmas = tmpTab[2,1: ]
alphas = tmpTab[3,1: ]
monoThresh = tmpTab[0,1: ]
monoFlag = (sampleProp >= monoThresh)


################
##CALCULATIONS##
################

if np.any(monoFlag):
    # note: if the sample were truly mono-Ancestry, then the results from this if clause and those from the else clause would be identical
    monoIndex = (np.nonzero(monoFlag)[0]).item()
    Z_score = (rawPRS - POPmeans[monoIndex]) / POPsigmas[monoIndex]
    alphaPar = alphas[monoIndex]
else:
    ##################
    ##WORK ON SIGMAS##
    ##################

    ## calculate the possible overall sigma
    SNPtmp1 = np.sum((EAFtab**2) * sampleProp, axis=1)
    SNPtmp2 = np.sum(EAFtab * sampleProp, axis=1)

    stDev = np.sum((POPsigmas**2) * sampleProp)
    stDev += 4 * np.sum((SNPtmp1 - (SNPtmp2**2)) * (logORs**2))
    stDev = np.sqrt(stDev)

    ## calculate the weighted alphas and PRS
    tmp1 = np.sum(alphas * sampleProp * (rawPRS - POPmeans) / POPsigmas)
    tmp2 = np.sum(alphas * sampleProp / POPsigmas)
    tmp3 = np.sum(alphas * sampleProp * POPsigmas)

    Z_score = (tmp1 / tmp2) / stDev
    alphaPar = tmp3 / stDev

print(f"alpha= {alphaPar:.6f}; PRS_Zscore= {Z_score:.6f}")
