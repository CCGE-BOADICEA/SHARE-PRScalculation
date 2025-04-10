"""
:authors: Andrew Lee, Lorenzo Ficorella
:organization: CCGE, DPHPC, University of Cambridge
:copyright:  2022 Cambridge University. All rights reserved.
:   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
:   without even the implied warranty of MERCHANTABILITY or FITNESS for a particular purpose.
:contact: lf364@medschl.cam.ac.uk
"""
import vcf2prs
import argparse
import sys
import csv


def main():
    """
    This function provides a command-line interface for the vcf2prs command.
    """
    resFormat1 = "raw PRS: {}\nalpha: {}\nz-score: {}"
    resFormat2 = "mean: {}\nstandard deviation: {}\nalpha: {}"

    description = """
    \nCommand line tool for calculating PRSs. It can perform 4 calculations:
    1. Calculate an individual's PRS from a VCF file, using the -g option.
    2. Convert a raw PRS to a z-score, using the -r option.
    3. Convert a z-score to a raw PRS, using the -z option.
    4. Display the characteristics of a PRS, using no option.
    \nExamples
    1. To calculate an individual's PRS:
    \t$ vcf2prs -g sample_VCF_files/sample_BCAC_313_genotypes.vcf -s Med\
     PRS_files/BCAC_313_PRS.prs
    2. To convert a raw PRS to a z-score:
    \t$ vcf2prs -r -0.3815 PRS_files/BCAC_313_PRS.prs
    3. To convert a z-score to a raw PRS:
    \t$ vcf2prs -z 2.80501 PRS_files/BCAC_313_PRS.prs
    4. To calculate the characteristics of a PRS:
    \t$ vcf2prs PRS_files/BCAC_313_PRS.prs
    \nIf ancestry proportions are provided, calculations 1 and 2 can be performed
    for mixed-ancestry individuals
    """
    # Set up a parser with the desired options
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        'prs_file_name',
        metavar='<PRS File name>',
        type=str,
        nargs=1,
        help='Name of PRS info file'
        )

    groupAncs = parser.add_mutually_exclusive_group()
    
    groupAncs.add_argument(
        '-a',
        '--ancestry',
        choices=[0,1,2,3,4],
        help=('Number corresponding to the ancestry to be considered.\n'
              '0= unspecified; 1= European; 2= African; '
              '3= East Asian; 4= South Asian.'),
        default=0,
        metavar='<ancestry code>',
        dest='ancestry_code',
        type=int,
        action='store'
        )

    groupAncs.add_argument(
        '-p',
        '--proportions',
        help=('Name of the file containing the ancestry proportions.\n'
              'Required for PGS calculations with mixed ancestry.'),
        default=None,
        metavar='<props file name>',
        dest='props_file_name',
        type=str,
        action='store'
        )

    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='%(prog)s {}'.format(vcf2prs.__version__)
        )

    groupCalc = parser.add_mutually_exclusive_group()

    groupCalc.add_argument(
        '-g',
        '--geno',
        help='Name of genotype file (VCF format).',
        default=None,
        metavar='<genotype file name>',
        dest='geno_file_name',
        type=str,
        action='store'
        )

    groupCalc.add_argument(
        '-z',
        '--z_Score',
        help='Value of the z-score.',
        default=None,
        metavar='<z-Score>',
        dest='z_Score',
        type=float,
        action='store'
        )

    groupCalc.add_argument(
        '-r',
        '--raw_PRS',
        help='Value of the raw PRS.',
        default=None,
        metavar='<raw PRS>',
        dest='raw_PRS',
        type=float,
        action='store'
        )

    parser.add_argument(
        '-s',
        '--sample',
        help=('Name of sample in the genotype file to be used.\nIt can only'
              ' be used in conjunction with the -g option.'),
        default=None,
        metavar='<sample name>',
        dest='sample_name',
        type=str,
        action='store'
        )

    parser.add_argument(
        '-o',
        '--output',
        help=('Name of the file where results should be stored.\n'
              'It is only used in conjunction with the -g option, if'
              ' the sample is not specified.'),
        default=None,
        metavar='<results file name>',
        dest='results_file_name',
        type=str,
        action='store'
        )

    # Parse the command line and extract the names of the snp file, vcf file and sample name
    prs_file = parser.parse_args().prs_file_name[0]
    ancestry = parser.parse_args().ancestry_code
    props_file = parser.parse_args().props_file_name
    geno_file = parser.parse_args().geno_file_name
    res_file = parser.parse_args().results_file_name
    sample = parser.parse_args().sample_name
    raw_PRS = parser.parse_args().raw_PRS
    z_Score = parser.parse_args().z_Score

    if props_file is not None:
        print('Ancestry proportions provided; options 3 and 4 disabled')
        ancestry=-1

    if (geno_file is None and props_file is None and sample is not None):
        parser.print_help()
        sys.exit(22)


    # Create an object of the vcf2prs class
    try:
        prs = vcf2prs.Prs(prs_file=prs_file, ancestry=ancestry)
    except vcf2prs.Vcf2PrsError as err:
        print(err, file=sys.stderr)
        sys.exit(22)

    #(if possible): calculate rawPRS and assign it to PRS object; if input by the user, assign it directly
    if geno_file is not None:
        prs.calculate_prs_from_vcf(geno_file, sample, 'GT', 'DS')
    elif raw_PRS is not None:
        prs.raw_PRS = [raw_PRS]

    if (ancestry<0):
        # Case 1, mixed ancestry -rawPRS in the proportions file
        prs.calculate_mixed_prs(props_file, sample)
    elif hasattr(prs, 'raw_PRS'):
        # Case 1 (other mixed ancestry), Case 1 (single ancestry), Case 2
        prs.calculate_z_from_raw(prs.raw_PRS)
        prs.alpha = [prs.alpha] * len(prs.raw_PRS)
    elif z_Score is not None:
        # Case 3 - Calculate the raw PRS from the z-score
        prs.calculate_raw_from_z([z_Score])
        prs.alpha = [prs.alpha] * len(prs.raw_PRS)


    #Print results
    if hasattr(prs, 'z_Score'):
        # Cases 1, 2, 3
        if len(prs.raw_PRS) == 1:
            results = resFormat1.format(prs.raw_PRS[0], prs.alpha[0], prs.z_Score[0])
            print(results)
        elif res_file is None:
            print("sample\t raw PRS\t alpha\t\t z-score")
            for (tName, tRaw, tAlpha, tZed) in zip(prs.genoNames, prs.raw_PRS, prs.alpha, prs.z_Score):
                print(f"{tName}\t{tRaw: #.6}\t{tAlpha:#.6}\t{tZed: #.6}")
        else:
            with open(res_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["sample", "raw_PRS", "alpha", "z-score"])
                for (tName, tRaw, tAlpha, tZed) in zip(prs.genoNames, prs.raw_PRS, prs.alpha, prs.z_Score):
                    writer.writerow([tName, tRaw, tAlpha, tZed])

    else:
        # Case 4 - Calculate the mean and sd of the PRS
        if len(prs.meanList) > 1:
            print("Parameters for the all ancestries in the file.\n mean\t st. dev\t alpha")
            for (tMean, tSD, tAlpha) in zip(prs.meanList, prs.sdList, prs.alphaList):
                print(f"{tMean}\t{tSD: #.6}\t{tAlpha:#.6}")
            print("\nParameters for the selected ancestry (default=first values in the list)")
        results = resFormat2.format(prs.mean, prs.sd, prs.alpha)
        print(results)

