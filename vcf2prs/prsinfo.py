"""
:authors: Andrew Lee, Lorenzo Ficorella
:organization: CCGE, DPHPC, University of Cambridge
:copyright:  2022 Cambridge University. All rights reserved.
:   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
:   without even the implied warranty of MERCHANTABILITY or FITNESS for a particular purpose.
:contact: lf364@medschl.cam.ac.uk
"""
import csv
import io
from vcf2prs.exception import Vcf2PrsError
from vcf2prs.snp import Snp


class PrsInfo(object):
    """
    Class to hold to PRS info.

    Attributes:
        alphaList (list of floats)
        alpha (float): the square root of the proportion of the overall
            polygenic variance explained by the PRS. A real number (0, 1).
        threshList (list of floats): the threshold of ancestry proportion 
            after which an individual is considered of single ancestry
        meanList (list of floats)
        mean (float): the mean polygenic load explained by the variants in
            the file.
        sdList (list of floats)
        sd (float): the standard deviation of the polygenic load of the
            variants in the file.
        lorList (list of floats): the log-odds ratios for the SNPs in the PGS model
        eafList (list of floats): the ancestry-specific effect allele 
            frequencies in the PGS model
        snps (dict(Snps)): a dictionary of the variants contained in the
            file, where the keys are the tuples (chromosome, position) (in
            standarised form).

    """
    __slots__ = '_alphaList','_alpha', '_threshList', '_thresh', \
                '_meanList', '_mean', '_sdList', '_sd', \
                '_lorList', '_eafList', '_snps'

    def __init__(self, prs_file, ancestry, columns={}):
        """
        Method to initialise the PRS model.

        Args:
            prs_file (str, io.StringIO or io.TextIOWrapper): A string
                containing the filepath to a PRS file or a file-like object
                containing the contents of a PRS file.
            ancestry (int, required): the number codifying for a specific ancestry
                (0 for population, -1 for mixed calculations)
            prs_columns (dictionary, optional): A dictionary containing the names
                of the columns in the PRS file with some or all of the keys:
                    - "chr": The title of the column containg the chromosomes
                    - "pos": The title of the column containg the positions
                    - "ref": The title of the column containg the reference alleles
                    - "eff": The title of the column containg the effect allleles
                    - "lor": The title of the column containg the log odds ratios
                    - "eaf": The title of the column containg the effect allele frequencies

        Raises:
            Vcf2PrsError: If the arguments are not a valid or if error when reading the PRS file.

        """

        titles = {'chr': 'Chromosome',
                  'pos': 'Position',
                  'ref': 'Reference_Allele',
                  'eff': 'Effect_Allele',
                  'lor': 'Log_Odds_Ratio',
                  'eaf': 'Effect_Allele_Frequency'}
        if   (ancestry==1) : titles['eaf'] = 'eaf_eur'
        elif (ancestry==2) : titles['eaf'] = 'eaf_afr'
        elif (ancestry==3) : titles['eaf'] = 'eaf_eas'
        elif (ancestry==4) : titles['eaf'] = 'eaf_sas'

        titleKeys = list(titles.keys())
        if not isinstance(columns, dict):
            raise Vcf2PrsError("Invalid columns passed to PrsInfo object.")
        titles.update(columns)
        self._read_prs_file(prs_file, ancestry, titles, titleKeys)
        return

    @property
    def alphaList(self):
        return self._alphaList
    @alphaList.setter
    def alphaList(self, alphaList):
        self._alphaList = alphaList
    @property
    def alpha(self):
        return self._alpha
    @alpha.setter
    def alpha(self, alpha):
        self._alpha = alpha

    @property
    def threshList(self):
        return self._threshList
    @threshList.setter
    def threshList(self, threshList):
        self._threshList = threshList

    @property
    def meanList(self):
        return self._meanList
    @meanList.setter
    def meanList(self, meanList):
        self._meanList = meanList
    @property
    def mean(self):
        return self._mean
    @mean.setter
    def mean(self, mean):
        self._mean = mean

    @property
    def sdList(self):
        return self._sdList
    @sdList.setter
    def sdList(self, sdList):
        self._sdList = sdList
    @property
    def sd(self):
        return self._sd
    @sd.setter
    def sd(self, sd):
        self._sd = sd

    @property
    def lorList(self):
        return self._lorList
    @lorList.setter
    def lorList(self, lorList):
        self._lorList = lorList
    @property
    def eafList(self):
        return self._eafList
    @eafList.setter
    def eafList(self, eafList):
        self._eafList = eafList

    @property
    def snps(self):
        return self._snps
    @snps.setter
    def snps(self, snps):
        if not isinstance(snps, dict):
            raise Vcf2PrsError("Invalid snps type passed to PrsInfo object. "
                               "The snps should be of type dict, but was "
                               f"passed type {type(snps)}.")
        self._snps = snps


    @staticmethod
    def extract_coeffList(line, chCoeff):
        """
        Method to extract the values of coefficients from the first lines of the file.

        Args:
            line (str): One of the uncommented line of the PRS file. E.g, line should be of one of the forms:
                1) 0.45
                2) alpha/mean/sd/threshold = 0.45
            Blank spaces are ignored.

        Returns:
            list: the value(s) in that line.

        Raises:
            Vcf2PrsError: If values can not be extracted.

        """
        if not isinstance(line, str):
            raise Vcf2PrsError('extract_coeffList passed an invalid type: '
                               f'"{line}" of type {type(line)}.')

        if chCoeff not in line:
            print(f'Unable to extract values for "{chCoeff}"; '
                   'if possible, they will be recalculated')
            coeffList = [999.999]

        else:
            stripped = line
            try:
                if '=' in stripped:
                    stripped = stripped.split('=', 1)[1]
                if '#' in stripped:
                    stripped = stripped.split('#', 1)[0]
                stripped = stripped.split(',')
                coeffList = [float(i) for i in stripped]
            except (TypeError, ValueError):
                print(f'Unable to extract values for "{chCoeff}"; '
                    'they will be recalculated')
                coeffList = [999.999]

        return coeffList

    def extract_coeff(coeffList,ancestry, lowLim=None, uppLim=None):
        """
        Method to extract the value of the coefficient from its list.

        Args:
            coeffList (list): the list of possible values

        Returns:
            float: the value of the coefficient.

        Raises:
            Vcf2PrsError: If value can not be extracted or is invalid.

        """
        try:
            coeff = coeffList[ancestry]
        except (TypeError, ValueError):
            print('Unable to extract value; if possible, it will be recalculated')
            coeff = 999.999
        except (IndexError):
            raise Vcf2PrsError('Selected ancestry not available for this set')

        # Check that the value is within limits lowLim < coeff < uppLim
        if lowLim is not None:
            if not (lowLim <= coeff):
                raise Vcf2PrsError(f'value should be greater than {lowLim}. The value, '
                                    f'{coeff}, is out of bounds.')
        if uppLim is not None:
            if not (uppLim >= coeff):
                raise Vcf2PrsError(f'value should be smaller than {uppLim}. The value, '
                                    f'{coeff}, is out of bounds.')

        return coeff

    def reset_found(self):
        for snp in self.snps.values():
            snp.found = False

    def _read_prs_file(self, prs_file, ancestry, columns, column_keys):
        """
        Method to read in the PRS information form the file.

        Args:
            prs_file (str, io.StringIO or io.TextIOWrapper): A string
                containing the filepath to a PRS file or a file-like object
                containing the contents of a PRS file.
            ancestry (int, required): the number codifying for a specific ancestry
                (0 for population, -1 for mixed calculations)
            columns (dictionary, optional): A dictionary containing the names
                of the columns in the PRS file with some or all of the keys:
                    - "chr": The title of the column containg the chromosomes
                    - "pos": The title of the column containg the positions
                    - "ref": The title of the column containg the reference
                             alleles
                    - "eff": The title of the column containg the effect
                             allleles
                    - "lor": The title of the column containg the log odds
                             ratios
                    - "eaf": The title of the column containg the effect allele
                            frequencies
            column_keys (list (str)): A list of the keys to the column names
                necessary for a snp object.

        Raises:
            Vcf2PrsError: If parameters can not be extracted or is invalid.

        """

        # Check that the prs file is of the correct type
        if not isinstance(prs_file, (str, io.StringIO, io.TextIOWrapper)):
            raise Vcf2PrsError("Invalid geno_file passed to PrsInfo object, "
                               f"of type {type(prs_file)}.")

        # Open the PRS file.
        if isinstance(prs_file, str):
            try:
                fsock = open(prs_file, 'r')
                prs_data = fsock
            except Exception as ex:
                raise Vcf2PrsError('Unable to open the PRS file '
                                   f'{type(ex).__name__}: {ex.args}')
        else:
            prs_data = prs_file

        # The first lines contains the coefficients for the PRS model in the file.
        # Open the rest of the file as a csv reader, and determine the
        # relevant columns from the headers.
        print(ancestry)
        tmpAnc = max(0, ancestry)
        try:
            tmp = prs_data.__next__(); 
            while tmp[0] == "#": tmp = prs_data.__next__()
            self.threshList = PrsInfo.extract_coeffList(tmp, 'threshold')
            if (self.threshList[0] != 999.999):
                tmp = prs_data.__next__()
                while tmp[0] == "#": tmp = prs_data.__next__()
            elif (ancestry==-1):
                raise Vcf2PrsError('No thresholds available for mixed calculations.')

            self.meanList = PrsInfo.extract_coeffList(tmp, 'mean')
            self.mean = PrsInfo.extract_coeff(self.meanList,tmpAnc)
            if self.mean != 999.999:
                tmp = prs_data.__next__()
                while tmp[0] == "#": tmp = prs_data.__next__()
            self.sdList = PrsInfo.extract_coeffList(tmp, 'sd')
            self.sd = PrsInfo.extract_coeff(self.sdList,tmpAnc,0)
            if self.sd != 999.999:
                tmp = prs_data.__next__()
                while tmp[0] == "#": tmp = prs_data.__next__()
            self.alphaList = PrsInfo.extract_coeffList(tmp, 'alpha')
            self.alpha = PrsInfo.extract_coeff(self.alphaList,tmpAnc,0,1)
            snp_data = csv.DictReader(prs_data)
        #if (ancestry<0):
        #    import numpy as np
        #    tmpTab = np.genfromtxt(prs_data, delimiter=',',skip_header=1,usecols = (4,6,7,8,9))
        #    lorList = tmpTab[:,0]
        #    eafList = tmpTab[:,1:]
        except (IOError, UnicodeDecodeError, StopIteration):
            if 'fsock' in locals():
                fsock.close()
            raise Vcf2PrsError('Unable to read the PRS file')

        # Check that the header is present in the file.
        if snp_data.fieldnames is None:
            if 'fsock' in locals():
                fsock.close()
            raise Vcf2PrsError('Unable to read the PRS file.')

        # Check that all the columns have been found.
        for key in column_keys:
            if columns[key] not in snp_data.fieldnames:
                if 'fsock' in locals():
                    fsock.close()
                raise Vcf2PrsError(f'The column "{columns[key]}" is not in the'
                                   ' PRS file.')

        try:
            # Loop through the lines of the file, and add each snp to the
            # dictionary, keyed by chromosome and position.
            self.snps = {}

            eafList = []
            for row in snp_data:
                # Use chromosome and position as the key for the snp object.
                c = Snp.convert_chromosome(row[columns['chr']])
                p = Snp.convert_position(row[columns['pos']])

                # Initialise the SNP
                self.snps[(c, p)] = Snp(**{key: row[columns[key]] for key in
                                           column_keys})
                if (ancestry<0):
                    tmpList = []
                    for col in ['eaf_eur', 'eaf_afr','eaf_eas','eaf_sas']:
                        tmpList.append(float(row[col]))
                    eafList.append(tmpList)
        except:
            raise Vcf2PrsError('The chosen PRS files does not contain '
                               'the EAFs for all required ancestries')
        finally:
            # Close the file
            if 'fsock' in locals():
                fsock.close()

        if len(self.snps) == 0:
            raise Vcf2PrsError('No variants found in the PRS file.')

        self.lorList = [s.lor for s in self.snps.values()]
        self.eafList = eafList
        
        # Calculate the mean and standard deviation for the set of variants.
        if self.mean == 999.999:
           self.mean = sum([s.mean for s in self.snps.values()])     #mostly useless; we generally read it from the PRS file
        if self.sd == 999.999:
           self.sd = (sum([s.var for s in self.snps.values()]))**0.5  #mostly useless; we generally read it from the PRS file
        return
