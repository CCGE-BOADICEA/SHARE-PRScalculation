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
        eafAncestry = ['unknown','eaf_eur', 'eaf_afr','eaf_eas','eaf_sas']
        if (ancestry>0):
            titles['eaf'] = eafAncestry[ancestry]

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
    def extract_coeffList(line, chCoeff, essFlag=0, lowLim=None, uppLim=None):
        """
        Method to extract the values of coefficients from the first lines of the file.

        Args:
            line (str): One of the uncommented line of the PRS file. E.g, line should be of one of the forms:
                1) 0.45 #alpha/mean/sd/threshold
                2) alpha/mean/sd/threshold = 0.45
                Blank spaces are ignored.
            chCoeff: name of the coefficient to be extracted
            essFlag (-1,0,1): are those coefficients essential or can be recalculated?
            lowLim, uppLim (float): lower and upper limits

        Returns:
            resList (list of floats): the value(s) in the coefficient list.
            XOR [888.888] for signalling that chCoeff wasn't found
            XOR [999.999] for signalling that it was found but something went wrong

        Raises:
            Vcf2PrsError: If values can not be extracted or are invalid.
        """

        if not isinstance(line, str):
            raise Vcf2PrsError('extract_coeffList was passed an invalid type: '
                               f'"{line}" of type {type(line)}.')
        if chCoeff not in line.lower():
            if (essFlag==1):
                raise Vcf2PrsError(f'extract_coeffList was unable to find "{chCoeff}".')
            elif (essFlag==0):
                print(f'extract_coeffList was unable to find "{chCoeff}". They will be recalculated')
            return [888.888]
        
        stripped = line
        if '=' in stripped:
            stripped = stripped.split('=', 1)[1]
        if '#' in stripped:
            stripped = stripped.split('#', 1)[0]
        stripped = stripped.split(',')

        coeffN = len(stripped)
        resList = [999.999] * coeffN
        for i in range(coeffN):
            try:
                coeff = float(stripped[i])
            except:
                if (essFlag==1):
                    raise Vcf2PrsError(f'the {i+1}th coefficient was not a float.')
                elif (essFlag==0):
                    print(f'the {i+1}th coefficient was not a float. It will be recalculated')
                continue
            if (lowLim is not None) and (coeff < lowLim):
                if (essFlag==1):
                    raise Vcf2PrsError(f'value "{coeff}" should be >= "{lowLim}".')
                elif (essFlag==0):
                    print(f'value "{coeff}" should be >= "{lowLim}". It will be recalculated')
                continue
            if (uppLim is not None) and (coeff > uppLim):
                if (essFlag==1):
                    raise Vcf2PrsError(f'value "{coeff}" should be <= "{uppLim}".')
                elif (essFlag==0):
                    print(f'value "{coeff}" should be <= "{uppLim}". It will be recalculated')
                continue
            resList[i] = coeff

        return resList

    def extract_coeff(coeffList, ancestry, essFlag=False):
        """
        Method to extract the value of the coefficient from its list.

        Args:
            coeffList (list): the list of possible values
            ancestry (integer): the specific ancestry chosen
            essFlag (-1,0,1): is this coefficient essential or can be recalculated?

        Returns:
            float: the value of the coefficient.

        Raises:
            Vcf2PrsError: If value can not be extracted or is invalid.
        """

        try:
            coeff = coeffList[ancestry]
        except (IndexError):
            raise Vcf2PrsError('Selected ancestry not available for this set')
        if not isinstance(coeff, float):
            if essFlag:
                raise Vcf2PrsError('extract_coeff was unable to extract value.')
            else:
                print('extract_coeff was unable to extract value. It will be recalculated')
                coeff = 999.999
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
        try:
            tmp = prs_data.__next__(); 
            while tmp[0] == "#": tmp = prs_data.__next__()
            essFlag = 1 if (ancestry==-1) else -1   #threshold is only relevant for mixed ancestry calculations
            self.threshList = PrsInfo.extract_coeffList(tmp, 'threshold', essFlag, 0, 1)
            if (self.threshList[0] != 888.888):
                tmp = prs_data.__next__()
                while tmp[0] == "#": tmp = prs_data.__next__()
            essFlag = 1 if (ancestry==-1) else 0   #mean and sd cannot be recalculated in mixed ancestry calculations
            self.meanList = PrsInfo.extract_coeffList(tmp, 'mean', essFlag, None, None)
            if self.meanList[0] != 888.888:
                tmp = prs_data.__next__()
                while tmp[0] == "#": tmp = prs_data.__next__()
            self.sdList = PrsInfo.extract_coeffList(tmp, 'sd', essFlag, 0, None)
            if self.sdList[0] != 888.888:
                tmp = prs_data.__next__()
                while tmp[0] == "#": tmp = prs_data.__next__()
            self.alphaList = PrsInfo.extract_coeffList(tmp, 'alpha', 1, 0, 1)
            snp_data = csv.DictReader(prs_data)
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
                    tmpList = [float(row[i]) for i in ['eaf_eur', 'eaf_afr','eaf_eas','eaf_sas']]
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
        
        # Extract or calculate alpha, mean and standard deviation for the set of variants.
        if (ancestry>=0):
            self.alpha = PrsInfo.extract_coeff(self.alphaList,ancestry,True)
            self.mean = PrsInfo.extract_coeff(self.meanList,ancestry,False)
            if (self.mean == 888.888) or (self.mean == 999.999):
                self.mean = sum([s.mean for s in self.snps.values()])      #mostly useless; we generally read it from the PRS file
            self.sd = PrsInfo.extract_coeff(self.sdList,ancestry,False)
            if (self.sd == 888.888) or (self.sd == 999.999):
                self.sd = (sum([s.var for s in self.snps.values()]))**0.5  #mostly useless; we generally read it from the PRS file
        return
