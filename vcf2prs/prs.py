"""
:authors: Andrew Lee, Lorenzo Ficorella
:organization: CCGE, DPHPC, University of Cambridge
:copyright:  2022 Cambridge University. All rights reserved.
:   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
:   without even the implied warranty of MERCHANTABILITY or FITNESS for a particular purpose.
:contact: lf364@medschl.cam.ac.uk
"""
import vcf
import io
import numpy as np
from vcf2prs.exception import Vcf2PrsError
from vcf2prs.snp import Snp
from vcf2prs.prsinfo import PrsInfo


class Prs(object):
    """
    Class to calculate PRS from a vcf file.

    Attributes:
        alphaList (list of floats)
        alpha (float): The alpha of the PRS.
        threshList (list of floats)
        meanList (list of floats)
        mean (float): The mean of the PRS.
        sdList (list of floats)
        sd (float): The sd of the PRS.
        raw_PRS (float): The value of the raw PRS.
        z_Score (float): The value of the z-Score PRS.
        lorList (list of floats)
        eafList (list of floats)
        genoNames (list of strings): sample names in
            the genotype file
    """
    __slots__ = "_alphaList", "_alpha", "_threshList", \
                "_meanList", "_mean", "_sdList", "_sd", \
                "_prs_info", "_raw_PRS", "_z_Score", \
                '_lorList', '_eafList', '_genoNames'

    def __init__(self, prs_file, ancestry=0, prs_columns={}):
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
            Vcf2PrsError: If an error occurs.

        """

        self._prs_info = PrsInfo(prs_file, ancestry, prs_columns )
        return

    @property
    def alphaList(self):
        return self._prs_info.alphaList
# setter not required here: it's only used, but it's set in prsinfo.py
#    @alphaList.setter
#    def alphaList(self, alphaList):
#        self._prs_info.alphaList = alphaList

    @property
    def alpha(self):
        return self._prs_info.alpha
    @alpha.setter
    def alpha(self, alpha):
        self._prs_info.alpha = alpha

    @property
    def threshList(self):
        return self._prs_info.threshList
# setter not required here: it's only used, but it's set in prsinfo.py
#    @threshList.setter
#    def threshList(self, threshList):
#        self._prs_info.threshList = threshList

    @property
    def meanList(self):
        return self._prs_info.meanList
# setter not required here: it's only used, but it's set in prsinfo.py
#    @meanList.setter
#    def meanList(self, meanList):
#        self._prs_info.meanList = meanList
    @property
    def mean(self):
        return self._prs_info.mean
# setter not required here: it's only used, but it's set in prsinfo.py
#    @mean.setter
#    def mean(self, mean):
#        self._prs_info.mean = mean

    @property
    def sdList(self):
        return self._prs_info.sdList
# setter not required here: it's only used, but it's set in prsinfo.py
#    @sdList.setter
#    def sdList(self, sdList):
#        self._prs_info.sdList = sdList
    @property
    def sd(self):
        return self._prs_info.sd
# setter not required here: it's only used, but it's set in prsinfo.py
#    @sd.setter
#    def sd(self, sd):
#        self._prs_info.sd = sd

    @property
    def lorList(self):
        return self._prs_info.lorList
# setter not required here: it's only used, but it's set in prsinfo.py
#    @lorList.setter
#    def lorList(self, lorList):
#        self._prs_info.lorList = lorList

    @property
    def eafList(self):
        return self._prs_info.eafList
# setter not required here: it's only used, but it's set in prsinfo.py
#    @eafList.setter
#    def eafList(self, eafList):
#        self._prs_info.eafList = eafList

    @property
    def genoNames(self):
        return self._genoNames
    @genoNames.setter
    def genoNames(self, genoNames):
        self._genoNames = genoNames

    @property
    def raw_PRS(self):
        return self._raw_PRS
    @raw_PRS.setter
    def raw_PRS(self, raw_PRS):
        if not isinstance(raw_PRS, list):
            raise Vcf2PrsError("Invalid raw_PRS passed to Prs object, "
                               f"of type {type(raw_PRS)}.")
        else:
            for tRaw in raw_PRS:
                if not isinstance(tRaw, float):
                    raise Vcf2PrsError("Invalid rawPRS contained in raw_PRS list, "
                                       f"of type {type(tRaw)}.")
        self._raw_PRS = raw_PRS

    @property
    def z_Score(self):
        return self._z_Score
    @z_Score.setter
    def z_Score(self, z_Score):
        if not isinstance(z_Score, list):
            raise Vcf2PrsError("Invalid z_Score passed to Prs object, "
                               f"of type {type(z_Score)}.")
        else:
            for tZed in z_Score:
                if not isinstance(tZed, float):
                    raise Vcf2PrsError("Invalid zScore contained in z_score list, "
                                       f"of type {type(tZed)}.")
        self._z_Score = z_Score

    @staticmethod
    def validate_record(record, snp, gt='GT', ds='DS'):
        """
        Method to validate the vcf record.

        Args:
            record (vcf.model._Record): a vcf record.
            snp (Snp): The reference variant from the prs file.
            gt (str, optional): the key for genotypes in the VCF record format
                statement.
            ds (str, optional): the key for dosages in the VCF record format
                statement.

        Raises:
            Vcf2PrsError: If the record is not valid.

        """
        # Check that the arguments are of the correct type
        if not isinstance(gt, str):
            raise Vcf2PrsError('Error: The method Vcf2Prs.validate_record has '
                               'been passed a genotype key, {0} of the '
                               'incorrect type, {1}. Should be of type str.'.
                               format(gt, type(gt)))
        if not isinstance(ds, str):
            raise Vcf2PrsError('Error: The method Vcf2Prs.validate_record has '
                               'been passed a dosage key, {0} of the '
                               'incorrect type, {1}. Should be of type str.'.
                               format(ds, type(ds)))

        # Check that the record is of the correct type
        if not isinstance(record, vcf.model._Record):
            raise Vcf2PrsError('Error: The method Vcf2Prs.validate_record has '
                               'been passed a record of the incorrect type, '
                               '{0}. Should be of type "vcf.model._Record".'.
                               format(type(record)))

        # Check the chromosome and position.
        Snp.convert_chromosome(record.CHROM)
        Snp.convert_position(record.POS)

        # Check that the reference allele is valid
        Snp.validate_allele_bases(record.REF)

        # Validate the alternative alleles
        if not isinstance(record.ALT, list):
            raise Vcf2PrsError('Error: In the method Vcf2Prs.validate_record '
                               'for the record on chromosome {0} at position '
                               '{1}, the alternative alleles '
                               'is not of type list, but of type {2}'.
                               format(record.CHROM, record.POS,
                                      type(record.ALT)))
        if record.ALT != [None]:
            for alt in record.ALT:
                if not isinstance(alt, vcf.model._Substitution):
                    raise Vcf2PrsError('Error: In the method '
                                       'Vcf2Prs.validate_record for the record'
                                       ' on chromosome {0} at position {1}'
                                       ', the alternative allele "{2}" is not '
                                       'of type "vcf.model._Substitution", but'
                                       ' of type {3}'.
                                       format(record.CHROM, record.POS, alt,
                                              type(alt)))
                Snp.validate_allele_bases(str(alt))

        # Check that the reference allele is not in the alternative alleles
        if str(record.REF).upper() in [str(i).upper() for i in record.ALT]:
            raise Vcf2PrsError('Error: In Vcf2Prs.validate_record for the '
                               'record on chromosome {0} at position {1}, '
                               'the reference allele, {2} is in the list '
                               'of alternative alleles {3}'.
                               format(record.CHROM, record.POS, record.REF,
                                      record.ALT))

        # Check that the SNP is correct
        if not isinstance(snp, Snp):
            raise Vcf2PrsError('Invalid reference variant passed to '
                               'extract_count.')

        # Check that the reference alleles are the same
        if str(record.REF).upper() != snp.ref.upper():
            raise Vcf2PrsError('Error: in the method Vcf2Prs.extract_count for'
                               ' VCF record on chromosome {0} at position {1},'
                               ' the reference allele for the PRS {2} does not'
                               ' match the VCF reference allele, {3}.'.
                               format(record.CHROM, record.POS, snp.ref,
                                      record.REF))

        # Check that the effect allele is in the list of alternative alleles
        vcf_alt_alleles = [str(i).upper() for i in record.ALT]
        if snp.eff.upper() not in vcf_alt_alleles:
            raise Vcf2PrsError('Error: in the method Vcf2Prs.extract_count for'
                               ' VCF record on chromosome {0} at position {1},'
                               ' the effect allele for the PRS {2} is not in '
                               'the VCF alternative alleles, {3}.'.
                               format(record.CHROM, record.POS, snp.eff,
                                      record.ALT))

        return

    @staticmethod
    def validate_sample(record, sample, gt='GT', ds='DS'):
        """
        Method to validate the vcf record.

        Args:
            record (vcf.model._Record): a vcf record.
            sample (str): the name of the sample
            gt (str, optional): the key for genotypes in the VCF record format
                statement.
            ds (str, optional): the key for dosages in the VCF record format
                statement.

        Raises:
            Vcf2PrsError: If the record is not valid.

        """
        # Check that the arguments are of the correct type
        if not isinstance(gt, str):
            raise Vcf2PrsError('Error: The method Vcf2Prs.validate_record has '
                               'been passed a genotype key, {0} of the '
                               'incorrect type, {1}. Should be of type str.'.
                               format(gt, type(gt)))
        if not isinstance(ds, str):
            raise Vcf2PrsError('Error: The method Vcf2Prs.validate_record has '
                               'been passed a dosage key, {0} of the '
                               'incorrect type, {1}. Should be of type str.'.
                               format(ds, type(ds)))

        # Check that the record is of the correct type
        if not isinstance(record, vcf.model._Record):
            raise Vcf2PrsError('Error: The method Vcf2Prs.validate_record has '
                               'been passed a record of the incorrect type, '
                               '{0}. Should be of type "vcf.model._Record".'.
                               format(type(record)))

        # Check that the sample is in the genotypes
        try:
            record.genotype(sample)
        except (KeyError, IndexError):
            raise Vcf2PrsError('Error: In the method Vcf2Prs.validate_record '
                               'for the record on chromosome {0} at position '
                               '{1}, the sample, {2} is not present in the '
                               'genotypes.'.
                               format(record.CHROM, record.POS, sample))

        # Check that the sample has some information
        if not hasattr(record.genotype(sample).data, gt):
            raise Vcf2PrsError('Error: In the method Vcf2Prs.validate_record '
                               'for the record on chromosome {0} at position '
                               '{1}, the sample, {2} is not valid'.
                               format(record.CHROM, record.POS, sample))

        # Check that the sample has either a called genotype or a dosage
        if not ((hasattr(record.genotype(sample).data, gt) and
                 record.genotype(sample).called) or
                hasattr(record.genotype(sample).data, ds)):
            raise Vcf2PrsError('Error: In the method Vcf2Prs.validate_record '
                               'for the record on chromosome {0} at position '
                               '{1}, the sample, {2} is not has neither a '
                               'genotype or a dosage.'.
                               format(record.CHROM, record.POS, sample))

        # If the sample has a genotype and it is called, validate it.
        if hasattr(record.genotype(sample).data, gt) and \
           record.genotype(sample).called:

            # Check that the genotype is of the correct type
            if not isinstance(record.genotype(sample)[gt], str):
                raise Vcf2PrsError('Error: In the method Vcf2Prs.validate'
                                   'Record for the record on chromosome {0} at'
                                   ' position {1}, the genotype, {2}, of '
                                   'sample, {3}, is not of type str, but of '
                                   'type {4}'.
                                   format(record.CHROM, record.POS,
                                          record.genotype(sample)[gt],
                                          sample,
                                          type(record.genotype(sample)[gt])))

            # Parse the genotype
            genotypes = record.genotype(sample)[gt].replace('/', '|')
            genotypes = genotypes.split('|')

            # Check that two genotypes are present.
            # In principle, for male on chromosme X or Y there could be only 1.
            # This may be the case for prostate cancer.
            if len(genotypes) != 2:
                raise Vcf2PrsError('Error: In the method Vcf2Prs.validate'
                                   'Record for the record on chromosome {0} at'
                                   ' position {1}, the sample, {2} has an '
                                   'invalid genotype, {3}.'.
                                   format(record.CHROM, record.POS, sample,
                                          record.genotype(sample)[gt]))

            # Determine the list of allowed genotypes.
            allowd_genotyes = ['0']
            if record.ALT != [None]:
                altN = len(record.ALT)
                allowd_genotyes = [str(i) for i in range(altN+1)]

            # Check that each genotype is in the allowed genotypes.
            for g in genotypes:
                if g not in allowd_genotyes:
                    raise Vcf2PrsError('Error: In the method Vcf2Prs.validate'
                                       'Record for the record on chromosome '
                                       '{0} at position {1}, the genotype for '
                                       'sample {2}, {3} is invalid.'.
                                       format(record.CHROM, record.POS, sample,
                                              record.genotype(sample)[gt]))
            genotype_dosage = genotypes.count('1')

        # If the sample has a dosage, validate it.
        if hasattr(record.genotype(sample).data, ds) and \
           record.genotype(sample)[ds] is not None:
            # Check that the dosage is of the correct type
            if not isinstance(record.genotype(sample)[ds], (float, int)):
                raise Vcf2PrsError('Error: In the method Vcf2Prs.validate'
                                   'Record for the record on chromosome {0} at'
                                   ' position {1}, the dosage for sample {2}, '
                                   '{3}, is not of type float, but of type '
                                   '{4}'.
                                   format(record.CHROM, record.POS, sample,
                                          record.genotype(sample)[ds],
                                          type(record.genotype(sample)[ds])))

            # Check that only one alternative allele is present
            if len(record.ALT) != 1:
                raise Vcf2PrsError('Error: In the method Vcf2Prs.validate'
                                   'Record for the record on chromosome {0} at'
                                   ' position {1}, trying to use dosages, but '
                                   'multiple alternative alleles are '
                                   'specified, {2}'.
                                   format(record.CHROM, record.POS,
                                          record.ALT))

            # Check that the dosage is within bounds
            if record.genotype(sample)[ds] < 0.0 or \
               record.genotype(sample)[ds] > 2.0:
                raise Vcf2PrsError('Error: In the method Vcf2Prs.validate'
                                   'Record for the record on chromosome {0} at'
                                   ' position {1}, for sample {2}, the dosage '
                                   '{3}, is out of bounds. Doages should be '
                                   'within the bounds 0<=Doagae<=2'.
                                   format(record.CHROM, record.POS, sample,
                                          record.genotype(sample)[ds]))

            # If no alternative allele is specified, the dosage must be 0.0
            if record.ALT == [None] and \
               abs(record.genotype(sample)[ds]) > 1e-4:
                raise Vcf2PrsError('Error: In the method Vcf2Prs.validate'
                                   'Record for the record on chromosome {0} at'
                                   ' position {1}, trying to use dosages, but '
                                   'no alternative allele is specified and the'
                                   ' dosage is not eqaul to 0'.
                                   format(record.CHROM, record.POS))

        # If the sample has both a genotype and a dosage, check that they are
        # consistent.
        if hasattr(record.genotype(sample).data, gt) and \
           record.genotype(sample).called and \
           hasattr(record.genotype(sample).data, ds) and \
           record.genotype(sample)[ds] is not None and \
           abs(float(genotype_dosage) - record.genotype(sample)[ds]) > 1e-4:
            raise Vcf2PrsError('Error: In the method Vcf2Prs.validate_record '
                               'for the record on chromosome {0} at position '
                               '{1}, for sample {2} the dosage {3}, and '
                               'genotype {4} are inconsistent.'.
                               format(record.CHROM, record.POS, sample,
                                      record.genotype(sample)[ds],
                                      record.genotype(sample)[gt]))

        return

    @staticmethod
    def extract_count(record, sample, snp, gt='GT', ds='DS'):
        """
        Method to calculate the number of copies of the effect allele in the
        genotype from the vcf file.

        Args:
            record (vcf.model._Record): a vcf record.
            sample (str): the name of the sample
            snp (Snp): The reference variant from the prs file.
            gt (str, optional): the key for genotypes in the VCF record format
                statement.
            ds (str, optional): the key for dosages in the VCF record format
                statement.

        Returns:
            float: The number of copies of the effect allele in the genotype.
                The number is [0,2].

        Raises:
            Vcf2PrsError: If the calculation is not possible.

        """

        Prs.validate_sample(record, sample, gt, ds)

        # Use genotyping information
        if (hasattr(record.genotype(sample).data, gt) and
                record.genotype(sample).called):

            # Determine the list of altenative alleles
            vcf_alt_alleles = [str(i).upper() for i in record.ALT]

            # Determine numeric symbol for the effect allele.
            eff_allele_symbol = str(vcf_alt_alleles.index(snp.eff.upper()) + 1)

            # Determine the list of numeric genotypes.
            genotypes = record.genotype(sample)[gt].replace('/', '|')
            genotypes = genotypes.split('|')

            # Count the number of effect alleles
            return float(genotypes.count(eff_allele_symbol))
        else:
            if record.genotype(sample)[ds] <= 1e-4:
                return 0.0
            else:
                return record.genotype(sample)[ds]

    def calculate_prs_from_vcf(self, geno_file, sample=None, gt='GT', ds='DS'):
        """
        Method to calculate the raw PRS from the genotypes from a VCF file.

        Args:
            geno_file (str, io.StringIO or io.TextIOWrapper): A string
                containing the filepath to a VCF file or a file-like object
                containing the contents of a VCF file.
            sample (str, optional): the name of the sample in the genotype
                file to be used to calculate the PRS.
            gt(str, optional): the key for genotypes in the VCF record format
                statement.
            ds(str, optional): the key for dosages in the VCF record format
                statement.

        Raises:
            Vcf2PrsError: If an error occurs.

        """

        # Check that the genotype file is of the correct type
        if not isinstance(geno_file, (str, io.StringIO, io.TextIOWrapper)):
            raise Vcf2PrsError("Invalid geno_file passed to Prs object, "
                               f"of type {type(geno_file)}.")

        # Check that the sample is of the correct type
        if sample is not None and not isinstance(sample, (str)):
            raise Vcf2PrsError("Invalid sample name passed to Prs object, "
                               f"of type {type(sample)}.")

        # Use the genotype file as a vcf Reader object.
        try:
            if isinstance(geno_file, str):
                fsock = open(geno_file, 'r')
                geno_data = vcf.Reader(fsock)
            else:
                geno_data = vcf.Reader(geno_file)
        except Exception as ex:
            if 'fsock' in locals():
                fsock.close()
            raise Vcf2PrsError('Unable to open the VCF file '
                               f'{type(ex).__name__}: {ex.args}')

        # Determine if the sample is in the file, and raise an exception if it is not.
        if len(geno_data.samples) < 1:
            geno_data._reader.close()
            raise Vcf2PrsError('Error: No samples found in the genotype file.')
        if (sample is None):
            self.genoNames = geno_data.samples
            genoN = len(geno_data.samples)
        else:
            if (sample not in geno_data.samples):
                geno_data._reader.close()
                raise Vcf2PrsError(f'Error: the sample "{sample}" is not '
                                   'in the genotype file.')
            if geno_data.samples.count(sample) != 1:
                geno_data._reader.close()
                raise Vcf2PrsError('Error: the sample "{sample}" is in '
                               'the genotype file more than once.')
            self.genoNames = [sample]
            genoN = 1


        # Reset the found property for each snp
        self._prs_info.reset_found()
        # Initialise the raw_PRS to 0, so that it can be incremented.
        self.raw_PRS = [0.0] * genoN
        try:
            # Loop through the variants in the file, determine if it is in the
            # dictionary of snps, and if so calculate the raw_PRS, otherwise
            # ignore.
            for record in geno_data:
                # Convert the chromosome and position to standatd form, so that
                # they can be used as a key.
                chromo = Snp.convert_chromosome(record.CHROM)
                pos = Snp.convert_position(record.POS)

                # Check if the record is in the set of snps for the PRS
                if (chromo, pos) in self._prs_info.snps:
                    snp = self._prs_info.snps[(chromo, pos)]

                    # Check if the variant has already been found, if so raise
                    # an exception.
                    if snp.found:
                        raise Vcf2PrsError('Error: More than one occurance of '
                                           f'a variant on chromosome {chromo} '
                                           f'and position position {pos} in '
                                           'the geno_file.')
                    # Mark the variant as found
                    snp.found = True

                    # Validate the record
                    Prs.validate_record(record, snp, gt, ds)

                    for sc in range(genoN):
                        try:
                            # Count the number of copies of the effect allele.
                            count = Prs.extract_count(record, self.genoNames[sc], snp, gt, ds)
                            # Increment the raw_PRS
                            self.raw_PRS[sc] += count * snp.lor
                        except:
                            self.raw_PRS[sc] = float('nan')

        except IndexError:
            raise Vcf2PrsError('Error: Invalid record in the geno_file.')
        finally:
            geno_data._reader.close()

        # Check that all variants have been found.
        for chrpos, snp in self._prs_info.snps.items():
            if not snp.found:
                raise Vcf2PrsError(('Error: The variant at position {1} on '
                                    'chromosome {0} was not found in the '
                                    'genotype file.').format(*chrpos))

        return

    def calculate_raw_from_z(self, z_Score):
        """
        Method to calculate the raw polygenic load from the z score.

        Args:
            z_Score (float): the z score PRS.
        """
        self.z_Score = z_Score
        ## Convert the Z score to a raw PRS.
        self.raw_PRS = [(tZed * self.sd + self.mean) for tZed in z_Score]
        return

    def calculate_z_from_raw(self, raw_PRS):
        """
        Method to calculate the Z-score from the raw polygenic load.

        Args:
            raw_PRS (float): the raw PRS.
        """
        self.raw_PRS = raw_PRS
        ## Convert the Z score to a raw PRS.
        self.z_Score = [((tRaw - self.mean) / self.sd) for tRaw in raw_PRS]
        return

    def calculate_mixed_prs(self, props_file, sample=None):
        """
        Method to calculate the alpha and Z-score from the raw polygenic load,
               for individuals of mixed ancestry.

        Args:
            props_file (str, io.StringIO or io.TextIOWrapper): A string
                containing the filepath to a csv file
                containing the ancestry proportions.
            sample (str, optional): the name of the sample in the proportions
                file to be used to calculate the PRS.

        Raises:
            Vcf2PrsError: If an error occurs.

        """

        # Check that the proportion file is of the correct type
        if not isinstance(props_file, (str, io.StringIO, io.TextIOWrapper)):
            raise Vcf2PrsError("Invalid props_file passed to Prs object, "
                               f"of type {type(props_file)}.")

        # The first line contains the headers; open the rest as a csv reader
        try:
            tmpTab = np.genfromtxt(props_file, dtype=float, delimiter=',',skip_header=1)
            if len(tmpTab.shape) == 1:
                propsN = 1
                colsN  = tmpTab.shape[0]
                tmpTab = np.array([tmpTab])
            else:
                propsN = tmpTab.shape[0]
                colsN  = tmpTab.shape[1]
            if colsN==6:          #robust, but has to be updated if other ancestries are considered
                propNames = np.genfromtxt(props_file, dtype=str, delimiter=',',skip_header=1, usecols=0)
                tmpTab = tmpTab[:,1:]
            else:
                propNames = np.array(['row'+str(i) for i in range(1,propsN+1)])
        except (IOError, ValueError, UnicodeDecodeError, StopIteration):
            raise Vcf2PrsError('Unable to read the props file')

        propPRS = tmpTab[:,0]
        propDistr = tmpTab[:,1:5]
        for idx,row in enumerate(propDistr):
            if any(row<0):
                raise Vcf2PrsError('Negative proportions!')
            if not( 0.99< np.sum(row) <1.01 ):
                raise Vcf2PrsError('Ancestry proportions do not add up to 1!')
            propDistr[idx,:] = row/np.sum(row)

        if sample is not None:
            propsN = 1
            # check that the sample is of the correct type
            if not (isinstance(sample, (str))):
                raise Vcf2PrsError("Invalid sample name passed to Prs object, "
                                    f"of type {type(sample)}.")
            # check it exists in the proportion file
            sampleInd = np.where(propNames == sample)[0]
            if len(sampleInd)==1:
                propPRS = propPRS[sampleInd]
                propDistr = propDistr[sampleInd,:]
            elif len(sampleInd)>1:
                raise Vcf2PrsError('Error: In the method Vcf2Prs.calculate_mixed_prs: '
                               'the sample, {0} is  present multiple times '
                               'in the proportions file.'. format(sample))
            else:
                raise Vcf2PrsError('Error: In the method Vcf2Prs.calculate_mixed_prs: '
                               'the sample, {0} is not present in the proportions file.'.
                               format(sample))

        if not hasattr(self, 'raw_PRS'):
            #if rawPRS wasn't otherwise available, retrieve it from the proportion file
            self.raw_PRS = propPRS.tolist()
            self.genoNames = propNames
        elif propsN != len(self.raw_PRS):
            # Check concordance between genotype and proportion files (if needed)
            raise Vcf2PrsError('Number of genotype samples and number of proportions do not match!')
        elif (propsN>1) and any(self.genoNames != propNames):
            # Check concordance between genotype and proportion files (if needed)
            raise Vcf2PrsError('Names of genotype samples and names of proportions do not match!')

        alphaList = np.asarray(self.alphaList[1:5], dtype=np.float32)
        threshList = np.asarray(self.threshList[1:5], dtype=np.float32)
        meanList = np.asarray(self.meanList[1:5], dtype=np.float32)
        sdList = np.asarray(self.sdList[1:5], dtype=np.float32)
        lorList = np.asarray(self.lorList, dtype=np.float32)
        eafList = np.asarray(self.eafList, dtype=np.float32)

        self.z_Score = [0.0] * propsN
        self.alpha = [0.0] * propsN
        for i in range(propsN):
            monoFlag = (propDistr[i,:] >= threshList)
            if np.any(monoFlag):
                # note: if the sample were truly mono-Ancestry, then the results from this if clause and those from the else clause would be identical
                monoIndex = (np.nonzero(monoFlag)[0]).item()
                self.z_Score[i] = (self.raw_PRS[i] - meanList[monoIndex]) / sdList[monoIndex]
                self.alpha[i] = alphaList[monoIndex]
            else:
                ## calculate the possible overall sigma
                SNPtmp1 = np.sum((eafList**2) * propDistr[i,:], axis=1)
                SNPtmp2 = np.sum(eafList * propDistr[i,:], axis=1)

                stDev = np.sum((sdList**2) * propDistr[i,:])
                stDev += 4 * np.sum((SNPtmp1 - (SNPtmp2**2)) * (lorList**2))
                stDev = np.sqrt(stDev)

                ## calculate the weighted alphaList and PRS
                tmp1 = np.sum(alphaList * propDistr[i,:] * (self.raw_PRS[i] - meanList) / sdList)
                tmp2 = np.sum(alphaList * propDistr[i,:] / sdList)
                tmp3 = np.sum(alphaList * propDistr[i,:] * sdList)

                self.z_Score[i] = (tmp1 / tmp2) / stDev
                self.alpha[i] = tmp3 / stDev

        return