"""
:authors: Andrew Lee, Lorenzo Ficorella
:organization: CCGE, DPHPC, University of Cambridge
:copyright:  2022 Cambridge University. All rights reserved.
:   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
:   without even the implied warranty of MERCHANTABILITY or FITNESS for a particular purpose.
:contact: lf364@medschl.cam.ac.uk
"""
from vcf2prs.exception import Vcf2PrsError
from math import exp, log


class Snp(object):
    """
    Class to hold information on the variants in a PRS.
    Calculates the mean and variance of the variant.

    Args:
        chr (str): the chromosome on which the variant is located. Must be in
            Snp.CHROMOSOME_NAMES
        pos (str/int): the position of the variant on its chromosome - must
            be >=0
        ref (str): the reference allele in bases - must be composed of the
            characters in Snp.BASES.
        eff (str): the effect allele in bases - must be composed of the
            characters in Snp.BASES.
        lor (str/float): the log odds ratio.
        eaf (str/float): the effect allele frequency - must be (0,1).

    Attributes:
        chr (str): the chromosome on which the variant is located.
        pos (int): the position of the variant on its chromosome.
        ref (str): the reference allele in bases.
        eff (str): the effect allele in bases.
        lor (float): the log odds ratio.
        eaf (float): the effect allele frequency
        mean (float): the mean log odds ratio.
        var (float): the variance of the log odds ratio.
        found (bool): if the variant has been found for the sample.
        CHROMOSOME_NAMES (list(str)): The standarised chromosome names.
        BASES (list(str)): The allowed bases: A,T,C,G - all alleles must be
            composed of these.

    Raises:
        Vcf2PrsError: If the arguments are not valid.

    """
    CHROMOSOME_NAMES = ['X', 'Y'] + [str(i) for i in range(1,24)]
    BASES = ['A', 'T', 'C', 'G']

    __slots__ = ('_chr', '_pos', '_ref', '_eff', '_lor', '_eaf',
                 '_mean', '_var', '_found')

    def __init__(self, chr=None, pos=None, ref=None, eff=None,
                 lor=None, eaf=None):

        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.eff = eff
        self.lor = lor
        self.eaf = eaf
        self.found = False
        self._calculate_mean_and_variance()

        return

    @property
    def chr(self):
        """The chromosome"""
        return self._chr

    @chr.setter
    def chr(self, chr):
        self._chr = Snp.convert_chromosome(chr)

    @property
    def pos(self):
        """The position"""
        return self._pos

    @pos.setter
    def pos(self, pos):
        self._pos = Snp.convert_position(pos)

    @property
    def ref(self):
        """The reference allele"""
        return self._ref

    @ref.setter
    def ref(self, ref):
        if not Snp.validate_allele_bases(ref):
            raise Vcf2PrsError(f'Invalid reference allele {ref} passed to Snp '
                               'object.')
        self._ref = ref.upper()
        if (hasattr(self, '_eff') and (self.eff is not None) and
                (self.ref is not None) and (self.eff == self.ref)):
            raise Vcf2PrsError(f'Invalid reference {self.ref} and effect '
                               f'{self.eff} allele pair.')

    @property
    def eff(self):
        """The effect allele"""
        return self._eff

    @eff.setter
    def eff(self, eff):
        if not Snp.validate_allele_bases(eff):
            raise Vcf2PrsError(f'Invalid effect allele {eff} passed to Snp '
                               'object.')
        self._eff = eff.upper()
        if (hasattr(self, '_ref') and (self.eff is not None) and
                (self.ref is not None) and (self.eff == self.ref)):
            raise Vcf2PrsError(f'Invalid reference {self.ref} and effect '
                               f'{self.eff} allele pair.')

    @property
    def lor(self):
        """The log odds ratio"""
        return self._lor

    @lor.setter
    def lor(self, lor):
        if isinstance(lor, bool):
            raise Vcf2PrsError("Invalid log odds ratio passed to a Snp object,"
                               "of type bool.")
        try:
            self._lor = float(lor)
        except (TypeError, ValueError):
            raise Vcf2PrsError("Invalid log odds ratio passed to a Snp object."
                               f" Unable to convert {lor} to a float.")

    @property
    def eaf(self):
        """The effect allele frequency"""
        return self._eaf

    @eaf.setter
    def eaf(self, eaf):
        if isinstance(eaf, bool):
            raise Vcf2PrsError("Invalid effect allele frequency passed to a"
                               "Snp object, of type bool.")
        try:
            self._eaf = float(eaf)
        except (TypeError, ValueError):
            raise Vcf2PrsError("Invalid effect allele frequency passed to a "
                               f"Snp object. Unable to convert {eaf} to a "
                               "float.")

        # Check that eaf is within limits 0 <= eaf <= 1
        if not (0.0 <= self.eaf <= 1.0):
            raise Vcf2PrsError(f"Invalid effect allele frequency ({self.eaf})"
                               " passed to a Snp object. It must be in the"
                               " interval [0, 1]")

    @property
    def mean(self):
        """The mean log odds ratio"""
        return self._mean

    @mean.setter
    def mean(self, mean):
        if not isinstance(mean, float):
            raise Vcf2PrsError("Invalid mean passed to Snp object.")
        self._mean = mean

    @property
    def var(self):
        """The variance of the log odds ratio"""
        return self._var

    @var.setter
    def var(self, var):
        if not isinstance(var, float):
            raise Vcf2PrsError("Invalid variance passed to Snp object.")
        if var < -0.000000000000001:
            raise Vcf2PrsError("Invalid variance (less than 0) passed to Snp"
                               "object.")
        self._var = var

    @property
    def found(self):
        """If the snp has been found for the sample."""
        return self._found

    @found.setter
    def found(self, found):
        if not isinstance(found, bool):
            raise Vcf2PrsError("Invalid found passed to Snp object.")
        self._found = found

    @staticmethod
    def convert_chromosome(chromo):
        """
        Function to convert the name of the chromosome to standard format, and
        check that it is valid.

        Args:
            chromo (str/int): The name of the chromosome.

        Returns:
            str: The element of Snp.CHROMOSOME_NAMES corresponding to the
                input i.e. a string containing the name of the chromosome in
                standard format.

        Raises:
            Vcf2PrsError: If the argument cannot be converted to the standard
                format.

        """

        if (not isinstance(chromo, (str, int))) or (type(chromo) is bool):
            raise Vcf2PrsError('Error: the method convert_chromosome has'
                               ' been passed a non-string/integer argument, '
                               '"{0}", of type {1}'.
                               format(chromo, type(chromo)))

        chromosome = str(chromo).lower()
        chromosome = chromosome.replace('chromosome', '')
        chromosome = chromosome.replace('chr', '')
        chromosome = chromosome.replace(' ', '')
        chromosome = chromosome.upper()

        if chromosome not in Snp.CHROMOSOME_NAMES:
            raise Vcf2PrsError('Error: the method convert_chromosome has '
                               'been passed a non-valid chromosome, {0}.'.
                               format(chromo))

        return chromosome

    @staticmethod
    def convert_position(pos):
        """
        Function to convert the chromosoidal position to standard
        format, and check that it is valid.

        Args:
            pos (str/int): The chromosoidal position.

        Returns:
            int: An integer containing the position.

        Raises:
            Vcf2PrsError: If the argument is not a valid position.

        """
        # Check that it is of the correct type
        if (not isinstance(pos, (str, int))) or (type(pos) is bool):
            raise Vcf2PrsError('Error: the method convert_position has '
                               'been passed a non-string/integer argument, '
                               '"{0}", of type {1}'.
                               format(pos, type(pos)))

        # Convert strings to int
        try:
            position = int(pos)
        except (TypeError, ValueError):
            raise Vcf2PrsError('Error: the method convert_position could not '
                               'convert the "{0}", to integer.'.format(pos))

        # Check that the position is not negative.
        if position < 0:
            raise Vcf2PrsError('Error: the method convert_position has been '
                               'passed a negative integer, {0}.'.format(pos))
        else:
            return position

    @staticmethod
    def validate_allele_bases(allele):
        """
        Function to verify that a string is composed of the characters in
        Snp.BASES.

        Args:
            allele (str): The allele

        Returns:
            bool: True if the string is composed of bases, False if not.

        Raises:
            Vcf2PrsError: If an error occurs.

        """
        if not isinstance(allele, str):
            raise Vcf2PrsError('Error: the method validate_allele_bases has '
                               'been passed a non-string argument, "{0}", of '
                               'type {1}'.format(allele, type(allele)))
        for a in allele:
            if a.upper() not in Snp.BASES:
                return False
        return True

    def _calculate_mean_and_variance(self):
        """
        Method to calculate the mean and variance.
        """
        # In principle, for a male on chromosme X or Y this may need to be
        # altered, which may be the case for prostate cancer.

        low = (1.0 - self.eaf) ** 2
        med = 2.0 * self.eaf * (1.0 - self.eaf)
        hi = self.eaf ** 2
        eor = exp(self.lor)
        self.var = log((low + med * (eor ** 2) + hi * (eor ** 4)) /
                       (low + med * eor + hi * (eor ** 2)) ** 2)
        self.mean = (med + 2.0 * hi) * self.lor
        return
