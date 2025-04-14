PRS from VCF file
=================


Description
-----------

This package provides a command line tool, ``vcf2prs``, that performs 4 
different calculations depending on the options supplied:

#. Calculate an individual's PRS from a VCF file.
#. Convert a raw PRS to a z-score.
#. Convert a z-score to a raw PRS.
#. Display the characteristics of a PRS.

A z-score is a standard normal PRS.
It can also work on mixed-ancestry individuals, if ancestry proportions are provided.

If you are only interested in calculating PRS for individuals of mixed-ancestry and you already have the raw PRS 
(and the ancestry proportions), then you can also use the R or python scripts in the standaloneScripts folder.
[Ficorella L, Yang X, Mavaddat N, Carver T, Hassan H, Dennis J, et al. Adapting the BOADICEA breast and ovarian cancer risk models for the ethnically diverse UK population. medRxiv. 2025:2025.02.14.25322307.]


Inputs
------

All calculations require a PRS file - a file containing information on the 
variants used to construct the PRS (detail given below).

Calculating an individual's PRS from their genotype requires:

-  A genotype file, i.e. a VCF file containing the genotypes of the individual of interest.
-  If there is more than one sample in the genotype file, users can either
   supply the name of the sample of interest or calculate PRS on all samples.

Converting a raw PRS to a z-score requires

-  A raw PRS (float number)

Converting a z-score to a raw PRS requires

-  A z-score (float number)

Calculating z-score for mixed-ancestry individuals requires

-  A proportions file, i.e. a .csv file containing the raw PRS (if available) and
   the proportions for each ancestry
-  If there is more than one sample in the genotype file, users can either
   supply the name of the sample of interest or calculate z-score on all samples.


The PRS file
~~~~~~~~~~~~

The PRS file contains the information on the variants used in the PRS.
It should be a plain text file (UTF-8).

The first lines should of the file contain the values of ``threshold``, ``mean``, ``sd`` and ``alpha``.
They must be provided in that order; however, for non-mixed calculations, only ``alpha`` is required.
``alpha`` is the square root of the proportion of the overall model polygenic variance explained by the PRS.
``threshold`` is the proportion above which an individual is considered as a single-ancestry one.
Values should be in one of the forms (blank spaces are ignored):

#. ``0.45``
#. ``alpha = 0.45``
#. ``0.497,0.497,0.196,0.326,0.329   #alpha``


The rest of the file should be a comma separated value (csv) file with
fields separated by commas, ",".
There should be a line containing the column titles, which are case 
sensitive, followed by one line for each variant.
It should contain the following fields:

-  'Chromosome' - the chromosome on which the variant resides
-  'Position' - the position of the variant
-  'Reference\_Allele' - the reference allele for the variant
-  'Effect\_Allele' - the effect allele for the variant
-  'Log\_Odds\_Ratio' - the log odds ratio for the variant
-  'Effect\_Allele\_Frequency'- the (population) allele frequency for the variant
-  'eaf_eur','eaf_afr','eaf_eas' and 'eaf_sas' - the (ancestry-specific) allele frequencies
    eur = European; afr = African; eas = East Asian; sas = South Asian

The file can contain extra columns, which will be ignored.
The order of the columns is not important.

Genotypes from the VCF file are matched to the PRS variants via
chromosome and position, so build should be consistent between the two
files.

Some PRS files are given in the directory ``PRSmodels_CanRisk``, e.g.:

-  PRSmodels\_CanRisk/BCAC\_313\_PRS.prs
-  PRSmodels\_CanRisk/BRIDGES\_306\_PRS.prs
-  PRSmodels\_CanRisk/PERSPECTIVE\_295\_PRS.prs
-  PRSmodels\_CanRisk/PRISM\_289\_PRS.prs
-  PRSmodels\_CanRisk/OCAC\_36\_PRS.prs

BCAC\_313\_PRS.prs, BRIDGES\_306\_PRS.prs, PERSPECTIVE\_295\_PRS.prs, and  
PRISM\_289\_PRS.prs are taken from analyses using the BCAC v11
dataset (Europeans), using the subset of individuals in the validation
and test sets combined, as per the analysis in Mavaddat et al. (2019).
Individuals with unknown age and age>=80 were excluded.
OCAC\_36\_PRS.prs is taken from Dareng et al. (in prep, 2021) with
frequencies from BCAC v11 dataset (Europeans), using the subset of
individuals in the validation and test sets combined, as per the
analysis in Mavaddat et al (2019).
Individuals with unknown age and age>=80 were excluded.

**NB: All alleles are assumed to be on the forward strand and use
positions from NCBI Build 37, hg19.**


The VCF file
~~~~~~~~~~~~

The VCF file should be in VCF format v4.0 or v4.1.
It can contain additional samples (which are not used), and may also
contain variants not used in the PRS (which are ignored).

The script can use either genotypes or dosages, or a mixture of both.
The specified sample should have a valid genotype or dosage for each
variant in the PRS file, if not an exception is raised, and no result is
returned.
If both a genotypes and a dosages are present and the genotype is
called, the dosage and genotype must be consistent (up to a tolerance of
1e-4).
Genotypes are included in the VCF file with the variable ``GT``, as
per VCFv4.1.

Genotypes from the VCF file are matched to the PRS variants via
chromosome and position, so build should be consistent between the two
files.

An example of a VCF file that uses genotyping is given in
``Examples/VCFfiles/sample_BCAC_313_genotypes.vcf``, where the corresponding
PRS file is ``PRS_files/BCAC_313_PRS.prs``.

Dosages
^^^^^^^

A sample's dosage for a variant is given by

.. code::

   DS = P(0/1) + 2 x P(1/1)

where ``P(0/1)`` is the probability of the sample having one copy of
the alternative allele, and ``P(1/1)`` is the probability of the sample
having two copies.
These probabilities can be imputed or derived from the sample's
genotype (i.e. 0, or 1).
Dosages should be presented as float (real) or integer values.

When using dosages the VCF file should contain the meta-information
line

.. code::

   ##FORMAT=<ID=DS,Number=1,Type=Float,Description="Alternate Allele Dosage">

and the variable ``DS`` should be specified in the record format
field.

When using dosages, there should only be one alternative allele for the
record, and it must be the same as the effect allele specified in the
PRS file.

An example of a VCF file that uses dosages is given in
``Examples/VCFfiles/sample_BCAC_313_dosages.vcf``, where the corresponding
PRS file is ``PRS_files/BCAC_313_PRS.prs``.

Examples
^^^^^^^^

Examples VCF files are given in the directory: ``Examples/VCFfiles``

-  'sample\_BCAC\_313\_genotypes.vcf' - A sample file for the BCAC 313
   variant PRS using genotypes only, with samples 'Low', 'Med', 'Mod'
   and 'High' corresponding to roughly 5%, 50%, 75% and 95% PRS
   percentiles respectively.
-  'sample\_BCAC\_313\_dosages.vcf' - A sample file for the BCAC 313
   variant PRS using a combination of genotypes and dosages, with
   samples 'Low', 'Med', 'Mod' and 'High' corresponding to roughly 5%,
   50%, 75% and 95% PRS percentiles respectively.
-  'sample\_BCAC\_313.vcf' - A sample file for the BCAC 313 variant PRS
   using a combination of genotypes and dosages, with samples '0.01',
   '0.025', '0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7',
   '0.8', '0.9', '0.95', '0.975', '0.99' corresponding to roughly 1%,
   2.5%, 5%, 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 95%, 97.5% and
   99% PRS percentiles respectively.
-  'sample\_OCAC\_36.vcf' - A sample file for the OCAC 36 variant PRS
   using a combination of genotypes and dosages, with samples '0.01',
   '0.025', '0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7',
   '0.8', '0.9', '0.95', '0.975', '0.99' corresponding to roughly 1%,
   2.5%, 5%, 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 95%, 97.5% and
   99% PRS percentiles respectively.

Tri-Allelic Variants
^^^^^^^^^^^^^^^^^^^^

Some variants may have more than two possible alleles.
This case can be accommodated by noting that the the VCF format allows
a single reference allele and multiple alternative alleles.
The reference allele should be taken as the same as that in the PRS
file, with all other alleles taken as the alternative alleles.
The effect allele in the PRS file should be among the alternative
alleles.
The script counts the number of copies of the effect allele (0, 1, or
1) in the sample's genotype.
The other alternative alleles are not counted, so in effect they will
be treated as if they are copies of the reference allele.

If dosages are required at such a variant, and assuming that the
dosage is the dosage of the effect allele, then the extraneous
alternative alleles should be removed from the VCF file.
Any genotypes that feature the extraneous alternative alleles should
be recoded as having the reference allele in its place.


Converting between a Raw PRS and a Z-score
------------------------------------------

A raw PRS can be converted to a z-score and vice versa using the mean
(mu) and standard deviation (sigma) of the raw PRS.
The raw PRS can be converted the a z-score using the formula:

.. code::

   z-score = (rawPRS - mu) / sigma

A z-score can be converted to a raw PRS using the formula:

.. code::

   rawPRS = z-score * sigma + mu

The mean and standard deviation are calculated in the PRS
characteristics.


Installation
-----------

These scripts require Python 3.

The repository is structured as a Python package, and the command line tool
and module can be instlalled using the ``pip`` Python package installer.
First download the repository from GitHub, either via SSH

.. code:: console

   $ git clone git@github.com:CCGE-BOADICEA/SHARE-PRScalculation.git

or via HTTPS

.. code:: console

   $ git clone https://github.com/CCGE-BOADICEA/SHARE-PRScalculation.git

Then install the package via

.. code:: console

   $ cd SHARE-PRScalculation
   $ pip install .

This will make the command available system-wide.


The package can be uninstalled using the command

.. code:: console

   $ pip uninstall -y vcf2prs

Requirements
~~~~~~~~~~~~

This package uses the PyVCF and numpy Python packages, which are listed in
``requirements.txt``.
They can be installed using the command

.. code:: console
   $ pip install -r requirements.txt

This package also requires the following packages: argparse, csv, io, math, sys.
Since they are usually already installed in most system, they are not currently included in the requirements file.


Command line Usage
------------------

The package provides the command ``vcf2prs``, designed to be used on the 
Unix/Linux command line. Further details on usage are given by the help 
function:

.. code:: console

   $ vcf2prs -h


Use as a Module
~~~~~~~~~~~~

In addition to providing command line tools the package also provides a Python
module ``vcf2prs``.
This provides the class ``Vcf2Prs``.
Once the package has been installed it can be imported via:

.. code:: python3

   import vcf2prs
