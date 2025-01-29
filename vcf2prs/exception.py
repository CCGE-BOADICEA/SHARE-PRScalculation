"""
:authors: Andrew Lee, Lorenzo Ficorella
:organization: CCGE, DPHPC, University of Cambridge
:copyright:  2022 Cambridge University. All rights reserved.
:   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
:   without even the implied warranty of MERCHANTABILITY or FITNESS for a particular purpose.
:contact: lf364@medschl.cam.ac.uk
"""


class Vcf2PrsError(Exception):
    """
    Vcf2PrsError is the exception that is raised in this module PRS from VCF.

    Args:
        message (str): Explanation of why the exception was raised.

    """

    def __init__(self, message):
        super().__init__(message)
