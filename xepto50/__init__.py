"""
******************
***   Xepto50  ***
******************

Compute efficacy, we adapted xepto from zepto, zeptomolar (zM) is equal to 10^-21 moles/liter .
"""
__author__ = 'Julhash Kazi'
__email__ = 'xepto@kazilab.se'
__description__ = 'Drug efficacy calculation'
__version__ = '0.0.2'
__url__ = 'https://www.kazilab.se'

from .cal import XCal
from .run import x50
from .__main__ import main
