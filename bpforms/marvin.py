""" Calculate major microspecies and draw molecules

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-02-27
:Copyright: 2020, Karr Lab
:License: MIT
"""

try:
    from wc_utils.util.chem import marvin
except ImportError:
    marvin = None


def get_major_micro_species(*args, **kwargs):
    if not marvin:
        raise ImportError("ChemAxon Marvin must be installed")
    return marvin.get_major_micro_species(*args, **kwargs)


def draw_molecule(*args, **kwargs):
    if not marvin:
        raise ImportError("ChemAxon Marvin must be installed")
    return marvin.draw_molecule(*args, **kwargs)
