""" Configuration

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2019-03-11
:Copyright: 2017-2019, Karr Lab
:License: MIT
"""

import configobj
import os
import pkg_resources
import wc_utils.config


def get_config(extra=None):
    """ Get configuration

    Args:
        extra (:obj:`dict`, optional): additional configuration to override

    Returns:
        :obj:`configobj.ConfigObj`: nested dictionary with the configuration settings loaded from the configuration source(s).
    """
    paths = wc_utils.config.ConfigPaths(
        default=pkg_resources.resource_filename('bpforms', 'config/core.default.cfg'),
        schema=pkg_resources.resource_filename('bpforms', 'config/core.schema.cfg'),
        user=(
            'bpforms.cfg',
            os.path.expanduser('~/.wc/bpforms.cfg'),
        ),
    )

    return wc_utils.config.ConfigManager(paths).get_config(extra=extra)
