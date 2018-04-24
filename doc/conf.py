"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.cbp


_g = globals()
_g.update(build_package_configs(
    project_name='cbp',
    version=lsst.cbp.version.__version__))
