.. py:currentmodule:: lsst.cbp

.. _lsst.cbp:

########
lsst.cbp
########

The ``cbp`` package provides code for the collimated beam projector (CBP).
The primary objects of interest are:

`lsst.cbp.CoordinateConverter`: compute the telescope and CBP pointing that will give you a desired beam arrangement, such as placing beam B at point P on the pupil and point D on a specified detector.

`lsst.cbp.computeHolePositions`: compute hole positions for a CBP mask.

To construct an `lsst.cbp.CoordinateConverter` you will need to learn about configuration:

.. toctree::

    configuration

In order to :ref:`configure <lsst.cbp.configuration>` and use an `lsst.cbp.CoordinateConverter` it may help to have some understanding of the coordinate systems involved:

.. toctree::

    coordinateFrames

Python API reference
====================

.. automodapi:: lsst.cbp
.. automodapi:: lsst.cbp.coordUtils
.. automodapi:: lsst.cbp.testUtils
