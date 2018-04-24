# This file is part of cbp.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["CoordinateConverterConfig"]

import numpy as np

from lsst.afw.geom import degrees


class CoordinateConverterConfig:
    """Configuration for the CoordinateConverter

    :ref:`Configuration <lsst.cbp.configuration>` for the
    lsst.cbp.CoordinateConverter.

    Parameters
    ----------
    telPupilOffset : `float`
        Offset of the telescope pupil from the center of the telescope
        (mm, + if closer to the sky)
    telPupilDiameter : `float`
        Diameter of telescope pupil (mm)
    telPupilObscurationDiameter : `float`
        Diameter of obscuration in center of telescope pupil (mm)
    telFocalPlaneDiameter : `float`
        Diameter of telescope flocal plane (mm)
    telFlipX : `bool`
        True if the x axis of the telescope focal plane
        is flipped with respect to the pupil frame
    telAzimuthOffset : `float`
        Azimuth `offset`_ (degrees); defaults to 0
    telAzimuthScale : `float`
        Azimuth `scale`_; defaults to 1; must be ±1
        (in order to handle wrap correctly).
    telAltitudeOffset : `float` (optional)
        Telescope altitude `offset`_ (degrees); defaults to 0
    telAltitudeScale : `float` (optional)
        Telescope altitude `scale`_; defaults to 1
    telAltitudeLimits : pair of `float`
        Telescope minimum, maximum allowed observed altitude (degrees)
    telRotOffset : `float` (optional)
        Telescope camera rotator offset (degrees); defaults to 0
    telRotScale : `float`
        Telescope camera rotator scale; defaults to 1; must be ±1
        (in order to handle wrap correctly)
    defaultDetector : `str`
        Name of default detector

    cbpPosition : triplet of `float`
        CBP x, y, z position of center of CBP relative to the center
        of the telescope, in the base frame (mm)
    cbpFocalLength : `float`
        Effective focal length of the CBP (mm);
        635 mm is an estimate for the LSST's CBP
    cbpFlipX : `bool`
        True if the x axis of the CBP focal plane is flipped
        with respect to the pupil frame?
    cbpAzimuthOffset : `float` (optional)
        CBP azimuth `offset`_ (degrees); defaults to 0
    cbpAzimuthScale : `float` (optional)
        CBP azimuth `scale`_; defaults to 1; must be ±1
        (in order to handle wrap correctly)
    cbpAltitudeOffset : `float` (optional)
        CBP altitude `offset`_ (degrees); defaults to 0
    cbpAltitudeScale : `float` (optional)
        CBP altitude `scale`_; defaults to 1
    cbpAltitudeLimits : pair of `float`
        CBP minimum, maximum allowed observed altitude (degrees)

    Raises
    ------
    AssertionError
        Raised if ``telAzimuthScale``, ``cbpAzimuthScale``
        and/or ``telRotScale`` is not ±1

    Notes
    -----

    .. _offset:

    .. _scale:

    **Offset and Scale:**

    Azimuth, altitude and rotator offset and scale define the mapping
    between :ref:`internal angle <lsst.cbp.internal_angles>`
    and :ref:`observed angle <lsst.cbp.observed_angles>` as follows::

        observed angle = internal angle * scale + offset

    .. _fields:

    **Fields:**

    telPupilOffset : `float`
        Offset of the telescope pupil from the center of the telescope
        (mm, + if closer to the sky)
    telPupilDiameter : `float`
        Diameter of telescope pupil (mm)
    telPupilObscurationDiameter : `float`
        Diameter of obscuration in center of telescope pupil (mm)
    telFocalPlaneDiameter : `float`
        Diameter of telescope flocal plane (mm)
    telFlipX : `bool`
        True if the x axis of the telescope focal plane
        is flipped with respect to the pupil frame
    telAzAltOffset : pair of `lsst.afw.geom.Angle`
        Telescope azimuth and altitude `offset`_ (degrees)
    telAzAltScale : pair of `float`
        Telescope azimuth and altitude `scale`_;
        azimuth scale is ±1
    telAltitudeLimits : pair of `lsst.afw.geom.Angle`
        Telescope minimum, maximum allowed observed altitude
    telRotOffset : `float` (optional)
        Telescope camera rotator offset (degrees)
    telRotScale : `float`
        Telescope camera rotator scale; ±1
    defaultDetector : `str`
        Name of default detector

    cbpFocalLength : `float`
        Effective focal length of the CBP (mm);
        635 mm is an estimate for the LSST's CBP
    cbpFlipX : `bool`
        True if the x axis of the CBP focal plane is flipped
        with respect to the pupil frame?
    cbpAzAltOffset : pair of `lsst.afw.geom.Angle`
        CBP azimuth and altitude `offset`_ (degrees)
    cbpAzAltScale : pair of `float`
        CBP azimuth and altitude `scale`_; azimuth scale is ±1
    cbpAltitudeLimits : pair of `lsst.afw.geom.Angle`
        CBP minimum, maximum allowed observed altitude
    """
    def __init__(self, *, telPupilOffset,
                 telPupilDiameter, telPupilObscurationDiameter, telFocalPlaneDiameter, telFlipX,
                 telAzimuthOffset=0*degrees, telAzimuthScale=1,
                 telAltitudeOffset=0*degrees, telAltitudeScale=1, telAltitudeLimits,
                 telRotOffset=0*degrees, telRotScale=1,
                 defaultDetector,
                 cbpPosition, cbpFocalLength, cbpFlipX,
                 cbpAzimuthOffset=0*degrees, cbpAzimuthScale=1,
                 cbpAltitudeOffset=0*degrees, cbpAltitudeScale=1, cbpAltitudeLimits):
        self.telPupilOffset = telPupilOffset
        self.telPupilDiameter = telPupilDiameter
        self.telPupilObscurationDiameter = telPupilObscurationDiameter
        self.telFocalPlaneDiameter = telFocalPlaneDiameter
        self.telFlipX = telFlipX
        self.telAzAltOffset = (telAzimuthOffset*degrees, telAltitudeOffset*degrees)
        assert abs(telAzimuthScale) == 1
        self.telAzAltScale = (telAzimuthScale, telAltitudeScale)
        assert len(telAltitudeLimits) == 2
        self.telAltitudeLimits = tuple(val*degrees for val in telAltitudeLimits)
        self.telRotOffset = telRotOffset*degrees
        assert abs(telRotScale) == 1
        self.telRotScale = telRotScale

        self.defaultDetector = defaultDetector

        self.cbpPosition = cbpPosition
        self.cbpFocalLength = cbpFocalLength
        self.cbpFlipX = cbpFlipX
        self.cbpAzAltOffset = (cbpAzimuthOffset*degrees, cbpAltitudeOffset*degrees)
        assert abs(cbpAzimuthScale) == 1
        self.cbpAzAltScale = (cbpAzimuthScale, cbpAltitudeScale)
        assert len(cbpAltitudeLimits) == 2
        self.cbpAltitudeLimits = tuple(val*degrees for val in cbpAltitudeLimits)

    @property
    def cbpPosition(self):
        """Get the position of the CBP
        """
        return self._cbpPosition

    @property
    def cbpDistance(self):
        """Get the distance from the telescope to the CBP (mm)"""
        return self._cbpDistance

    @cbpPosition.setter
    def cbpPosition(self, cbpPosition):
        """Set the position of the CBP

        Parameters
        ----------
        cbpPosition : triplet of `float`
            CBP x, y, z position of center of CBP relative to the center
            of the telescope, in the base frame (mm)
        """
        assert len(cbpPosition) == 3
        self._cbpPosition = np.array(cbpPosition)
        self._cbpDistance = np.linalg.norm(self.cbpPosition)
