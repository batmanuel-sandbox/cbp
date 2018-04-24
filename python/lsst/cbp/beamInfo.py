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

__all__ = ["BeamInfo"]

import numpy as np

from lsst.afw.geom import Point2D
from lsst.afw.cameraGeom import PIXELS, FOCAL_PLANE, FIELD_ANGLE


class BeamInfo:
    """Information about a beam at the telescope

    Parameters
    ----------
    cameraGeom : `lsst.afw.cameraGeom.Camera`
        Camera geometry
    name : `str`
        Beam name
    holePos : pair of `float`
        Hole position on CBP mask (x, y mm)
    isOnPupil : `bool`
        See `fields`_ below
    isOnFocalPlane : `bool`
        See `fields`_ below
    focalPlanePos : pair of `float`
        See `fields`_ below
    pupilPos : pair of `float`
        See `fields`_ below
    focalFieldAngle : pair of `float`
        See `fields`_ below
    pupilFieldAngle : pair of `float`
        See `fields`_ below

    Notes
    -----
    .. _fields:

    **Fields:**

    name : `str`
        Name of beam
    holePos : pair of `float`
        Position of hole on CBP mask (x, y mm)
    isOnPupil : `bool`
        Beam is likely on the pupil and not obscured by the secondary
    isOnFocalPlane : `bool`
        Beam is likely on the focal plane;
        this is independent of ``isOnPupil``; both must be true
        for the beam to be detected
    isOnDetector : `bool`
        Beam is probably on a detector; this is unreliable
    focalPlanePos : `lsst.afw.geom.Point2D`
        Telescope :ref:`focal plane <lsst.cbp.focal_plane>` position
        of beam (x, y mm)
    focalFieldAngle : `lsst.afw.geom.Point2D`
        :ref:`Focal plane field angle <lsst.cbp.focal_plane_field_angle>`
        of beam (x, y rad)
    pupilFieldAngle : `lsst.afw.geom.Point2D`
        :ref:`Pupil field angle <lsst.cbp.pupil_field_angle>`
        of beam (x, y rad)
    pupilPos : `lsst.afw.geom.Point2D`
        Telescope :ref:`pupil plane position <lsst.cbp.pupil_position>`
        of beam (x, y mm)

    """
    def __init__(self, cameraGeom, name, holePos, isOnPupil, isOnFocalPlane, focalPlanePos, pupilPos,
                 focalFieldAngle, pupilFieldAngle):
        self._cameraGeom = cameraGeom
        self.name = name
        self.holePos = holePos
        self.isOnPupil = isOnPupil
        self.isOnFocalPlane = isOnFocalPlane
        self.focalPlanePos = Point2D(*focalPlanePos)
        self.focalFieldAngle = Point2D(*focalFieldAngle)
        self.pupilFieldAngle = Point2D(*pupilFieldAngle)
        self.pupilPos = Point2D(*pupilPos)
        fieldAngleToFocalPlane = self._cameraGeom.getTransform(FIELD_ANGLE, FOCAL_PLANE)
        self.focalPlanePos = fieldAngleToFocalPlane.applyForward(self.focalFieldAngle)
        self._isOnDetector = None
        self._detector = None

    @property
    def isOnDetector(self):
        """True if the beam is on a detector
        """
        if self._isOnDetector is None:
            self._findDetector()
        return self._isOnDetector

    @property
    def detectorPos(self):
        """The position of the spot on the detector,
        or (nan, nan) if isOnDetector is False
        """
        if not self.isOnDetector:
            return Point2D(np.nan, np.nan)
        pixelSys = self._detector.makeCameraSys(PIXELS)
        focalPlaneToPixels = self._cameraGeom.getTransform(FOCAL_PLANE, pixelSys)
        return focalPlaneToPixels.applyForward(self.focalPlanePos)

    @property
    def detectorName(self):
        """The name of the detector that the beam falls on,
        or None if isOnDetector is False
        """
        if not self.isOnDetector:
            return None
        return self._detector.getName()

    def _findDetector(self):
        detectors = self._cameraGeom.findDetectors(self.focalFieldAngle, FIELD_ANGLE)
        if not detectors:
            self._isOnDetector = False
        else:
            self._detector = detectors[0]
            self._isOnDetector = True
