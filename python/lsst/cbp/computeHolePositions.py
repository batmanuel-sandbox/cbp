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

__all__ = ["computeHolePositions"]

import math

from lsst.afw.geom import Point2D
from lsst.afw.cameraGeom import PIXELS, FIELD_ANGLE
from .coordUtils import fieldAngleToVector, vectorToFieldAngle


def computeHolePositions(detectorNames, detectorPositions, cameraGeom, cbpFlipX, cbpFocalLength):
    """Compute hole positions for a CBP mask

    Given the desired locations of one or more spots on on each detector,
    and assuming the telescope and CBP are pointing directly at each other,
    compute hole positions for a CBP mask.

    Parameters
    ----------
    detectorNames : `iterable` of `str`, or None,
        List of detector names
    detectorPositions : `iterable` of pair of `float`
        Detector x, y positions (pixels)
        Note that the center of a 1kx1k detector is (499.5, 499.5)
    cameraGeom : `lsst.afw.cameraGeom.Camera`
        Camera geometry
    cbpFlipX : `bool`
        Is the CBP focal plane flipped?

    Returns
    -------
    holePositions : `list` of `tuple` of pair of `float`
        CBP mask x, y hole positions mm
    """
    holePositions = []
    pixelPosList = [Point2D(*val) for val in detectorPositions]
    if detectorNames is None:
        detectorNames = list(cameraGeom.getNameIter())
    for detectorName in detectorNames:
        detector = cameraGeom[detectorName]
        pixelSys = detector.makeCameraSys(PIXELS)
        pixelsToFieldAngle = cameraGeom.getTransform(pixelSys, FIELD_ANGLE)
        telFieldAngleList = pixelsToFieldAngle.applyForward(pixelPosList)
        for telFieldAngle in telFieldAngleList:
            # compute hole positions for when the telescope and CBP are pointing at each other;
            # thus CBP pupil vector = telescope pupil vector with z negated.
            # This is simply a shortcut for setting telAzAlt and cbpAzAlt and then transforming
            # vectors from tel pupil to base and then to CBP pupil.
            telVector = fieldAngleToVector(telFieldAngle, False)
            cbpVector = (telVector[0], telVector[1], -telVector[2])
            cbpFieldAngle = vectorToFieldAngle(cbpVector, cbpFlipX)
            holePos = tuple(math.tan(ang) * cbpFocalLength for ang in cbpFieldAngle)
            holePositions.append(holePos)
    return holePositions
