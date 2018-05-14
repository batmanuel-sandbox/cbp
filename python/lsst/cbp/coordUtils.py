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

__all__ = ["fieldAngleToVector", "vectorToFieldAngle", "pupilPositionToVector",
           "computeShiftedPlanePos", "convertVectorFromBaseToPupil", "convertVectorFromPupilToBase",
           "computeAzAltFromBasePupil", "getFlippedPos", "rotate2d"]

import math

import numpy as np

from lsst.sphgeom import Vector3d
from lsst.afw.geom import SpherePoint, degrees, radians

# error data from computeAzAltFromBasePupil;
# error is determined by calling convertVectorFromPupilToBase
# using the az/alt computed by computeAzAltFromBasePupil
# and measuring the separation between that base vector and the input
_RecordError = False  # record errors?
_ErrorLimitArcsec = None  # only record errors larger than this
_ErrorList = None  # a list of tuples: (errorArcsec, vectorBase, vectorPupil)


def startRecordingErrors(errorLimit=1e-11*radians):
    """Start recording numeric errors in computeAzAltFromBasePupil
    and reset the error list.

    Parameters
    ----------
    errorLimit : `lsst.afw.geom.Angle`
        Only errors larger than this limit are recorded
    """
    global _RecordError, _ErrorLimitArcsec, _ErrorList
    _RecordError = True
    _ErrorLimitArcsec = errorLimit.asArcseconds()
    _ErrorList = []


def stopRecordingErrors():
    """Stop recording numeric errors in computeAzAltFromBasePupil
    """
    global _RecordError
    _RecordError = False


def getRecordedErrors():
    """Get recorded numeric errors in computeAzAltFromBasePupil,
    sorted by increasing error

    Returns
    -------
    errorList : `list` of tuples
        Recorded error as a list of tuples containing:
        - |error| in radians
        - vectorBase
        - vectorPupil
    """
    global _ErrorList
    if _ErrorList is None:
        raise RuntimeError("Errors were never recorded")
    return sorted(_ErrorList, key=lambda elt: elt[0])


RefSpherePoint = SpherePoint(0, 0, radians)


def getFlippedPos(xyPos, flipX):
    """Get a 2-dimensional position with the x axis properly flipped

    Parameters
    ----------
    xyPos : pair of `float`
        Position to rotate
    flipX : `bool`
        True if the x axis should be flipped (negated)

    Returns
    -------
    xyresult : pair of `float`
        ``xyPos`` with the x axis flipped or not, according to ``flipX``
    """
    return (-xyPos[0], xyPos[1]) if flipX else xyPos


def fieldAngleToVector(xyrad, flipX):
    """Convert a pupil field angle to a pupil unit vector

    Parameters
    ----------
    xyrad : `tuple` of 2 `float`
        x,y :ref:`pupil field angle <lsst.cbp.pupil_field_angle>` (radians)

    Returns
    -------
    vector : a numpy array of 3 `float`
        A unit vector in the :ref:`pupil frame <lsst.cbp.pupil_frame>`
    """
    assert len(xyrad) == 2
    xyradRightHanded = getFlippedPos(xyrad, flipX=flipX)
    amount = math.hypot(*xyradRightHanded)*radians
    bearing = math.atan2(xyradRightHanded[1], xyradRightHanded[0])*radians
    return np.array(RefSpherePoint.offset(bearing=bearing, amount=amount).getVector())


def vectorToFieldAngle(vec, flipX):
    """Convert a vector to a pupil field angle

    Parameters
    ----------
    vec : `iterable` of 3 `float`
        3-dimensional vector in the
        :ref:`pupil frame <lsst.cbp.pupil_frame>`;
        the magnitude is ignored, but must be large enough
        to compute an accurate unit vector
    flipX : bool
        Set True if the x axis of the focal plane is flipped
        with respect to the pupil

    Returns
    -------
    fieldAngle : `tuple` of 2 `float`
        x,y :ref:`pupil field angle <lsst.cbp.pupil_field_angle>` (radians)
    """
    sp = SpherePoint(Vector3d(*vec))
    amountRad = RefSpherePoint.separation(sp).asRadians()
    bearingRad = RefSpherePoint.bearingTo(sp).asRadians()
    xyrad = (amountRad*math.cos(bearingRad),
             amountRad*math.sin(bearingRad))
    return getFlippedPos(xyrad, flipX=flipX)


def pupilPositionToVector(xyPos, flipX):
    """Convert a pupil plane position to a 3D vector

    Parameters
    ----------
    xyPos : sequence of 2 `float`
        :ref:`pupil plane position <lsst.cbp.pupil_position>` (mm)
    flipX : `bool`
        True if the x axis of the position flipped (negated)

    Returns
    -------
    vector : sequence of 3 `float`
        3-dimensional vector in the
        :ref:`pupil frame <lsst.cbp.pupil_frame>`:
        x = 0, y = plane position x, z = plane position y
    """
    xyPosRightHanded = getFlippedPos(xyPos, flipX=flipX)
    return np.array([0, xyPosRightHanded[0], xyPosRightHanded[1]])


def computeShiftedPlanePos(planePos, fieldAngle, shift):
    """Compute the plane position of a vector on a plane
    shifted along the optical axis

    Parameters
    ----------
    planePos : pair of `float`
        Plane position at which the vector intersects the plane (x, y)
    fieldAngle : pair of `float`
        Field angle of vector (x, y radians)
    shift : `float`
        Amount by which the new plane is shifted along the optical axis.
        If ``shift`` and both components of ``fieldAngle``
        are positive then both axes of the shifted plane position
        will be larger (more positive) than ``planePos``.

    Returns
    -------
    shiftedPlanePos : tuple of `float`
        Plane position at which vector intersects the shifted plane (x, y)

    Notes
    -----
    `flipX` is not an input because for this computation it does not matter
    if the x axis is flipped: fieldAngle and planePos are either both
    flipped or not, and that cancels out the effect.
    """
    # unit vector y = plane x, unit vector z = plane y
    # (ignoring flipped x axis, which cancels out)
    unitVector = fieldAngleToVector(fieldAngle, False)
    dxy = [shift*unitVector[i]/unitVector[0] for i in (1, 2)]
    return tuple(planePos[i] + dxy[i] for i in range(2))


def convertVectorFromBaseToPupil(vectorBase, pupilAzAlt):
    """Given a vector in base coordinates and the pupil pointing,
    compute the vector in pupil coordinates

    Parameters
    ----------
    vectorBase: sequence of 3 `float`
        3-dimensional vector in the :ref:`base frame <lsst.cbp.base_frame>`
    pupilAzAlt: `lsst.afw.geom.SpherePoint`
        Pointing of the pupil frame as
        :ref:`internal azimuth, altitude <lsst.cbp.internal_angles>`

    Returns
    -------
    vectorPupil: `np.array` of 3 `float`
        3-dimensional vector in the
        :ref:`pupil frame <lsst.cbp.pupil_frame>`

    Note
    ----
    This can also be implemented as the following Euler angle
    rotation matrix, where angle 1 = azimuth, angle 2 = altitude,
    sx = sine(anglex) and cx = cosine(anglex):

    c1*c2    -s1    -c1*s2
    c2*s1     c1     s1s2
    s1        0      c2
    """
    vectorMag = np.linalg.norm(vectorBase)
    vectorSpBase = SpherePoint(Vector3d(*vectorBase))

    telBase = SpherePoint(0, 0, radians)

    # rotate vector around base z axis by -daz
    daz = pupilAzAlt[0] - telBase[0]
    zaxis = SpherePoint(Vector3d(0, 0, 1))
    vectorSpRot1 = vectorSpBase.rotated(axis=zaxis, amount=-daz)

    # rotate vector around pupil -y axis by -dalt
    dalt = pupilAzAlt[1] - telBase[1]
    negYAxis = SpherePoint(Vector3d(0, -1, 0))
    vectorSpPupil = vectorSpRot1.rotated(axis=negYAxis, amount=-dalt)
    return np.array(vectorSpPupil.getVector()) * vectorMag


def convertVectorFromPupilToBase(vectorPupil, pupilAzAlt):
    """Given a vector in pupil coordinates and the pupil pointing,
    compute the vector in base coords

    Parameters
    ----------
    vectorPupil: sequence of 3 `float`
        3-dimesional vector in the
        :ref:`pupil frame <lsst.cbp.pupil_frame>`
    pupilAzAlt: `lsst.afw.geom.SpherePoint`
        Pointing of the pupil frame as
        :ref:`internal azimuth, altitude <lsst.cbp.internal_angles>`

    Returns
    -------
    vectorBase: `np.array` of 3 `float`
        3-dimensional vector in the :ref:`base frame <lsst.cbp.base_frame>`
    """
    vectorMag = np.linalg.norm(vectorPupil)
    vectorSpPupil = SpherePoint(Vector3d(*vectorPupil))

    telBase = SpherePoint(0, 0, radians)

    # rotate vector around pupil -y axis by dalt
    dalt = pupilAzAlt[1] - telBase[1]
    negYAxis = SpherePoint(Vector3d(0, -1, 0))
    vectorSpRot1 = vectorSpPupil.rotated(axis=negYAxis, amount=dalt)

    # rotate that around base z axis by daz
    daz = pupilAzAlt[0] - telBase[0]
    zaxis = SpherePoint(Vector3d(0, 0, 1))
    vectorSpBase = vectorSpRot1.rotated(axis=zaxis, amount=daz)
    return np.array(vectorSpBase.getVector()) * vectorMag


def computeAzAltFromBasePupil(vectorBase, vectorPupil):
    """Compute az/alt from a vector to in the base frame
    and the same vector in the pupil frame

    Parameters
    ----------
    vectorBase : `iterable` of three `float`
        3-dimensional vector in the
        :ref:`base frame <lsst.cbp.base_frame>`
    vectorPupil : `iterable` of `float`
        The same vector in the :ref:`pupil frame <lsst.cbp.pupil_frame>`.
        This vector should be within 45 degrees or so of the optical axis
        for accurate results.

    Returns
    -------
    pupilAzAlt : `lsst.afw.geom.SpherePoint`
        Pointing of the pupil frame as
        :ref:`internal azimuth, altitude <lsst.cbp.internal_angles>`

    Notes
    -----
    The magnitude of each vector is ignored, except that a reasonable
    magnitude is required in order to compute an accurate unit vector.
    """
    spBase = SpherePoint(Vector3d(*vectorBase))
    spPupil = SpherePoint(Vector3d(*vectorPupil))

    # compute telescope altitude
    # Use an equation derived from the following, which is from Euler angles:
    # base z = -sin(negAlt) pupil x + cos(negAlt) pupil z
    #        = sin(alt) pupil x + cos(alt) pupil z

    xb, yb, zb = spBase.getVector()
    xp, yp, zp = spPupil.getVector()

    factor = 1 / math.fsum((xp**2, zp**2))
    addend1 = xp * zb
    addend2 = zp * math.sqrt(math.fsum((xp**2, zp**2, -zb**2)))
    sinAltToTry = [val for val in (
        factor * (addend1 + addend2),
        factor * (addend1 - addend2),
    ) if abs(val) <= 1.0]
    if not sinAltToTry:
        raise RuntimeError("No solution found; required altitude out of range")

    zerror = 1
    goodSinAlt = None
    for i, sinAlt in enumerate(sinAltToTry):
        cosAlt = math.cos(math.asin(sinAlt))
        predzb = sinAlt*xp + cosAlt*zp
        # print("case", i, "sinAlt=", sinAlt, "predzb=", predzb, "zb=", zb, "delta=", predzb-zb)
        if abs(predzb - zb) < zerror:
            zerror = abs(predzb - zb)
            goodSinAlt = sinAlt
    if goodSinAlt is None:
        raise RuntimeError("No solution found")
    alt = math.asin(goodSinAlt)*radians

    # Consider the spherical triangle connecting the telescope pointing (pupil frame x axis),
    # the vector, and zenith (the base frame z axis). The length of all sides is known.
    # Solve for the angle between the sides at the zenith;
    # that angle is the difference in azimuth between the
    # telescope pointing and the azimuth of the base vector.
    # A is the side connecting the vector to the pupil frame x axis
    # (since 0, 0 is a unit vector pointing along pupil frame x);
    # B is the side connecting telescope pointing to the zenith
    # C is the side connecting the vector to the zenith
    sideA = SpherePoint(0, 0, radians).separation(spPupil).asRadians()
    sideB = (90*degrees - alt).asRadians()
    sideC = (90*degrees - spBase[1]).asRadians()

    # sideA can be small or zero so use a half angle formula
    # sides B and C will always be well away from 0 and 180 degrees
    semiPerimeter = 0.5*math.fsum((sideA, sideB, sideC))
    sinHalfAngleA = math.sqrt(math.sin(semiPerimeter - sideB) * math.sin(semiPerimeter - sideC) /
                              (math.sin(sideB) * math.sin(sideC)))
    daz = 2*math.asin(sinHalfAngleA)*radians
    if spPupil[0].wrapCtr() > 0:
        daz = -daz
    az = spBase[0] + daz
    global _RecordError
    if _RecordError:  # to study sources of numerical imprecision
        global _ErrorLimitArcsec, _ErrorList
        sp = SpherePoint(az, alt)
        vectorBaseRT = convertVectorFromPupilToBase(vectorPupil, sp)
        errorArcsec = SpherePoint(Vector3d(vectorBaseRT)).separation(
            SpherePoint(Vector3d(*vectorBase))).asArcseconds()
        if errorArcsec > _ErrorLimitArcsec:
            _ErrorList.append((errorArcsec, vectorBase, vectorPupil))
    return SpherePoint(az, alt)


def rotate2d(pos, angle):
    """
    Rotate a 2-dimensional position by a given angle.

    Parameters
    ----------
    pos : pair of `float`
        Position to rotate
    angle : `lsst.afw.geom.Angle`
        Amount of rotation

    Examples
    --------

    ``rotate2d((1, 2), 90*lsst.afw.geom.degrees, False)`` returns `(-2, 1)`
    """
    angRad = angle.asRadians()
    sinAng = math.sin(angRad)
    cosAng = math.cos(angRad)
    return (
        cosAng*pos[0] - sinAng*pos[1],
        sinAng*pos[0] + cosAng*pos[1]
    )
