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

__all__ = ["CoordinateConverter"]

import math

import numpy as np
import scipy.optimize

from . import coordUtils
from lsst.afw.geom import Extent2D, Point2D, SpherePoint, radians
from lsst.afw.cameraGeom import PIXELS, FOCAL_PLANE, FIELD_ANGLE
from .beamInfo import BeamInfo

# Record errors from setFocalFieldAngle?
_RecordErrors = False
# list of errors from the root finder: error (abs value, radians), pupilPos, focalFieldAngle, beam
_ErrorList = None


def startRecordingErrors():
    """Start recording setFocalFieldAngle errors and reset accumulators
    """
    global _RecordErrors, _ErrorList
    _RecordErrors = True
    _ErrorList = []


def stopRecordingErrors():
    """Stop recording errors"""
    global _RecordErrors
    _RecordErrors = False


def getRecordedErrors():
    """Return setFocalFieldAngle errors

    Returns
    -------
    errorLists : `tuple` of three `list`

        - rootFinderErrors: each element is a tuple:

            - |error| in radians
            - pupilPos
            - focalFieldAngle
            - beam name

        - minimizerErrors: each element is a tuple
          that is the same as root finder errors

        - minimizerFailures: each element is a tuple:

            - pupilPos
            - focalFieldAngle
            - beam name
    """
    return _ErrorList


class CoordinateConverter:
    """Coordinate conversions for the collimated beam projector (CBP)

    This object supports the following tasks:

    - Given a desired CBP "arrangement" (e.g. a particular beam
      should fall on a particular spot on the pupil
      and be focused to spot on a particular position of a detector),
      compute the necessary telescope and CBP pointing
      and the information about where each beam is directed.
    - Given the current telescope and CBP pointing and a desired
      offset to the resulting CBP arrangement, compute the new
      telescope and CBP pointing.

    See `Notes` for a summary of how to use this object

    Parameters
    ----------
    config : `lsst.cbp.CoordinateConverterConfig`
        Telescope and CBP :ref:`configuration <lsst.cbp.configuration>`
    maskInfo : `lsst.cbp.MaskInfo`
        Information about the CBP mask
    camera : `lsst.afw.cameraGeom.Camera`
        Camera geometry

    Notes
    -----

    **How to Use This Object**

    Call a ``set`` method, such as `setDetectorPos` to specify
    a desired CBP arrangement, or an offset method, such as
    `offsetDetectorPos`, to adjust the current arrangement.
    These methods update the telescope and CBP pointing,
    and thus the beam information. That is all they do;
    they return nothing and save none of the specifics of your request.

    Read the computed telescope and CBP pointing as attributes
    `telAzAltObserved`, `telRotObserved` and `cbpAzAltObserved`
    and move the telescope and CBP accordingly.

    There are two ways to get information about beams:

    - To get information for all beams, iterate on this object
      to get one `lsst.cbp.BeamInfo` per beam. Also the length
      of this object is the number of beams.
    - To get information for a specific beam, use ``[beamNameOrIndex]``;
      see `__getitem__` for details.
      Also attribute `beamNames` provides an iterable of beam names.

    Whenever the telescope, camera rotator or CBP are moved,
    you must update the appropriate attribute(s) accordingly.
    Otherwise the beam information will be incorrect
    and offset commands will not work as expected.
    After each such update you can read the new beam information
    as usual.

    You are free to update configuration information at any time,
    by setting the appropriate attribute(s).
    The two items you are most likely to update are:

    - `maskInfo`: set this when you change the mask
    - `config.telRotOffset`: set this if you want to correct
      the orientation of the spot pattern on the focal plane.

    After each such update you can read the new beam info and telescope,
    and CBP position to see the effect of the change.

    """
    def __init__(self, config, maskInfo, cameraGeom):
        self.config = config
        self.maskInfo = maskInfo
        self.cameraGeom = cameraGeom
        self._fieldAngleToFocalPlane = cameraGeom.getTransform(FIELD_ANGLE, FOCAL_PLANE)
        # amount to add to default hole position to compute telescope rotator angle
        self._holeDelta = Extent2D(1, 0)
        self._telAzAlt = SpherePoint(np.nan, np.nan, radians)
        self._telRot = np.nan*radians
        self._cbpAzAlt = SpherePoint(np.nan, np.nan, radians)

    def setFocalFieldAngle(self, pupilPos, focalFieldAngle=None, beam=None):
        """Set the focal plane field angle of a beam

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilPos : pair of `float`
            Position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm)
        focalFieldAngle : pair of `float` (optional)
            :ref:`Focal plane field angle <lsst.cbp.focal_plane_field_angle>`
            of the specified beam (x, y rad); defaults to (0, 0)
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to ``self.maskInfo.defaultBeam``
        """
        beam = self.maskInfo.asHoleName(beam)
        if focalFieldAngle is None:
            focalFieldAngle = Point2D(0, 0)
        else:
            focalFieldAngle = Point2D(*focalFieldAngle)

        # if the field angle is too small matter, and too small
        # to determine its orientation, treat it as 0,0 (no iteration required)
        if math.hypot(*focalFieldAngle) < 1e-10:
            self.setPupilFieldAngle(pupilPos, focalFieldAngle, beam)
            return

        # minimize the field angle error as a function of the orientation
        # of the field angle; start with rotation = 0 so pupil field angle
        # equals focal plane field angle

        # record initial conditions, in case the miminizer fails to converge
        telAzAlt = self._telAzAlt
        telRot = self._telRot
        cbpAzAlt = self._cbpAzAlt
        try:
            class TrialFunctor:
                """Functor to compute error in focal plane orientation

                Parameters
                ----------
                cco : `lsst.cbp.CoordinateConverter`
                    A coordinate converter object
                pupilPos : pair of `float`
                    Position of the specified beam on the
                    :ref:`telescope pupil <lsst.cbp.pupil_position>`
                    (x, y mm)
                focalFieldAngle : pair of `float` (optional)
                    :ref:`Focal plane field angle <lsst.cbp.focal_plane_field_angle>`
                    of the specified beam (x, y rad); defaults to (0, 0)
                beam : `int` or `str` (optional)
                    Name or index of beam; defaults to ``self.maskInfo.defaultBeam``
                """
                def __init__(self, cco, pupilPos, focalFieldAngle, beam):
                    self.cco = cco
                    self.pupilPos = pupilPos
                    # desired focal plane field angle
                    self.focalFieldAngle = focalFieldAngle
                    self.beam = beam
                    self.focalFieldAngleOrientation = math.atan2(focalFieldAngle[1],
                                                                 focalFieldAngle[0])*radians

                def __call__(self, rotAngleRadArr):
                    """Compute the error in focal plane orientation
                    (in radians) at the specified camera rotation angle

                    Parameters
                    ----------
                    rotAngleRadArr : sequence of one `float`
                        The internal camera rotator angle, in radians.
                        It is passed in as a sequence of one element,
                        as required by scipy.optimize.
                    """
                    rotAngleRad = rotAngleRadArr[0]
                    pupilFieldAngle = coordUtils.rotate2d(pos=self.focalFieldAngle, angle=rotAngleRad*radians)
                    self.cco.setPupilFieldAngle(pupilPos=self.pupilPos,
                                                pupilFieldAngle=pupilFieldAngle,
                                                beam=self.beam)
                    measuredFocalFieldAngle = self.cco[beam].focalFieldAngle
                    measuredFocalFieldAngleOrientation = math.atan2(measuredFocalFieldAngle[1],
                                                                    measuredFocalFieldAngle[0])*radians
                    return self.focalFieldAngleOrientation.separation(
                        measuredFocalFieldAngleOrientation).asRadians()

            funcToFindRoot = TrialFunctor(cco=self, pupilPos=pupilPos,
                                          focalFieldAngle=focalFieldAngle, beam=beam)

            # Fit the rotator angle using a root finder
            with np.errstate(divide="ignore", invalid="ignore"):
                iterResult = scipy.optimize.root(fun=funcToFindRoot, x0=np.array([0.0], dtype=float),
                                                 options=dict(fatol=1e-11), method="broyden1")

            if not iterResult.success:
                raise RuntimeError("Iteration failed to converge")
            # call the function again to make sure the final value found is the one that is used
            err = funcToFindRoot(iterResult.x)

            global _RecordErrors, _ErrorList
            if _RecordErrors:
                _ErrorList.append((abs(err), pupilPos, focalFieldAngle, beam))

        except Exception:
            self._telAzAlt = telAzAlt
            self._telRot = telRot
            self._cbpAzAlt = cbpAzAlt
            raise

    def setFocalPlanePos(self, pupilPos, focalPlanePos=None, beam=None):
        """Set the position of a spot on the focal plane

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilPos : pair of `float`
            Position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm)
        focalPlanePos : pair of `float` (optional)
            :ref:`Focal plane position <lsst.cbp.focal_plane>` of the spot
            formed by the specified beam (x, y mm); defaults to (0, 0)
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam
        """
        beam = self.maskInfo.asHoleName(beam)
        if focalPlanePos is None:
            focalPlanePos = Point2D(0, 0)
        else:
            focalPlanePos = Point2D(*focalPlanePos)
        focalFieldAngle = self._fieldAngleToFocalPlane.applyInverse(focalPlanePos)
        self.setFocalFieldAngle(pupilPos=pupilPos, focalFieldAngle=focalFieldAngle, beam=beam)

    def setDetectorPos(self, pupilPos, detectorPos=None, detector=None, beam=None):
        """Set the position of a spot on a detector

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilPos : pair of `float`
            Position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm)
        detectorPos : pair of `float` (optional)
            Position of the spot formed by the specified beam
            on the specified detector (x, y pixels);
            defaults to the center of the detector
        detector : `str` (optional
            Name of detector; defaults to self.config.defaultDetector
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam
        """
        beam = self.maskInfo.asHoleName(beam)
        if detector is None:
            detector = self.config.defaultDetector
        detectorInfo = self.cameraGeom[detector]
        if detectorPos is None:
            detectorPos = detectorInfo.getCenter(PIXELS)
        else:
            detectorPos = Point2D(*detectorPos)
        pixelSys = detectorInfo.makeCameraSys(PIXELS)
        pixelsToFieldAngle = self.cameraGeom.getTransform(pixelSys, FIELD_ANGLE)
        focalFieldAngle = pixelsToFieldAngle.applyForward(detectorPos)
        self.setFocalFieldAngle(pupilPos=pupilPos, focalFieldAngle=focalFieldAngle, beam=beam)

    def offsetDetectorPos(self, pupilOffset=None,
                          detectorOffset=None, beam=None):
        """Offset the detector position and/or pupil position of a beam

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilOffset : pair of `float` (optional)
            Offset of the position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm)
        detectorOffset : pair of `float` (optional)
            Offset of the position of the specified spot
            on the detector it is presently on (x, y pixels);
            defaults to (0, 0)
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam
        """
        beamInfo = self[beam]
        if not beamInfo.isOnDetector:
            raise RuntimeError("This beam is not on a detector")
        if pupilOffset is None:
            pupilOffset = (0, 0)
        pupilOffset = Extent2D(*pupilOffset)
        if detectorOffset is None:
            detectorOffset = (0, 0)
        detectorOffset = Extent2D(*detectorOffset)
        newPupilPos = beamInfo.pupilPos + pupilOffset
        newDetectorPos = beamInfo.detectorPos + detectorOffset
        self.setDetectorPos(pupilPos=newPupilPos,
                            detectorPos=newDetectorPos,
                            detector=beamInfo.detectorName,
                            beam=beam)

    def offsetFocalPlanePos(self, pupilOffset=None,
                            focalPlaneOffset=None, beam=None):
        """Offset the position focal plane position and/or
        pupil position of a beam

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilOffset : pair of `float` (optional)
            Offset of the position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm)
        focalPlaneOffset : pair of `float` (optional)
            Offset of the position of the specified spot
            on the :ref:`focal plane <lsst.cbp.focal_plane>` (x, y mm);
            defaults to (0, 0)
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam
        """
        beamInfo = self[beam]
        if pupilOffset is None:
            pupilOffset = (0, 0)
        pupilOffset = Extent2D(*pupilOffset)
        if focalPlaneOffset is None:
            focalPlaneOffset = (0, 0)
        focalPlaneOffset = Extent2D(*focalPlaneOffset)
        newPupilPos = beamInfo.pupilPos + pupilOffset
        newFocalPlanePos = beamInfo.focalPlanePos + focalPlaneOffset
        self.setFocalPlanePos(pupilPos=newPupilPos,
                              focalPlanePos=newFocalPlanePos,
                              beam=beam)

    def offsetFocalFieldAngle(self, pupilOffset=None,
                              focalFieldAngleOffset=None, beam=None):
        """Offset the focal plane field angle and/or pupil position
        of a beam

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilOffset : pair of `float` (optional)
            Offset of the position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm)
        focalFieldAngleOffset : pair of `float` (optional)
            Offset of the
            :ref:`focal plane field angle <lsst.cbp.focal_plane_field_angle>`
            of the specified beam (x, y mm); defaults to (0, 0)
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam
        """
        beamInfo = self[beam]
        if pupilOffset is None:
            pupilOffset = (0, 0)
        pupilOffset = Extent2D(*pupilOffset)
        if focalFieldAngleOffset is None:
            focalFieldAngleOffset = (0, 0)
        focalFieldAngleOffset = Extent2D(*focalFieldAngleOffset)
        newPupilPos = beamInfo.pupilPos + pupilOffset
        newFocalFieldAngle = beamInfo.focalFieldAngle + focalFieldAngleOffset
        self.setFocalFieldAngle(pupilPos=newPupilPos,
                                focalFieldAngle=newFocalFieldAngle,
                                beam=beam)

    def __getitem__(self, beam):
        """Dict-like access to beam information

        Get a BeamInfo for a beam specified by integer index or name
        """
        beam = self.maskInfo.asHoleName(beam)
        return self.getBeamInfo(beam=beam)

    def __iter__(self):
        """Iterator over beam information

        Do not modify this object during iteration, e.g. by modifying
        attributes or calling set or offset methods.
        """
        for beam in self.beamNames:
            yield self[beam]

    def __len__(self):
        """The number of beams"""
        return self.maskInfo.numHoles

    @property
    def cbpInBounds(self):
        """Return True if the observed CBP altitude is in bounds
        """
        alt = self.cbpAzAltObserved[1]
        return self.config.cbpAltitudeLimits[0] <= alt <= self.config.cbpAltitudeLimits[1]

    @property
    def telInBounds(self):
        alt = self.telAzAltObserved[1]
        return self.config.telAltitudeLimits[0] <= alt <= self.config.telAltitudeLimits[1]

    @property
    def beamNames(self):
        """An `iterable` of beam names, in index order"""
        return self.maskInfo.holeNames

    @property
    def cbpAzAltObserved(self):
        """The observed az/alt pointing of the CBP
        as an `lsst.afw.geom.SpherePoint`
        """
        return SpherePoint(*[self._internalToObserved(
            internal=self._cbpAzAlt[i],
            offset=self.config.cbpAzAltOffset[i],
            scale=self.config.cbpAzAltScale[i]) for i in range(2)])

    @cbpAzAltObserved.setter
    def cbpAzAltObserved(self, cbpObs):
        """Set the observed az/alt pointing of the CBP

        Parameters
        ----------
        cbpObs : `lsst.afw.geom.SpherePoint`
            Observed az/alt of the CBP
        """
        self._cbpAzAlt = SpherePoint(*[self._observedToInternal(
            observed=cbpObs[i],
            offset=self.config.cbpAzAltOffset[i],
            scale=self.config.cbpAzAltScale[i]) for i in range(2)])

    @property
    def telAzAltObserved(self):
        """The observed az/alt pointing of the telescope
        as an `lsst.afw.geom.SpherePoint`
        """
        return SpherePoint(*[self._internalToObserved(
            internal=self._telAzAlt[i],
            offset=self.config.telAzAltOffset[i],
            scale=self.config.telAzAltScale[i]) for i in range(2)])

    @telAzAltObserved.setter
    def telAzAltObserved(self, telObs):
        """Set the observed az/alt pointing of the telescope

        Parameters
        ----------
        telObs : `lsst.afw.geom.SpherePoint`
            Observed az/alt of the telescope
        """
        self._telAzAlt = SpherePoint(*[self._observedToInternal(
            observed=telObs[i],
            offset=self.config.telAzAltOffset[i],
            scale=self.config.telAzAltScale[i]) for i in range(2)])

    @property
    def telRotObserved(self):
        """The observed angle of the telescope camera rotator
        as an `lsst.afw.geom.Angle`
        """
        return self._internalToObserved(
            internal=self._telRot,
            offset=self.config.telRotOffset,
            scale=self.config.telRotScale)

    @telRotObserved.setter
    def telRotObserved(self, telRotObserved):
        """Set the observed angle of the telescope camera rotator

        Parametres
        ----------
        telRotObserved : `lsst.afw.geom.Angle`
            The observed angle of the telescope camera rotator
        """
        self._telRot = self._observedToInternal(
            observed=telRotObserved,
            offset=self.config.telRotOffset,
            scale=self.config.telRotScale)

    @property
    def cbpAzAltInternal(self):
        """Return the internal az/alt pointing of the CBP

        Primarily intended for testing
        """
        return self._cbpAzAlt

    @property
    def telAzAltInternal(self):
        """Return the internal az/alt pointing of the telescope

        Primarily intended for testing
        """
        return self._telAzAlt

    @property
    def telRotInternal(self):
        return self._telRot
        """Return the internal angle of the telescope camera rotator

        Primarily intended for testing
        """

    def setPupilFieldAngle(self, pupilPos, pupilFieldAngle=None, beam=None):
        """Set the pupil field angle and pupil position of a beam

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        This method is primarily intended for internal use,
        to support the other set methods. It is public so it can be
        unit-tested.

        Parameters
        ----------
        pupilPos : pair of `float`
            Position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm)
        pupilFieldAngle : pair of `float` (optional)
            Pupil field angle of specified beam (x, y rad);
            defaults to (0, 0)
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam
        """
        beam = self.maskInfo.asHoleName(beam)
        if pupilFieldAngle is None:
            pupilFieldAngle = Point2D(0, 0)
        else:
            pupilFieldAngle = Point2D(*pupilFieldAngle)
        beamPosAtCtr = coordUtils.computeShiftedPlanePos(pupilPos, pupilFieldAngle,
                                                         -self.config.telPupilOffset)
        beamVectorInCtrPupil = self._computeBeamVectorInCtrPupilFrame(
            beamPosAtCtr=beamPosAtCtr, pupilFieldAngle=pupilFieldAngle)
        cbpVectorInCtrPupil = self._computeCbpVectorInCtrPupilFrame(
            beamPosAtCtr=beamPosAtCtr,
            beamVectorInCtrPupil=beamVectorInCtrPupil)

        telAzAlt = coordUtils.computeAzAltFromBasePupil(
            vectorBase=self.config.cbpPosition,
            vectorPupil=cbpVectorInCtrPupil)

        beamVectorBase = coordUtils.convertVectorFromPupilToBase(
            vectorPupil=beamVectorInCtrPupil,
            pupilAzAlt=telAzAlt,
        )
        beamFieldAngleCbp = self._getBeamCbpFieldAngle(beam)
        beamUnitVectorCbpPupil = coordUtils.fieldAngleToVector(beamFieldAngleCbp, self.config.cbpFlipX)
        cbpAzAlt = coordUtils.computeAzAltFromBasePupil(
            vectorBase=-beamVectorBase,
            vectorPupil=beamUnitVectorCbpPupil,
        )

        self._telRot = self._computeCameraRotatorAngle(telAzAlt=telAzAlt, cbpAzAlt=cbpAzAlt)
        self._telAzAlt = telAzAlt
        self._cbpAzAlt = cbpAzAlt

    def _computeCameraRotatorAngle(self, telAzAlt, cbpAzAlt):
        """Compute the internal camera rotator angle needed a given
        telescope and CBP pointing

        Parameters
        ----------
        telAzAlt : `lsst.afw.geom.SpherePoint`
            Telescope internal azimuth and altitude
        cbpAzAlt : `lsst.afw.geom.SpherePoint`
            CBP internal azimuth and altitude
        """
        # compute focal plane position, ignoring self._telRot,
        # for two holes separated by x in the CBP equidistant from the center;
        # compute the angle that would make the spots line up with the x axis in the focal plane
        ctrHolePos = Point2D(0, 0)
        holeDelta = Extent2D(*coordUtils.getFlippedPos(self._holeDelta, flipX=self.config.cbpFlipX))
        holePos1 = ctrHolePos - holeDelta
        holePos2 = ctrHolePos + holeDelta
        pupilUnitVector1 = self._computeTelPupilUnitVectorFromHolePos(holePos1, telAzAlt=telAzAlt,
                                                                      cbpAzAlt=cbpAzAlt)
        pupilUnitVector2 = self._computeTelPupilUnitVectorFromHolePos(holePos2, telAzAlt=telAzAlt,
                                                                      cbpAzAlt=cbpAzAlt)
        # Rotation is done in a right-handed system, regardless of telFlipX
        pupilFieldAngle1 = coordUtils.vectorToFieldAngle(pupilUnitVector1, flipX=False)
        pupilFieldAngle2 = coordUtils.vectorToFieldAngle(pupilUnitVector2, flipX=False)
        focalPlane1 = self._fieldAngleToFocalPlane.applyForward(Point2D(*pupilFieldAngle1))
        focalPlane2 = self._fieldAngleToFocalPlane.applyForward(Point2D(*pupilFieldAngle2))
        deltaFocalPlane = np.subtract(focalPlane2, focalPlane1)
        return -math.atan2(deltaFocalPlane[1], deltaFocalPlane[0])*radians

    def _getBeamCbpFieldAngle(self, beam):
        """Return the field angle of the specified beam with respect to
        the CBP pupil frame

        Parameters
        ----------
        beam : `int`, `str` or None
            Name or index of beam; if None then
            ``self.maskInfo.defaultBeam``

        Returns
        -------
        fieldAngle : a pair of floats, in radians
        """
        holePos = self.maskInfo.getHolePos(beam)
        return tuple(math.atan(pos / self.config.cbpFocalLength) for pos in holePos)

    def _computeBeamVectorInCtrPupilFrame(self, beamPosAtCtr, pupilFieldAngle):
        """Compute the beam vector to the CBP in the centered pupil frame

        beamPosAtCtr: position of beam on centered pupil (x,y mm)
        pupilFieldAngle: incident angle of beam on pupil (x,y rad)
        cbpDistance: distance between telescope and CBP centers (mm)
        """
        # compute beam position on "center pupil": a plane parallel to pupil at center of telescope
        beamPosVec = coordUtils.pupilPositionToVector(beamPosAtCtr, self.config.telFlipX)
        beamUnitVec = coordUtils.fieldAngleToVector(pupilFieldAngle, self.config.telFlipX)
        abyz = beamPosVec[1]*beamUnitVec[1] + beamPosVec[2]*beamUnitVec[2]
        beamPosMag = np.linalg.norm(beamPosAtCtr)
        cbpDistance = self.config.cbpDistance
        beamLength = -abyz + math.sqrt(math.fsum((cbpDistance**2, abyz**2, -beamPosMag**2)))
        return beamLength*np.array(beamUnitVec)

    def _computeCbpVectorInCtrPupilFrame(self, beamPosAtCtr, beamVectorInCtrPupil):
        """Compute the vector from telescope to CBP in the telescope's
        centered pupil frame

        beamPosAtCtr: position of beam on centered pupil (x,y mm)
        beamAngle: incident angle of beam on pupil (x,y rad)
        """
        beamPosVec = coordUtils.pupilPositionToVector(beamPosAtCtr, self.config.telFlipX)
        return beamVectorInCtrPupil + beamPosVec

    def _rotateFocalPlaneToPupil(self, focalPlane):
        """Rotate a position or field angle in telescope focal plane
        orientation to one in telescope pupil orientation
        """
        return coordUtils.rotate2d(pos=focalPlane, angle=-self._telRot)

    def _rotatePupilToFocalPlane(self, pupil):
        """Rotate a position or field angle in telescope pupil orientation
        to one in telescope focal plane orientation
        """
        return coordUtils.rotate2d(pos=pupil, angle=self._telRot)

    def getBeamInfo(self, beam, *, holePos=None):
        """Get beam information for a beam from a specified CBP
        beam or hole position.

        You may specify a hole position. This can be useful for unit
        tests and "what if" scenarios.

        Parameters
        ----------
        beam : `str`
            Beam name; defaults to ""
        holePos : pair of `float` (optional)
            Hole position on CBP mask (x, y mm);
            defaults to the actual hole position of the named beam
        """
        if holePos is None:
            holePos = self.maskInfo.getHolePos(beam)

        # compute focal plane field angle of the beam
        beamPupilUnitVector = self._computeTelPupilUnitVectorFromHolePos(holePos, telAzAlt=self._telAzAlt,
                                                                         cbpAzAlt=self._cbpAzAlt)
        pupilFieldAngle = coordUtils.vectorToFieldAngle(beamPupilUnitVector, self.config.telFlipX)
        focalFieldAngle = self._rotatePupilToFocalPlane(pupilFieldAngle)

        # compute focal plane position of the beam
        focalPlanePos = self._fieldAngleToFocalPlane.applyForward(Point2D(*pupilFieldAngle))
        isOnFocalPlane = math.hypot(*focalPlanePos) < self.config.telFocalPlaneDiameter

        # vector from telescope centered pupil to actual pupil
        pupilOffset = np.array((self.config.telPupilOffset, 0, 0), dtype=float)

        # compute the pupil position of the beam on the actual telescope pupil
        # (first compute on the centered pupil, then subtract the pupil offset)
        cbpPositionPupil = coordUtils.convertVectorFromBaseToPupil(
            vectorBase=self.config.cbpPosition,
            pupilAzAlt=self._telAzAlt,
        ) - pupilOffset
        pupilNormalVector = np.array((1, 0, 0), dtype=float)
        pupilDistance = np.dot(cbpPositionPupil, pupilNormalVector) / \
            np.dot(beamPupilUnitVector, pupilNormalVector)
        pupilPosVector = cbpPositionPupil - beamPupilUnitVector*pupilDistance
        # the x component should be zero; y, z -> plane position ±x, y
        pupilPos = coordUtils.getFlippedPos((pupilPosVector[1], pupilPosVector[2]),
                                            flipX=self.config.telFlipX)
        pupilPosDiameter = math.hypot(*pupilPos)
        isOnPupil = self.config.telPupilObscurationDiameter <= pupilPosDiameter and \
            pupilPosDiameter <= self.config.telPupilDiameter
        return BeamInfo(
            cameraGeom=self.cameraGeom,
            name=beam,
            holePos=holePos,
            isOnPupil=isOnPupil,
            isOnFocalPlane=isOnFocalPlane,
            focalPlanePos=focalPlanePos,
            pupilPos=pupilPos,
            focalFieldAngle=focalFieldAngle,
            pupilFieldAngle=pupilFieldAngle,
        )

    def _computeCbpPupilUnitVectorFromHolePos(self, holePos):
        beamCbpFieldAngle = [math.atan(pos/self.config.cbpFocalLength) for pos in holePos]
        return coordUtils.fieldAngleToVector(beamCbpFieldAngle, self.config.cbpFlipX)

    def _computeTelPupilUnitVectorFromHolePos(self, holePos, telAzAlt, cbpAzAlt):
        """Compute the telescope pupil unit vector of a beam
        from its CBP hole position

        Parameters
        ----------
        telAzAlt : `lsst.afw.geom.SpherePoint`
            Telescope internal azimuth and altitude
        cbpAzAlt : `lsst.afw.geom.SpherePoint`
            CBP internal azimuth and altitude
        """
        beamCbpPupilUnitVector = self._computeCbpPupilUnitVectorFromHolePos(holePos)
        beamCbpBaseUnitVector = coordUtils.convertVectorFromPupilToBase(
            vectorPupil=beamCbpPupilUnitVector,
            pupilAzAlt=cbpAzAlt,
        )
        beamTelBaseUnitVector = -beamCbpBaseUnitVector
        return coordUtils.convertVectorFromBaseToPupil(
            vectorBase=beamTelBaseUnitVector,
            pupilAzAlt=telAzAlt,
        )

    def _observedToInternal(self, observed, offset, scale):
        """Convert an observed angle into an internal angle

        Returns observed/scale - offset
        """
        return (observed - offset)/scale

    def _internalToObserved(self, internal, offset, scale):
        """Convert an internal angle into an observed angle

        Returns (internal + offset) * scale
        """
        return internal*scale + offset
