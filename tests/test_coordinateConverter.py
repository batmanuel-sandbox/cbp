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

import itertools
import math
import unittest

import numpy as np

from lsst.sphgeom import Vector3d
from lsst.afw.cameraGeom import FIELD_ANGLE, FOCAL_PLANE
from lsst.afw.geom import Box2D, Point2D, SpherePoint, arcseconds, degrees, radians
import lsst.cbp.coordUtils as coordUtils
from lsst.cbp.testUtils import SampleCoordinateConverter
import lsst.cbp.coordinateConverter
import lsst.utils.tests

# set True to record and report numeric error in convertVectorFromPupilToBase
# and in the iterative component of setFocalFieldAngle
ReportRecordedErrors = False


class CoordConverterTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        if ReportRecordedErrors:
            lsst.cbp.coordUtils.startRecordingErrors()
            lsst.cbp.coordinateConverter.startRecordingErrors()
        self.scc = SampleCoordinateConverter()
        self.cco = self.scc.coordinateConverter
        # a list of all detector names plus None for the default detector
        self.detectorNames = itertools.chain([None], self.cco.cameraGeom.getNameIter())
        # set values for maximum error that are "plenty good enough"
        self.maxPupilPosErr = 1e-4  # mm
        self.maxDetectorPosErr = 1e-4  # pixels
        self.maxFocalPlanePosErr = self.scc.pixelSizeMm * self.maxDetectorPosErr  # mm
        self.maxFieldAngleErrRad = (0.001*arcseconds).asRadians()

    def tearDown(self):
        if ReportRecordedErrors:
            utilsErrList = coordUtils.getRecordedErrors()
            coordUtils.stopRecordingErrors()
            if len(utilsErrList) > 0:
                print("\nWarning: recorded %s significant errors in computeAzAltFromBasePupil;"
                      " the worst 5 are:" % len(utilsErrList))
                for elt in utilsErrList[-5:]:
                    print("error=%0.5f arcsec, vectorBase=%s, vectorPupil=%s" % elt)

            errorList = lsst.cbp.coordinateConverter.getRecordedErrors()
            if len(errorList) > 0:
                errorArray = np.array([val[0] for val in errorList])
                print("\nroot finder mean error=%0.3g, sdev=%0.3g, min=%0.3g, max=%0.3g" %
                      (errorArray.mean(), errorArray.std(),
                       errorArray.min(), errorArray.max()))
        del self.cco
        del self.scc

    def testOffsetDetectorPos(self):
        pupilPos = (3000, 4000)
        detectorPos = (400, 550)
        for pupilOffset, detectorOffset, detector, beam in itertools.product(
            (None, (0, 0), (500, -300)),
            (None, (0, 0), (-50, 350)),
            self.detectorNames,
            self.cco.beamNames,
        ):
            with self.subTest(pupilOffset=pupilOffset, detectorPos=detectorPos, detector=detector, beam=beam):
                self.cco.setDetectorPos(pupilPos=pupilPos, detectorPos=detectorPos, detector=detector,
                                        beam=beam)
                initialBeamInfo = self.cco[beam]
                self.assertPairsAlmostEqual(initialBeamInfo.pupilPos, pupilPos, self.maxPupilPosErr)
                self.assertPairsAlmostEqual(initialBeamInfo.detectorPos, detectorPos,
                                            maxDiff=self.maxDetectorPosErr)
                self.cco.offsetDetectorPos(pupilOffset=pupilOffset, detectorOffset=detectorOffset,
                                           beam=beam)
                finalBeamInfo = self.cco[beam]
                if pupilOffset is None:
                    desiredPupilPos = pupilPos
                else:
                    desiredPupilPos = np.add(pupilPos, pupilOffset)
                if detectorOffset is None:
                    desiredDetectorPos = detectorPos
                else:
                    desiredDetectorPos = np.add(detectorPos, detectorOffset)
                self.assertPairsAlmostEqual(finalBeamInfo.pupilPos, desiredPupilPos, self.maxPupilPosErr)
                self.assertPairsAlmostEqual(finalBeamInfo.detectorPos, desiredDetectorPos,
                                            maxDiff=self.maxDetectorPosErr)

    def testOffsetFocalFieldAngle(self):
        pupilPos = (3000, 4000)
        focalFieldAngle = (-0.03, 0.04)
        for pupilOffset, focalFieldAngleOffset, beam in itertools.product(
            (None, (0, 0), (500, -300)),
            (None, (0, 0), (0.02, -0.03)),
            self.cco.beamNames,
        ):
            with self.subTest(pupilOffset=pupilOffset, focalFieldAngleOffset=focalFieldAngleOffset,
                              beam=beam):
                self.cco.setFocalFieldAngle(pupilPos=pupilPos, focalFieldAngle=focalFieldAngle, beam=beam)
                initialBeamInfo = self.cco[beam]
                self.assertPairsAlmostEqual(initialBeamInfo.pupilPos, pupilPos, self.maxPupilPosErr)
                self.assertPairsAlmostEqual(initialBeamInfo.focalFieldAngle, focalFieldAngle)

                self.cco.offsetFocalFieldAngle(pupilOffset=pupilOffset,
                                               focalFieldAngleOffset=focalFieldAngleOffset, beam=beam)
                finalBeamInfo = self.cco[beam]
                if pupilOffset is None:
                    desiredPupilPos = pupilPos
                else:
                    desiredPupilPos = np.add(pupilPos, pupilOffset)
                if focalFieldAngleOffset is None:
                    desiredFocalFieldAngle = focalFieldAngle
                else:
                    desiredFocalFieldAngle = np.add(focalFieldAngle, focalFieldAngleOffset)
                self.assertPairsAlmostEqual(finalBeamInfo.pupilPos, desiredPupilPos, self.maxPupilPosErr)
                self.assertPairsAlmostEqual(finalBeamInfo.focalFieldAngle, desiredFocalFieldAngle,
                                            maxDiff=self.maxFieldAngleErrRad)

    def testOffsetFocalPlanePos(self):
        pupilPos = (3000, 4000)
        focalPlanePos = (-400, 550)
        for pupilOffset, focalPlaneOffset, beam in itertools.product(
            (None, (0, 0), (500, -300)),
            (None, (0, 0), (22, -350)),
            self.cco.beamNames,
        ):
            with self.subTest(pupilOffset=pupilOffset, focalPlanePos=focalPlanePos, beam=beam):
                self.cco.setFocalPlanePos(pupilPos=pupilPos, focalPlanePos=focalPlanePos, beam=beam)
                initialBeamInfo = self.cco[beam]
                self.assertPairsAlmostEqual(initialBeamInfo.pupilPos, pupilPos, self.maxPupilPosErr)
                self.assertPairsAlmostEqual(initialBeamInfo.focalPlanePos, focalPlanePos,
                                            maxDiff=self.maxFocalPlanePosErr)
                self.cco.offsetFocalPlanePos(pupilOffset=pupilOffset, focalPlaneOffset=focalPlaneOffset,
                                             beam=beam)
                finalBeamInfo = self.cco[beam]
                if pupilOffset is None:
                    desiredPupilPos = pupilPos
                else:
                    desiredPupilPos = np.add(pupilPos, pupilOffset)
                if focalPlaneOffset is None:
                    desiredFocalPlanePos = focalPlanePos
                else:
                    desiredFocalPlanePos = np.add(focalPlanePos, focalPlaneOffset)
                self.assertPairsAlmostEqual(finalBeamInfo.pupilPos, desiredPupilPos, self.maxPupilPosErr)
                self.assertPairsAlmostEqual(finalBeamInfo.focalPlanePos, desiredFocalPlanePos,
                                            maxDiff=self.maxFocalPlanePosErr)

    def testSampleBasics(self):
        """Test basic elements of the sample coordinate converter
        """
        self.assertEqual(len(self.cco.cameraGeom), len(self.scc.detectorFracPosList))
        self.assertEqual(list(self.cco.cameraGeom.getNameIter()),
                         ["D" + str(i) for i in range(len(self.scc.detectorFracPosList))])
        self.assertEqual(len(self.cco), 6)
        self.assertEqual(tuple(self.cco.beamNames), ("0", "1", "2", "3", "4", "5"))
        self.assertEqual(self.cco.maskInfo.numHoles, 6)
        self.assertPairsAlmostEqual(self.cco.maskInfo.getHolePos("0"), (0, 0))
        for beamInfo, beam in zip(self.cco, self.cco.beamNames):
            self.assertEqual(beamInfo.name, beam)

    def testSetPupilFieldAngleTrivial(self):
        """Test setPupilFieldAngle for the trivial case of hole 0
        aimed perpendicular to the center of the pupil
        """
        self.cco.setPupilFieldAngle(pupilPos=(0, 0))

        # the telescope should be pointed at the center of the CBP and vice-versa
        # NOTE: It would be nice to get better than the 0.0028" that I measure
        self.assertSpherePointsAlmostEqual(self.cco.telAzAltInternal,
                                           SpherePoint(Vector3d(*self.cco.config.cbpPosition)),
                                           maxSep=0.01*arcseconds)
        self.assertSpherePointsAlmostEqual(self.cco.cbpAzAltInternal,
                                           SpherePoint(Vector3d(*(-self.cco.config.cbpPosition))),
                                           maxSep=0.01*arcseconds)

        self.assertAnglesAlmostEqual(self.cco.telRotInternal, 0*degrees)

        # beam 0 should be pointed to the center of the pupil, normal to the pupil
        # and land on the center of the focal plane and the center of detector D0
        beamInfo0 = self.cco["0"]
        self.assertEqual(beamInfo0.name, "0")
        self.assertPairsAlmostEqual(beamInfo0.holePos, (0, 0))
        self.assertFalse(beamInfo0.isOnPupil)  # blocked by the central obscuration
        self.assertTrue(beamInfo0.isOnFocalPlane)
        self.assertTrue(beamInfo0.isOnDetector)
        self.assertPairsAlmostEqual(beamInfo0.focalPlanePos, (0, 0))
        self.assertPairsAlmostEqual(beamInfo0.focalFieldAngle, (0, 0))
        self.assertPairsAlmostEqual(beamInfo0.pupilFieldAngle, (0, 0))
        self.assertPairsAlmostEqual(beamInfo0.pupilPos, (0, 0), self.maxPupilPosErr)
        self.assertEqual(beamInfo0.detectorName, "D0")
        bboxd = Box2D(self.cco.cameraGeom["D0"].getBBox())
        detectorCtrPos = bboxd.getCenter()
        detector34Pos = bboxd.getMin() + bboxd.getDimensions()*0.75
        self.assertPairsAlmostEqual(beamInfo0.detectorPos, detectorCtrPos)
        self.assertTrue(beamInfo0.isOnDetector)

        # beam 1 should on detector D0, 3/4 of the way from LL to UR
        beamInfo1 = self.cco["1"]
        self.assertTrue(beamInfo1.isOnDetector)
        self.assertEqual(beamInfo1.detectorName, "D0")
        self.assertPairsAlmostEqual(beamInfo1.detectorPos, detector34Pos)

        # beam 2 should land on the center of detector D1
        beamInfo2 = self.cco["2"]
        self.assertTrue(beamInfo2.isOnDetector)
        self.assertEqual(beamInfo2.detectorName, "D1")
        self.assertPairsAlmostEqual(beamInfo2.detectorPos, detectorCtrPos)

        # beam 3 should on detector D1, 3/4 of the way from LL to UR
        beamInfo3 = self.cco["3"]
        self.assertTrue(beamInfo3.isOnDetector)
        self.assertEqual(beamInfo3.detectorName, "D1")
        self.assertPairsAlmostEqual(beamInfo3.detectorPos, detector34Pos)

        # beam 4 should land on the center of detector D2
        beamInfo4 = self.cco["4"]
        self.assertTrue(beamInfo4.isOnDetector)
        self.assertEqual(beamInfo4.detectorName, "D2")
        self.assertPairsAlmostEqual(beamInfo4.detectorPos, detectorCtrPos, maxDiff=1e-4)

        # beam 5 should on detector D2, 3/4 of the way from LL to UR
        # The measured error is 5e-7 pixels, which is small enough not to worry.
        # I strongly suspect it is due to inaccuracy in the inverse of the
        # field angle to focal plane transform.
        beamInfo5 = self.cco["5"]
        self.assertTrue(beamInfo5.isOnDetector)
        self.assertEqual(beamInfo5.detectorName, "D2")
        self.assertPairsAlmostEqual(beamInfo5.detectorPos, detector34Pos, maxDiff=1e-4)

    def testSetFocalFieldAngle(self):
        fieldAngleToFocalPlane = self.cco.cameraGeom.getTransform(FIELD_ANGLE, FOCAL_PLANE)
        for focalFieldAngle in ((0, 0), (0, 0.05), (-0.05, -0.03)):
            desiredFocalPlanePos = fieldAngleToFocalPlane.applyForward(Point2D(*focalFieldAngle))
            for pupilPos, beam in itertools.product(
                ((0, 0), (0, 5000), (-5000, 0), (5000, -5000)),
                self.cco.beamNames,
            ):
                with self.subTest(focalFieldAngle=focalFieldAngle, pupilPos=pupilPos, beam=beam):
                    self.cco.setFocalFieldAngle(pupilPos=pupilPos, focalFieldAngle=focalFieldAngle, beam=beam)
                    beamInfo = self.cco[beam]
                    self.assertPairsAlmostEqual(beamInfo.pupilPos, pupilPos, self.maxPupilPosErr)
                    self.assertPairsAlmostEqual(beamInfo.focalFieldAngle, focalFieldAngle,
                                                maxDiff=self.maxFieldAngleErrRad)
                    self.assertPairsAlmostEqual(beamInfo.focalPlanePos, desiredFocalPlanePos,
                                                maxDiff=self.maxFocalPlanePosErr)
                    self.checkOrientation()

    def testSetDetectorPos(self):
        for detectorPos, pupilPos, detector in itertools.product(
            ((0, 0), (500, 1000), (25, 1850)),
            ((0, 0), (0, 5000), (-5000, 0), (5000, -5000)),
            self.detectorNames,
        ):
            with self.subTest(detectorPos=detectorPos, pupilPos=pupilPos, detector=detector):
                if detector is None:
                    desiredDetector = self.cco.config.defaultDetector
                else:
                    desiredDetector = detector
                for pupilPos in (
                ):
                    for beam in range(4):
                        self.cco.setDetectorPos(pupilPos=pupilPos, detectorPos=detectorPos, detector=detector,
                                                beam=beam)
                        beamInfo = self.cco[beam]
                        self.assertTrue(beamInfo.isOnDetector)
                        self.assertEqual(beamInfo.detectorName, desiredDetector)
                        self.assertPairsAlmostEqual(beamInfo.pupilPos, pupilPos,
                                                    maxDiff=self.maxPupilPosErr)
                        self.assertPairsAlmostEqual(beamInfo.detectorPos, detectorPos,
                                                    maxDiff=self.maxDetectorPosErr)
                        self.checkOrientation()

    def testSetFocalPlanePos(self):
        maxErrMm = self.scc.pixelSizeMm * 0.0001
        for focalPlanePos, pupilPos, beam in itertools.product(
            ((0, 0), (500, 200), (25, 850)),
            ((0, 0), (0, 5000), (-5000, 0), (5000, -5000)),
            self.cco.beamNames,
        ):
            with self.subTest(focalPlanePos=focalPlanePos, pupilPos=pupilPos, beam=beam):
                self.cco.setFocalPlanePos(pupilPos=pupilPos, focalPlanePos=focalPlanePos, beam=beam)
                beamInfo = self.cco[beam]
                self.assertTrue(beamInfo.isOnFocalPlane)
                errMm = math.hypot(*np.subtract(beamInfo.focalPlanePos, focalPlanePos))
                self.assertLess(errMm, maxErrMm)
                self.assertPairsAlmostEqual(beamInfo.pupilPos, pupilPos, self.maxPupilPosErr)

    def testSetPupilFieldAngleZero(self):
        """Test setPupilFieldAngle for zero field angle and various points on the pupil
        """
        for pupilPos in ((0, 5000), (-5000, 0), (5000, -5000)):
            with self.subTest(pupilPos=pupilPos):
                self.cco.setPupilFieldAngle(pupilPos=pupilPos)

                # the telescope should be pointing in the opposite direction of the CBP
                telDir = self.cco.telAzAltInternal.getVector()
                cbpDir = self.cco.cbpAzAltInternal.getVector()
                negativeCbpDir = -np.array(cbpDir, dtype=float)
                np.testing.assert_allclose(telDir, negativeCbpDir, atol=1e-15)

                maxDetectorPosError = 0.0001

                # beam 0 should be pointed to the center of the pupil, normal to the pupil
                # and land on the center of the focal plane and the center of detector D0
                beamInfo0 = self.cco["0"]
                self.assertEqual(beamInfo0.name, "0")
                self.assertPairsAlmostEqual(beamInfo0.holePos, (0, 0))
                self.assertTrue(beamInfo0.isOnPupil)
                self.assertTrue(beamInfo0.isOnFocalPlane)
                self.assertTrue(beamInfo0.isOnDetector)
                self.assertPairsAlmostEqual(beamInfo0.focalPlanePos, (0, 0), maxDiff=maxDetectorPosError)
                self.assertPairsAlmostEqual(beamInfo0.focalFieldAngle, (0, 0))
                self.assertPairsAlmostEqual(beamInfo0.pupilFieldAngle, (0, 0))
                self.assertPairsAlmostEqual(beamInfo0.pupilPos, pupilPos, self.maxPupilPosErr)
                self.assertEqual(beamInfo0.detectorName, "D0")
                bboxd = Box2D(self.cco.cameraGeom["D0"].getBBox())
                detectorCtrPos = bboxd.getCenter()
                detector34Pos = bboxd.getMin() + bboxd.getDimensions()*0.75
                self.assertPairsAlmostEqual(beamInfo0.detectorPos, detectorCtrPos,
                                            maxDiff=maxDetectorPosError)
                self.assertTrue(beamInfo0.isOnDetector)

                # beam 1 should on detector D0, 3/4 of the way from LL to UR
                beamInfo1 = self.cco["1"]
                self.assertTrue(beamInfo1.isOnDetector)
                self.assertEqual(beamInfo1.detectorName, "D0")
                self.assertPairsAlmostEqual(beamInfo1.detectorPos, detector34Pos,
                                            maxDiff=maxDetectorPosError)

                # beam 2 should land on the center of detector D1
                beamInfo2 = self.cco["2"]
                self.assertTrue(beamInfo2.isOnDetector)
                self.assertEqual(beamInfo2.detectorName, "D1")
                self.assertPairsAlmostEqual(beamInfo2.detectorPos, detectorCtrPos,
                                            maxDiff=maxDetectorPosError)

                # beam 3 should on detector D1, 3/4 of the way from LL to UR
                beamInfo3 = self.cco["3"]
                self.assertTrue(beamInfo3.isOnDetector)
                self.assertEqual(beamInfo3.detectorName, "D1")
                self.assertPairsAlmostEqual(beamInfo3.detectorPos, detector34Pos,
                                            maxDiff=maxDetectorPosError)

                # beam 4 should land on the center of detector D2
                beamInfo4 = self.cco["4"]
                self.assertTrue(beamInfo4.isOnDetector)
                self.assertEqual(beamInfo4.detectorName, "D2")
                self.assertPairsAlmostEqual(beamInfo4.detectorPos, detectorCtrPos,
                                            maxDiff=maxDetectorPosError)

                # beam 5 should on detector D2, 3/4 of the way from LL to UR
                beamInfo5 = self.cco["5"]
                self.assertTrue(beamInfo5.isOnDetector)
                self.assertEqual(beamInfo5.detectorName, "D2")
                self.assertPairsAlmostEqual(beamInfo5.detectorPos, detector34Pos,
                                            maxDiff=maxDetectorPosError)

    def testSetPupilFieldAngle(self):
        for pupilFieldAngle, pupilPos, beam in itertools.product(
            ((0, 0), (0, 0.05), (-0.05, -0.03)),
            ((0, 0), (0, 5000), (-5000, 0), (5000, -5000)),
            self.cco.beamNames,
        ):
            with self.subTest(pupilFieldAngle=pupilFieldAngle, pupilPos=pupilPos, beam=beam):
                self.cco.setPupilFieldAngle(pupilPos=pupilPos, pupilFieldAngle=pupilFieldAngle, beam=beam)
                beamInfo = self.cco[beam]
                self.assertPairsAlmostEqual(beamInfo.pupilPos, pupilPos, self.maxPupilPosErr)
                self.assertPairsAlmostEqual(beamInfo.pupilFieldAngle, pupilFieldAngle, maxDiff=1e-5)
                self.checkOrientation()

    def testAngleProperties(self):
        """Test cbpAzAltObserved, telAzAltObserved, telRotObserved
        and their Internal equivalents
        """
        # set non-trivial scale and offset for all axes;
        # azimuth and rotator scales must be ±1 in order to handle wrap correctly;
        # altitude scales should be nearly 1 and altitudeoffsets should be small
        # to avoid hitting limits
        self.cco.config.cbpAzAltScale = (-1, 0.98)
        self.cco.config.cbpAzAltOffset = (33.2*degrees, -2.67*degrees)
        self.cco.config.telAzAltScale = (1, 0.95)
        self.cco.config.telAzAltOffset = (-31.5*degrees, 2.3*degrees)
        self.cco.config.telRotScale = -1
        self.cco.config.telRotOffset = 222.2*degrees
        for cbpAzAltObserved in (
            SpherePoint(1, 2, degrees),
            SpherePoint(-45, 75.2, degrees),
        ):
            self.cco.cbpAzAltObserved = cbpAzAltObserved
            self.assertSpherePointsAlmostEqual(self.cco.cbpAzAltObserved, cbpAzAltObserved)

            # observed angle = internal angle * scale + offset
            # internal angle = (observed angle - offset) / scale
            predictedCbpAzAltInternal = SpherePoint(
                *[(cbpAzAltObserved[i] - self.cco.config.cbpAzAltOffset[i])/self.cco.config.cbpAzAltScale[i]
                  for i in range(2)])
            self.assertSpherePointsAlmostEqual(self.cco.cbpAzAltInternal, predictedCbpAzAltInternal)

        for telAzAltObserved in (
            SpherePoint(-3, 5, degrees),
            SpherePoint(37, -.2, degrees),
        ):
            self.cco.telAzAltObserved = telAzAltObserved
            self.assertSpherePointsAlmostEqual(self.cco.telAzAltObserved, telAzAltObserved)
            predictedTelAzAltInternal = SpherePoint(
                *[(telAzAltObserved[i] - self.cco.config.telAzAltOffset[i])/self.cco.config.telAzAltScale[i]
                  for i in range(2)])
            self.assertSpherePointsAlmostEqual(self.cco.telAzAltInternal, predictedTelAzAltInternal)

        for telRotObserved in (0*degrees, -32*degrees, 167*degrees):
            self.cco.rotAzAltObserved = telRotObserved
            self.assertAnglesAlmostEqual(self.cco.rotAzAltObserved, telRotObserved)
            predictedRotInternal = telRotObserved/self.cco.config.telRotScale - self.cco.config.telRotOffset
            self.assertAnglesAlmostEqual(self.cco.telRotInternal, predictedRotInternal)

    def testInBounds(self):
        """Test the telInBounds and cbpInBounds properties
        """
        # shrink the upper altitude limits so there's somewhere to go
        self.cco.config.cbpAltitudeLimits = (-88*degrees, 88*degrees)
        self.cco.config.telAltitudeLimits = (5*degrees, 85*degrees)
        self.cco.setFocalPlanePos((0, 0))
        self.assertTrue(self.cco.cbpInBounds)
        self.assertTrue(self.cco.telInBounds)
        originalCbpAzAltObserved = self.cco.cbpAzAltObserved
        cbpAltLim = self.cco.config.cbpAltitudeLimits
        for cbpAltObserved in cbpAltLim:
            self.cco.cbpAzAltObserved = SpherePoint(originalCbpAzAltObserved[0], cbpAltObserved)
            self.assertTrue(self.cco.cbpInBounds)
            self.assertTrue(self.cco.telInBounds)
        for badCbpAltObserved in (
            cbpAltLim[0] - 1e-15*radians,
            cbpAltLim[0] - 2*degrees,
            cbpAltLim[1] + 1e-15*radians,
            cbpAltLim[1] + 2*degrees,
        ):
            self.cco.cbpAzAltObserved = SpherePoint(originalCbpAzAltObserved[0], badCbpAltObserved)
            self.assertFalse(self.cco.cbpInBounds)
            self.assertTrue(self.cco.telInBounds)
        self.cco.cbpAzAltObserved = originalCbpAzAltObserved

        originalTelAzAltObserved = self.cco.telAzAltObserved
        telAltLim = self.cco.config.telAltitudeLimits
        for telAltObserved in telAltLim:
            self.cco.telAzAltObserved = SpherePoint(originalTelAzAltObserved[0], telAltObserved)
            self.assertTrue(self.cco.telInBounds)
            self.assertTrue(self.cco.cbpInBounds)
        for badTelAltObserved in (
            telAltLim[0] - 1e-15*radians,
            telAltLim[0] - 2*degrees,
            telAltLim[1] + 1e-15*radians,
            telAltLim[1] + 2*degrees,
        ):
            self.cco.telAzAltObserved = SpherePoint(originalTelAzAltObserved[0], badTelAltObserved)
            self.assertFalse(self.cco.telInBounds)
            self.assertTrue(self.cco.cbpInBounds)

    def testHolePositionsFlipX(self):
        for telFlipX, cbpFlipX in itertools.product((False, True), (False, True)):
            with self.subTest(telFlipX=telFlipX, cbpFlipX=cbpFlipX):
                scc = SampleCoordinateConverter(telFlipX=telFlipX, cbpFlipX=cbpFlipX)
                holePositions = [scc.maskInfo.getHolePos(name) for name in scc.maskInfo.holeNames]
                if not (telFlipX or cbpFlipX):
                    unflippedHolePositions = holePositions[:]
                    flippedHolePositions = [(-pos[0], pos[1]) for pos in unflippedHolePositions]
                elif cbpFlipX:
                    self.assertPairListsAlmostEqual(holePositions, flippedHolePositions)
                else:
                    self.assertPairListsAlmostEqual(holePositions, unflippedHolePositions)

    def testSetPupilFieldAngleFlipX(self):
        """Test setPupilFieldAngle with varying flipX for telescope and CBP
        """
        beam = 3  # pick a beam that is a bit off center
        for telFlipX, cbpFlipX in itertools.product((False, True), (False, True)):
            with self.subTest(telFlipX=telFlipX, cbpFlipX=cbpFlipX):
                # flip pupilPos and fieldAngle so that everything else is identical
                # e.g. all 3-D vectors; that makes debugging easier if the test fails
                unflippedPupilPos = (5000, -5000)
                unflippedPupilFieldAngle = (0.1, 0.1)
                pupilPos = coordUtils.getFlippedPos(unflippedPupilPos, flipX=telFlipX)
                pupilFieldAngle = coordUtils.getFlippedPos(unflippedPupilFieldAngle, flipX=telFlipX)

                scc = SampleCoordinateConverter(telFlipX=telFlipX, cbpFlipX=cbpFlipX)
                cco = scc.coordinateConverter
                cco.setPupilFieldAngle(pupilPos=pupilPos, pupilFieldAngle=pupilFieldAngle, beam=beam)
                beamInfo = cco[beam]
                self.assertPairsAlmostEqual(beamInfo.pupilFieldAngle, pupilFieldAngle)
                self.assertPairsAlmostEqual(beamInfo.pupilPos, pupilPos, self.maxPupilPosErr)

    def testSetDetectorPosFlipX(self):
        """Test setDetectorPos with varying flipX for telescope and CBP
        """
        pupilPos = (5000, -5000)
        detectorName = "D2"  # pick a detector well away from the center of the focal plane
        detectorPos = (750, 250)
        beam = 3  # pick a beam that is a bit off center
        for telFlipX, cbpFlipX in itertools.product((False, True), (False, True)):
            with self.subTest(telFlipX=telFlipX, cbpFlipX=cbpFlipX):
                scc = SampleCoordinateConverter(telFlipX=telFlipX, cbpFlipX=cbpFlipX)
                cco = scc.coordinateConverter
                cco.setDetectorPos(pupilPos=pupilPos, detectorPos=detectorPos, detector=detectorName,
                                   beam=beam)
                beamInfo = cco[beam]
                self.assertTrue(beamInfo.isOnDetector)
                self.assertEqual(beamInfo.detectorName, detectorName)
                self.checkOrientation()
                self.assertPairsAlmostEqual(beamInfo.detectorPos, detectorPos)
                self.assertPairsAlmostEqual(beamInfo.pupilPos, pupilPos, self.maxPupilPosErr)

    def checkOrientation(self, maxDiff=50*arcseconds):
        """Check that the orientation of the focal plane is correct

        The definition of correct orientation is that two points to the
        left and right of the CBP center (by a buried delta)
        line up in the focal plane. We'll just use +/-1 for our delta
        """
        ctrHolePos = Point2D(0, 0)
        dx = -1 if self.cco.config.cbpFlipX else 1
        holePos1 = ctrHolePos[0] + dx, ctrHolePos[1]
        holePos2 = ctrHolePos[0] + dx, ctrHolePos[1]
        beamInfo1 = self.cco.getBeamInfo(beam="virtual1", holePos=holePos1)
        beamInfo2 = self.cco.getBeamInfo(beam="virtual2", holePos=holePos2)
        deltaFocalPlane = np.subtract(beamInfo2.focalPlanePos, beamInfo1.focalPlanePos)
        orientError = -math.atan2(deltaFocalPlane[1], deltaFocalPlane[0])*radians
        self.assertAnglesAlmostEqual(orientError, 0*radians)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
