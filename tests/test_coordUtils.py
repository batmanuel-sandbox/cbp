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
from lsst.afw.geom import SpherePoint, degrees, radians
from lsst.cbp import coordUtils
import lsst.utils.tests

RadPerDeg = math.pi / 180

# set True to record and report numeric error in convertVectorFromPupilToBase
ReportRecordedErrors = False


class CoordUtilsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        if ReportRecordedErrors:
            coordUtils.startRecordingErrors()

    def tearDown(self):
        if ReportRecordedErrors:
            errList = coordUtils.getRecordedErrors()
            coordUtils.stopRecordingErrors()
            if len(errList) > 0:
                print("\nWarning: recorded %s significant errors in computeAzAltFromBasePupil;"
                      " the worst 5 are:" % len(errList))
                for elt in errList[-5:]:
                    print("error=%0.5f arcsec, vectorBase=%s, vectorPupil=%s" % elt)

    def testGetFlippedPos(self):
        floatList = (0, -5.1, 4.3)
        for x, y, flipX in itertools.product(floatList, floatList, (False, True)):
            with self.subTest(x=x, y=y, flipX=flipX):
                if flipX:
                    desiredResult = (-x, y)
                else:
                    desiredResult = (x, y)
                result = coordUtils.getFlippedPos(xyPos=(x, y), flipX=flipX)
                self.assertEqual(result, desiredResult)

    def testFieldAngleToVector(self):
        sp00 = SpherePoint(0, 0, degrees)
        degList = (-90, -89.9, -20, 0, 10, 89.9, 90)
        for xdeg, ydeg, flipX in itertools.product(degList, degList, (False, True)):
            with self.subTest(xdeg=xdeg, ydeg=ydeg, flipX=flipX):
                xrad = xdeg * RadPerDeg
                signx = -1 if flipX else 1
                testOrientation = xdeg != 0 or ydeg != 0
                yrad = ydeg * RadPerDeg
                fieldAngle = (xrad, yrad)
                vector = coordUtils.fieldAngleToVector(fieldAngle, flipX)
                self.assertAlmostEqual(np.linalg.norm(vector), 1)
                if testOrientation:
                    # orientation should match
                    orientationFromFieldAngle = math.atan2(yrad, signx*xrad)*radians
                    # field angle x = vector y, field angle y = vector z
                    orientationFromVector = math.atan2(vector[2], vector[1])*radians
                    self.assertAnglesAlmostEqual(orientationFromVector, orientationFromFieldAngle)

                # now test as spherical geometry
                sp = SpherePoint(Vector3d(*vector))
                separation = sp00.separation(sp)
                predictedSeparation = math.hypot(xrad, yrad)*radians
                self.assertAnglesAlmostEqual(predictedSeparation, separation)
                if testOrientation:
                    bearing = sp00.bearingTo(sp)
                    self.assertAnglesAlmostEqual(orientationFromFieldAngle, bearing)

                # test round trip through vectorToFieldAngle
                fieldAngleFromVector = coordUtils.vectorToFieldAngle(vector, flipX)
                np.testing.assert_allclose(fieldAngleFromVector, fieldAngle, atol=1e-15)

    def testVectorToFieldAngle(self):
        # note: more sophisticated cases are tested by testFieldAngleToVector
        for flipX in (False, True):
            signx = -1 if flipX else 1
            for magMultiplier in (0.001, 1, 1000):
                for vector, predictedFieldAngleDeg in (
                    ((1, 0, 0), (0, 0)),
                    ((0, 1, 0), (signx*90, 0)),
                    ((0, 0, 1), (0, 90)),
                ):
                    predictedFieldAngle = [val*RadPerDeg for val in predictedFieldAngleDeg]
                    fieldAngle = coordUtils.vectorToFieldAngle(vector, flipX)
                    np.testing.assert_allclose(predictedFieldAngle, fieldAngle, atol=1e-15)

    def testComputeShiftedPlanePosZeroFieldAngle(self):
        """Test computeShiftedPlanePos with zero field angle

        This should result in no change
        """
        zeroFieldAngle = (0, 0)
        for planePos in (
            (-1000, -2000),
            (0, 0),
            (5000, 4000),
        ):
            for shift in (-500, 0, 500):
                shiftedPlanePos = coordUtils.computeShiftedPlanePos(planePos, zeroFieldAngle, shift)
                self.assertPairsAlmostEqual(planePos, shiftedPlanePos)

    def testComputeShiftedPlanePosZeroShift(self):
        """Test computeShiftedPlanePos with zero shift

        This should result in no change
        """
        zeroShift = 0
        for planePos in (
            (-1000, -2000),
            (0, 0),
            (5000, 4000),
        ):
            for fieldAngle in (
                (0.5, 0.5),
                (0, 1),
                (-0.5, 0.3),
            ):
                shiftedPlanePos = coordUtils.computeShiftedPlanePos(planePos, fieldAngle, zeroShift)
                self.assertPairsAlmostEqual(planePos, shiftedPlanePos)

    def testComputeShiftedPlanePos(self):
        """Test computeShiftedPlanePos for the general case
        """
        # in the general case the increase in x and y equals the shift
        # times y/x, z/x of a vector equivalent to the field angle
        for ratios in (
            (0.0, 0.0),
            (0.5, 0.5),
            (-0.23, 0.75),
            (0.3, 0.1),
        ):
            vector = (1, ratios[0], ratios[1])
            fieldAngle = coordUtils.vectorToFieldAngle(vector, False)
            for planePos in (
                (-1000, -2000),
                (0, 0),
                (5000, 4000),
            ):
                for shift3 in (-550, 0, 375):
                    predictedShiftedPlanePos = [planePos[i] + shift3*ratios[i] for i in range(2)]
                    shiftedPlanePos = coordUtils.computeShiftedPlanePos(planePos, fieldAngle, shift3)
                    self.assertPairsAlmostEqual(shiftedPlanePos, predictedShiftedPlanePos)

    def testConvertVectorFromPupilToBase(self):
        """Test convertVectorFromPupilToBase and convertVectorFromBaseToPupil
        """
        cos30 = math.cos(30 * math.pi / 180)
        sin30 = math.sin(30 * math.pi / 180)
        for magMultiplier in (0.001, 1, 1000):
            for azAltDeg, vectorPupil, predictedVectorBase in (
                # at az=0, alt=0: base = pupil
                ((0, 0), (1, 0, 0), (1, 0, 0)),
                ((0, 0), (0, 1, 0), (0, 1, 0)),
                ((0, 0), (0, 0, -1), (0, 0, -1)),
                ((0, 0), (1, -1, 1), (1, -1, 1)),
                # at az=90, alt=0: base x = pupil -y, base y = pupil x, base z = pupil z
                ((90, 0), (1, 0, 0), (0, 1, 0)),
                ((90, 0), (0, 1, 0), (-1, 0, 0)),
                ((90, 0), (0, 0, -1), (0, 0, -1)),
                ((90, 0), (1, -1, 1), (1, 1, 1)),
                # at az=0, alt=90: base x = - pupil z, base y = pupil y, base z = pupil x
                ((0, 90), (1, 0, 0), (0, 0, 1)),
                ((0, 90), (0, 1, 0), (0, 1, 0)),
                ((0, 90), (0, 0, -1), (1, 0, 0)),
                ((0, 90), (1, -1, 1), (-1, -1, 1)),
                # at az=90, alt=90: base x = -pupil y, base y = - pupil z, base z = pupil x
                ((90, 90), (1, 0, 0), (0, 0, 1)),
                ((90, 90), (0, 1, 0), (-1, 0, 0)),
                ((90, 90), (0, 0, -1), (0, 1, 0)),
                ((90, 90), (1, -1, 1), (1, -1, 1)),
                # at az=0, alt=45:
                # base x = cos(30) * pupil x - sin(30) * pupil z
                # base y = pupil y
                # base z = sin(30) * pupil x + cos(30) * pupil z
                ((0, 30), (1, 0, 0), (cos30, 0, sin30)),
                ((0, 30), (0, 1, 0), (0, 1, 0)),
                ((0, 30), (0, 0, -1), (sin30, 0, -cos30)),
                ((0, 30), (1, -1, 1), (cos30 - sin30, -1, sin30 + cos30)),
                # at az=30, alt=0:
                # base x = cos(30) * pupil x - sin(30) * pupil y
                # base y = sin(30) * pupil x + cos(30) * pupil y
                # base z = pupil z
                ((30, 0), (1, 0, 0), (cos30, sin30, 0)),
                ((30, 0), (0, 1, 0), (-sin30, cos30, 0)),
                ((30, 0), (0, 0, -1), (0, 0, -1)),
                ((30, 0), (1, -1, 1), (cos30 + sin30, sin30 - cos30, 1)),
            ):
                vectorPupil = np.array(vectorPupil) * magMultiplier
                predictedVectorBase = np.array(predictedVectorBase) * magMultiplier
                pupilAzAlt = SpherePoint(*azAltDeg, degrees)
                vectorBase = coordUtils.convertVectorFromPupilToBase(vectorPupil=vectorPupil,
                                                                     pupilAzAlt=pupilAzAlt)
                atol = max(magMultiplier, 1) * 1e-15
                msg = "azAltDeg=%s, vectorPupil=%s" % (azAltDeg, vectorPupil)
                np.testing.assert_allclose(vectorBase, predictedVectorBase, atol=atol,
                                           err_msg=msg, verbose=True)

                vectorPupilRoundTrip = coordUtils.convertVectorFromBaseToPupil(vectorBase=vectorBase,
                                                                               pupilAzAlt=pupilAzAlt)
                np.testing.assert_allclose(vectorPupil, vectorPupilRoundTrip, atol=atol,
                                           err_msg=msg, verbose=True)

    def testComputeAzAltFromPupilBaseWithBaseEqualsPupil(self):
        """Test computeAzAltFromBasePupil with baseVector=pupilVector,
        so the telescope will to internal az, alt=0
        """
        zeroSp = SpherePoint(0, 0, radians)
        for vector, pupilMagFactor, baseMagFactor in itertools.product(
            ((1, 0, 0), (0, -1, 0), (0, -0.5, 0.5), (0.5, 0, 0.5), (1, 0.7, -0.8)),
            (1, 1000),
            (1, 1000),
        ):
            with self.subTest(vector=vector, pupilMagFactor=pupilMagFactor, baseMagFactor=baseMagFactor):
                vectorPupil = np.array(vector, dtype=float) * pupilMagFactor
                vectorBase = np.array(vector, dtype=float) * baseMagFactor
                obs = coordUtils.computeAzAltFromBasePupil(vectorPupil=vectorPupil,
                                                           vectorBase=vectorBase)
                sep = zeroSp.separation(obs).asRadians()
                self.assertLess(sep, 1e-14)

    def testComputeAzAltFromPupilBaseWithVectorPupil100(self):
        """Test computeAzAltFromBasePupil with vectorPupil = (1, 0, 0),
        so internal az/alt points along vectorBase
        """
        vectorPupil = (1, 0, 0)
        for vectorBase, pupilMagFactor, baseMagFactor in itertools.product(
            ((1, 0, 0), (0, -1, 0), (0, -0.5, 0.5), (0.5, 0, 0.5), (1, 0.7, -0.8)),
            (1, 1000),
            (1, 1000),
        ):
            with self.subTest(vectorBase=vectorBase, pupilMagFactor=pupilMagFactor,
                              baseMagFactor=baseMagFactor):
                predictedPupilAzalt = SpherePoint(Vector3d(*vectorBase))
                vectorPupilScaled = np.array(vectorPupil, dtype=float) * pupilMagFactor
                vectorBaseScaled = np.array(vectorBase, dtype=float) * baseMagFactor
                pupilAzAlt = coordUtils.computeAzAltFromBasePupil(vectorPupil=vectorPupilScaled,
                                                                  vectorBase=vectorBaseScaled)
                sep = pupilAzAlt.separation(predictedPupilAzalt)
                if sep.asRadians() > 1e-14:
                    print("Warning: sep=%0.5f asec for vectorPupilScaled=%s, vectorBaseScaled=%s" %
                          (sep.asArcseconds(), vectorPupilScaled, vectorBaseScaled))
                # the worst error I see is 0.0026" for vectorBase=(1, 0.7, -0.8)
                # that's bad enough to scary, but is acceptable
                self.assertLess(sep.asArcseconds(), 0.01)

    def testComputeAzAltFromPupilBase(self):
        """Test computeAzAltFromBasePupil with general values
        """
        # transform the pupil vector back to the base vector
        # using the computed internal az/alt position
        for vectorPupil, vectorBase, pupilMagFactor, baseMagFactor in itertools.product(
            ((1, 0, 0), (2, 1, 0), (2, 0, 1), (2, 0.7, -0.8)),
            ((1, 0, 0), (0, 1, 0), (1, -0.7, 0.8)),
            (1, 1000),
            (1, 1000),
        ):
            with self.subTest(vectorPupil=vectorPupil, vectorBase=vectorBase, pupilMagFactor=pupilMagFactor,
                              baseMagFactor=baseMagFactor):
                vectorPupilScaled = np.array(vectorPupil, dtype=float) * pupilMagFactor
                pupilMag = np.linalg.norm(vectorPupilScaled)
                vectorBaseScaled = np.array(vectorBase, dtype=float) * baseMagFactor
                pupilAzAlt = coordUtils.computeAzAltFromBasePupil(vectorPupil=vectorPupilScaled,
                                                                  vectorBase=vectorBaseScaled)
                # check the round trip; note that the magnitude of the returned vector
                # will equal the magnitude of the input vector
                vectorBaseRoundTrip = coordUtils.convertVectorFromPupilToBase(
                    vectorPupil=vectorPupilScaled,
                    pupilAzAlt=pupilAzAlt)
                vectorBaseRoundTripMag = np.linalg.norm(vectorBaseRoundTrip)
                self.assertAlmostEqual(vectorBaseRoundTripMag, pupilMag, delta=1e-15*pupilMag)
                spBase = SpherePoint(Vector3d(*vectorBase))
                spBaseRoundTrip = SpherePoint(Vector3d(*vectorBaseRoundTrip))
                sep = spBase.separation(spBaseRoundTrip)
                self.assertLess(sep.asRadians(), 2e-15)

    def testRotate2d(self):
        for pos, angleDeg, expectedPos in (
            ((1, 2), 0, (1, 2)),
            ((1, 0), -30, (math.cos(RadPerDeg*30), -0.5)),
            ((0, 1), -30, (0.5, math.cos(RadPerDeg*30))),
            ((1, 2), 90, (-2, 1)),
        ):
            with self.subTest(pos=pos, angleDeg=angleDeg, expectedPos=expectedPos):
                angle = angleDeg*degrees
                rotatedPos = coordUtils.rotate2d(pos, angle)
                self.assertPairsAlmostEqual(rotatedPos, expectedPos)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
