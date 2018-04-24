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

import unittest

from lsst.afw.geom import degrees
from lsst.cbp import computeHolePositions, MaskInfo
import lsst.utils.tests
from lsst.cbp.testUtils import SampleCoordinateConverter


class ComputeHolePositionsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.scc = SampleCoordinateConverter()
        self.cco = self.scc.coordinateConverter
        # point the telescope and CCB right at each other
        self.cco.setFocalPlanePos(pupilPos=(0, 0))

    def tearDown(self):
        del self.scc
        del self.cco

    def testSetup(self):
        self.assertAnglesAlmostEqual(self.cco.telAzAltInternal.separation(self.cco.cbpAzAltInternal),
                                     180*degrees)
        self.assertAnglesAlmostEqual(self.cco.telRotInternal, 0*degrees)

    def testOneDetector(self):
        detectorPositions = ((100, 250), (250, 100), (500, 1500))
        for detector in self.cco.cameraGeom:
            detectorName = detector.getName()
            holePositions = computeHolePositions(
                detectorNames=[detectorName],
                detectorPositions=detectorPositions,
                cameraGeom=self.cco.cameraGeom,
                cbpFlipX=self.cco.config.cbpFlipX,
                cbpFocalLength=self.cco.config.cbpFocalLength,
            )
            self.assertEqual(len(holePositions), len(detectorPositions))
            self.cco.maskInfo = MaskInfo(
                name="OneDetector",
                holePositions=holePositions,
                defaultHole=0,
            )
            self.assertEqual(len(self.cco), len(detectorPositions))
            for i, beamInfo in enumerate(self.cco):
                self.assertTrue(beamInfo.isOnDetector)
                self.assertEqual(beamInfo.detectorName, detectorName)
                self.assertPairsAlmostEqual(beamInfo.detectorPos, detectorPositions[i])

    def testOneHoleAllDetectors(self):
        detectorPosition = (449.2, 732.1)
        detectorNames = set(self.cco.cameraGeom.getNameIter())
        numDetectors = len(self.cco.cameraGeom)
        for cbpFlipX in (False, True):
            self.cco.config.cbpFlipX = cbpFlipX
            holePositions = computeHolePositions(
                detectorNames=None,
                detectorPositions=[detectorPosition],
                cameraGeom=self.cco.cameraGeom,
                cbpFlipX=self.cco.config.cbpFlipX,
                cbpFocalLength=self.cco.config.cbpFocalLength,
            )
            if not cbpFlipX:
                flippedHolePositions = [(-val[0], val[1]) for val in holePositions]
            else:
                self.assertPairListsAlmostEqual(holePositions, flippedHolePositions)
            self.assertEqual(len(holePositions), numDetectors)
            self.cco.maskInfo = MaskInfo(
                name="AllDetectors",
                holePositions=holePositions,
                defaultHole=0,
            )
            self.assertEqual(len(self.cco), numDetectors)
            for beamInfo in self.cco:
                self.assertTrue(beamInfo.isOnDetector)
                self.assertIn(beamInfo.detectorName, detectorNames)
                self.assertPairsAlmostEqual(beamInfo.detectorPos, detectorPosition)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
