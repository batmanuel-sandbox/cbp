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

__all__ = ["SampleCoordinateConverter"]

import numpy as np

import lsst.afw.cameraGeom as cameraGeom
from lsst.afw.cameraGeom.cameraFactory import makeDetector
import lsst.afw.geom as afwGeom
from lsst.afw.table import AmpInfoCatalog, AmpInfoTable, LL
from .coordinateConverterConfig import CoordinateConverterConfig
from .coordinateConverter import CoordinateConverter
from .computeHolePositions import computeHolePositions
from .maskInfo import MaskInfo


class SampleCoordinateConverter:
    """An object containing a CoordinateConverter and the information
    used to create it.

    Parameters
    ----------
    detectorFracPosList : `iterable` of pair of `float` (optional)
        Position of the center of each detector, as a fraction of
        the width and height of the detector.
        The first element must have value (0, 0).
        See the field of the same name for more information.
        Defaults to::

            (
                (0, 0),
                (1.01, 0),  # 1.01: leave a 1% gap
                (-4, 7),    # a corner detector in the LSST camera
            )

    holeFracPosList : `iterable` of pair of `float` (optional)
        Positions of holes on a given detector,
        as a fraction of the distance from lower left corner
        to upper right corner. Thus (0.5, 0.5) is centered
        on the detector.
        Defaults to `((0, 0), (0.75, 0.75))`

    Notes
    -----
    **Fields:**

    detectorWidthPix : `int`
        Width of each detector, in pixels
    detectorHeightPix : `int`
        Height of each detector, in pixels
    pixelSizeMm : `float`
        Width = height of each pixel, in mm
    plateScale : `lsst.afw.geom.Angle`
        Plate scale: in angle on the sky per mm on the focal plane
    detectorFracPosList : `iterable` of pair of `float`
        Position of the center of each detector, as a fraction of
        the width and height of the detector. For instance
        (0, 0) is a detector centered on the focal plane
        and (1, 0) is adjacent to a centered detector,
        in the direction of increasing focal plane x.
    holeFracPosList : `iterable` of pair of `float`
        Positions of holes on a given detector,
        as a fraction of the distance from lower left corner
        to upper right corner. Thus (0.5, 0.5) is centered
        on the detector.
    cameraGeom : `lsst.afw.cameraGeom.Camera`
        Camera geometry. There will be one detector per entry in
        detectorFracPosList with names "D0", "D1", ...
        Detector "D0" is centered on the focal plane.
    config : `lsst.cbp.CoordinateConverterConfig`
        Basic configuration for ``coordinateConverter``
    maskInfo : `lsst.cbp.MaskInfo`
        CBP mask information
    coordinateConverter : `lsst.cbp.CoordinateConverter`
        The test coordinate converter
    """

    def __init__(self, detectorFracPosList=None, holeFracPosList=None,
                 telFlipX=False, cbpFlipX=False):
        # these value are close to LSST but the detectors are
        # slightly rectangular as that can expose some errors
        self.detectorWidthPix = 4000
        self.detectorHeightPix = 3998
        self.pixelSizeMm = 0.01
        self.plateScale = 20 * afwGeom.arcseconds
        if holeFracPosList is None:
            holeFracPosList = ((0.5, 0.5), (0.75, 0.75))
        self.holeFracPosList = holeFracPosList
        if detectorFracPosList is None:
            detectorFracPosList = (
                (0, 0),
                (1.01, 0),  # 1.01: leave a 1% gap
                (-4, 7),    # a corner detector in the LSST camera
            )
        self.detectorFracPosList = detectorFracPosList

        self.cameraGeom = self.makeCameraGeom()
        self.config = self.makeCoordinateConverterConfig(
            telFlipX=telFlipX,
            cbpFlipX=cbpFlipX,
        )
        self.maskInfo = self.makeMaskInfo()
        self.coordinateConverter = CoordinateConverter(
            config=self.config,
            maskInfo=self.maskInfo,
            cameraGeom=self.cameraGeom,
        )

    def makeCoordinateConverterConfig(self, telFlipX, cbpFlipX):
        """Make a CoordinateConverterConfig"""
        return CoordinateConverterConfig(
            telPupilOffset=101,
            telPupilDiameter=8500,
            telPupilObscurationDiameter=1000,
            telFocalPlaneDiameter=3000,
            telFlipX=telFlipX,
            telAzimuthOffset=-180,  # offset=-180, scale=-1 for az=0 N, 90 E
            telAzimuthScale=-1,
            telAltitudeOffset=0,
            telAltitudeScale=1,
            telAltitudeLimits=(0, 89),
            telRotOffset=0,
            telRotScale=-1,
            defaultDetector="D0",
            cbpPosition=(10000, 3000, 5000),
            cbpFocalLength=635,  # nominal value for LSST CBP
            cbpFlipX=cbpFlipX,
            cbpAzimuthOffset=-180,
            cbpAzimuthScale=-1,
            cbpAltitudeOffset=0,
            cbpAltitudeScale=1,
            cbpAltitudeLimits=(-70, 70),
        )

    def makeMaskInfo(self):
        """Make a MaskInfo

        The mask will have one hole per entry in self.holeFracPosList
        per detector.

        self.cameraGeom and self.config must be set before calling this
        method.
        """
        detectorBBox = afwGeom.Box2I(afwGeom.Point2I(0, 0),
                                     afwGeom.Extent2I(self.detectorWidthPix, self.detectorHeightPix))
        bboxd = afwGeom.Box2D(detectorBBox)
        bboxdMin = np.array(bboxd.getMin())
        bboxdDim = np.array(bboxd.getDimensions())
        detectorPositions = [bboxdMin + bboxdDim*np.array(fracPos) for fracPos in self.holeFracPosList]
        holePositions = computeHolePositions(
            detectorNames=None,
            detectorPositions=detectorPositions,
            cameraGeom=self.cameraGeom,
            cbpFlipX=self.config.cbpFlipX,
            cbpFocalLength=self.config.cbpFocalLength,
        )
        return MaskInfo(
            name="test",
            holePositions=holePositions,
            defaultHole=0,
        )

    def makeCameraGeom(self):
        """Make camera geometry

        There is one field per entry in self.detectorFracPosList
        with specifications set byself.detectorWidthPix,
        self.detectorHeightPix, and self.pixelSizeMm.

        The plate scale is set by self.plateScale
        and the amount of optical distortion is fixed
        """
        radialCoeff = np.array([0.0, 1.0, 0.0, 0.925]) / self.plateScale.asRadians()
        fieldAngleToFocalPlane = afwGeom.makeRadialTransform(radialCoeff)
        focalPlaneToFieldAngle = fieldAngleToFocalPlane.getInverse()
        cameraTransformMap = cameraGeom.TransformMap(cameraGeom.FOCAL_PLANE,
                                                     {cameraGeom.FIELD_ANGLE: focalPlaneToFieldAngle})
        detectorList = self._makeDetectorList(focalPlaneToFieldAngle)
        return cameraGeom.Camera("test", detectorList, cameraTransformMap)

    def _makeDetectorList(self, focalPlaneToFieldAngle):
        """Make a list of detectors

        Parameters
        ----------
        focalPlaneToFieldAngle : `lsst.afw.geom.TransformPoint2ToPoint2`
            A transform from FOCAL_PLANE to FIELD_ANGLE coordinates
            in the forward direction

        Returns
        -------
        A list of detectors, each an `lsst.afw.cameraGeom.Detector`
        """
        detectorList = []
        for i, fpPos in enumerate(self.detectorFracPosList):
            detectorConfig = self._makeDetectorConfig(id=i, fpPos=fpPos)
            ampInfoCatalog = self._makeAmpInfoCatalog()
            detector = makeDetector(detectorConfig, ampInfoCatalog, focalPlaneToFieldAngle)
            detectorList.append(detector)
        return detectorList

    def _makeDetectorConfig(self, id, fpPos):
        """Make a DetectorConfig for a detector

        Parameters
        ----------
        fpPos : pair of `float`
            Focal plane position of detector, in units of detector
            width/height. For example:

            - (0, 0) is a detector centered on the focal plane
            - (1, 0) is adjacent to a centered detector,
              in the direction of increasing focal plane x
        """
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0),
                             afwGeom.Extent2I(self.detectorWidthPix, self.detectorHeightPix))
        ctr = afwGeom.Box2D(bbox).getCenter()
        pixelSizeMm = 0.01
        config = cameraGeom.DetectorConfig()
        config.name = "D{}".format(id)
        config.id = id
        config.serial = '0000011'
        config.detectorType = 0
        config.bbox_x0 = bbox.getMinX()
        config.bbox_x1 = bbox.getMaxX()
        config.bbox_y0 = bbox.getMinY()
        config.bbox_y1 = bbox.getMaxY()
        config.pixelSize_x = pixelSizeMm
        config.pixelSize_y = pixelSizeMm
        config.transformDict.nativeSys = 'Pixels'
        config.transformDict.transforms = None
        config.refpos_x = ctr[0]
        config.refpos_y = ctr[1]
        config.offset_x = fpPos[0] * pixelSizeMm * self.detectorWidthPix
        config.offset_y = fpPos[1] * pixelSizeMm * self.detectorHeightPix
        config.transposeDetector = False
        config.pitchDeg = 0.0
        config.yawDeg = 0.0
        config.rollDeg = 0.0
        return config

    def _makeAmpInfoCatalog(self):
        """Construct an amplifier info catalog
        """
        xExtent = 1000
        yExtent = 2000
        readNoise = 3.975  # amplifier read noise, in e-
        saturationLevel = 65535
        linearityType = cameraGeom.NullLinearityType
        linearityCoeffs = [0, 0, 0, 0]

        schema = AmpInfoTable.makeMinimalSchema()

        self.ampInfoDict = {}
        ampCatalog = AmpInfoCatalog(schema)
        # amplifier gain (e-/ADU) and read noiuse (ADU/pixel) from lsstSim raw data
        # note that obs_test amp <ampX><ampY> = lsstSim amp C<ampY>,<ampX> (axes are swapped)
        gain = 1.8
        readNoise = 3.9
        record = ampCatalog.addNew()
        record.setName("0")
        ampBBox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(xExtent, yExtent))
        record.setBBox(ampBBox)
        readCorner = LL
        record.setRawBBox(ampBBox)
        record.setRawDataBBox(ampBBox)
        record.setRawHorizontalOverscanBBox(afwGeom.Box2I(afwGeom.Point2I(0, 0),
                                                          afwGeom.Extent2I(0, yExtent)))
        record.setRawXYOffset(afwGeom.Extent2I(0, 0))
        record.setReadoutCorner(readCorner)
        record.setGain(gain)
        record.setReadNoise(readNoise)
        record.setSaturation(saturationLevel)
        record.setSuspectLevel(float("nan"))
        record.setLinearityCoeffs([float(val) for val in linearityCoeffs])
        record.setLinearityType(linearityType)
        record.setHasRawInfo(True)
        record.setRawFlipX(False)
        record.setRawFlipY(False)
        record.setRawVerticalOverscanBBox(afwGeom.Box2I())
        record.setRawPrescanBBox(afwGeom.Box2I())
        return ampCatalog
