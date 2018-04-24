.. _lsst.cbp.configuration:

################################################
How to Configure an lsst.cbp.CoordinateConverter
################################################

You will need to provide three items of configuration to construct an lsst.cbp.CoordinateConverter:

- `basic configuration`_ as an `lsst.cbp.CoordinateConverterConfig`
- `mask information`_ as an `lsst.cbp.MaskInfo`
- `camera geometry`_ as an `lsst.afw.cameraGeom.Camera`

Basic Configuration
-------------------

Basic configuration, an instance of `lsst.cbp.CoordinateConverterConfig`, includes the position of the collimated beam projector relative to the telescope and information about the axes of the telescope and CBP. Fields that require careful thought are described here.

The first step is to pick a :ref:`base coordinate frame <lsst.cbp.base_frame>`, which is tied to the observatory. The z axis must be vertical and the x and y axes can point anywhere you like, as long as the frame is right-handed. One good choice is: x points south and y points east.

Determine the position of the CBP with respect to the telescope in the :ref:`base frame <lsst.cbp.base_frame>`: ``cbpPosition``. This is measured from the center of the telescope to the center of the CBP, where "center" is that point where the azimuth and altitude axes of the telescope or CBP intersect.

Determine telescope azimuth offset and scale: ``telAzimuthOffset``. ``telAzimuthScale``. Mentally point the telescope along the base frame x axis. The :ref:`observed azimuth <lsst.cbp.observed_angles>` is `telAzimuthOffset``. Now mentally slew the telescope an addition +90° in :ref:`observed azimuth <lsst.cbp.observed_angles>`. If the telescope is pointing along the base frame y axis then ``telAzimuthScale`` is +1; if along base frame -y then ``telAzimuthScale`` is -`.
For example if you choose your base x axis pointing south (and thus base y pointing east) and your control system uses azimuth = 0° north and 90° east (a common convention), then ``telAzimuthOffset`` = 180° and ``telAzimuthScale`` = -1.

Determine CBP azimuth and offset scale, ``cbpAzimuthOffset``, and ``cbpAzimuthScale``, in the same way.

Determine the parity of the telescope plane: ``telFlipX``. Mentally point the telescope along x (:ref:`internal azimuth, altitude <lsst.cbp.internal_angles>` = (0°, 0°)) and rotate the camera such that a point at positive :ref:`focal plane y <lsst.cbp.focal_plane>` is emitted in the :ref:`base frame x, z plane <lsst.cbp.base_frame>`, with the ray angled towards base +z. Now a point along :ref:`focal plane x <lsst.cbp.focal_plane>` has a principle ray in the :ref:`base frame x, y plane <lsst.cbp.base_frame>`. If that ray is tilted towards base +y then focal plane is not flipped. If the ray is tilted towards base -y then the focal plane is flipped.

Determine the parity of the CBP focal plane ``cbpFlipX`` in the same way. For the CBP, focal plane x, y is mask hole position x, y.

Determine the telescope camera rotator offset and scale: ``telRotatorScale``, ``telRotatorOffset``. Mentally rotate the camera until a point along focal plane +y is emitted in the :ref:`base frame x, z plane <lsst.cbp.base_frame>`, tilting towards the base z axis. The current observed rotator angle is ``telRotatorOffset``. Now rotate the camera an additional +90°. The principle ray from a point along focal plane +y should now be in the :ref:`base frame x, y plane <lsst.cbp.base_frame>`. If the principle ray points along base -x then ``telRotatorScale`` is -1; if it tilts towards base +x then ``telRotatorScale`` is +1.

Mask Information
----------------

Information about holes in the CBP mask, as an `lsst.cbp.MaskInfo`. You must provide the position of each hole on the CBP mask and the index or name of a default hole (typically the hole closest to the center of the mask). In the case of large or irregularly shaped holes you must pick a suitable point in the hole for the position.

Camera Geometry
---------------

Camera geometry is an instances of `lsst.afw.cameraGeom.Camera`. The contents of interest to the CBP software includes a transformation between `focal plane position <lsst.cbp.focal_plane>` and `focal plane field angle <lsst.cbp.focal_plane_field_angle>`, a mapping from detector position (in pixels) to focal plane position for each detector, and the dimensions of each detector.

Defining camera geometry is out of scope for this document. However, every "obs" package defines the appropriate camera geometry. To obtain it you will need to open a data butler on a data repository; you can then get the camera geometry using `butler.get("camera")`.
