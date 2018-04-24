.. _lsst.cbp.coordinate_frames:

#################
Coordinate Frames
#################

In order to :ref:`configure <lsst.cbp.configuration>` and use the lsst.cbp.CoordinateConverter object effectively
it is important to understand the following coordinate frames.

.. _lsst.cbp.observed_angles:

Observed Azimuth, Altitude and Rotator
======================================
Azimuth, altitude and rotator angle of the telescope or CBP axes, using the conventions of the telescope, but in an ideal frame in which imperfections such as tilt and non-perpendicularity of the axes are ignored.
The transformation from "observed" az/alt to commands sent to the axis actuators consists of applying a pointing model, and is left to other software.

In order to accommodate different azimuth and rotator conventions, while simplifying the math, all internal computations are performed using :ref:`internal angles <lsst.cbp.internal_angles>`.
Internal values are mapped to observed values using an offset and scale for each axis,
which is specified in `lsst.cbp.CoordinateConverterConfig`.

.. _lsst.cbp.internal_angles:

Internal Azimuth, Altitude and Rotator
======================================
Internal azimuth and altitude are the pointing of the telescope or CBP `pupil frame`_, as the longitude and latitude axes of `lsst.afw.geom.SpherePoint`.
This specifies the zero point and sign conventions in terms of the axes of the `pupil frame`_, which see for more information.
In order to accommodate the azimuth and altitude convention used by the telescope and CBP axis controllers, this software also supports :ref:`observed azimuth and altitude <lsst.cbp.observed_angles>`, which differ from the internal values by a configurable offset and scale.

Internal rotation angle of the telescope's camera rotator is defined as the orientation of the `focal plane`_ x,y axes relative to the pupil y,z axes. Thus:

- At rotator angle zero the direction of increasing `focal plane`_ y is along increasing azimuth and `±x <lsst.cbp.flipped_x_axis>` is along the direction of increasing altitude.
- At rotator angle 90° the direction of increasing `focal plane`_ y is along the direction of decreasing altitude and `±x <lsst.cbp.flipped_x_axis>` is along the direction of increasing azimuth.

In order to accommodate the rotator convention used by the telescope's camera rotator controller, this software also supports :ref:`internal rotator angle <lsst.cbp.internal_angles>`, which differs from :ref:`observed rotator angle <lsst.cbp.observed_angles>` by a configurable offset and scale.

.. _lsst.cbp.base_frame:

Base Frame
==========
The base frame is a 3-dimensional coordinate frame that is fixed with respect to the observatory.
It is used to express the position of the CBP with respect to the telescope.
The z axis of the base frame must point up.
You are free to choose a convenient orientation for the x and y axes, so long as the coordinate system is right handed.
It is common to align x and y axes with cardinal points, e.g. x points south and y points east.

The base frame matches the `pupil frame`_ when the telescope or CBP is pointed to :ref:`internal azimuth and altitude <lsst.cbp.internal_angles>` = 0°.
Thus:
- x is the optical axis at azimuth=0°, altitude=0°
- y is the optical axis azimuth=90°, altitude=0°
- z is the optical axis altitude=90°
and in this way the base frame determines the offset and scale for telescope and CBP :ref:`observed azimuth and altitude <lsst.cbp.observed_angles>` with respect to :ref:`internal azimuth and altitude <lsst.cbp.internal_angles>`.

.. _lsst.cbp.pupil_frame:

Pupil Frame
===========
The pupil frame is simply the `base frame`_ rotated by the :ref:`internal azimuth and altitude <lsst.cbp.internal_angles>` of the telescope or CBP in the standard fashion: first rotate the `base frame`_ about base z by azimuth, then rotate that about the rotated y axis by -altitude.
At internal :ref:`internal azimuth and altitude <lsst.cbp.internal_angles>` = 0 the `pupil frame`_ matches the `base frame`_.

- x is the optical axis (+ points away from the origin)
- y at altitude = 0 it is the direction of increasing azimuth; pupil y is also `pupil field angle`_ ±x and `pupil position`_ ±x (-x if the x axis is `flipped <lsst.cbp.flipped_x_axis>`, +x if not)
- z is the direction of increasing altitude; pupil z is also field angle y and pupil position y

Note that the origin of the `pupil frame`_ will typically be offset along the optical axis from the optical pupil.
This is true of most telescopes but not the CBP (which, by design, has its optical pupil at the center).

.. _lsst.cbp.focal_plane:

Focal Plane
===========
The focal plane is a 2-dimensional plane approximation to the actual focal surface (which typically has some curvature).
The :ref:`internal rotation angle <lsst.cbp.internal_angles>` is the angle of the  ±x,y focal plane axes with respect to the y,z `pupil frame`_.

.. _lsst.cbp.flipped_x_axis:

If -x then the x axis of the focal plane and all other 2-dimensional plane positions (`pupil position`_, `focal plane field angle`_ and `pupil field angle`_) is said to be "flipped".
Determining this parity for the telescope and CBP is part of :ref:`configuration <lsst.cbp.configuration>`.

.. _lsst.cbp.pupil_position:

Pupil Position
==============
A 2-dimensional plane approximation to the primary mirror of the telescope.
This is used to specify the position of a beam on the telescope pupil.

The pupil position plane is the y,z plane of the `pupil frame`_:

- `pupil position`_ `±x <lsst.cbp.flipped_x_axis>` is along `pupil frame`_ y
- `pupil position`_ y is along `pupil frame`_ z

The `pupil frame`_ plane is typically offset along the optical axis from the optical pupil.
This is true of most telescopes but not the CBP (which, by design, has its optical pupil at the center).

If the `focal plane`_ x axis is flipped then the x axis of all other 2-dimensional plane coordinates are flipped, including this one.

.. _lsst.cbp.pupil_field_angle:

Pupil Field Angle
=================
The angle of incidence of a ray on the pupil, expressed in x,y radians.
The two components of the field angle define a great circle arc:
- arc length = hypot(x, y)
- bearing = atan2(y, `±x <lsst.cbp.flipped_x_axis>`) with 0 along `pupil frame`_ y and 90° along `pupil frame`_ z
The incident ray is the pupil x axis offset by this great circle arc.

.. _lsst.cbp.focal_plane_field_angle:

Focal Plane Field Angle
=======================
`Pupil field angle`_ with the components expressed in `focal plane`_ x,y instead of pupil x,y.
Thus this is a rotation of `pupil field angle`.
Thus y is always the pupil z axis and `±x <lsst.cbp.flipped_x_axis>` is always the pupil z.

Note that the camera geometry includes a transform from `focal plane`_ position to `focal plane field angle`_ (typically a 3rd order radial polynomial).
"""
