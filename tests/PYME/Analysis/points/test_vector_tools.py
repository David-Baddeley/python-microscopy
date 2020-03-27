
from PYME.Analysis.points import vector_tools
import numpy as np

# Make a shell for coordinate testing
from skimage import morphology
r = 49
ball = morphology.ball(r, dtype=int)
SPHERICAL_SHELL = morphology.binary_dilation(ball) - ball
X, Y, Z = np.where(SPHERICAL_SHELL)
X_C, Y_C, Z_C = X - r, Y - r, Z - r
R = np.sqrt(X_C**2 + Y_C**2 + Z_C**2)

def test_cart2sph():
    azi, zen, r = vector_tools.cart2sph(X_C, Y_C, Z_C)
    np.testing.assert_almost_equal(r, R)

def test_cartesian_to_spherical():
    azi, zen, r = vector_tools.cartesian_to_spherical(X_C, Y_C, Z_C)
    np.testing.assert_almost_equal(r, R)

def test_spherical_to_cartesian():
    azi, zen, r = vector_tools.cartesian_to_spherical(X_C, Y_C, Z_C)
    x, y, z = vector_tools.spherical_to_cartesian(azi, zen, r)
    np.testing.assert_array_almost_equal(X_C, x)
    np.testing.assert_array_almost_equal(Y_C, y)
    np.testing.assert_array_almost_equal(Z_C, z)

def test_find_principal_axes():
    x = np.arange(1000)
    y = np.zeros_like(x)
    z = y

    standard_deviations, principal_axes = vector_tools.find_principal_axes(x, y, z, sample_fraction=1.)

    np.testing.assert_array_equal([0., 0.], standard_deviations[1:])
    np.testing.assert_almost_equal(np.std(x), standard_deviations[0], decimal=0)

def test_unity_scaled_projection():
    x = np.arange(1000)
    y = np.zeros_like(x)
    z = np.zeros_like(x)

    scaling_axes = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    scaling_factors = [1., 1., 1.]
    xs, ys, zs = vector_tools.scaled_projection(x, y, z, scaling_factors, scaling_axes)
    np.testing.assert_array_equal([x, y, z], [xs, ys, zs])

def test_scalar_projection():
    x = np.arange(1000)
    y = np.zeros_like(x)
    z = np.zeros_like(x)

    scaling_axes = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    scaling_factors = [5., 5., 5.]
    xs, ys, zs = vector_tools.scaled_projection(x, y, z, scaling_factors, scaling_axes)
    np.testing.assert_array_equal([5*x, y, z], [xs, ys, zs])

def test_rotated_projection():
    x = np.arange(1000)
    y = np.zeros_like(x)
    z = np.zeros_like(x)

    scaling_axes = np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., 1.]])
    scaling_factors = [5., 5., 5.]
    xs, ys, zs = vector_tools.scaled_projection(x, y, z, scaling_factors, scaling_axes)
    np.testing.assert_array_equal([y, 5*x, z], [xs, ys, zs])


def test_direction_to_n_nearest_points():
    test_n = 5
    x0, y0, z0 = 0., 0., 0.
    x = np.arange(2 * test_n)
    y = x[:]
    z = 4 * test_n * np.ones_like(x)
    z[:test_n] = x[:test_n]

    az, zen, r, v = vector_tools.direction_to_nearest_n_points(x, y, z, x0, y0, z0, test_n)
    np.testing.assert_almost_equal(az, np.pi / 4)
    np.testing.assert_almost_equal(r, np.sqrt(3) * 2)
    np.testing.assert_almost_equal(zen, np.arccos(2./r))

    # try shifting
    x0, y0, z0 = 1., 1., 1.
    az, zen, r, v = vector_tools.direction_to_nearest_n_points(x, y, z, x0, y0, z0, test_n)
    np.testing.assert_almost_equal(az, np.pi / 4)
    np.testing.assert_almost_equal(r, np.sqrt(3))
    np.testing.assert_almost_equal(zen, np.arccos(1. / r))

    # try multiple query points
    x0 = np.arange(2.)
    y0, z0 = x0[:], x0[:]
    az, zen, r, v = vector_tools.direction_to_nearest_n_points(x, y, z, x0, y0, z0, test_n)
    np.testing.assert_array_almost_equal(az, (np.pi / 4) * np.ones_like(x0))
    np.testing.assert_almost_equal(r, np.array([2, 1]) * np.sqrt(3))
    np.testing.assert_almost_equal(zen, np.arccos(np.array([2, 1]) / r))

def test_find_points_within_cylinder():
    x, y = np.arange(10), np.arange(10)
    # r = np.sqrt((x ** 2)[:, None] + (y ** 2)[None, :])
    xx, yy, zz = np.meshgrid(x, y, np.arange(10))
    # rr = np.broadcast_to(r, zz.shape)
    v0 = np.array([0., 0., 1.])
    v1, v2 = np.array([1., 0., 0.]), np.array([0., 1., 0.]),
    radius = 2.5
    length = 2.
    x0, y0, z0 = 4.5, 4.5, 4.5
    inside, axial_distance = vector_tools.find_points_within_cylinder(xx.ravel(), yy.ravel(), zz.ravel(), x0, y0, z0, radius, length, v0, v1, v2)
    for pi, point in enumerate(zip(xx.ravel(), yy.ravel(), zz.ravel())):
        if inside[pi]:
            assert (np.sqrt((point[0] - x0) ** 2 + (point[1] - y0) ** 2) <= radius)
            assert (point[2] <= z0 + length and point[2] >= z0)
        else:
            outside_radially = np.sqrt((point[0] - x0) ** 2 + (point[1] - y0) ** 2) > radius
            outside_axially = (point[2] > z0 + length or point[2] < z0)
            assert (outside_radially or outside_axially)

def test_distance_to_image_mask():
    from PYME.IO import tabular
    from PYME.IO.image import ImageStack
    from PYME.IO.MetaDataHandler import CachingMDHandler
    x_p = np.arange(4)
    radius = 3
    x = np.arange(10)
    x0, y0, z0 = np.array([0, 0, 0])  # np.array([4.5, 4.5, 4.5])
    points = tabular.DictSource({'x': x_p + x0, 'y': x_p + y0, 'z': x_p + z0})
    points.mdh = CachingMDHandler({'voxelsize.x': 1, 'voxelsize.y': 1})
    # xx, yy, zz = np.meshgrid(x, y, z)

    r = np.sqrt((x[:, None, None] - x0) ** 2 + (x[None, :, None] - y0) ** 2 + (x[None, None, :] - z0) ** 2)
    mask = ImageStack(r <= radius)

    distances = vector_tools.distance_to_image_mask(mask, points)
    np.testing.assert_almost_equal(np.zeros(2), distances[:2])
    assert(distances[2] == 1)
    assert(distances[-1] == np.sqrt((3 - 2)**2 + (3 - 2)**2 + (3 - 1)**2))
