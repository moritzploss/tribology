
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from Constants import PltOpts, CMaps


def plt_surf(m, n, pixel_width, z):
    """
    function to plot output of artificial surf. only used for debugging
    """
    x = np.linspace(0, (m - 1) * pixel_width, m)
    y = np.linspace(0, (n - 1) * pixel_width, n)
    x_grid, y_grid = np.meshgrid(x, y, sparse=False, indexing='ij')

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(x_grid, y_grid, np.transpose(z),
                    rstride=1,
                    cstride=1,
                    cmap=CMaps.default.value,
                    linewidth=0.5,
                    antialiased=True)
    plt.savefig('{}/test-surf.png'.format(os.getcwd()),
                bbox_inches='tight',
                dpi=PltOpts.res.value)
    plt.close()


def artificial_surf(sigma, h, lx, m, n, qr, eng=None):
    """
    Calls matlab function "artificial_surf.m" developed by Mona Mahboob Kanafi.
    See file "artificial_surf.m" for copyright notice
    
    parameters
    sigma: standard deviation, i.e.root - mean - square roughness Rq(m)
    H: Hurst exponent(roughness exponent), 0 <= H <= 1
    It relates to the fractal dimension of a surface by D = 3 - H.
    Lx: length of topography in x direction (m)
    m: number of pixels in x
    n: number of pixels in y

    example:
    [z , PixelWidth, PSD] = artificial_surf(0.5e-3, 0.8, 0.1, 512 , 512);
    generates surface z with sigma = 0.5 mm, H = 0.8 Hurst exponent
    [i.e. a surface with fractal dimension D = 2.2], Lx = 10 cm [the size of
    topography in x direction], m = n = 512 [results in a square surface,
    i.e. Lx = Ly = 10 cm]. With these parameters simply the Pixel Width of
    the final topography is Lx/(m-1) = 1.96e-4.
    """

    z, pixel_width, qx, qy, cq, q, c = \
        eng.artificial_surf(float(sigma), float(h), float(lx), float(m),
                            float(n), float(qr), nargout=7)
    [n, m] = z.size
    return np.asarray(z)*1000, pixel_width, qx, qy, cq, q, c, n, m


if __name__ == "__main__":
    pass
