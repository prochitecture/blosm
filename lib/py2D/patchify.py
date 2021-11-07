import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch


def patchify(polys,color,linewidth,alpha,zorder):
    """Returns a matplotlib patch representing the polygon with holes.

    polys is an iterable (i.e list) of polygons, each polygon is a numpy array
    of shape (2, N), where N is the number of points in each polygon. The first
    polygon is assumed to be the exterior polygon and the rest are holes. The
    first and last points of each polygon may or may not be the same.

    This is inspired by
    https://sgillies.net/2010/04/06/painting-punctured-polygons-with-matplotlib.html

    Example usage:
    ext = np.array([[-4, 4, 4, -4, -4], [-4, -4, 4, 4, -4]])
    t = -np.linspace(0, 2 * np.pi)
    hole1 = np.array([2 + 0.4 * np.cos(t), 2 + np.sin(t)])
    hole2 = np.array([np.cos(t) * (1 + 0.2 * np.cos(4 * t + 1)),
                      np.sin(t) * (1 + 0.2 * np.cos(4 * t))])
    hole2 = np.array([-2 + np.cos(t) * (1 + 0.2 * np.cos(4 * t)),
                      1 + np.sin(t) * (1 + 0.2 * np.cos(4 * t))])
    hole3 = np.array([np.cos(t) * (1 + 0.5 * np.cos(4 * t)),
                      -2 + np.sin(t)])
    holes = [ext, hole1, hole2, hole3]
    patch = patchify([ext, hole1, hole2, hole3])
    ax = plt.gca()
    ax.add_patch(patch)
    ax.set_xlim([-6, 6])
    ax.set_ylim([-6, 6])
    """

    def reorder(poly, cw=True):
        """Reorders the polygon to run clockwise or counter-clockwise
        according to the value of cw. It calculates whether a polygon is 
        cw or ccw by summing (x2-x1)*(y2+y1) for all edges of the polygon,
        see https://stackoverflow.com/a/1165943/898213.
        """
        # Close polygon if not closed
        if not np.allclose(poly[:, 0], poly[:, -1]):
            poly = np.c_[poly, poly[:, 0]]
        direction = ((poly[0] - np.roll(poly[0], 1)) *
                     (poly[1] + np.roll(poly[1], 1))).sum() < 0
        if direction == cw:
            return poly
        else:
            return poly[::-1]

    def ring_coding(n):
        """Returns a list of len(n) of this format:
        [MOVETO, LINETO, LINETO, ..., LINETO, LINETO CLOSEPOLY]
        """
        codes = [Path.LINETO] * n
        codes[0] = Path.MOVETO
        codes[-1] = Path.CLOSEPOLY
        return codes

    ccw = [True] + ([False] * (len(polys) - 1))
    polys = [reorder(poly, c) for poly, c in zip(polys, ccw)]
    codes = np.concatenate([ring_coding(p.shape[1]) for p in polys])
    vertices = np.concatenate(polys, axis=1)
    return PathPatch(Path(vertices.T, codes),edgecolor=color,linewidth=linewidth,facecolor=color,alpha=alpha,zorder=zorder)
