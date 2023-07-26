import numpy as np

def poly_winding_number(pnts, poly, return_points=False):
    """
    Return points in polygon using a winding number algorithm in numpy.

    Parameters
    ----------
    pnts : Nx2 array
        Points represented as an x,y array.
    poly : Nx2 array
        Polygon consisting of at least 4 points oriented in a clockwise manner.
    return_points : boolean
        If True, returns mask and points coordinates. If False, return only mask. Default: False.

    Returns
    -------
    The mask (and points if selected) within or on the boundary of the geometry.

    References
    ----------
    `<https://github.com/congma/polygon-inclusion/blob/master/
    polygon_inclusion.py>`_.  inspiration for this numpy version
    """
    x0, y0 = poly[:-1].T  # polygon `from` coordinates
    x1, y1 = poly[1:].T   # polygon `to` coordinates
    x, y = pnts.T         # point coordinates
    y_y0 = y[:, None] - y0
    x_x0 = x[:, None] - x0
    diff_ = (x1 - x0) * y_y0 - (y1 - y0) * x_x0  # diff => einsum in original
    chk1 = (y_y0 >= 0.0)
    chk2 = np.less(y[:, None], y1)  # pnts[:, 1][:, None], poly[1:, 1])
    chk3 = np.sign(diff_).astype('int')
    pos = (chk1 & chk2 & (chk3 > 0)).sum(axis=1, dtype=int)
    neg = (~chk1 & ~chk2 & (chk3 < 0)).sum(axis=1, dtype=int)
    wn = pos - neg
    if return_points:
        out_ = pnts[np.nonzero(wn)]
        return wn.astype('bool'), out_
    return wn.astype('bool')