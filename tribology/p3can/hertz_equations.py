import math


def hertz_displ(e1, e2, ny1, ny2, r1x, r1y, r2x, r2y, f):
    """Normal elastic displacement for arbitrary bodies according to hertz
    contact theory"""
    apb = 0.5 * (1 / r1x + 1 / r1y + 1 / r2x + 1 / r2y)
    bma = 0.5 * math.sqrt((1 / r1x - 1 / r1y) ** 2 + (1 / r2x - 1 / r2y) ** 2)
    ra = 1 / (apb - bma)
    rb = 1 / (apb + bma)
    rc = math.sqrt(ra * rb)

    f1 = 1 - math.pow(math.pow(ra / rb, 0.0602) - 1, 1.456)
    f2 = 1 - math.pow(math.pow(ra / rb, 0.0684) - 1, 1.531)
    ec = 1 / ((1 - ny1 ** 2) / e1 + (1 - ny2 ** 2) / e2)

    c = math.pow(3 * f * rc / (4 * ec), 1 / 3) * f1
    e = 1 - math.pow(rb / ra, 4 / 3)
    a = c * math.pow(1 - e ** 2, 1 / 4)
    b = c * math.pow(1 - e ** 2, 1 / 4)
    delta = a * b / rc * (f2 / f1)
    return delta


def get_diam_approx(x_axis, x_profile):
    """Approximate hertz contact radii based on actual contact body geometry"""
    x1 = x_axis[0]
    x2 = x_axis[math.floor(len(x_axis) / 2)]
    x3 = x_axis[-1]
    y1 = x_profile[0]
    y2 = x_profile[math.floor(len(x_axis) / 2)]
    y3 = x_profile[-1]

    a = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    b = math.sqrt((x2 - x3) ** 2 + (y2 - y3) ** 2)
    c = math.sqrt((x3 - x1) ** 2 + (y3 - y1) ** 2)
    s = (a + b + c) / 2
    d = math.sqrt(s * (s - a) * (s - b) * (s - c))
    try:
        diam = 2 * a * b * c / (4 * d)
    except ZeroDivisionError:
        diam = 40000
    return diam


if __name__ == "__main__":
    hertz_displ(210000, 210000, 0.3, 0.3, 12.7 / 2, 12.7 / 2, 41039, 41039, 50)
