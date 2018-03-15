
def secant_method(x, fx):
    """Applies secant method to find root x of function f(x). If not enough
    (x, f(x)) value pairs are known to apply secant method, a new x value is
    guessed by slightly changing the initial x value"""
    if fx[-1] != 0:
        if len(x) > 1 and fx[-1] != 0 and abs(fx[-1]) != abs(fx[-2]):
            x0 = x[-2]
            x1 = x[-1]
            fx0 = fx[-2]
            fx1 = fx[-1]
            m = (fx1 - fx0) / (x1 - x0)
            return x1 + (-fx1 / m)
        else:
            return x[0] * 0.9 + 0.0001
    else:
        return x[-1]
