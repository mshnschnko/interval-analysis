# from intervalpy import Interval
# from intervalpy import Interval
from Interval import *


def min_delta_search(interval: Interval, delta: int | float) -> int | float:

    def f(interval: Interval, delta: int | float) -> int | float:
        if abs(interval.left) <= abs(interval.right):
            return interval.left + delta
        else:
            return interval.right - delta

    min_delta = delta
    a = 0
    b = delta
    while abs(f(interval, min_delta)) > 0.01:
        min_delta = (a + b) / 2
        if f(interval, a) * f(interval, min_delta) < 0:
            b = min_delta
        else:
            a = min_delta
    return min_delta


if __name__ == "__main__":
    print("Enter A_11, A_12, A_21, A_22: ", end='')
    a, b, c, d = [float(n) for n in input().split()]
    if (a < 0 or b < 0 or c < 0 or d < 0):
        raise ValueError("Values must be non-negative.")
    delta = min(a, b, c, d)
    print(delta)
    A = [[Interval(a - delta, a + delta), Interval(b - delta, b + delta)],
         [Interval(c - delta, c + delta), Interval(d - delta, d + delta)]]
    
    det_A = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    
    print(det_A)

    if 0 in Interval(det_A.left + delta, det_A.right - delta):
        min_delta = min_delta_search(det_A, delta)
        print(min_delta)

    # print(Interval(det_A.left + delta, det_A.right - delta))