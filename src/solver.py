from geometry import (
    Point, Polygon, 
    initQ0, getCenter, 
    computeF, findExtremeDirections, 
    clipPolygon, checkZeroInConvexHull
)
from typing import Tuple


def solve(A: Polygon, B: Polygon, max_iter: int = 1000, eps: float = 1e-11) -> Tuple[float, Point, str]:
    Q = initQ0(A, B)
    x = getCenter(Q)

    best_F = float('inf')
    best_shift = Point(0, 0)

    for _ in range(max_iter):

        function_value = computeF(x, A, B)

        if (function_value.value < best_F):
            best_F = function_value.value
            best_shift = x.copy()

        extreme_directions = findExtremeDirections(function_value.directions)

        Q = clipPolygon(Q, x, extreme_directions.min_angle)

        if extreme_directions.min_angle != extreme_directions.max_angle:
            Q = clipPolygon(Q, x, extreme_directions.max_angle)

        if (not Q):
            return best_F, best_shift, "Q empty"

        new_x = getCenter(Q)

        x_change = (new_x - x).norm()
        x = new_x

        zero_in_hull = checkZeroInConvexHull(
            function_value.directions, 
            extreme_directions
        )

        if zero_in_hull: 
            return best_F, best_shift, "zero in subgradient hull"
            
        if x_change < eps: 
            return best_F, best_shift, "argument stopped changing"
    
    return best_F, best_shift, "reach max iter"