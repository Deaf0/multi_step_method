import geometry
from typing import Tuple


def solve(A: geometry.Polygon, B: geometry.Polygon, max_iter: int = 1000, eps: float = 1e-11) -> Tuple[float, geometry.Point, str]:
    Q = geometry.initQ0(A, B)
    x = geometry.getCenter(Q)

    best_F = float('inf')
    best_shift = geometry.Point(0, 0)

    for _ in range(max_iter):

        function_value = geometry.compute_F(x, A, B)

        if (function_value.value < best_F):
            best_F = function_value.value
            best_shift = x.copy()

        extreme_directions = geometry.findExtremeDirections(function_value.directions)

        Q = geometry.clipPolygon(Q, x, extreme_directions.min_angle)

        if extreme_directions.min_angle != extreme_directions.max_angle:
            Q = geometry.clipPolygon(Q, x, extreme_directions.max_angle)

        if (not Q):
            return best_F, best_shift, "Q empty"

        new_x = geometry.getCenter(Q)

        x_change = (new_x - x).norm()
        x = new_x

        zero_in_hull = geometry.isZeroInConvexHull(
            function_value.directions, 
            extreme_directions
        )

        if zero_in_hull: 
            return best_F, best_shift, "zero in subgradient hull"
            
        if x_change < eps: 
            return best_F, best_shift, "argument stopped changing"
    
    return best_F, best_shift, "reach max iter"


if __name__ == "__main__":

    A = [
        geometry.Point(0,0),
        geometry.Point(6,0),
        geometry.Point(3,5)
    ]

    B = [
        geometry.Point(8,2),
        geometry.Point(10,2),
        geometry.Point(10,4),
        geometry.Point(8,4)
    ]

    hausdorf_dist, best_shift, stop_reason = solve(A, B)

    print(f"F(x) = {hausdorf_dist}, x = ({best_shift.x}, {best_shift.y}), reason for stoppage: {stop_reason}")

