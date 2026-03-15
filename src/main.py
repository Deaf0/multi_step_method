import geometry


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

Q = geometry.initQ0(A, B)

x = geometry.getCenter(Q)

best_F = float('inf')
best_shift = geometry.Point(0, 0)

iter = 0
MAX_ITER = 1000

while (iter < MAX_ITER):
    function_value = geometry.compute_F(x, A, B)
    if (function_value.value < best_F):
        best_F = function_value.value
        best_shift = geometry.Point(x.x, x.y)

    print(f"Q = {Q}")
    print(f"x = {x}, F(x) = {function_value.value}, dF(x) = {function_value.directions}\n")
    print(f"best_F: {best_F}, best_shift: {best_shift}")

    extreme_directions = geometry.findExtremeDirections(function_value.directions)
    print(f"extreme_directions: {extreme_directions}")
    Q = geometry.clipPolygon(Q, x, extreme_directions.min_angle)
    if extreme_directions.min_angle != extreme_directions.max_angle:
        Q = geometry.clipPolygon(Q, x, extreme_directions.max_angle)

    if (not Q):
        print("Q became empty")
        print(f"Best shift: ({best_shift.x}, {best_shift.y})")
        break

    new_x = geometry.getCenter(Q)
    x_change = (new_x - x).norm()
    x = new_x

    zero_in_hull = geometry.isZeroInConvexHull(function_value.directions, extreme_directions)

    if zero_in_hull or x_change < 1e-11:
        print("STOP: ")
        if zero_in_hull: 
            print("0 in subgradient hull")
        if x_change < 1e-11: 
            print("argument stopped changing")
        print(f"Best shift: ({best_shift.x}, {best_shift.y})")
        print(f"F(x) = {function_value.value}")
        break
    
    iter += 1
