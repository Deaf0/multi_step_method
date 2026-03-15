import math
from typing import Optional
from typing import List, Tuple, Optional
from dataclasses import dataclass, field


class Point:
    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = x
        self.y = y

    def __eq__(self, other):
        if not isinstance(other, Point):
            return False
        
        eps = 1e-10
        return (abs(self.x - other.x) < eps and
                abs(self.y - other.y) < eps)

    def __add__(self, other: 'Point') -> 'Point':
        return Point(self.x + other.x, self.y + other.y)
    
    def __sub__(self, other: 'Point') -> 'Point':
        return Point(self.x - other.x, self.y - other.y)
    
    def __mul__(self, scalar: float) -> 'Point':
        return Point(self.x * scalar, self.y * scalar)
    
    def __rmul__(self, scalar: float) -> 'Point':
        return self.__mul__(scalar)
    
    def dot(self, other: 'Point') -> float:
        return self.x * other.x + self.y * other.y
    
    def norm(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y)
    
    def __repr__(self) -> str:
        return f"Point({self.x}, {self.y})"


Polygon = List[Point]


@dataclass
class BoundingBox:
    min: Point
    max: Point
    
    def __init__(self, poly: Polygon):
        self.min = Point(float('inf'), float('inf'))
        self.max = Point(-float('inf'), -float('inf'))
        
        for p in poly:
            self.min.x = min(self.min.x, p.x)
            self.min.y = min(self.min.y, p.y)
            self.max.x = max(self.max.x, p.x)
            self.max.y = max(self.max.y, p.y)


@dataclass
class FunctionValue:
    value: float = 0.0
    directions: List[Point] = field(default_factory=list)


@dataclass
class ExtremeDirections:
    min_angle: Point
    max_angle: Point
    max_gap: float = None

    def __eq__(self, other):
        if not isinstance(other, ExtremeDirections):
            return False
        
        return (self.min_angle.x == other.min_angle.x and
                self.min_angle.y == other.min_angle.y and
                self.max_angle.x == other.max_angle.x and
                self.max_angle.y == other.max_angle.y)


def initQ0(A: Polygon, B: Polygon) -> Polygon:
    box_a = BoundingBox(A)
    box_b = BoundingBox(B)
    
    p1 = Point(
        box_a.min.x - box_b.max.x,
        box_a.min.y - box_b.max.y
    )

    p2 = Point(
        box_a.max.x - box_b.min.x,
        box_a.max.y - box_b.min.y
    )
    
    min_x = min(p1.x, p2.x)
    max_x = max(p1.x, p2.x)
    min_y = min(p1.y, p2.y)
    max_y = max(p1.y, p2.y)

    Q0 = [
        Point(min_x, min_y),
        Point(max_x, min_y),
        Point(max_x, max_y),
        Point(min_x, max_y)
    ]

    return Q0


def getCenter(poly: Polygon) -> Point:
    center = Point(0, 0)
    for p in poly:
        center = center + p
    return center * (1.0 / len(poly))


def add_direction(dirs: List[Point], g: Point, eps=1e-10):
    for v in dirs: 
        if (v - g).norm() < eps:
            return
    dirs.append(g)


def compute_F(x: Point, A: Polygon, B: Polygon) -> FunctionValue: 
    result = FunctionValue(-float('inf')) 
    eps = 1e-10 
    
    for a in A: 
        shortest_distance = float('inf')
        closest_points = []
        
        for b in B:
            diff = a - x - b
            diff_norm = diff.norm()

            if diff_norm < shortest_distance - eps:
                shortest_distance = diff_norm
                closest_points = [b]
            elif abs(diff_norm - shortest_distance) < eps:
                closest_points.append(b)

        if shortest_distance > result.value + eps:
            result.value = shortest_distance
            result.directions = []

            for b in closest_points:
                diff = a - x - b
                diff_norm = diff.norm()
                if diff_norm > eps:
                    l = diff * (-1.0 / diff_norm)
                    add_direction(result.directions, l)

        elif abs(shortest_distance - result.value) < eps:
            for b in closest_points:
                diff = a - x - b
                diff_norm = diff.norm()
                if diff_norm > eps:
                    l = diff * (-1.0 / diff_norm)
                    add_direction(result.directions, l)

    for b in B: 
        shortest_distance = float('inf')
        closest_points = []
        
        for a in A:
            diff = b + x - a
            diff_norm = diff.norm()

            if diff_norm < shortest_distance - eps:
                shortest_distance = diff_norm
                closest_points = [a]
            elif abs(diff_norm - shortest_distance) < eps:
                closest_points.append(a)

        if shortest_distance > result.value + eps:
            result.value = shortest_distance
            result.directions = []

            for a in closest_points:
                diff = b + x - a
                diff_norm = diff.norm()
                if diff_norm > eps:
                    l = diff * (1.0 / diff_norm)
                    add_direction(result.directions, l)

        elif abs(shortest_distance - result.value) < eps:
            for a in closest_points:
                diff = b + x - a
                diff_norm = diff.norm()
                if diff_norm > eps:
                    l = diff * (1.0 / diff_norm)
                    add_direction(result.directions, l)

    return result


def findExtremeDirections(subgrads):
    if not subgrads:
        raise RuntimeError("no subgradients")
    
    if len(subgrads) == 1:
        return ExtremeDirections(subgrads[0], subgrads[0])

    angles = []

    for v in subgrads:
        a = (math.atan2(v.y, v.x) + 2*math.pi) % (2*math.pi) 
        angles.append((a, v))

    angles.sort(key=lambda x: x[0])

    n = len(angles)

    max_gap = -1
    max_i = 0

    for i in range(n):
        a1 = angles[i][0]
        a2 = angles[(i+1) % n][0]

        if i == n-1:
            a2 += 2*math.pi

        gap = a2 - a1

        if gap > max_gap:
            max_gap = gap
            max_i = i

    min_vec = angles[(max_i+1) % n][1]
    max_vec = angles[max_i][1]

    return ExtremeDirections(min_vec, max_vec, max_gap)


def signedDistance(origin_point: Point, n: Point, destination_point: Point) -> float:
    return n.dot(destination_point - origin_point)


def clipPolygon(poly: Polygon, origin_point: Point, n: Point) -> Polygon:
    result = []

    size = len(poly)
    for i in range(size):
        curr = poly[i]
        next_point = poly[(i + 1) % size]
        
        curr_dist = signedDistance(origin_point, n, curr)
        next_dist = signedDistance(origin_point, n, next_point)
        
        if curr_dist <= 0:
            result.append(curr)
        
        if curr_dist * next_dist < 0:
            t = curr_dist / (curr_dist - next_dist)
            intersection = Point(
                curr.x + t * (next_point.x - curr.x),
                curr.y + t * (next_point.y - curr.y)
            )
            result.append(intersection)
    
    return result


def isZeroInConvexHull(subgrads: List[Point], extreme_directions: ExtremeDirections) -> bool:
    n = len(subgrads)
    eps = 1e-10

    if n == 0:
        return False
    
    if len(subgrads) == 1:
        return False
    
    if len(subgrads) == 2:  
        return (subgrads[0].dot(subgrads[1]) + 1) < eps
    
    return extreme_directions.max_gap <= math.pi + eps


if __name__ == '__main__':

    # A = [
    #     Point(0,0),
    #     Point(6,0),
    #     Point(3,5)
    # ]

    # B = [
    #     Point(8,2),
    #     Point(10,2),
    #     Point(10,4),
    #     Point(8,4)
    # ]

    # x = Point(-6, -1)

    # function_value = compute_F(x, A, B) 
    # print(f"function_value: {function_value.value}")
    


    