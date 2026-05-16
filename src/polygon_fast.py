import os
from dataclasses import dataclass
from typing import List

import numpy as np
from scipy.spatial import cKDTree

from geometry import FunctionValue, Point, Polygon, addDirection


def default_workers() -> int:
    count = os.cpu_count() or 1
    return max(1, min(count, 32))


def to_array(polygon: Polygon) -> np.ndarray:
    return np.asarray([(p.x, p.y) for p in polygon], dtype=np.float64)


@dataclass
class PolygonPair:
    A: Polygon
    B: Polygon
    A_arr: np.ndarray
    B_arr: np.ndarray
    tree_A: cKDTree
    tree_B: cKDTree

    @classmethod
    def from_polygons(cls, A: Polygon, B: Polygon) -> "PolygonPair":
        A_arr = to_array(A)
        B_arr = to_array(B)
        return cls(A, B, A_arr, B_arr, cKDTree(A_arr), cKDTree(B_arr))


def _collect_directions_a_to_b(
    result: FunctionValue,
    A_arr: np.ndarray,
    B_arr: np.ndarray,
    tree_B: cKDTree,
    offset: np.ndarray,
    eps: float,
) -> None:
    shifted = A_arr - offset
    dists, _ = tree_B.query(shifted)
    max_dist = float(np.max(dists))

    if max_dist > result.value + eps:
        result.value = max_dist
        result.directions.clear()

    if abs(max_dist - result.value) > eps:
        return

    near_max = np.flatnonzero(dists >= max_dist - eps)
    for i in near_max:
        target = shifted[i]
        r = float(dists[i]) + eps
        for j in tree_B.query_ball_point(target, r=r):
            diff = target - B_arr[j]
            diff_norm = float(np.hypot(diff[0], diff[1]))
            if diff_norm > eps:
                addDirection(result.directions, Point(-diff[0] / diff_norm, -diff[1] / diff_norm))
            else:
                addDirection(result.directions, Point(1.0, 0.0))
                addDirection(result.directions, Point(-1.0, 0.0))


def _collect_directions_b_to_a(
    result: FunctionValue,
    A_arr: np.ndarray,
    B_arr: np.ndarray,
    tree_A: cKDTree,
    offset: np.ndarray,
    eps: float,
) -> None:
    shifted = B_arr + offset
    dists, _ = tree_A.query(shifted)
    max_dist = float(np.max(dists))

    if max_dist > result.value + eps:
        result.value = max_dist
        result.directions.clear()

    if abs(max_dist - result.value) > eps:
        return

    near_max = np.flatnonzero(dists >= max_dist - eps)
    for i in near_max:
        target = shifted[i]
        r = float(dists[i]) + eps
        for j in tree_A.query_ball_point(target, r=r):
            diff = target - A_arr[j]
            diff_norm = float(np.hypot(diff[0], diff[1]))
            if diff_norm > eps:
                addDirection(result.directions, Point(diff[0] / diff_norm, diff[1] / diff_norm))
            else:
                addDirection(result.directions, Point(1.0, 0.0))
                addDirection(result.directions, Point(-1.0, 0.0))


def computeF_fast(x: Point, pair: PolygonPair, eps: float = 1e-10) -> FunctionValue:
    result = FunctionValue()
    offset = np.array([x.x, x.y], dtype=np.float64)
    _collect_directions_a_to_b(result, pair.A_arr, pair.B_arr, pair.tree_B, offset, eps)
    _collect_directions_b_to_a(result, pair.A_arr, pair.B_arr, pair.tree_A, offset, eps)
    return result
