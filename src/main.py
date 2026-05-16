import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

from geometry import Point
from polygon_fast import PolygonPair, default_workers
from solver import solve


def get_polygon_from_file(file: Path) -> List[Point]:
    polygon = []
    with file.open("r", encoding="utf-8") as f:
        next(f, None)

        for line in f:
            line = line.strip()
            if not line:
                continue

            x, y = map(float, line.split())
            polygon.append(Point(x, y))

    return polygon


@dataclass
class CaseResult:
    path_to_convex: str
    path_to_nonconvex: str
    hausdorf_dist: float


def _process_subdir(subdir: Path) -> Optional[CaseResult]:
    convex_polygon: List[Point] = []
    nonconvex_polygon: List[Point] = []
    path_to_convex = ""
    path_to_nonconvex = ""

    for file in subdir.iterdir():
        if file.name == f"{subdir.name}_polygon_convex.txt":
            convex_polygon = get_polygon_from_file(file)
            path_to_convex = f"{subdir.name}/{file.name}"
        elif file.name == f"{subdir.name}_polygon_nonconvex.txt":
            nonconvex_polygon = get_polygon_from_file(file)
            path_to_nonconvex = f"{subdir.name}/{file.name}"

    if not convex_polygon or not nonconvex_polygon:
        return None

    pair = PolygonPair.from_polygons(convex_polygon, nonconvex_polygon)
    hausdorf_dist, _, _ = solve(convex_polygon, nonconvex_polygon, pair=pair)
    return CaseResult(path_to_convex, path_to_nonconvex, hausdorf_dist)


def run_batch(data_path: Path, workers: Optional[int]) -> None:
    df = pd.read_csv(data_path / "metadata.csv", sep=";", index_col="path")

    if "hausDistMultiStep" not in df.columns:
        df["hausDistMultiStep"] = None

    subdirs = [p for p in data_path.iterdir() if p.is_dir()]
    n_workers = workers if workers is not None else default_workers()

    if n_workers <= 1 or len(subdirs) <= 1:
        results = [_process_subdir(subdir) for subdir in subdirs]
    else:
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            results = list(
                pool.map(
                    _process_subdir,
                    subdirs,
                    chunksize=max(1, len(subdirs) // (n_workers * 4)),
                )
            )

    for item in results:
        if item is None:
            continue
        if item.path_to_convex in df.index:
            df.at[item.path_to_convex, "hausDistMultiStep"] = item.hausdorf_dist
        if item.path_to_nonconvex in df.index:
            df.at[item.path_to_nonconvex, "hausDistMultiStep"] = item.hausdorf_dist

    df.to_csv(data_path / "metadata.csv", sep=";")


def main() -> None:
    parser = argparse.ArgumentParser(description="УЛЛ multi-step: расстояние Хаусдорфа")
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="число потоков для пакетной обработки (по умолчанию — число ядер, до 32)",
    )
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent.parent
    data_path = base_dir / "data"
    run_batch(data_path, args.workers)


if __name__ == "__main__":
    main()
