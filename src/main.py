import re
import pandas as pd
from geometry import Point
from solver import solve
from pathlib import Path
from typing import List


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

base_dir = Path(__file__).resolve().parent.parent
data_path = base_dir / 'data'

df = pd.read_csv('data/metadata.csv', sep=';', index_col='path')

if 'hausDistMultiStep' not in df.columns:
    df['hausDistMultiStep'] = None

for subdir in data_path.iterdir():
    if not subdir.is_dir():
        continue

    convex_polygon = []
    nonconvex_polygon = []
    path_to_convex = ''
    path_to_nonconvex = ''
    hausdorf_dist = None

    for file in subdir.iterdir():
        if file.name == f"{subdir.name}_polygon_convex.txt":
            convex_polygon = get_polygon_from_file(file)
            path_to_convex = f"{subdir.name}/{file.name}"
        elif file.name == f"{subdir.name}_polygon_nonconvex.txt":
            nonconvex_polygon = get_polygon_from_file(file)
            path_to_nonconvex = f"{subdir.name}/{file.name}"
    if convex_polygon and nonconvex_polygon:
        hausdorf_dist, _, _ = solve(convex_polygon, nonconvex_polygon)

        if path_to_convex in df.index:
            df.at[path_to_convex, 'hausDistMultiStep'] = hausdorf_dist
            path_to_convex = ''

        if path_to_nonconvex in df.index:
            df.at[path_to_nonconvex, 'hausDistMultiStep'] = hausdorf_dist
            path_to_nonconvex = ''

df.to_csv(data_path / 'metadata.csv', sep=';')


