"""hilbert：希尔伯特曲线的编号映射、子方块区间、区域分解与邻格查询内核。"""

from .core import (
    MAX_ORDER,
    HilbertError,
    adjacent,
    d_to_xy,
    index_bits,
    max_index,
    neighbors,
    quadrant_range,
    region_ranges,
    side_length,
    xy_to_d,
)

__all__ = [
    "MAX_ORDER",
    "HilbertError",
    "adjacent",
    "d_to_xy",
    "index_bits",
    "max_index",
    "neighbors",
    "quadrant_range",
    "region_ranges",
    "side_length",
    "xy_to_d",
]
