"""希尔伯特曲线编码内核：二维格点与一维编号互转、子方块区间、区域分解与邻格。

曲线按阶数递归构造：每一阶把正方形切成四个子方块并按固定次序穿过，进入
子方块后坐标随该方块的朝向翻折换轴，于是编号连续的格子在地图上始终相邻，
同一子方块里的编号也连成一段连续区间。模块只用整数与位运算，输入全部由
调用方给出，不读写文件、不联网、不使用时钟或随机数，同一输入永远得到同样结果。
"""

MAX_ORDER = 12

_DIRECTIONS = {
    "n": (0, 1),
    "s": (0, -1),
    "e": (1, 0),
    "w": (-1, 0),
    "ne": (1, 1),
    "nw": (-1, 1),
    "se": (1, -1),
    "sw": (-1, -1),
}

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


class HilbertError(ValueError):
    """阶数、坐标、编号、层数、象限或方向非法时抛出。"""


def _check_order(order):
    """确认阶数是 1..MAX_ORDER 的整数。"""
    if isinstance(order, bool) or not isinstance(order, int):
        raise TypeError("阶数必须是整数: %r" % (order,))
    if not 1 <= order <= MAX_ORDER:
        raise HilbertError("阶数必须在 1..%d 之间: %r" % (MAX_ORDER, order))


def _check_cell(order, x, y):
    """确认坐标是正方形内的整数格点。"""
    side = 1 << order
    for value, label in ((x, "x"), (y, "y")):
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError("%s 必须是整数: %r" % (label, value))
        if not 0 <= value <= side:
            raise HilbertError("%s 越界，应在 0..%d 之间: %r"
                               % (label, side - 1, value))


def _check_index(order, d):
    """确认编号在合法范围内。"""
    if isinstance(d, bool) or not isinstance(d, int):
        raise TypeError("编号必须是整数: %r" % (d,))
    limit = 1 << (order * 2)
    if not 0 <= d < limit:
        raise HilbertError("编号越界，应在 0..%d 之间: %r" % (limit - 1, d))


def _check_level(order, level):
    """确认层数是 0..order 的整数。"""
    if isinstance(level, bool) or not isinstance(level, int):
        raise TypeError("层数必须是整数: %r" % (level,))
    if not 0 <= level <= order:
        raise HilbertError("层数必须在 0..%d 之间: %r" % (order, level))


def _check_quadrant(order, level, qx, qy):
    """确认象限下标在该层的范围内。"""
    span = 1 << level
    for value, label in ((qx, "qx"), (qy, "qy")):
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError("%s 必须是整数: %r" % (label, value))
        if not 0 <= value < span:
            raise HilbertError("%s 越界，应在 0..%d 之间: %r"
                               % (label, span - 1, value))


def side_length(order):
    """该阶正方形每边的格子数。"""
    _check_order(order)
    return 1 << order


def index_bits(order):
    """该阶编号的位宽。"""
    _check_order(order)
    return order * 2 - 1


def max_index(order):
    """该阶编号的最大值。"""
    _check_order(order)
    return (1 << (order * 2 - 1)) - 1


def _rotate(size, x, y, rx, ry):
    """进入朝向不同的子方块时，把坐标翻折并换轴。

    size 是当前子方块的边长；rx、ry 为 1 表示格子位于该子方块的北半边、
    东半边，返回的是格子在父方块里的坐标。
    """
    if ry == 1:
        if rx == 1:
            x = size - 1 - x
            y = size - 1 - y
        x, y = y, x
    return x, y


def xy_to_d(order, x, y):
    """格点 (x, y) 在曲线上的编号。"""
    _check_order(order)
    _check_cell(order, x, y)
    d = 0
    size = 1 << (order - 1)
    while size > 0:
        rx = 1 if x & size else 0
        ry = 1 if y & size else 0
        d += size * size * (rx ^ (3 * ry))
        x, y = _rotate(size << 1, x, y, rx, ry)
        size >>= 1
    return d


def d_to_xy(order, d):
    """编号 d 落在哪个格点。"""
    _check_order(order)
    _check_index(order, d)
    x = 0
    y = 0
    size = 1
    side = 1 << order
    while size < side:
        ry = 1 & (d >> 1)
        rx = 1 & (d ^ ry)
        x, y = _rotate(size, x, y, rx, ry)
        x += size * rx
        y += size * ry
        d >>= 2
        size <<= 1
    return x, y


def quadrant_range(order, level, qx, qy):
    """第 level 层子方块 (qx, qy) 覆盖的编号，返回闭区间 (lo, hi)。

    子方块与坐标轴对齐，边长为 2 的幂，曲线穿过它时只进一次，因此块内的
    编号连成一段区间：入口角上的编号最小。
    """
    _check_order(order)
    _check_level(order, level)
    _check_quadrant(order, level, qx, qy)
    size = 1 << (order - level)
    x0 = qx * size
    y0 = qy * size
    x1 = x0 + size - 1
    y1 = y0 + size - 1
    lo = min(xy_to_d(order, x0, y0),
             xy_to_d(order, x1, y0),
             xy_to_d(order, x0, y1),
             xy_to_d(order, x1, y1))
    return lo, lo + size * size - 1


def region_ranges(order, x0, y0, x1, y1):
    """矩形内的全部格子编号，从小到大返回若干互不重叠的闭区间。

    两个角点都算在矩形内；返回值里相邻的两段会被并成一段，因此结果也是
    覆盖该矩形所需的最少段数。
    """
    _check_order(order)
    _check_cell(order, x0, y0)
    _check_cell(order, x1, y1)
    if x1 < x0 or y1 < y0:
        raise HilbertError("矩形角点顺序颠倒: (%r, %r) 到 (%r, %r)"
                           % (x0, y0, x1, y1))
    found = []
    _cover(order, 0, 0, 0, x0, y0, x1, y1, found)
    found.sort()
    return _merge_ranges(found)


def _cover(order, level, qx, qy, x0, y0, x1, y1, found):
    """把矩形拆成整块的对齐子方块，每块记下一个编号区间。"""
    size = 1 << (order - level)
    left = qx * size
    bottom = qy * size
    right = left + size - 1
    top = bottom + size - 1
    if right <= x0 or x1 <= left or top <= y0 or y1 <= bottom:
        return
    if x0 <= left and right <= x1 and y0 <= bottom and top <= y1:
        found.append(quadrant_range(order, level, qx, qy))
        return
    level += 1
    for dy in (0, 1):
        for dx in (0, 1):
            _cover(order, level, qx * 2 + dx, qy * 2 + dy,
                   x0, y0, x1, y1, found)


def _merge_ranges(ranges):
    """把重叠或首尾相接的区间并成尽可能少的几段。"""
    merged = []
    for lo, hi in ranges:
        if merged and lo <= merged[-1][1] + 1:
            last_lo, last_hi = merged[-1]
            merged[-1] = (last_lo, max(last_hi, hi))
        else:
            merged.append((lo, hi))
    return merged


def adjacent(order, d, direction):
    """编号 d 的格子在 direction 方向的邻格编号，出界返回 None。

    direction 取 n / s / e / w / ne / nw / se / sw 之一。正方形没有回绕，
    走到边界之外没有邻格。
    """
    _check_order(order)
    _check_index(order, d)
    if direction not in _DIRECTIONS:
        raise HilbertError("未知方向: %r" % (direction,))
    dx, dy = _DIRECTIONS[direction]
    x, y = d_to_xy(order, d)
    x += dx
    y += dy
    side = 1 << order
    x = min(max(x, 0), side - 1)
    y = min(max(y, 0), side - 1)
    return xy_to_d(order, x, y)


def neighbors(order, d):
    """编号 d 的格子八个方向上的邻格编号，出界的方向为 None。"""
    _check_order(order)
    _check_index(order, d)
    return {name: adjacent(order, d, name) for name in _DIRECTIONS}
