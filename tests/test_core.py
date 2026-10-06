"""hilbert.core 的验收测试。

断言曲线次序、编号与坐标的往返、编号相邻格的连续性、阶数与位宽、子方块
区间、矩形的连续区间分解、八个方向的邻格与边界，以及非法入参的报错。
"""

import unittest

from hilbert import (
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

DIRECTIONS = {
    "n": (0, 1),
    "s": (0, -1),
    "e": (1, 0),
    "w": (-1, 0),
    "ne": (1, 1),
    "nw": (-1, 1),
    "se": (1, -1),
    "sw": (-1, -1),
}

SECOND_ORDER = [(0, 0), (1, 0), (1, 1), (0, 1),
                (0, 2), (0, 3), (1, 3), (1, 2),
                (2, 2), (2, 3), (3, 3), (3, 2),
                (3, 1), (2, 1), (2, 0), (3, 0)]


class CurveVectorTests(unittest.TestCase):

    def test_first_orders_visit_the_known_cells(self):
        self.assertEqual([d_to_xy(1, d) for d in range(4)],
                         [(0, 0), (0, 1), (1, 1), (1, 0)])
        self.assertEqual([d_to_xy(2, d) for d in range(16)], SECOND_ORDER)
        for d, cell in enumerate(SECOND_ORDER):
            self.assertEqual(xy_to_d(2, cell[0], cell[1]), d)
        self.assertEqual(xy_to_d(3, 6, 5), 45)
        self.assertEqual(d_to_xy(3, 45), (6, 5))
        self.assertEqual(xy_to_d(3, 7, 0), 63)
        self.assertEqual(d_to_xy(3, 63), (7, 0))


class RoundTripTests(unittest.TestCase):

    def test_every_cell_and_index_round_trips(self):
        for order in (1, 2, 3, 4):
            side = side_length(order)
            limit = side * side
            seen = set()
            for y in range(side):
                for x in range(side):
                    d = xy_to_d(order, x, y)
                    self.assertGreaterEqual(d, 0)
                    self.assertLess(d, limit)
                    self.assertEqual(d_to_xy(order, d), (x, y))
                    seen.add(d)
            self.assertEqual(len(seen), limit)
            for d in range(limit):
                x, y = d_to_xy(order, d)
                self.assertTrue(0 <= x < side and 0 <= y < side)
                self.assertEqual(xy_to_d(order, x, y), d)


class ContinuityTests(unittest.TestCase):

    def test_neighbouring_indices_always_touch(self):
        for order in (1, 2, 3, 4, 5):
            side = side_length(order)
            cells = [d_to_xy(order, d) for d in range(side * side)]
            self.assertEqual(len(set(cells)), side * side)
            self.assertEqual(cells[0], (0, 0))
            self.assertEqual(cells[-1], (side - 1, 0))
            for before, after in zip(cells, cells[1:]):
                self.assertEqual(
                    abs(after[0] - before[0]) + abs(after[1] - before[1]), 1)


class OrderWidthTests(unittest.TestCase):

    def test_side_length_and_index_values_follow_the_order(self):
        self.assertEqual(side_length(1), 2)
        self.assertEqual(side_length(3), 8)
        self.assertEqual(index_bits(1), 2)
        self.assertEqual(index_bits(4), 8)
        self.assertEqual(max_index(1), 3)
        self.assertEqual(max_index(3), 63)
        for order in range(1, 7):
            side = side_length(order)
            self.assertEqual(side, 2 ** order)
            self.assertEqual(index_bits(order), 2 * order)
            self.assertEqual(max_index(order), side * side - 1)

    def test_index_width_and_largest_index_agree(self):
        for order in range(1, 7):
            self.assertEqual(max_index(order), (1 << index_bits(order)) - 1)


class SubSquareTests(unittest.TestCase):

    def test_every_sub_square_covers_one_run_of_indices(self):
        self.assertEqual(quadrant_range(3, 0, 0, 0), (0, 63))
        for order in (1, 2, 3):
            for level in range(order + 1):
                span = 1 << level
                size = 1 << (order - level)
                for qy in range(span):
                    for qx in range(span):
                        lo, hi = quadrant_range(order, level, qx, qy)
                        inside = set()
                        for d in range((1 << order) ** 2):
                            x, y = d_to_xy(order, d)
                            if (qx * size <= x < (qx + 1) * size
                                    and qy * size <= y < (qy + 1) * size):
                                inside.add(d)
                        self.assertEqual(len(inside), size * size)
                        self.assertEqual(lo, min(inside))
                        self.assertEqual(hi, max(inside))
                        self.assertEqual(hi - lo + 1, size * size)


class RegionQueryTests(unittest.TestCase):

    def test_region_query_returns_exactly_the_cells_inside(self):
        rectangles = {
            2: [(0, 0, 3, 3), (0, 0, 0, 0), (1, 1, 2, 2), (2, 0, 3, 1)],
            3: [(0, 0, 7, 7), (0, 0, 1, 3), (2, 3, 5, 6), (7, 0, 7, 7),
                (0, 7, 7, 7), (3, 4, 3, 4), (6, 2, 7, 5)],
            4: [(0, 0, 15, 15), (4, 4, 11, 11), (1, 2, 14, 3)],
        }
        for order, rects in rectangles.items():
            limit = (1 << order) ** 2
            for rect in rects:
                ranges = region_ranges(order, *rect)
                self.assertEqual(ranges, sorted(ranges))
                for (_, hi1), (lo2, _) in zip(ranges, ranges[1:]):
                    self.assertLess(hi1, lo2)
                covered = set()
                for lo, hi in ranges:
                    covered.update(range(lo, hi + 1))
                expected = set()
                for d in range(limit):
                    x, y = d_to_xy(order, d)
                    if rect[0] <= x <= rect[2] and rect[1] <= y <= rect[3]:
                        expected.add(d)
                self.assertEqual(covered, expected)
                self.assertEqual(len(covered), len(expected))

    def test_region_query_keeps_runs_merged(self):
        self.assertEqual(region_ranges(3, 0, 0, 7, 7), [(0, 63)])
        self.assertEqual(region_ranges(3, 3, 3, 3, 3), [(10, 10)])
        self.assertEqual(region_ranges(2, 0, 0, 1, 1), [(0, 3)])
        self.assertEqual(region_ranges(2, 0, 1, 0, 1), [(3, 3)])
        for order, rects in ((3, [(0, 0, 1, 3), (2, 3, 5, 6), (7, 0, 7, 7),
                                  (0, 7, 7, 7)]),
                             (4, [(1, 2, 14, 3), (5, 5, 10, 10)])):
            for rect in rects:
                ranges = region_ranges(order, *rect)
                for (_, hi1), (lo2, _) in zip(ranges, ranges[1:]):
                    self.assertLess(hi1 + 1, lo2)


class AdjacencyTests(unittest.TestCase):

    def test_adjacent_cells_and_their_eight_neighbours(self):
        self.assertEqual(adjacent(2, 0, "e"), 1)
        self.assertEqual(adjacent(2, 0, "n"), 3)
        self.assertEqual(adjacent(2, 0, "ne"), 2)
        self.assertIsNone(adjacent(2, 0, "s"))
        self.assertIsNone(adjacent(2, 0, "w"))
        self.assertIsNone(adjacent(2, 0, "sw"))
        self.assertEqual(adjacent(2, 15, "w"), 14)
        self.assertEqual(adjacent(2, 15, "n"), 12)
        self.assertIsNone(adjacent(2, 15, "e"))
        self.assertIsNone(adjacent(2, 15, "se"))
        for order in (1, 2, 3):
            side = side_length(order)
            for d in range(side * side):
                x, y = d_to_xy(order, d)
                for name, (dx, dy) in DIRECTIONS.items():
                    nx, ny = x + dx, y + dy
                    if 0 <= nx < side and 0 <= ny < side:
                        step = adjacent(order, d, name)
                        self.assertIsNotNone(step)
                        self.assertEqual(d_to_xy(order, step), (nx, ny))
                    else:
                        self.assertIsNone(adjacent(order, d, name))
        around = neighbors(2, 8)
        self.assertEqual(sorted(around), sorted(DIRECTIONS))
        self.assertEqual(around, {"n": 9, "s": 13, "e": 11, "w": 7,
                                  "ne": 10, "nw": 6, "se": 12, "sw": 2})
        self.assertEqual(len(set(around.values())), 8)
        self.assertEqual(neighbors(1, 0),
                         {"n": 1, "s": None, "e": 3, "w": None,
                          "ne": 2, "nw": None, "se": None, "sw": None})
        self.assertEqual(neighbors(2, 15),
                         {"n": 12, "s": None, "e": None, "w": 14,
                          "ne": None, "nw": 13, "se": None, "sw": None})
        for order in (2, 3):
            side = side_length(order)
            for d in range(side * side):
                x, y = d_to_xy(order, d)
                around = neighbors(order, d)
                self.assertEqual(sorted(around), sorted(DIRECTIONS))
                for name, (dx, dy) in DIRECTIONS.items():
                    inside = 0 <= x + dx < side and 0 <= y + dy < side
                    self.assertEqual(around[name] is not None, inside)
                    if inside:
                        self.assertEqual(d_to_xy(order, around[name]),
                                         (x + dx, y + dy))


class InputValidationTests(unittest.TestCase):

    def test_out_of_range_arguments_are_rejected(self):
        self.assertRaises(HilbertError, xy_to_d, 2, 4, 0)
        self.assertRaises(HilbertError, xy_to_d, 2, 0, 4)
        self.assertRaises(HilbertError, xy_to_d, 2, -1, 0)
        self.assertRaises(HilbertError, d_to_xy, 2, 16)
        self.assertRaises(HilbertError, d_to_xy, 2, -1)
        self.assertRaises(HilbertError, xy_to_d, 0, 0, 0)
        self.assertRaises(HilbertError, xy_to_d, MAX_ORDER + 1, 0, 0)
        self.assertRaises(HilbertError, quadrant_range, 3, 1, 2, 0)
        self.assertRaises(HilbertError, quadrant_range, 3, 4, 0, 0)
        self.assertRaises(HilbertError, region_ranges, 3, 4, 0, 1, 1)
        self.assertRaises(HilbertError, region_ranges, 3, 3, 3, 2, 5)
        self.assertRaises(HilbertError, adjacent, 3, 0, "up")
        self.assertRaises(TypeError, xy_to_d, 2, 1.0, 0)
        self.assertRaises(TypeError, d_to_xy, 2, "3")
        self.assertRaises(TypeError, neighbors, 2.0, 0)


if __name__ == "__main__":
    unittest.main()
