# hilbert

只依赖 Python 标准库的希尔伯特曲线编码内核：二维格点与一维编号的互转、阶数对应的
边长与编号位宽、与坐标轴对齐的子方块覆盖的连续编号区间、矩形查询拆成的若干段连续
区间，以及一个格子在八个方向上的邻格与边界行为。所有输入由调用方给出，不使用真实
时钟、网络或随机数，同一输入永远得到同样的结果。

- `hilbert/core.py` — 编码内核：编号映射与逆映射、子方块区间、区域分解、邻格。
- `tests/test_core.py` — 验收用例。

## 运行测试

在项目根目录执行：

```
python3 -m unittest discover -s tests -v
```

Windows 上如果没有 `python3`，可用：

```
python -m unittest discover -s tests -v
```
