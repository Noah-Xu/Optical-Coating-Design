# Coating Design

## Knowledge Map

00_基础物理/
01_麦克斯韦与波动/
02_波矢与相位/
03_界面与菲涅尔/
04_导纳与阻抗/
05_多层膜理论/
06_减反射膜设计/
07_宽带与多层优化/
08_带通与截止滤光膜/
09_分光膜与偏振/
10_工程与设计方法/

## Description

Thin film optics / coating design knowledge base.

## TMM Python 工具

仓库新增了 `tmm_reflectance.py`，可用于根据膜层结构计算反射率（Transfer Matrix Method, TMM）。

### 功能

- 输入多层膜结构：基底 `nk`、每层 `nkd`。
- 支持周期膜快捷输入：`repeat.times + repeat.layers`。
- 输出两类曲线：
  - 固定波长，扫描入射角（`angle_scan`）；
  - 固定入射角，扫描波长（`wavelength_scan`）。
- 支持偏振：`s` / `p` / `avg`（默认平均偏振）。
- 输出 CSV；可选输出 PNG 曲线图。

### 运行方式

```bash
python tmm_reflectance.py --print-template
python tmm_reflectance.py -i examples/tmm_input_angle_scan.json -o angle_scan.csv --plot angle_scan.png
python tmm_reflectance.py -i examples/tmm_input_wavelength_scan.json -o wavelength_scan.csv --plot wavelength_scan.png
```

### 输入示例

见：

- `examples/tmm_input_angle_scan.json`
- `examples/tmm_input_wavelength_scan.json`
