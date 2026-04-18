#!/usr/bin/env python3
"""TMM 计算多层膜反射率。"""

from __future__ import annotations

import argparse
import csv
import json
import math
import cmath
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class Layer:
    n: float
    k: float
    d_nm: float

    @property
    def n_complex(self) -> complex:
        return complex(self.n, -self.k)


def parse_nk(value: list[float] | tuple[float, float]) -> complex:
    if len(value) != 2:
        raise ValueError(f"nk 必须是长度为 2 的数组 [n, k]，收到: {value}")
    n, k = value
    return complex(float(n), -float(k))


def parse_nkd(value: dict[str, Any] | list[float] | tuple[float, float, float]) -> Layer:
    if isinstance(value, dict):
        return Layer(n=float(value["n"]), k=float(value.get("k", 0.0)), d_nm=float(value["d_nm"]))
    if len(value) != 3:
        raise ValueError(f"nkd 必须是长度为 3 的数组 [n, k, d_nm]，收到: {value}")
    n, k, d_nm = value
    return Layer(n=float(n), k=float(k), d_nm=float(d_nm))


def expand_layers(raw_layers: list[Any]) -> list[Layer]:
    layers: list[Layer] = []
    for item in raw_layers:
        if isinstance(item, dict) and "repeat" in item:
            repeat = item["repeat"]
            times = int(repeat["times"])
            block_layers = [parse_nkd(x) for x in repeat["layers"]]
            for _ in range(times):
                layers.extend(block_layers)
        else:
            layers.append(parse_nkd(item))
    return layers


def layer_admittance(nc: complex, cos_theta: complex, pol: str) -> complex:
    if pol == "s":
        return nc * cos_theta
    if pol == "p":
        return cos_theta / nc
    raise ValueError("pol 必须是 's' 或 'p'")


def reflectance_tmm(
    wavelength_nm: float,
    angle_deg: float,
    ambient_nk: complex,
    substrate_nk: complex,
    layers: list[Layer],
    pol: str = "avg",
) -> float:
    theta0 = math.radians(angle_deg)
    sin_theta0 = math.sin(theta0)

    def calc_for_pol(one_pol: str) -> float:
        n0, ns = ambient_nk, substrate_nk

        cos0 = cmath.sqrt(1 - (sin_theta0 / n0) ** 2)
        coss = cmath.sqrt(1 - (n0 * sin_theta0 / ns) ** 2)

        m11, m12, m21, m22 = 1 + 0j, 0 + 0j, 0 + 0j, 1 + 0j

        for layer in layers:
            nj = layer.n_complex
            cosj = cmath.sqrt(1 - (n0 * sin_theta0 / nj) ** 2)
            qj = layer_admittance(nj, cosj, one_pol)
            delta = 2 * math.pi * nj * cosj * layer.d_nm / wavelength_nm

            cd, sd = cmath.cos(delta), cmath.sin(delta)
            a11, a12 = cd, 1j * sd / qj
            a21, a22 = 1j * qj * sd, cd

            nm11 = m11 * a11 + m12 * a21
            nm12 = m11 * a12 + m12 * a22
            nm21 = m21 * a11 + m22 * a21
            nm22 = m21 * a12 + m22 * a22
            m11, m12, m21, m22 = nm11, nm12, nm21, nm22

        q0 = layer_admittance(n0, cos0, one_pol)
        qs = layer_admittance(ns, coss, one_pol)

        num = q0 * m11 + q0 * qs * m12 - m21 - qs * m22
        den = q0 * m11 + q0 * qs * m12 + m21 + qs * m22
        r = num / den
        return float(abs(r) ** 2)

    if pol == "avg":
        return 0.5 * (calc_for_pol("s") + calc_for_pol("p"))
    return calc_for_pol(pol)


def frange(start: float, stop: float, step: float) -> list[float]:
    if step <= 0:
        raise ValueError("step 必须 > 0")
    values: list[float] = []
    x = start
    eps = step * 1e-9
    while x <= stop + eps:
        values.append(x)
        x += step
    return values


def run_config(config: dict[str, Any]) -> tuple[str, list[float], list[float]]:
    ambient_nk = parse_nk(config.get("ambient_nk", [1.0, 0.0]))
    substrate_nk = parse_nk(config["substrate_nk"])
    layers = expand_layers(config["layers"])
    pol = config.get("polarization", "avg")

    mode = config["mode"]
    if mode == "angle_scan":
        wavelength_nm = float(config["wavelength_nm"])
        x = frange(float(config["angle_deg_start"]), float(config["angle_deg_stop"]), float(config["angle_deg_step"]))
        y = [reflectance_tmm(wavelength_nm=wavelength_nm, angle_deg=a, ambient_nk=ambient_nk, substrate_nk=substrate_nk, layers=layers, pol=pol) for a in x]
        return "angle_deg", x, y

    if mode == "wavelength_scan":
        angle_deg = float(config["angle_deg"])
        x = frange(float(config["wavelength_nm_start"]), float(config["wavelength_nm_stop"]), float(config["wavelength_nm_step"]))
        y = [reflectance_tmm(wavelength_nm=wl, angle_deg=angle_deg, ambient_nk=ambient_nk, substrate_nk=substrate_nk, layers=layers, pol=pol) for wl in x]
        return "wavelength_nm", x, y

    raise ValueError("mode 必须是 'angle_scan' 或 'wavelength_scan'")


def save_csv(path: Path, x_name: str, x: list[float], y: list[float]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([x_name, "reflectance"])
        for xi, yi in zip(x, y):
            writer.writerow([f"{xi:.8g}", f"{yi:.10g}"])


def maybe_plot(path: Path, x_name: str, x: list[float], y: list[float]) -> bool:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"[WARN] 无法导入 matplotlib，跳过绘图: {exc}")
        return False

    plt.figure(figsize=(7, 4.2))
    plt.plot(x, y, lw=2)
    plt.xlabel(x_name)
    plt.ylabel("reflectance")
    plt.grid(True, ls="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    return True


def print_template() -> None:
    template = {
        "ambient_nk": [1.0, 0.0],
        "substrate_nk": [1.52, 0.0],
        "polarization": "avg",
        "mode": "wavelength_scan",
        "angle_deg": 0,
        "wavelength_nm_start": 400,
        "wavelength_nm_stop": 800,
        "wavelength_nm_step": 2,
        "layers": [{"repeat": {"times": 5, "layers": [[2.1, 0.0, 55], [1.45, 0.0, 80]]}}, [1.38, 0.0, 105]],
    }
    print(json.dumps(template, ensure_ascii=False, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description="TMM 多层膜反射率计算")
    parser.add_argument("-i", "--input", type=Path, help="输入 JSON 文件")
    parser.add_argument("-o", "--output", type=Path, default=Path("reflectance.csv"), help="输出 CSV")
    parser.add_argument("--plot", type=Path, help="输出曲线图 PNG（可选）")
    parser.add_argument("--print-template", action="store_true", help="打印输入模板 JSON")
    args = parser.parse_args()

    if args.print_template:
        print_template()
        return

    if not args.input:
        raise SystemExit("请指定 --input 配置文件，或使用 --print-template")

    with args.input.open("r", encoding="utf-8") as f:
        config = json.load(f)

    x_name, x, y = run_config(config)
    save_csv(args.output, x_name, x, y)
    print(f"已写出 CSV: {args.output}")

    if args.plot:
        if maybe_plot(args.plot, x_name, x, y):
            print(f"已写出图像: {args.plot}")


if __name__ == "__main__":
    main()
