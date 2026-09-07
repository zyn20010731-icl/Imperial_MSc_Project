# Step 3（Y4H5）最终代码指南

版本：2026-09-07，正式重算版本 `run_20260907_formal_v2`  
代码包：`Step3_Y4H5_final_0907`  
适用对象：Imperial College CX3，27 个正式工况

本文件是 Step 3 的最终交接说明，覆盖数据、计算、输出、验证和下载。它取代早期“阶段版”和试跑说明。第四章污染影响的分组计算不属于本任务；本版本只计算 27 个正式 case。

---

## 1. 计算范围与 case 清单

Step2 的 `cases_28.csv` 是 case 顺序和原始文件名的唯一来源。程序读取其中 28 行，排除唯一的：

```text
40D / SRG111_SC39d
```

因此正式计算共 27 个 case：

| 距离 | case 数量 |
|---|---:|
| 5D | 7 |
| 10D | 7 |
| 30D | 7 |
| 40D | 6 |

`SFG312_SC237` 和 `SFG312_SC237d` 在报告中统一写成 `SFG312_SC237d`。原始输入路径仍使用 CSV 中的 `source` 名称，以保证能找到文件。

每个 case 的 phi 标签由 `PhiThEven` 计算为：

```text
phiTag = phi% × 1000 四舍五入后的四位整数，例如 phi0.237 → phi0237
```

---

## 2. 目录和只读约束

### 2.1 CX3 目录

```bash
MASTER=$HOME/MSc_Project/MSc_2026_final
CODE=$MASTER/01_CODE/Step3_Y4H5_from_scratch_0906_v3
RUN_ID=run_20260907_formal_v2
OUT=$MASTER/06_OUTPUTS_Y4H5/$RUN_ID
```

主要输入：

```text
$MASTER/01_CODE/Step2_PLIF_interface_2.0/cases_28.csv
$MASTER/05_OUTPUTS/<distance>/<source>/Y_interface_detection_method1_even_<phiTag>_final_0905/
    interface_holes_<distance>_<source>_<phiTag>_0905.mat
    Step2_PLIF_interface_2.0_DONE_0905.txt
$MASTER/02_INPUTS/STITCHED_PIV/<distance>/<source>/stitched_PIV.nc
$MASTER/02_INPUTS/SPATIAL_CALIBRATION/tform_normal_to_PLIF_matrix.mat
$MASTER/02_INPUTS/SPATIAL_CALIBRATION/tform_tilted_to_PLIF_matrix.mat
```

代码只向 `$MASTER/06_OUTPUTS_Y4H5/$RUN_ID/` 写入。Step1、Step2、原始 PLIF、PIV、标定矩阵以及旧的 `RDS` 文件均按只读方式打开，不能修改、覆盖、重存或改变权限。

每个 case 的 `runtime/identity.json` 保存输入文件的大小、修改时间和代码 SHA-256。计算结束后再次检查这些信息；如果输入在运行期间变化，程序报错。输出 manifest 和 `DONE` 文件只在独立验证通过后生成。

`_work/eta/*.pkl`、`_work/flux/*.pkl` 是断点续跑缓存，不是论文数据；可在确认正式输出已保存后删除。它们不会写回 Step1/Step2。

### 2.2 输出根目录

27 个 case 的结果分别位于：

```text
06_OUTPUTS_Y4H5/run_20260907_formal_v2/<distance>/<canonical-case>/
```

汇总结果位于输出根目录，不在一个名为 `ALL_CASES` 的根目录下：

```text
OUT/5D/ALL_CASES/
OUT/10D/ALL_CASES/
OUT/30D/ALL_CASES/
OUT/40D/ALL_CASES/
OUT/BY_CASE/
OUT/ALL27_HOLE_SIZE_PDF/
OUT/TABLES_0906/
OUT/Y_all_cases_summary_0906.csv
OUT/run_summary_0906.json
```

这是下载时最容易出错的地方：`ALL_CASES` 是每个距离目录的子目录，不能写成 `OUT/ALL_CASES/`。

---

## 3. 固定物理参数和单位

所有 27 个 case 使用同一组参数：

| 符号 | 数值 | 单位 | 用途 |
|---|---:|---|---|
| \(\rho\) | 1000 | kg m⁻³ | 通量 |
| \(\nu\) | \(1.0\times10^{-6}\) | m² s⁻¹ | 耗散率、Kolmogorov 尺度 |
| \(d\) | 0.010 | m | 圆柱直径、归一化 |
| \(U_\infty\) | 0.38 | m s⁻¹ | 无量纲参考速度 |
| PLIF 标定 | 23.92 | pixel mm⁻¹ | 两个方向相同 |
| \(p\) | 23920 | pixel m⁻¹ | 像素物理尺寸 |
| \(\Delta\) | \(1/p\) | m pixel⁻¹ | 面积权重 \(\Delta^2\) |
| \(L\) | 0.030 | m | 固定流向归一化长度 |
| 对称因子 | 2 | — | 半场积分扩展到两侧 |

最终通量均为“单位流向长度”的量：

```text
质量通量       kg/(m s)
流向动量通量   N/m
动能通量       W/m
```

有量纲结果的无量纲参考量为：

\[
q_m^*=\frac{q_m}{\rho U_\infty d},\qquad
q_p^*=\frac{q_p}{\rho U_\infty^2d},\qquad
q_k^*=\frac{q_k}{\tfrac12\rho U_\infty^3d}.
\]

---

## 4. 第一阶段：空间校准和映射

### 4.1 PIV 网格

从 `stitched_PIV.nc` 读取 `gridxFull`、`gridyFull` 和完整 2501 帧速度。要求：

- 数组为规则的 `(time, y, x)`；
- `gridxFull` 沿流向递增；
- `gridyFull` 在文件中实际递减；
- 速度变量为 `vel_xFulls/vel_yFulls` 或 `velxFulls/velyFulls`。

27 个 case 统一优先使用现成的 `smoothed` 速度变量。代码自动识别上面两种变量命名；不再为三个特殊 case 另行拼接 filtered 文件。

### 4.2 相机选择和双线性插值

对每个 PLIF 像素查询点，先用原坐标做几何变换。判断相机归属时使用最近整数像素：

\[
q_r=\lfloor r+0.5\rfloor,\qquad q_c=\lfloor c+0.5\rfloor .
\]

相机归属只用于选择 normal 或 tilted 标定矩阵；几何变换本身仍使用原始浮点坐标。允许恰好落在 PIV 最外层矢量中心的点。PIV 域外不外推，也不换相机重试。

速度采样采用双线性插值：只有权重非零的角点需要有限；只要一个非零权重角点是 NaN，该点即无效。固定映射只计算一次，保存为 `_shared/map.npz` 和 `Y_PLIF_TO_PIV_BILINEAR_MAP_0906.h5`。

通量统计使用有效且位于 \(Y>0\) 工作域的映射点。界面耗散率采样使用整条可映射界面，不额外限制 \(Y>0\)。

---

## 5. 第二阶段：速度、脉动、耗散率和 Kolmogorov 长度

### 5.1 3×3 中值滤波

对每一帧的瞬时 `U/V` 做 3×3 中值滤波：

- 边缘使用最近像素复制；
- 中心像素原来是 NaN 时，输出保持 NaN；
- 中心有效时，在窗口内有限值中取中位数；
- 窗口全为 NaN 时输出 NaN；
- 不再按上游 flag 额外删除有限速度。

结果写入 case 自己的 `_shared/Y_FILTERED_VELOCITY_CACHE_0906.h5`，输入 NetCDF 仍只读。

### 5.2 全时间平均场

用完整 2501 帧 PIV 序列逐格点平均有限样本，得到：

\[
\overline U(x,y)=\langle U_f(x,y,t)\rangle_t,qquad
\overline V(x,y)=\langle V_f(x,y,t)\rangle_t .
\]

这个二维平均场不受 PLIF 坏帧、unknown 或 invalid 标签影响。unknown/invalid 只在映射后的 PLIF 统计中排除。

### 5.3 脉动和梯度

\[
u'=U_f-\overline U,\qquad v'=V_f-\overline V .
\]

耗散率只使用脉动速度梯度。内部网格用中心差分，边界用一阶单边差分，并且沿真实坐标 \(x,y\) 求导；不能把递减的 `gridyFull` 先反转后再改变符号。

对所有有效界面微元的共同有效梯度样本池化：

\[
\varepsilon_I=\nu\left[
-\left\langle\left(\frac{\partial u'}{\partial x}\right)^2\right\rangle
+2\left\langle\left(\frac{\partial u'}{\partial y}\right)^2\right\rangle
+2\left\langle\left(\frac{\partial v'}{\partial x}\right)^2\right\rangle
+8\left\langle\left(\frac{\partial v'}{\partial y}\right)^2\right\rangle
\right],
\]

\[
\eta_I=\left(\frac{\nu^3}{\varepsilon_I}\right)^{1/4}.
\]

若一帧没有任何四个梯度同时有效的界面采样点，该帧只从 \(\eta\) 统计中排除；在它的 hole/flux 数据仍满足条件时，仍可参加第三至第五阶段。

18\(\eta_I\) 是 hole 的 engulfment 距离阈值。论文中相同量的耗散率系数和 \(\eta=(\nu^3/\varepsilon)^{1/4}\) 用于本实现的物理量核对；5D、10D、30D、40D 的共同使用是本项目的既定分析定义。

---

## 6. 第三阶段：hole 分类和 Actual flux

### 6.1 标签和有效像素

Step2 已经给出界面、wake、hole、unknown、invalid 的索引和边界。Step3 不重新分割 PLIF，也不重写 Step2 标签。

- invalid 像素不参加任何计算和统计；
- unknown 像素不参加速度或均值统计，但计入几何概率分母；
- 未标记且有效的像素属于观测 FoV，计入概率分母，E/W 指示量为 0；
- 有几何标签但 U/V 缺失的点仍计入概率分母；
- 速度核要求 U 和 V 同时有限；
- 完整 wake 存在但不与 \(Y>0\) 有效 PIV 相交时，该帧无效。

无界面、无 wake、\(Y>0\) 工作域没有有效速度的帧也排除。没有 engulfment 不构成失败，engulfed 通量可以为 0。

### 6.2 Hole 判据

每个 hole 先在完整 PLIF 图像上计算面积、形心、边界、界面最近距离和尺寸；不能先裁剪到 PIV 域再重算这些量。

设 hole 形心到界面中点树的最近距离为 \(d_{I,h}\)：

\[
\mathrm{engulfed}_h
=\mathbb 1\left[d_{I,h}\le18\eta_I\right].
\]

分辨率门槛使用完整 hole 面积与 2×2 PIV 单元面积比较：

\[
A_h\ge 4\,\mathrm{median}(\Delta x)\,\mathrm{median}(|\Delta y|).
\]

详细表同时保留 All 和 Resolved 计数。通量使用通过该分辨率门槛且落入有效 PIV 映射域的 engulfed hole 像素；internal/wake 使用 Step2 的完整 wake 集合。

### 6.3 瞬时通量核

所有通量均使用 3×3 滤波后的瞬时速度，不使用平均速度构造非线性核：

\[
K_m=U_f,qquad
K_p=U_f^2,qquad
K_k=\frac12(U_f^2+V_f^2)U_f .
\]

每个像素贡献面积为 \(\Delta^2\)。第 \(t\) 帧先把所有 engulfed 或 wake 像素逐个相加，再除以固定流向长度并乘对称因子：

\[
q_{\alpha,s}^{\mathrm{Actual}}(t)
=\frac{2\rho\Delta^2}{L}
\sum_{i\in s(t)}K_{\alpha,i}(t),
\quad s\in\{E,W\}.
\]

最后对共同有效帧取平均：

\[
q_{\alpha,s}^{\mathrm{Actual}}
=\frac1{|K|}\sum_{t\in K}q_{\alpha,s}^{\mathrm{Actual}}(t).
\]

这就是后面 y 和 \(\xi\) 近似的参考值。逐 hole 的原始量也保存，但原始量没有除以 \(L\) 或乘对称因子。

---

## 7. 第四阶段：laboratory-coordinate 的 y 近似

### 7.1 统计顺序

对每个有效映射像素用物理 y 坐标分箱。y 的首箱从 0 开始。所有共同有效帧、所有有效观测像素放在一起池化；不能先对每帧做概率乘积再平均，也不能先平均速度再平方。

池化后计算：

\[
\overline U(y)=\langle U_f\rangle,
\qquad
\overline{U^2}(y)=\langle U_f^2\rangle,
\qquad
\overline{K_k}(y)=\left\langle\frac12(U_f^2+V_f^2)U_f\right\rangle,
\]

\[
P_E(y)=\frac{N_E(y)}{N_P(y)},qquad
P_W(y)=\frac{N_W(y)}{N_P(y)}.
\]

其中 \(N_P\) 是有效几何观测像素数，包含 unknown；速度均值的样本数可以较少。\(P_E\le P_W\le1\) 是统计约束，unknown 不被误写为 background。

剖面密度为：

\[
J_{m,s}(y)=\rho\,\overline U(y)P_s(y),
\]
\[
J_{p,s}(y)=\rho\,\overline{U^2}(y)P_s(y),
\]
\[
J_{k,s}(y)=\rho\,\overline{K_k}(y)P_s(y).
\]

可观测部分统计出的剖面代表完整条带，不再乘观测覆盖比例。于是：

\[
q_{\alpha,s}^{(y)}
=2\Delta y\sum_jJ_{\alpha,s}(y_j).
\]

这里的 2 是对称扩展，\(\Delta y\) 为 y 分箱宽度。y 图为了显示对称性，将已统计的正侧剖面镜像画到负 y；镜像只改变显示，不再把图上的负侧重复积分一次。

`u_bar` 图直接画 \(\overline U\)，纵轴单位为 m/s，不能除以 \(U_\infty\)。

---

## 8. 第五阶段：local-interface-coordinate 的 \(\xi\) 近似

### 8.1 符号和分箱

每个有效像素到当帧界面中点树的最近距离换算为物理距离。wake 侧取正号，界面另一侧取负号：

\[
\xi_i=\begin{cases}
+d_{I,i},&i\in\mathrm{wake},\\
-d_{I,i},&i\notin\mathrm{wake}.
\end{cases}
\]

\(\xi=0\) 单独保留一个箱；正、负侧不能跨零插值。\(\xi_{\max}\) 取本 case 所有有效帧、有效观测像素中实际出现的最大正 \(\xi\)。不使用 `18 eta` 截断，积分上限就是这个 case 的观测 \(\xi_{\max}\)。

### 8.2 池化剖面和像素数权重

和 y 方法一样，先把所有共同有效帧的瞬时核和概率池化成 \(\xi\) 剖面。代码另外保存每帧每个 \(\xi\) 箱中的像素数 \(n_\xi(t,j)\)。因此积分虽然写成 \(\xi\) 积分，实际含义是对每一个像素积分：

\[
q_{\alpha,s}^{(\xi)}(t)
=\frac{2\Delta^2}{L}
\sum_j n_\xi(t,j)J_{\alpha,s}(\xi_j).
\]

再对帧平均：

\[
q_{\alpha,s}^{(\xi)}
=\frac1{|K|}\sum_{t\in K}q_{\alpha,s}^{(\xi)}(t).
\]

所以不同箱的权重由该帧实际出现的像素数量决定；不能改成每个箱统一一个 \(\Delta\xi\) 权重，也不能把每个箱当成相同面积。\(L=30\) mm 只负责单位流向长度归一化。

---

## 9. 第六阶段：输出文件

### 9.1 每个 case 的正式输出

每个 case 目录下的主要文件如下：

```text
results/01_eta/
    Y_eta_case_0906.h5
    Y_eta_frame_terms_0906.csv
    Y_eta_summary_0906.txt

results/02_holes/
    Y_holes_case_0906.h5
    Y_holes_detailed_case_0906.csv
    Y_holes_annotation_case_0906.h5
    Y_hole_counts_by_frame_case_0906.csv
    Y_hole_characteristic_size_PDF_data_case_0906.csv
    Y_hole_characteristic_size_PDF_case_0906.png/.pdf
    Y_holes_case_summary_0906.txt

results/03_actual_flux/
    Y_actual_flux_case_0906.h5
    Y_actual_flux_frames_case_0906.csv
    Y_actual_flux_resolved_engulfed_holes_case_0906.csv
    Y_actual_flux_case_summary_0906.txt

results/04_y_approximation/
    Y_y_approximation_case_0906.h5
    Y_profiles_y_case_0906.csv
    Y_ubar_y_JFM_case_0906.png/.pdf
    Y_PE_y_JFM_case_0906.png/.pdf
    Y_Pwake_y_JFM_case_0906.png/.pdf
    Y_mass_flux_density_y_JFM_case_0906.png/.pdf
    Y_momentum_flux_density_y_JFM_case_0906.png/.pdf
    Y_kinetic_energy_flux_density_y_JFM_case_0906.png/.pdf
    Y_y_approximation_case_summary_0906.txt

results/05_xi_approximation/
    Y_mass_flux_approximations_case_0906.h5
    Y_profiles_xi_case_0906.csv
    Y_xi_fluxes_by_frame_case_0906.csv
    与 04_y_approximation 对应的六类 PNG/PDF 图
    Y_mass_flux_approximations_case_summary_0906.txt

results/06_final/
    Y_FINAL_case_0906.h5
    Y_FINAL_case_summary_0906.json
    Y_FINAL_case_summary_0906.txt
    Y_stage4_diagnostics_case_0906.h5
    diagnostic_4panel/*.png
```

每个正式 case 还应有：

```text
Step3_Y4H5_outputs_0906.tsv
Step3_Y4H5_DONE_0906.txt
runtime/identity.json
runtime/status_0906.json
runtime/validation_0906.json
```

### 9.2 五张 27 行汇总表

`OUT/TABLES_0906/` 是论文表格的正式来源。每张表均有 27 行数据（另加一行表头），并同时输出 CSV、TXT、Overleaf-ready TEX：

| 表 | 内容 | 文件主名 |
|---|---|---|
| A1 | Actual：E/W 有量纲、无量纲、E/W contribution | `Y_flux_table_A1_actual_0906` |
| A2 | y approximation：同上 | `Y_flux_table_A2_yBased_0906` |
| A3 | \(\xi\) approximation：同上 | `Y_flux_table_A3_xiBased_0906` |
| A4 | y、\(\xi\) 相对 Actual 的百分比差异 | `Y_flux_table_A4_percentage_differences_0906` |
| A5 | engulfment contribution：Actual 及 y/\(\xi\) 差异 | `Y_flux_table_A5_contribution_comparison_0906` |

另外：

```text
Y_complete_flux_tables_A1_A5_0906.txt
Y_flux_all_methods_0906.csv/.txt
Y_case_physics_summary_0906.csv/.txt
Y_TABLES_0906_manifest.txt
```

A1–A3 的 contribution 在文件中存为比例，例如 0.35；论文若写百分数需乘 100。A4 已经是百分数。A5 的 CSV 同时给出比例差和 percentage-point 差，不能把比例再乘错一次。

### 9.3 图形规范

普通剖面图同时输出 PNG 和 PDF；PDF 用于论文排版，线条、文字和坐标轴是矢量。四联诊断图保留 PNG，因为背景包含原始 PLIF 和栅格 enstrophy 图像，使用位图更合适。

四联图固定为横向 1×4：

- 每张图旋转显示 90°（逆时针视觉方向）；
- 四张图都使用同一完整 stitched-PIV 覆盖总范围，不拼接 PLIF 背景；
- `(a)` 原始 PLIF 亮度；
- `(b)` 原始亮度归一化后的 `log10`，叠加 interface、engulfed、internal、unknown；
- `(c)` 在 (b) 上叠加滤波后的瞬时速度场；
- `(d)` 同样的线条叠加 `ln(omega_z^2)` enstrophy 背景；
- interface 金黄色，engulfed 绿色，internal 黑色，unknown 灰色；
- 每张子图有自己的 colorbar；
- 不写子图标题，只在图外左上方写 `(a)`–`(d)`；
- 横轴 `x/d`，纵轴 `y/d`；
- 四边都有黑色边框，但只有左、下边框带向外刻度线；上、右边框无刻度线；
- 多 case 曲线使用多种 JFM 风格颜色；图内不写标题。

hole PDF 使用对数坐标：横轴 \(\sqrt{A_h}/\eta_I\)，纵轴 \(P(\sqrt{A_h}/\eta_I)\)，engulfed 用不同颜色的三角形，internal 用方形，并完整写出 `−3 Slope`、`−5 Slope` 参考线标签。

---

## 10. 运行、续跑和验证

### 10.1 本地代码检查

在代码包根目录：

```bash
python3 -m py_compile step3_y4h5/*.py
../.venv-step3/bin/python -m unittest -v tests/test_science.py
```

### 10.2 CX3 提交

```bash
cd "$HOME/MSc_Project/MSc_2026_final/01_CODE/Step3_Y4H5_from_scratch_0906_v3"
bash submit_formal_v4.sh
```

脚本提交 `STEP3_CASE=1..27` 共 27 个 PBS 作业。每个作业申请 1 node、8 CPU、32 GB、4 小时，并把 NumPy/BLAS 线程限制为 1，避免过度抢占。每个作业顺序执行：

1. 单元测试；
2. 本 case 的完整计算；
3. `validate --mark-done`。

查询：

```bash
/opt/pbs/bin/qstat -u iz425
find "$HOME/MSc_Project/MSc_2026_final/06_OUTPUTS_Y4H5/run_20260907_formal_v2" \\
  -name status_0906.json -print
```

只有 27 个 `Step3_Y4H5_DONE_0906.txt` 都存在并通过 manifest 校验后，才能执行汇总：

```bash
qsub run_combine_v4.pbs
```

汇总作业读取已验证的最终 H5，不重新读取 Step1/Step2，也不重算科学量；它生成四个距离的 `ALL_CASES`、`BY_CASE`、27-case hole PDF 和 `TABLES_0906`。

### 10.3 出错处理

- `runtime/status_0906.json` 为 `FAILED` 时，先看该 case 的 `runtime/errors_0906.log` 和 PBS 输出；不能把异常当作普通缺测静默跳过。
- 已写出的 `_work` checkpoint 支持同一 run ID 续跑；科学输入或代码指纹改变时必须使用新的 run ID。
- `DONE` 不代表“部分完成”；它只在全部必需文件、数值恒等式、图像和输入指纹检查通过后产生。
- `run_summary_0906.json` 只有 combine 成功后才出现；不能用单个 case 的 summary 代替它。

---

## 11. Windows PowerShell 下载

目标目录固定为：

```text
C:\Users\iz425\MSc_project\Step3
```

建议使用 SFTP 的绝对远程路径，避免 `cd` 到错误层级。下面的命令只下载论文需要的结果，不下载 `_work/*.pkl` 和大体积 filtered cache：

```powershell
$localRoot = 'C:\Users\iz425\MSc_project\Step3'
$remoteRoot = '/rds/general/user/iz425/home/MSc_Project/MSc_2026_final/06_OUTPUTS_Y4H5/run_20260907_formal_v2'

New-Item -ItemType Directory -Force $localRoot | Out-Null

sftp iz425@login.cx3.hpc.imperial.ac.uk
```

进入 `sftp>` 后：

```text
lcd C:/Users/iz425/MSc_project/Step3
cd /rds/general/user/iz425/home/MSc_Project/MSc_2026_final/06_OUTPUTS_Y4H5/run_20260907_formal_v2

get Y_all_cases_summary_0906.csv
get run_summary_0906.json
get TABLES_0906/Y_complete_flux_tables_A1_A5_0906.txt
get TABLES_0906/Y_TABLES_0906_manifest.txt
get TABLES_0906/Y_flux_table_A1_actual_0906.csv
get TABLES_0906/Y_flux_table_A1_actual_0906.txt
get TABLES_0906/Y_flux_table_A1_actual_0906.tex
get TABLES_0906/Y_flux_table_A2_yBased_0906.csv
get TABLES_0906/Y_flux_table_A2_yBased_0906.txt
get TABLES_0906/Y_flux_table_A2_yBased_0906.tex
get TABLES_0906/Y_flux_table_A3_xiBased_0906.csv
get TABLES_0906/Y_flux_table_A3_xiBased_0906.txt
get TABLES_0906/Y_flux_table_A3_xiBased_0906.tex
get TABLES_0906/Y_flux_table_A4_percentage_differences_0906.csv
get TABLES_0906/Y_flux_table_A4_percentage_differences_0906.txt
get TABLES_0906/Y_flux_table_A4_percentage_differences_0906.tex
get TABLES_0906/Y_flux_table_A5_contribution_comparison_0906.csv
get TABLES_0906/Y_flux_table_A5_contribution_comparison_0906.txt
get TABLES_0906/Y_flux_table_A5_contribution_comparison_0906.tex
```

需要全部图和每个 case 的详细 CSV 时，再分别下载：

```text
get -r 5D/ALL_CASES
get -r 10D/ALL_CASES
get -r 30D/ALL_CASES
get -r 40D/ALL_CASES
get -r BY_CASE
get -r ALL27_HOLE_SIZE_PDF
```

这里的 `ALL_CASES` 必须带距离前缀；`cd` 到输出根目录后直接执行 `get -r ALL_CASES` 会报 `not found`。

下载后在 PowerShell 检查 27 个正式 case：

```powershell
(Get-ChildItem $localRoot -Filter 'Step3_Y4H5_DONE_0906.txt' -Recurse).Count
Get-Content "$localRoot\run_summary_0906.json"
(Get-Content "$localRoot\TABLES_0906\Y_flux_table_A1_actual_0906.csv").Count
```

预期分别为 27、包含 `ALL27_AGGREGATES_WRITTEN` 的 JSON，以及 28 行（表头 + 27 个 case）。

---

## 12. 最终核对清单

提交论文前逐项确认：

1. 27 个 case 都有 `DONE`、manifest 和 validation 文件；
2. 输入 MAT、NetCDF、标定矩阵和 Step2 输出的修改时间/大小未变化；
3. 每个 case 的 `Final/actual` 等于逐帧 Actual 的平均；
4. 每个 case 的 `Final/xiBased` 等于逐帧 \(\xi\) 通量的平均；
5. `eta_I=(nu^3/epsilon_I)^(1/4)`；
6. `u_bar` 的图和 CSV 使用 m/s 的直接 \(\overline U\)，没有除以 \(U_\infty\)；
7. A1–A5 每张表都有 27 个 case 行；
8. A4 是相对 Actual 的百分比，A5 的比例和 percentage-point 字段没有重复乘 100；
9. `TABLES_0906/Y_TABLES_0906_manifest.txt` 中列出的文件都存在；
10. 四个距离的 `ALL_CASES`、`BY_CASE`、`ALL27_HOLE_SIZE_PDF` 和根目录 summary 已下载；
11. 不把 `.pkl` 缓存、filtered velocity cache 或 `_work` 当作论文结果；
12. 普通 PDF 用于论文，四联诊断图使用 PNG。

最终代码包中包含本指南、Python 源码、测试和 PBS 脚本；它不包含任何 Step1/Step2 原始数据或 RDS 文件。
