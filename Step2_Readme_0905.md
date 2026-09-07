# Step2 完整交接文档：算法、参数、路径、HPC 执行与报错规避

文档日期：2026-09-05
用途：把 Step2（界面/洞检测）从算法逻辑到 HPC 实际运行的全部信息整理在一份文档里，供后续对接 Step3 时查阅，避免遗漏任何一个环节。
本文档整合了以下既有项目文档的内容，并补充了本轮（2026-09-05）HPC 并行执行层重构与真机调试过程中新产生的信息：`Step2_v5_final_algorithm_summary.md`、`Step2_v2_label_based_algorithm_spec.md`（v4）、`Step1_Step2_Step3_输出位置总览.md`、`Step2_v2_HPC_rebuild_reference.md`。

---

## 一、Step2 在整个流程里的位置

```
Step1（28 cases，已全部跑完）→ 背景剖面 + 中心线归一化剖面 + 好帧列表
        │
        ▼
Step2（28 cases，本文档主题）→ 界面(interface) + 孔洞(hole) 检测
        │
        ▼
Step3（27 cases，排除 40D/SRG111_SC39d）→ Y4H5 engulfment 三重通量计算
```

Step2 对每个 case 的每一个"好帧"独立处理：读取该 case 自己的 PLIF 原始图 + Step1 算好的背景/归一化剖面，做二值化和连通域分类，判定 interface（界面）、hole（孔洞）、unknown（视场外不确定区域）等标签，画诊断图和高分辨率图，最后把整个 case 的逐帧结果合并成一个 `.mat` 文件。

---

## 二、算法逻辑（v5 最终版，2026-09-04 用户确认）

**重要提醒（详见第三节）**：本节内容基于项目文档 `Step2_v5_final_algorithm_summary.md`。已直接搜索本地保存的、与 CX3 上正在跑的代码逐字节一致的 `Step2_PLIF_interface_2_0_chunk.m` / `merge.m` 源码，确认当前生产代码确实实现的就是本节描述的 v5 label-based 算法（`wake1ID`/`wake1PixelIdx` 等 v5 专属字段全部存在，旧版 `interface_rc`/`envelope_rc` 等字段完全不存在），细节见第三节。

### 2.1 术语与坐标约定

- 图像尺寸固定：`nRows = 1600`（行），`nCols = 2560`（列）。
- 四条边缘：上边缘 `row=1`；下边缘 `row=nRows`；左边缘 `col=1`；右边缘 `col=nCols`。
- wake 中心线：固定列坐标 `yCL_col`（来自 Step1 归一化剖面文件 `Snorm.yCL_col`，取不到时退回 `514.45`）。
- wake 中心线带宽：`centreBandHalfWidth = 25`（像素）。
- 小 hole 面积阈值：`minHoleAreaPx = 4`（像素），固定常数。
- wake1.0 面积阈值（v5 新增）：`200`（像素），固定常数，与 `minHoleAreaPx` 是两个独立常数。
- 二值化阈值 `phi_th`：来自 `cases_28.csv` 的 `PhiThEven` 列，逐 case 不同。
- **连通规则（v5 相对 v4 的关键改动）**：0 侧连通块划分严格 4 连通；1 侧连通块划分改为 8 连通（含对角相邻）。"缝"（interface/hole boundary/unknown boundary 依附的对象）仍只存在于 4 连通意义下共享一条边的两个像素之间——8 连通只影响"两个 1 侧像素算不算同一块"，不影响缝本身怎么定义、怎么画。

### 2.2 预处理与二值化（沿用旧版，未改动）

```
I_clean  = max(I_raw − B_profile_even, 0)
phiStar  = I_clean / max(Phi_ref_even_forNorm, epsDen)
```

非有限值（NaN 或 ±Inf）像素记为 invalid，直接丢弃，不参与任何连通性判断。

```
mask1（1 侧） = phiStar >= phi_th
mask0（0 侧） = phiStar <  phi_th   （仅对有效像素）
```

### 2.3 连通域标注

- 0 侧：`bwconncomp(zeroMask, 4)`，严格 4 连通。
- 1 侧：`bwconncomp(oneMask, 8)`，8 连通，直接作用于原始 1 侧像素，不再做形态学闭运算（v4 版本用过 `imclose` 焊合噪声缺口，v5 已去掉这一步，直接用 8 连通一次性解决）。

### 2.4 0 侧连通块分类：background / unknown / 候选

1. **background**：含有列坐标 `= nCols`（接触右边缘）像素的块，一旦判定不再撤销。一帧没有任何块接触右边缘则整帧无效。
2. **unknown**：剩下的块中，含有列坐标=1、行坐标=1、或行坐标=nRows 像素的块。
3. **候选池**：既非 background 也非 unknown 的块，进入候选池，等第 2.6 节判定。

### 2.5 1 侧连通块分类：wake1.0

需**同时满足两个条件**：
1. 穿过中心线：至少一个像素列坐标落在 `[yCL_col−25, yCL_col+25]`。
2. 面积达标（v5 新增）：连通块像素总数 > 200。

只满足①不满足②：退回"未标记"状态，不参与 interface 判断，但仍有机会在 2.7 节被判定为 wake2.0。

### 2.6 hole 判定（v5 相对 v4 放宽）

对每个候选 0 侧块，只要它 4 连通相邻的 1 侧像素中**至少有一个**带 wake1.0 标签，就判定为 hole（不再要求全部相邻像素都是 wake1.0，这是相对 v4"严格完全包围"的实质性放宽）。

### 2.7 wake2.0（v5 新增）

hole 判定完成后，对仍"未标记"的 1 侧块再判一次：如果它所有相邻的 0 侧像素都属于 2.6 节判定出的完整 hole 集合（面积清理之前）、且至少存在一个这样的相邻像素，判定为 wake2.0，计入最终 wake。wake2.0 采用严格"完全包围"规则。

### 2.8 wake 标签的最终构成

由三部分合并：① 2.5 节的 wake1.0（须同时满足穿过中心线+面积>200）；② 2.6 节选中的全部 hole 像素（不分面积大小）；③ 2.7 节的 wake2.0。

### 2.9 hole 面积清理（沿用旧版）

对 2.6 节选中的每个 hole，像素数 < `minHoleAreaPx=4` 则失去 hole 标签（仍保留 wake 标签）。

### 2.10 三种边界线标签定义

- **interface**：background 像素与"属于 wake1.0 的原 1 侧像素"之间 4 连通意义下的缝，按 wake1.0 逐块编号分组保存。
- **hole boundary**：每个正式 hole 内部 0 侧像素与相邻 1 侧像素之间的缝，按 hole 编号分组。
- **unknown boundary**：每个 unknown 块内部 0 侧像素与相邻 1 侧像素之间的缝，按 unknown 编号分组。

### 2.11 画线算法（v5 整体替换旧版"缝隙中点+距离阈值"方案）

- **真实边**：每条缝画它真实存在的、单位长度的像素边界线段本身，端点坐标都是"整数±0.5"，不压缩成中点。
- **顶点邻接拼接**：按真实端点坐标建立"顶点—缝"邻接表，顺着表把首尾相接的缝串成折线。一个顶点 2 条缝：直接连接；1 条缝：真正的开放端点；3/4 条缝（棋盘格分支点）：按"共享同一个所有权像素"配对，各自独立延伸成折线，不会被强行拼成穿过顶点的整线。
- **拐角平滑**：普通顶点（只连2条缝）按半径 `cornerRadiusPx`（默认 0.15）切圆角；分支顶点保持尖角。
- 这个方案完全不需要距离阈值判断（旧版依赖人为选定的 √2 或 1 这样的阈值），孤立缝如实画成一条独立线段（旧版画"×"）。

### 2.12 可视化与出图约定

- interface → 荧光黄 `[1.00,0.90,0.00]`；hole boundary → 红色 `[0.85,0.05,0.05]`；unknown boundary → 明亮灰色 `[0.65,0.65,0.65]`。
- 底色：纯白到蓝单色渐变（白→`[0.02,0.20,0.55]`，256 级）。
- 坐标轴边框：下、左边缘保留 MATLAB 默认外向刻度线；上、右边缘手动绘制纯黑平滑边框。
- 不写整图大标题；四联图内部 (a)/(b)/(c)/(d) 小标题保留。
- 每张图同时导出矢量 PDF（给 Overleaf/LaTeX）和 1200 dpi 高清 PNG（给 PowerPoint）两种格式。

### 2.13 正式管线产出范围

只产出：四联诊断图（(a) 原始灰度图 (b) 彩色 phiStar (c) 二值 mask (d) phiStar 底图叠加三种线）+ 高分辨率单联图（原始灰度图叠加三种线）。debug 脚本额外产出的整帧/局部放大标签色图不属于正式输出。

### 2.14 已知限制（v5 文档原文记录，尚待更大范围真实数据验证）

- 1 侧改纯 8 连通后理论上不再能焊合 2 像素及以上宽度的真实缺口，只在合成网格单元测试验证过，尚未在全部 28 个 case、更多帧上做大范围真实预览确认。
- hole 判定放宽后，理论上可能让真正处于孤立染料斑内部、但恰好也挨到一点真实 wake 边缘的候选块被误判为 hole，目前只在 frame 2880 一帧上做过目视确认。
- wake2.0 目前只影响"最终 wake 语义集合"，不产生独立可视化线条或统计字段；如后续 Step3 需要单独统计，需要在移植时补充相应输出字段。

---

## 三、【已确认】当前生产代码输出的是 v5 新版（label-based，按 wake1.0/wake2.0 分组）字段结构

项目里原本有两份关于 Step2 输出字段的记录，彼此不一致，因为对应两个不同时期的代码：`Step1_Step2_Step3_输出位置总览.md`（09-01，依据的是 v1.1 旧包，interface 是有序单曲线：`interface_rc`、`envelope_rc`、`interfaceYmin/Ymax_byRow`）vs `Step2_v2_label_based_algorithm_spec.md` 第12/13节（v4/v5 label-based 重新设计，interface 是按 wake1.0 分组的无序像素/线段集合）。

**已直接对本地保存的、与 CX3 上正在跑的完全一致（已用 diff 逐字节核对过）的 `Step2_PLIF_interface_2_0_chunk.m` / `merge.m` 源码做字段名搜索，确认结论如下，不是猜测或待验证：**

`ResultEven` 逐帧字段（直接从源码 `grep "^\s*R\.\w* = "` 得到的完整清单）：

```
frameNo, frameName, status, isValidFrame,
backgroundPixelIdx, backgroundRC, nBackgroundPixels,
unknownID, unknownPixelIdx, unknownRC, unknownBoundaryRC, unknownArea_px, nUnknown,
wake1ID, wake1PixelIdx, wake1RC, wake1Area_px, nWake1,
wake2ID, wake2PixelIdx, wake2RC, wake2Area_px, nWake2,
holeID, holePixelIdx, holeRC, holeBoundaryRC, holeArea_px, nHoles, nSmallHolesFiltered,
interfaceRC, nInterfaceSegments,
wakePixelIdx, wakeRC, nWakePixels, nWakeDyePixels,
invalidPixelIdx, invalidRC, nInvalidPixels
```

**完全没有出现**旧版的 `interface_rc`（下划线）、`envelope_rc`、`interfaceYmin/Ymax_byRow`、`interfaceSpansTopBottom`、`anyBoundaryTouchesTop/Bottom`、`topIntersection_rc`、`bottomIntersection_rc`、`holes_rc`、`holeMeanPhiStar` 这些字段——一个都不存在于当前代码里。

`Meta` 里也确认了一批 v5 专属参数字段（`grep "^\s*Meta\.\w* = "` 得到）：`connectivity`、`centreBandHalfWidth`、`cornerRadiusPx`、`minHoleAreaPx`、`minWakeOneAreaPx`（对应 200 像素 wake1.0 面积门槛）、`holeMeanPhiStarRecheck`、`interfaceCol`/`holeBoundaryCol`/`unknownBoundaryCol`（三种线的配色）、`areaBasedCleanup`/`areaBasedCleanupNote`、`instantaneousSpatialSmoothing`、`reflectionRemoval`，均与 v5 文档描述的算法参数一一对应。

**结论：Step3 现有代码（依赖有序单曲线 `interface_rc` 做弧长采样）无法直接读取当前这一版 Step2 输出，`interfaceRC` 是按 `wake1ID` 分组的无序像素/线段集合，交接给 Step3 前必须先设计"分组像素集合 → 有序曲线（供弧长采样）"的桥接逻辑**——这是真实存在、需要新写代码解决的问题，不是一个还需要去验证的假设。v4 文档第16节已经点名这是遗留课题，截至本文档撰写时刻仍未有对应代码落地。

**每个字段具体的数据类型、形状（比如 `Nx2`/`Nx4`/cell array）、含义，以及 `Meta` 每个字段的含义，见第六节 6.2（逐字段全表，含"局部序号不是跨帧追踪 ID"这类关键陷阱说明）。**

---

## 四、输入：Step1 的输出文件与 `cases_28.csv`

### 4.1 Step2 需要读取的 3 个 Step1 输出文件

| 文件 | 关键字段 | 用途 |
|---|---|---|
| `PLIF_black_stripe_check_widthGT10_PLIF/PLIF_black_stripe_check_widthGT10_results_PLIF.mat` | `frameNoEvenUse`（好帧的 PLIF 帧号列表，严格递增，都是偶数，范围100–5100）、`pivFrameNoEvenUse`（对应 PIV 帧号 `=(frameNoEvenUse-98)/2`）、`badEvenFrameNos`、`badEvenFileNames`、`excludedPIVFrameNos` | 决定这个 case 一共要处理哪些"好帧" |
| `PLIF_background_profile_col1980_2020/PLIF_background_profile_col1980_2020_even_goodFrames.mat` | `B_profile_even`（长度必须等于 1600） | `I_clean = raw - B_profile_even` |
| `PLIF_normalisation_profile_even_final/PLIF_normalisation_profile_rawMean_bgCol1980_2020_even_goodFrames.mat` | `Phi_ref_even_forNorm`（长度1600，全部正有限值）、`yCL_col`（缺失回退514.45） | 归一化分母；wake1.0 中心线判定 |

三个文件如果各自也存了 `frameNoEvenUse`/`pivFrameNoEvenUse`，chunk.m 会额外校验彼此一致，不一致直接报错——这是防止 Step1 三个输出文件互相对不上的保护，已锁定，不需要改。

### 4.2 `cases_28.csv` 结构

列：`Index`（1–28，严格顺序）、`Distance`（如 `5D`/`10D`/`30D`/`40D`）、`CaseName`（如 `SRG111_SC39d`）、`PhiThEven`（该 case 使用的二值化阈值：5D 距离全部 7 个 case 统一固定 **0.350**；10D/30D/40D 各 case 独立拟合值，范围 0.546–0.679）。文件放在 `01_CODE/Step2_PLIF_interface_2.0/cases_28.csv`，与代码同目录。**这个文件本轮没有被改动过，也不应该被覆盖。**

### 4.3 `phiTag` 构造规则

```matlab
phiTag = sprintf('phi%04d', round(phi_th_even * 1000));
% 例：0.350 -> 'phi0350'，0.556 -> 'phi0556'
```

---

## 五、路径与目录结构总览

`MASTER = $HOME/MSc_Project/MSc_2026_final`（已通过 symlink 指向 Ephemeral 实际存储，代码里的路径字符串不用改）。

```
MASTER/
├── 01_CODE/Step2_PLIF_interface_2.0/
│   ├── cases_28.csv
│   ├── Step2_PLIF_interface_2_0_chunk.m      ← 锁定文件（本轮加了 outputTag 参数）
│   ├── Step2_PLIF_interface_2_0_merge.m      ← 锁定文件（本轮加了 outputTag 参数）
│   ├── s2case_driver.m                        ← 本轮新增：case 级并行驱动
│   ├── run_s2case_single.pbs                  ← 本轮新增：单 case 独立提交脚本
│   └── (旧 case.m/case_serial.m — 已废弃，不属于新架构)
├── 02_INPUTS/PLIF/<distance>/<case>/           原始 PLIF tif 序列（只读）
├── 05_OUTPUTS/<distance>/<case>/
│   ├── PLIF_black_stripe_check_widthGT10_PLIF/         (Step1 产物)
│   ├── PLIF_background_profile_col1980_2020/           (Step1 产物)
│   ├── PLIF_normalisation_profile_even_final/          (Step1 产物)
│   └── Y_interface_detection_method1_even_<phiTag>_final_<outputTag>/   (Step2 最终产物)
│       ├── Diagnostic_4panel_<outputTag>/
│       ├── High_resolution_panel_<outputTag>/
│       ├── interface_holes_<distance>_<case>_<phiTag>_<outputTag>.mat
│       ├── interface_check_<distance>_<case>_<phiTag>_<outputTag>.txt
│       └── Step2_PLIF_interface_2.0_DONE_<outputTag>.txt   ← 完成判据
├── 04_LOGS/Step2_PLIF_interface_2.0/<runId>/cases/     ← 每个 case 的详细运行日志
└── 06_MANIFESTS/Step2_PLIF_interface_2.0/<runId>/case_<NN>/
    ├── chunk_<NNN>.mat            ← chunk级中间产物（临时，merge成功后不是最终数据源）
    ├── chunk_<NNN>_DONE.txt
    └── STAGE_case.txt             ← 驱动函数的阶段打点（新架构自己加的，不是老版本产物）
```

**`<outputTag>` 是本轮新引入的执行层参数**（详见第七节），当前这一轮实际跑的值是 `0905`。旧的 `0904` 输出树仍然完整保留在原地，两者互不覆盖。

---

## 六、Step2 最终输出内容

### 6.1 图片

- **四联诊断图**：`Diagnostic_4panel_<outputTag>/`，每个"好帧"一张，1200 dpi，(a) 原始 PLIF 灰度图 (b) 彩色 phiStar (c) 二值 mask (d) phiStar 底图叠加三种线。
- **高分辨率单联图**：`High_resolution_panel_<outputTag>/`，仅第 1 个好帧+之后每隔 100 个好帧一张，独立高分辨率 Panel(f)，真洞用 50% 透明填色。

### 6.2 核心结果 `.mat`：`interface_holes_<distance>_<case>_<phiTag>_<outputTag>.mat`

文件里只存两个顶层变量：`ResultEven` 和 `Meta`（源码里唯一的 `save(...)` 语句就是 `save(temporaryMat, 'ResultEven', 'Meta', '-v7.3')`，没有第三个变量）。以下是**逐字段**的数据类型、形状、含义——直接从 `Step2_PLIF_interface_2_0_chunk.m` 的 `makeEmptyResultV2()`/`extractLabelsV2()`/`diagToResultRecord()` 和 `Step2_PLIF_interface_2_0_merge.m` 的 `Meta.*` 赋值代码逐行核对得到，不是从字段名猜的。

#### 6.2.1 总体约定（对所有字段都适用，先看这一段再看下面的表）

- **坐标系**：`row` = 图像行号（1..`nRows`=1600，竖直方向），`col` = 图像列号（1..`nCols`=2560，水平方向）。这句话就是源码里 `Meta.coordinateConvention` 字面存的内容，直接引用而非转述。
- 每一个 `*RC` 字段都是 `Nx2` 的 `uint32` 矩阵，每一行 `[row col]`。
- 每一个 `*PixelIdx` 字段是 `uint32` 的**线性索引**（MATLAB `ind2sub([1600 2560], idx)` 意义下的列优先线性索引），和同名 `*RC` 字段是同一批像素的两种表示（一个解码过，一个没解码），信息完全重复——Step3 只需要用 `*RC`，不需要自己解 `*PixelIdx`。
- 四个"整体类"字段（`background*`、`wake*`（不含 `wake1`/`wake2`）、`invalid*`）是**单一像素集合**：不分连通块，同类像素全部合并成一个 `Nx1`/`Nx2` 数组。
- 四类"分块类"字段（`wake1*`、`wake2*`、`hole*`、`unknown*`）是 **cell array**：每一帧内部按 `bwconncomp` 的连通块顺序编号 1..n，`xxxID`、`xxxPixelIdx`、`xxxRC`、`xxxArea_px`（以及 `holeBoundaryRC`/`unknownBoundaryRC`）这几个字段的下标 k 互相对应，即"第 k 个块"的 ID、像素、面积、边界是同一个 k。
- **关键陷阱**：`wake1ID`/`wake2ID`/`holeID`/`unknownID` 只是"这一帧内部第几个连通块"的**局部序号**（本质就是 `1:n`），**不是跨帧的追踪 ID**——同一个物理团块在相邻两帧里的编号完全可能不同（甚至块的个数本身也逐帧变化）。如果 Step3 需要跨帧追踪同一个 wake1.0/hole 团块，这个追踪逻辑必须自己另外写，当前 Step2 输出完全没有提供。
- `interfaceRC{k}` 对应的是 `wake1ID(k)` / `wake1RC{k}` 这同一个 wake1.0 块（下标顺序严格一致，`nWake1` 个 cell），不是独立的一套编号；同理 `holeBoundaryRC{k}` 对应 `holeID(k)`，`unknownBoundaryRC{k}` 对应 `unknownID(k)`。
- 所有 `*BoundaryRC` 字段和 `interfaceRC` 都是 `Nx4` `uint32`（`[rowA colA rowB colB]`）——**每一行是一对相邻像素（一条"缝隙/crack"），不是一条已经排好序的曲线**。同一个 cell 内部行与行之间没有顺序关系，相邻两行不保证在空间上首尾相连。这正是第三节结论里说的"Step3 现成的弧长采样代码不能直接吃"的根源：要变成有序曲线，必须自己写"缝隙点集合 → 排序连线"的桥接代码。

#### 6.2.2 `ResultEven`（`Nx1` struct array，N = 这个 case 实际处理的"好帧"数，按 `frameNo` 升序排列，一帧一个元素）

| 字段 | 类型 / 形状 | 含义 |
|---|---|---|
| `frameNo` | `double` 标量 | 这一帧的 PLIF"偶数帧号"（来自 Step1 的 `frameNoEvenUse`，如 100,102,...,5100） |
| `frameName` | `char` | 这一帧对应的原始 PLIF `.tif` 文件名（直接取自 `dir()` 结果的 `.name`） |
| `status` | `char`，取值只有 `'ok'` 或 `'invalid_no_background'` | `'ok'`=这一帧找到了合法的 background 区域，正常完成全部分类；`'invalid_no_background'`=0 侧连通域里没有任何一块摸到图像右边缘（`col=nCols`），无法确定 background，这一帧除 `invalidPixelIdx/invalidRC/nInvalidPixels/nWakeDyePixels` 外其余字段全部保持空（模板默认值）。模板里还有第三个默认值 `'not_processed'`，但只要 merge 顺利完成（本项目 28 个 case 全部如此），这个值不会出现在最终文件里——它只代表"这个位置从未被任何 chunk 写入过"，是留给 Step3 做防御性检查用的。 |
| `isValidFrame` | `logical` 标量 | 等价于 `status=='ok'` |
| `backgroundPixelIdx` | `uint32` `Nbgx1` | background 区域（图像最右侧、代表"环境流体"的那一大片 0 侧区域）全部像素的线性索引，不分连通块合并成一个集合 |
| `backgroundRC` | `uint32` `Nbgx2` `[row col]` | 同上，解码成行列坐标 |
| `nBackgroundPixels` | `double` 标量 | = `size(backgroundRC,1)` |
| `wakePixelIdx` | `uint32` `Nwakex1` | **最终 wake 区域**＝wake1.0 ∪ wake2.0 ∪ 全部 hole（含被面积过滤掉"正式 hole 资格"的小 hole）的像素并集，合并成一个集合 |
| `wakeRC` | `uint32` `Nwakex2` `[row col]` | 同上 |
| `nWakePixels` | `double` 标量 | wake 区域总面积（像素数，含 hole 内部面积） |
| `nWakeDyePixels` | `double` 标量 | 只数 wake1.0 ∪ wake2.0（**不含 hole**）的像素数——即真正被荧光染料覆盖的面积，和 `nWakePixels` 的差就是 hole 内部面积 |
| `nWake1` | `double` 标量 | wake1.0 连通块个数（1 侧、8-连通、跨越中心线带、面积 > `minWakeOneAreaPx`=200px 的块） |
| `wake1ID` | `double` `nWake1x1` | 块的局部序号 `1:nWake1`（见 6.2.1 的"关键陷阱"） |
| `wake1PixelIdx` | `cell{nWake1x1}`，每格 `uint32` 列向量 | 第 k 个 wake1.0 块的像素线性索引 |
| `wake1RC` | `cell{nWake1x1}`，每格 `uint32` `Nkx2` | 第 k 个 wake1.0 块的 `[row col]` |
| `wake1Area_px` | `double` `nWake1x1` | 第 k 个 wake1.0 块的像素面积 |
| `nWake2` | `double` 标量 | wake2.0 连通块个数（1 侧未达到 wake1.0 资格、但被"完全被 hole 包围"规则提升为 wake 的块，v5 新增分类，v4 没有） |
| `wake2ID`/`wake2PixelIdx`/`wake2RC`/`wake2Area_px` | 同 wake1 系列，逐块含义相同 | 针对 wake2.0 块 |
| `nHoles` | `double` 标量 | **正式** hole 个数（面积 ≥ `minHoleAreaPx`=4px，且已通过面积清理过滤） |
| `holeID` | `double` `nHolesx1` | 局部序号 |
| `holePixelIdx`/`holeRC` | `cell{nHolesx1}` | 第 k 个 hole（0 侧、被 wake1.0 像素包围的"洞"）内部的像素 |
| `holeArea_px` | `double` `nHolesx1` | 第 k 个 hole 的面积 |
| `holeBoundaryRC` | `cell{nHolesx1}`，每格 `uint32` `Mkx4` `[rowA colA rowB colB]` | 第 k 个 hole 与相邻 1 侧像素之间的"缝隙对"列表（A=hole 内像素，B=1 侧像素），无序 |
| `nSmallHolesFiltered` | `double` 标量 | 摸到过 wake1.0、但面积 < `minHoleAreaPx` 因而**丢失正式 hole 资格**的候选块个数——这些块的像素仍然计入 `wakePixelIdx`/`nWakePixels`，但不出现在 `holeID`/`holeRC`/`holeArea_px`/`holeBoundaryRC`/`nHoles` 里 |
| `nUnknown` | `double` 标量 | unknown 连通块个数（0 侧、既不摸右边缘也摸到左边缘/上边缘/下边缘的块——因为被图像边界切断，真实归类无法判断，永远不参与 hole/wake 判定） |
| `unknownID`/`unknownPixelIdx`/`unknownRC`/`unknownArea_px` | 同 hole 系列结构 | 针对 unknown 块 |
| `unknownBoundaryRC` | `cell{nUnknownx1}`，`uint32` `Mkx4` | 第 k 个 unknown 块与相邻 1 侧像素的缝隙对列表 |
| `interfaceRC` | `cell{nWake1x1}`，每格 `uint32` `Mkx4` `[rowA colA rowB colB]` | **核心科学产出**：background 与第 k 个 wake1.0 块之间的缝隙对列表（A=background 像素，B=该 wake1.0 块像素；只统计 wake1.0，从不含 wake2.0）；下标 k 与 `wake1ID(k)`/`wake1RC{k}` 严格对应 |
| `nInterfaceSegments` | `double` 标量 | 全部 `interfaceRC` cell 合并后的缝隙对总条数（`sum(cellfun(@(c) size(c,1), interfaceRC))`），**不是** wake1.0 块数（那是 `nWake1`） |
| `invalidPixelIdx`/`invalidRC` | `uint32`，`Ninvx1`/`Ninvx2` | `phiStar` 非有限值（NaN/Inf，例如黑条遮挡、除零）的像素，分类前就被剔除，不参与任何判定 |
| `nInvalidPixels` | `double` 标量 | 同上计数 |

#### 6.2.3 `Meta`（`1x1` struct，整个 case 一份，不分帧）

| 字段 | 类型 / 形状 | 含义 |
|---|---|---|
| `version` | `char` | 固定 `'Step2_PLIF_interface_2.0'`，包标识 |
| `caseIndex` | `double` | 1–28 |
| `distance`/`caseName`/`caseLabel` | `char` | 距离名/case 名/两者拼成的显示字符串 |
| `coordinateConvention` | `char`（长说明文本） | 就是 6.2.1 引用的那段坐标/字段格式说明，源码原文存在 `.mat` 里，Step3 可以直接读出来当"字段说明书"用 |
| `connectivity` | `char` | `'v5: 0-side 4-connectivity, 1-side 8-connectivity...'` |
| `algorithmVersion` | `char` | 含算法版本号和这一轮 `outputTag` 的说明字符串 |
| `phi_th_even` | `double` | 这个 case 实际用的二值化阈值 |
| `phiTag` | `char` | 如 `'phi0350'` |
| `frameNoEvenUse` | `double` 列向量，长度=`nGoodEvenFramesProcessed` | 这个 case 全部"好帧"的 PLIF 帧号列表（与 `ResultEven` 逐元素一一对应） |
| `pivFrameNoEvenUse` | `double` 列向量，同长度 | 对应的 PIV 帧号，`=(frameNoEvenUse-98)/2` |
| `badEvenFrameNos` | `double` 向量 | 被 Step1 标记为坏帧、排除在外的 PLIF 帧号 |
| `badEvenFileNames` | `cell`/字符串数组 | 对应坏帧的文件名 |
| `excludedPIVFrameNos` | `double` 向量 | 被排除的 PIV 帧号 |
| `nNominalEvenFrames` | `double` | 2501（理论总帧数，含好坏帧） |
| `nGoodEvenFramesProcessed` | `double` | = `numel(ResultEven)` |
| `nBadEvenFramesExcluded` | `double` | = `numel(badEvenFrameNos)` |
| `firstGoodEvenFrame`/`lastGoodEvenFrame` | `double` | `frameNoEvenUse` 的 min/max |
| `nStatusOK`/`nStatusInvalidNoBackground` | `double` | `ResultEven.status` 里两种取值各自的帧数统计，两者之和应等于 `nGoodEvenFramesProcessed` |
| `nRows`/`nCols` | `double` | 1600 / 2560，整个项目固定的图像尺寸 |
| `yCL_col` | `double` | 中心线所在列号（Step1 归一化文件里给的，缺失则回退 514.45） |
| `centreBandHalfWidth` | `double` | 25，中心线带半宽（列数），wake1.0 判定用 |
| `minHoleAreaPx` | `double` | 4，hole 正式资格面积门槛 |
| `minWakeOneAreaPx` | `double` | 200，wake1.0 资格面积门槛 |
| `cornerRadiusPx` | `double` | 0.15，画图用的圆角参数，不影响分类 |
| `plifDir`/`caseOutputRoot`/`blackStripeFile`/`backgroundProfileFile`/`normalisationProfileFile`/`outputDir`/`diagnosticDir`/`highResolutionDir` | `char`（绝对路径） | 这一次运行实际读取/写入用的全部路径，留作可追溯凭证 |
| `diagnosticEveryNFrames`/`highResolutionEveryNFrames` | `double` | 20 / 100，诊断图/高分辨率图的抽样间隔 |
| `diagnosticFrameNoPlanned`/`highResolutionFrameNoPlanned` | `double` 向量 | 计划出图的帧号列表（merge 阶段拿它去核对图片是否真的都生成了） |
| `B_profile_even`/`Phi_ref_even_forNorm` | `double` 列向量，长度 1600 | 直接拷贝自 Step1 的背景/归一化 profile，存一份在这里是为了可追溯（不用回头再翻 Step1 文件） |
| `diagnosticResolutionDPI`/`highResolutionPanelDPI` | `double` | 均为 1200 |
| `interfaceCol`/`holeBoundaryCol`/`unknownBoundaryCol` | `double` `1x3`（RGB，0–1） | 三种边界线画图颜色：黄 `[1 0.9 0]`/红 `[0.85 0.05 0.05]`/灰 `[0.65 0.65 0.65]` |
| `reflectionRemoval`/`instantaneousSpatialSmoothing`/`holeMeanPhiStarRecheck` | `logical` | 均为 `false`——三个可选处理步骤本轮均未启用 |
| `areaBasedCleanup` | `logical` | `true`——hole 面积清理确实启用了 |
| `areaBasedCleanupNote` | `char` | 面积清理规则的文字说明 |
| `frameCheckTxt` | `char` | 对应 `interface_check_..._<outputTag>.txt` 的路径 |
| `elapsedMin`/`sumChunkElapsedMin` | `double` | merge 自身耗时（分钟）/ 51 个 chunk 耗时之和 |
| `framesPerChunk`/`nChunksPerCase` | `double` | 50 / 51 |
| `chunkManifestDirectory` | `char` | 这个 case 的 `06_MANIFESTS/.../case_NN` 路径 |
| `description` | `char`（长说明文本） | 算法逻辑的完整文字总结，内容与第二节一致 |
| `requestedWorkers` | `double`，**固定值 10** | ⚠️**已确认是锁定文件里的硬编码历史遗留值，和这一轮实际运行配置无关**——`merge.m` 函数签名根本没有 workers 参数，这个 10 是写死的字面量，不随本轮实际调用改变。本轮实际每个 case 用的是 **8** 个 worker（`s2case_driver(caseIndex, runId, 8, outputTag)`），真实值只记在 PBS 日志和 `s2case_driver.m` 调用参数里，`Meta.requestedWorkers` 不可信，Step3 不要用它做任何判断 |
| `executionArchitecture` | `char`，**同样是硬编码历史文本** | 字面写着"PBS case array, at most 8 cases; 10 MATLAB process workers per case"——这描述的是最初设计草稿的执行方案，**不是本轮实际用的方案**（本轮是 27 个独立 qsub 任务，不是 array，见第七节 7.4）。这个字段和 `requestedWorkers` 一样，只是锁定文件里留下的旧文本，读它时必须知道它可能与实际执行方式不符，以实际 PBS 提交记录（第七节/第十二节）为准 |

**给 Step3 的额外提醒**：`Meta` 里除了 `requestedWorkers`/`executionArchitecture` 这两个已知的"文本与实际不符"字段外，其余字段（路径、阈值、frameNo 列表等）都是当次运行现读现算/现拷贝得到的真实值，可以信任。

### 6.3 其它产出

- `interface_check_<distance>_<case>_<phiTag>_<outputTag>.txt`：逐帧检测结果的文本记录。
- `Step2_PLIF_interface_2.0_DONE_<outputTag>.txt`：完成判据，merge.m 在重新生成前会先删除旧的同名文件再重新生成——**这意味着 merge 没有"已完成就跳过"的保护，重跑 merge 一定会整个重做一次**，即使这个 case 早就成功过。

### 6.4 v7.3 (.mat) 格式提醒

所有 `.mat` 输出都用 `-v7.3`（HDF5）格式保存。Python 用 `scipy.io.loadmat` 读不了，需要 `h5py`/`mat73`/`pymatreader`；MATLAB 直接 `load()` 即可。

---

## 七、代码架构：锁定文件 vs 执行层

### 7.1 硬约束（本轮已与用户确认并遵守）

用户原话："计算，画图，输出内容和格式这些跟算法有关的都保持完全不变，但结构可跟着运行优化的需要做调整"——`Step2_PLIF_interface_2_0_chunk.m` 和 `Step2_PLIF_interface_2_0_merge.m` 两个文件的**输入参数含义、内部算法、输出内容/格式**不允许改变；**执行层（谁在什么时候、以什么方式调用这两个文件）允许重写**。

### 7.2 本轮对锁定文件的唯一改动：`outputTag` 参数化

两个文件的函数签名都从两参数变成三参数（第三个可选，默认 `'0904'`，保证向后兼容）：

```matlab
Step2_PLIF_interface_2_0_chunk(taskIndex, runId, outputTag)
Step2_PLIF_interface_2_0_merge(caseIndex, runId, outputTag)
```

原来每处硬编码的 `"0904"` 字符串（输出目录名、诊断图子目录名、DONE 标记文件名、`interface_check_*.txt` 文件名、`interface_holes_*.mat` 文件名、`Meta.algorithmVersion` 字符串）全部替换成 `outputTag` 变量。**已用 `diff -u` 逐字节核对原始上传 zip 与补丁后文件，确认只有这些命名相关的 hunk 发生变化，计算/画图逻辑一个字节没动。**

补丁脚本 `apply_outputTag_patch.py`（放在 `01_CODE/Step2_PLIF_interface_2.0/` 同目录）用"精确字符串匹配、必须恰好命中1次否则整体中止不写入任何内容"的方式实现，避免对未通读全文的大文件做手工重打字带来的转录风险。已在 CX3 上成功运行过一次（`chunk.m` 5处替换、`merge.m` 7处替换）。

### 7.3 执行层：`s2case_driver.m`（本轮新写，替代旧的 case.m/case_serial.m）

单 case 驱动函数，职责：

1. 扫描该 case 的 51 个 chunk，跳过已有 `chunk_NNN.mat` + `chunk_NNN_DONE.txt` 的（断点续算）。
2. 剩余 chunk 用 `parcluster('Processes')` + `parpool`（`JobStorageLocation` 建在 `$TMPDIR` 下，node-local）起进程池，`parfor` 并行调用锁定的 `Step2_PLIF_interface_2_0_chunk(taskIndex, runId, outputTag)`。
3. 验证全部 51 个 chunk 真的都完成。
4. 调用锁定的 `Step2_PLIF_interface_2_0_merge(caseIndex, runId, outputTag)`。
5. 全程用 `writeDriverStage` 往 `STAGE_case.txt` 打点（`driver_entered`→`chunk_scan_done`→`pool_starting`→`pool_started`→`parfor_completed`→`about_to_merge`→`complete`，或 `FAILED: ...`），用 `fopen/fprintf/fclose` 而非标准输出，且只由驱动进程自己写（不在 parfor worker 里写），避免并发写冲突，也保证进程被静默杀死时仍能定位卡在哪一步。

28 个 case 完全独立：各自的 `cases_28.csv` 行、各自的 PLIF 输入、各自的输出文件夹、各自的 `06_MANIFESTS/.../case_<NN>/` 检查点目录。

### 7.4 PBS 提交方式的演变（本轮实测得出的最终结论）

**最初方案**（PBS array `-J 2-28%8`）**在这次实测中不可用**：即使单个 case 的资源申请（8核256GB）在只提交一个非 array job 时能秒排上，把同样的资源申请放进 array（`%8` 并发上限）后，长时间卡在 `Not Running: Placement set is too small: Insufficient amount of resource`，即便 `pbsnodes` 能看到多个完全空闲、资源规格匹配的节点。真机测试后确认：**问题出在 array 本身，不是资源数字的问题**（用同样的 8核256GB 提交一个不带 `-J` 的独立测试 job，立刻排上运行）。

**最终方案**：放弃 array，把 27 个 case（2–28）当成 **27 个完全独立的 qsub 提交**，`CASE_INDEX` 通过 `qsub -v CASE_INDEX=$i` 传入（脚本内部用 `CASE_INDEX="${CASE_INDEX:?CASE_INDEX is not set}"` 读取，替代原来的 `PBS_ARRAY_INDEX`）：

```bash
CODEDIR="$HOME/MSc_Project/MSc_2026_final/01_CODE/Step2_PLIF_interface_2.0"
cd "$CODEDIR"
for i in $(seq 2 28); do
    qsub -N "S2IF_C$(printf '%02d' $i)" -v STEP2_RUN_ID="run_20260905_parallel_v1",CASE_INDEX=$i run_s2case_single.pbs
done
```

这样提交后没有了 array 自带的 `%N` 并发阀门，PBS 按普通 job 的方式各自独立判断能不能排，实测立刻有多个 case 转为 R（运行中），其余正常排队（Q，非"资源不足"报错）。

**代价/需要留意的点**：失去了 array 自带的并发节流保护——如果集群资源短期内非常宽裕，27 个 case 可能同时冲上去，对共享文件系统（尤其 merge 阶段单个 case 写 20+GB 文件）、MATLAB 并发 license 数量造成压力。如果观察到 license 报错或 merge 阶段明显变慢卡住，需要手动分批提交，而不是回头用 array（已证实这台机器上 array 对这种大资源请求不可靠）。

### 7.5 当前实际运行配置（截至本文档撰写时刻）

- `RUN_ID = run_20260905_parallel_v1`
- `outputTag = '0905'`
- PBS 资源：`select=1:ncpus=8:mem=256gb`，`walltime=08:00:00`
- 每 case 内部 8 个 MATLAB worker（`s2case_driver(caseIndex, runId, 8, '0905')`）
- case 1：已用非独立-job 方式（`validate_case1_parallel.pbs`）单独验证成功（STAGE=complete，PBS Exit_status=0，最终 23GB `.mat`、DONE 标记齐全，无残留临时文件）
- case 2–28：用上面 7.4 节的独立提交方式在跑

---

## 八、HPC 安全调用规则（PBS → MATLAB `-batch`，本轮及历史调试共同验证）

1. **不要在 shell 里把 `matlab -batch "..."` 命令字符串再包一层 shell 级 `try/catch`**——驱动函数自己已有 try/catch 并写 STAGE 文件，外层不需要重复包裹。
2. **日志重定向用单一 `{ 命令块 } > "$LOGFILE" 2>&1`**，不要拆成"先 `>` 写头部、再 `>>` 追加"两段式。
3. **`set +e` 贯穿全脚本**，不要用 `set -euo pipefail` 再局部切换。
4. **不要引入 `$HOME` 覆盖**（`FAKE_HOME` 一类技巧）——已直接测试排除，不是必需项。
5. **不要引入任何跨进程的 gate/锁目录机制**——历史上"HPC-STABLE"复杂脚本引入过这类机制，没有实锤证明必要，反而显著提高复杂度和排查难度。
6. **`MATLAB_PREFDIR` 固定指向 `$TMPDIR` 下的子目录**（node-local，PBS 分配的每任务独立临时目录），绝不指向 RDS/Ephemeral 共享路径——这是 v1.1 时代最重要的教训之一，MATLAB 的启动文件/锁文件/进程池控制数据放在网络文件系统上会导致近乎瞬间失败。
7. **`module purge` → `module load tools/prod`（允许静默失败）→ `module load MATLAB/2024b`**，这个顺序历次验证都成功，维持不变。
8. **PBS array 的元素范围不能是单元素**（`-J 1-1` 非法，会报 `illegal -J value`），至少要 `-J 1-2%1` 这种区间。
9. **`-batch` 后面的 MATLAB 命令字符串必须是单一逻辑行**，写成带字面换行的多行字符串会报 "No MATLAB command specified for -batch command line argument"，要折叠成用逗号/分号连接的单行。
10. **本轮新增教训——大资源请求的 array job 在这台集群上不可靠**（详见 7.4 节），单 case 独立 qsub 提交是目前唯一验证过对大资源请求可靠的方式。

---

## 九、报错与坑：历史遗留 + 本轮实测（完整清单）

### 9.1 v1.1 时代历史教训（书面回顾，8点，原始日志未留存但结论已验证有效）

1. 早期出现过 1428 个任务（28 case × 51 chunk）的调度失败问题。
2. 阈值选取修正：5D 距离统一固定用 0.350，不用 Step1 各自的 Method3 拟合值；10D/30D/40D 才用各自拟合值。
3. `/tmp` 配额耗尽：MATLAB 需要 node-local 临时文件系统存放启动文件/锁文件/进程池数据，但绝不能把大体积 PNG/chunk MAT payload 也放进 `/tmp`——只存运行时小文件，正式产物直接写 RDS/Ephemeral 持久化路径。
4. MATLAB 在 RDS（网络文件系统）上启动会近乎瞬间失败（CPU 时间趋近于0）——本轮曾怀疑 HOME 覆盖/gate 机制导致类似问题，直接测试后排除。
5. Home 目录配额耗尽，促成整个项目迁移到 Ephemeral 存储。
6. 通过 symlink 完成迁移，`$HOME` 路径仍是软链接指向 Ephemeral，代码里路径字符串不用改——本轮日志确认这个迁移依然有效。
7. chunk 级 checkpoint 续算机制（`chunk_NNN.mat` + `chunk_NNN_DONE.txt` 同时存在才跳过）——已原样沿用进新架构。
8. v1.1 曾完整验证成功过一次。

### 9.2 本轮（09-05）执行层重构时的调试问题清单

- **PBS array + 大资源请求排不上号**（详见 7.4 节）：即使单个非 array job 秒排上，同样资源规格放进 array 就长期卡在 `Placement set is too small`，即使 `pbsnodes` 显示有完全空闲、规格匹配的节点。改成 27 个独立 qsub job 后立刻解决。**根因未 100% 确认**，怀疑与该集群调度器对 array 的 `%N` 并发上限做整体资源预留判断有关，但没有做逐项单变量隔离测试证实。
- **walltime 卡在 merge 阶段**：第一次验证 case 1 用 1 小时 walltime，在 merge 阶段（序列化写 ~22-23GB `.mat` 文件）被 PBS 杀掉（`walltime 3627 exceeded limit 3600`，只超了27秒）。教训：merge 阶段需要给足够 walltime 余量，不要卡着预估值设置。
- **真实内存峰值 vs 申请内存**：case 1 完整跑完后，`qstat -xf` 显示 `resources_used.mem = 147495776kb ≈ 140.7GB`——这是目前唯一一个真实测得的峰值数据。本轮申请的内存（一路从256GB试到192GB又试回256GB）最终定在256GB，**主要是为了绕开 array 调度怪异行为、不是因为科学计算真的需要256GB**——如果之后要精简资源申请，140.7GB+合理余量（比如180-200GB）是有真实数据支撑的更精确参考值,但注意这是仅从 case 1（5D/NoGrid）测得的，其它 27 个 case 的帧数/复杂度不同，峰值可能有出入，缩减前建议先看看已完成 case 的峰值分布再定。
- **`qalter` 对正在运行（R状态）的 job 增加 walltime 不生效**：`qalter -lwalltime=05:00:00` 对一个 R 状态的 job 执行后，`qstat -f` 显示 walltime 没有变化——本轮遇到过，未深究是权限限制还是调度器策略，绕过方式是直接撤掉重新提交（对还在排队、Time Use=0 的 job 有效；对已经在跑、已经产生真实计算进度的 job，撤掉重提会丢失进度，需要谨慎，好在本项目有 chunk 级断点续算，就算撤了重投也不会丢已完成的 chunk）。
- **"error" 关键字误报**：用 `grep -i error` 扫描日志会大量命中锁定代码里本来就有的、无害的常规诊断输出——例如 `Error policy: any safety-check failure stops this case immediately`（很可能是每处理一个 chunk 就打印一次的策略说明文字，不是真出错了）、以及 MATLAB 自带的带 HTML 超链接格式的函数名引用（`<a href="matlab:...errorDocCallback...">`），这些在已确认 100% 成功的 case 1 日志里同样大量出现。**正确的报错检测方式**：看 `MATLAB exit code: N` 是否为 0（脚本自己在每次调用后打印），或者用更精确的关键词如 `^Error using`、`Undefined function`、`Out of memory`、`nonzero exit`，不要用裸的 `error` 子串匹配。
- **heredoc 粘贴可能丢行**：至少一次实测中，一次约80行的 `cat > file << 'EOF' ... EOF` 粘贴后 `wc -l` 显示只有72行（应为79行），中间几行在粘贴过程中被截断粘连。**每次 heredoc 粘贴创建/覆盖脚本文件后，必须用 `wc -l` 核对行数与预期完全一致，再继续下一步**，不要假设粘贴一定完整。
- **merge.m 没有"已完成即跳过"保护**：不像 chunk 级有断点续算，merge 每次被调用都会完整重做一次拼接+保存+校验，即使这个 case 早就成功完成过。批量重跑/重新提交整个流程时要注意，不要对已经成功的 case 无意义地重新触发 merge（会浪费时间重写一次20+GB文件，虽然不会产生错误结果，但纯属浪费）。
- **不要对锁定文件做整体重打字**：本轮对 chunk.m/merge.m 只读过局部内容（没有逐行通读全文件），任何改动都用"精确字符串匹配、必须恰好命中1次否则全部中止不写入"的补丁方式完成（Edit 工具的手工编辑 + 独立 Python 补丁脚本两种方式互相校验，结果字节级一致），避免凭记忆重建大文件内容带来的转录风险。
- **旧有可疑复杂度（未使用，仅记录以防未来误踩）**：早期"HPC-STABLE"版本脚本引入过 gate/锁目录互斥机制（`_MATLAB_2024b_START_GATE`，心跳文件+重试退避）、把 matlab 放后台执行再轮询状态等做法，虽然没有实锤证明是失败直接原因，但显著提高复杂度，新架构完全没有采用，也不建议以后重新引入。
- **`.mat` 输出里 `Meta.requestedWorkers`（写死的10）和 `Meta.executionArchitecture`（写死的"array, 8 cases, 10 workers"文字）是锁定文件里的历史遗留硬编码，和本轮实际配置（8 worker/case、27 个独立 job、非 array）不符**——详见第六节 6.2.3 表格最后两行。读 `Meta` 做任何"这次到底怎么跑的"判断时，不要信这两个字段,以第七节 7.5 记录的真实配置为准。

---

## 十、给 Step3 交接的关键提示

1. **第三节已确认：当前 Step2 输出是 v5 label-based 字段结构，Step3 现有代码无法直接读取**——对接前必须先设计"分组像素集合→有序曲线（供弧长采样）"的桥接逻辑，这是最优先、最重要的一步，不是可以跳过再看的事。
2. **RUN_ID 不要硬编码**：本轮用的是 `run_20260905_parallel_v1`，输出 tag 是 `0905`；如果之后重新跑，两者都会变，读取代码要能发现/传入这两个值，不要写死。
3. **完成判据要按 outputTag 拼出实际文件名再检查**：`Step2_PLIF_interface_2.0_DONE_<outputTag>.txt` 存在，才代表这个 case 真正完成；本轮之前遗留的 `_0904` 后缀输出树也还在原地，不要和 `_0905` 的搞混。
4. **`phiTag` 现算，不要写死**：`sprintf('phi%04d', round(threshold*1000))`。
5. **v7.3 (.mat) 读取**：Python 用 `h5py`/`mat73`，MATLAB 直接 `load()`。
6. **28 个 case 里第 28 个（`40D/SRG111_SC39d`）Step3 明确排除**，Step3 只处理 27 个 case，用 `cases_27.csv`（比 `cases_28.csv` 多 `ExpectedNX/NY/UVar/VVar` 列）。跨步骤核对 case 时按 `(Distance, CaseName)` 字符串匹配，不要直接假设两个 CSV 的 `Index` 数字通用。
7. **v5 算法相对 v1.1/v4 的实质性变化**（第三节已确认最终数据确实是 v5 风格）：interface 不再是有序单曲线，而是按 wake1.0 编号分组的无序像素/线段集合（`interfaceRC`）；新增 wake1.0/wake2.0 的区分；hole 判定标准放宽。这直接影响 Step3 现有的"沿 interface 弧长采样"逻辑，必须先规划桥接方案才能对接，不是等对接时才会发现的潜在问题。

---

## 十一、当前状态快照（2026-09-05 全部完成后更新）

- Step1：28 个 case 全部在 CX3 跑完（历史状态，本轮未重跑）。
- Step2：`run_20260905_parallel_v1` / `outputTag=0905`：**28 个 case 全部成功完成**（`select=1:ncpus=8:mem=256gb`，`walltime=08:00:00`，8 workers/case；case 1 用单独的验证脚本先跑通，case 2–28 用 27 个独立 PBS job 提交，详见 7.4/7.5 节）。终极核查六项全部通过：28 个 case 的 `STAGE_case.txt` 全部为 `complete` 且 51/51 chunk 齐全；`05_OUTPUTS` 下找到正好 28 个 `*_final_0905` 输出文件夹；每个文件夹 DONE 标记、最终 `interface_holes_*.mat`（22GB–38GB，随 distance 递增，5D 最小、40D 最大，与预期一致）齐全、无残留临时文件；日志里精确报错关键词零命中；`MATLAB exit code` 全部为 0；PBS 层面 28 个 job（1 个 validate + 27 个独立 case）`Exit_status` 全部为 0。
- Step3：v8 版本，只有 case 1（`run_20260901_003025`，5D/NoGrid）严格校验通过，其余 26 个 case 和距离级汇总尚未提交/跑完（详见项目文档 `Step3_v8_HPC_review_and_fixes.md`）——这份 Step3 结果是基于**旧的 Step2 输出**（v1.1 时代的字段结构）跑出来的。第三节已确认新一轮 Step2（v5 label-based）字段结构与旧版有实质差异（`interfaceRC` 分组无序集合 vs 旧版 `interface_rc` 有序单曲线），Step3 现有代码**不能**直接对接新一轮 Step2 输出，case 1 已跑出的 Step3 结果也不会自动对新的 Step2 输出成立——这不是需要"重新验证"的疑虑，而是已确认存在、需要新写桥接代码才能解决的差异。

---

## 十二、附录：本轮 28 个 case 的实际输出路径完整清单

以下路径均已在终极核查中逐一确认：DONE 标记存在、最终 `.mat` 文件存在（大小已列出）、无残留临时文件。根目录前缀均为 `$MASTER/05_OUTPUTS/`，`$MASTER = $HOME/MSc_Project/MSc_2026_final`。

| Distance | CaseName | phiTag | 最终 `.mat` 大小 | 完整目录（`$MASTER/05_OUTPUTS/` 之后的部分） |
|---|---|---|---|---|
| 5D | NoGrid | phi0350 | 23G | `5D/NoGrid/Y_interface_detection_method1_even_phi0350_final_0905/` |
| 5D | SFG312_SC11d | phi0350 | 23G | `5D/SFG312_SC11d/Y_interface_detection_method1_even_phi0350_final_0905/` |
| 5D | SFG312_SC237 | phi0350 | 24G | `5D/SFG312_SC237/Y_interface_detection_method1_even_phi0350_final_0905/` |
| 5D | SFG354_SC15d | phi0350 | 25G | `5D/SFG354_SC15d/Y_interface_detection_method1_even_phi0350_final_0905/` |
| 5D | SRG111_SC22d | phi0350 | 25G | `5D/SRG111_SC22d/Y_interface_detection_method1_even_phi0350_final_0905/` |
| 5D | SRG111_SC39d | phi0350 | 25G | `5D/SRG111_SC39d/Y_interface_detection_method1_even_phi0350_final_0905/`（Step3 明确排除的 case） |
| 5D | SRG38_SC20d | phi0350 | 25G | `5D/SRG38_SC20d/Y_interface_detection_method1_even_phi0350_final_0905/` |
| 10D | NoGrid | phi0556 | 22G | `10D/NoGrid/Y_interface_detection_method1_even_phi0556_final_0905/` |
| 10D | SFG312_SC11d | phi0623 | 22G | `10D/SFG312_SC11d/Y_interface_detection_method1_even_phi0623_final_0905/` |
| 10D | SFG312_SC237d | phi0584 | 23G | `10D/SFG312_SC237d/Y_interface_detection_method1_even_phi0584_final_0905/` |
| 10D | SFG354_SC15d | phi0581 | 23G | `10D/SFG354_SC15d/Y_interface_detection_method1_even_phi0581_final_0905/` |
| 10D | SRG111_SC22d | phi0666 | 23G | `10D/SRG111_SC22d/Y_interface_detection_method1_even_phi0666_final_0905/` |
| 10D | SRG111_SC39d | phi0679 | 23G | `10D/SRG111_SC39d/Y_interface_detection_method1_even_phi0679_final_0905/` |
| 10D | SRG38_SC20d | phi0570 | 23G | `10D/SRG38_SC20d/Y_interface_detection_method1_even_phi0570_final_0905/` |
| 30D | NoGrid | phi0638 | 27G | `30D/NoGrid/Y_interface_detection_method1_even_phi0638_final_0905/` |
| 30D | SFG312_SC11d | phi0625 | 30G | `30D/SFG312_SC11d/Y_interface_detection_method1_even_phi0625_final_0905/` |
| 30D | SFG312_SC237 | phi0608 | 28G | `30D/SFG312_SC237/Y_interface_detection_method1_even_phi0608_final_0905/` |
| 30D | SFG354_SC15d | phi0592 | 29G | `30D/SFG354_SC15d/Y_interface_detection_method1_even_phi0592_final_0905/` |
| 30D | SRG111_SC22d | phi0559 | 35G | `30D/SRG111_SC22d/Y_interface_detection_method1_even_phi0559_final_0905/` |
| 30D | SRG111_SC39d | phi0582 | 35G | `30D/SRG111_SC39d/Y_interface_detection_method1_even_phi0582_final_0905/` |
| 30D | SRG38_SC20d | phi0614 | 31G | `30D/SRG38_SC20d/Y_interface_detection_method1_even_phi0614_final_0905/` |
| 40D | NoGrid | phi0637 | 29G | `40D/NoGrid/Y_interface_detection_method1_even_phi0637_final_0905/` |
| 40D | SFG312_SC11d | phi0632 | 33G | `40D/SFG312_SC11d/Y_interface_detection_method1_even_phi0632_final_0905/` |
| 40D | SFG312_SC237d | phi0596 | 31G | `40D/SFG312_SC237d/Y_interface_detection_method1_even_phi0596_final_0905/` |
| 40D | SFG354_SC15d | phi0597 | 35G | `40D/SFG354_SC15d/Y_interface_detection_method1_even_phi0597_final_0905/` |
| 40D | SRG111_SC22d | phi0546 | 37G | `40D/SRG111_SC22d/Y_interface_detection_method1_even_phi0546_final_0905/` |
| 40D | SRG111_SC39d | phi0566 | 38G | `40D/SRG111_SC39d/Y_interface_detection_method1_even_phi0566_final_0905/` |
| 40D | SRG38_SC20d | phi0599 | 32G | `40D/SRG38_SC20d/Y_interface_detection_method1_even_phi0599_final_0905/` |

共 28 行，与 `cases_28.csv` 的 28 个 case 一一对应（按字符串排序展示，不代表 `cases_28.csv` 里的 `Index` 顺序）。每个目录内部结构一致：`Diagnostic_4panel_0905/`、`High_resolution_panel_0905/`、`interface_check_<distance>_<case>_<phiTag>_0905.txt`、`interface_holes_<distance>_<case>_<phiTag>_0905.mat`（上表大小列）、`Step2_PLIF_interface_2.0_DONE_0905.txt`。

注意：`5D` 距离 7 个 case 的 `phiTag` 全部是 `phi0350`（对应文档第4.3节提到的固定阈值 0.350），`10D`/`30D`/`40D` 各 case 的 `phiTag` 各不相同——这与项目文档 `Step1_Step2_Step3_输出位置总览.md` 里记录的阈值规则完全吻合，是这次真机数据对既有文档的一次交叉验证。
