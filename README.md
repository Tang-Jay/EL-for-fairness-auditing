# EL for Fairness Auditing

**仓库地址：** [https://github.com/Tang-Jay/EL-for-fairness-auditing](https://github.com/Tang-Jay/EL-for-fairness-auditing)

基于**经验似然（Empirical Likelihood, EL）**与**欧式经验似然（EEL）**的公平性审计复现与实验代码。本仓库包含论文中表格与图形的完整 R 实现，涵盖覆盖率比较、运行时间、功效分析以及 COMPAS 数据的公平性置信区间与 FFR 检验。

## 项目简介

- **EL / EEL**：在分组损失/公平性约束下构造置信区间与假设检验，无需 Bootstrap 重抽样。
- **公平性审计**：对 COMPAS 等预测系统做组间 PPV（阳性预测值）差异的置信区间估计，以及对 FFR（False Flagging Rate）的 EL 检验与多重比较（BH 程序）。

对比方法包括：Bootstrap、EL、EEL、T-test。

## 环境要求

- **R**（建议 4.0+）
- 常用包：`ggplot2`, `dplyr`, `readr`（Fig4–6 绘图与数据处理）
- 可选：`Cairo`（导出高质量图形）、`emplik`（Tab3 中 EL 检验）、`sysfonts`, `showtext`（Fig4–6 中文字体）

在 R 中安装依赖：

```r
install.packages(c("ggplot2", "dplyr", "readr", "Cairo", "emplik"))
```

## 项目结构

```
EL-for-fairness-auditing/
├── Tab1 and Fig1/     # Table 1 覆盖率 + Figure 1 QQ 图
├── Tab2/              # Table 2 运行时间（EL/EEL vs Bootstrap）
├── Tab3/              # Table 3 FFR 的 EL p 值与 BH 程序
├── Figs2-3/           # 功效比较图（Power curves）
├── Fig4-6/            # COMPAS 公平性审计图（PPV 差异置信区间等）
└── README.md
```

### 各目录说明

| 目录 | 内容 | 主要脚本 |
|------|------|----------|
| **Tab1 and Fig1** | Bootstrap / EL / EEL 覆盖率表；EL、EEL 与卡方分位数的 QQ 图 | `Table1.R`, `Figure1.R`, `Generate-EL-Data.R`, `Generate-Bt-Data.R` |
| **Tab2** | Model1 下 EL、EEL 与 Candes 方法的运行时间对比表 | `Table2.R`, `Generate-Model1-EL-Data.R`, `Generate-Model1-Bt-Data.R` |
| **Tab3** | FFR 的 EL 检验与 Benjamini–Hochberg 多重检验 | `Table3.R`（使用 `emplik`包） |
| **Figs2-3** | 功效比较（T04/T14/T24/T34 等统计量） | `Figs2-3.R`, `Generate-Data.R` |
| **Fig4-6** | COMPAS 数据：非裔/白人、性别-年龄等子群的 PPV 差异 90%/95% 置信区间 | `Fig4.R`, `Fig5a.R`, `Fig5b.R`, `Fig6a.R`, `Fig6b.R`, `Generate-CI-AfricanVSCaucasian.R` 等 |


各子目录中通常包含：

- `GlambdaChen.R`：Chen 型经验似然 Lagrange 乘子求解（核心 EL 计算）
- `data/`：该目录下脚本所需的 CSV 数据（部分需先运行“生成数据”脚本得到）

## 使用方法

### 1. Table 1（覆盖率）

在 `Tab1 and Fig1/` 下：

1. 生成 Bootstrap 覆盖率数据（按需调整 `n`, `m`, `modelname`）：
   ```r
   source("Generate-Bt-Data.R")
   ```
2. 生成 EL/EEL 数据：
   ```r
   source("Generate-EL-Data.R")
   ```
3. 合并并输出 Table 1：
   ```r
   source("Table1.R")
   ```
   或调用 `table1_merge(modelname, ns, ms, a)`。

### 2. Figure 1（QQ 图）

在 `Tab1 and Fig1/` 下，确保已有对应 `Model1`/`Model2` 的 EL、EEL 数据后：

```r
source("Figure1.R")
```

### 3. Table 2（运行时间）

在 `Tab2/` 下需先有 `data/` 中 `Time-Model1-EL-EEL-*.csv` 与 `Time-Model1-Candes-*.csv`，然后：

```r
source("Table2.R")
# 或指定 n：table2_merge(400)
```

### 4. Table 3（FFR）

在 `Tab3/` 下：

```r
source("Table3.R")
```

依赖 `emplik::el.test` 做单组 EL 检验。

### 5. Figs2–3（功效图）

在 `Figs2-3/` 下先生成功效数据（见 `Generate-Data.R`），再运行：

```r
source("Figs2-3.R")
```

### 6. Fig4–6（COMPAS 公平性审计）

在 `Fig4-6/` 下：

1. 将 COMPAS 数据放在 `data/compas_data.csv`（可参考 [ProPublica COMPAS](https://github.com/propublica/compas-analysis)）。
2. 生成审计用置信区间数据（如非裔 vs 白人、性别-年龄）：
   ```r
   source("Generate-CI-AfricanVSCaucasian.R")
   # 以及 Generate-Fig4a-Data.R / Generate-Fig4b-Data.R 等（按需）
   ```
3. 画图：
   ```r
   source("Fig4.R")   # 或 Fig5a.R, Fig5b.R, Fig6a.R, Fig6b.R
   ```

图中为各子群的 PPV 差异（disparity）及 90%/95% 置信区间。



## 数据说明

- **COMPAS**：`Fig4-6/data/compas_data.csv` 来自 ProPublica 的 COMPAS 两年再犯数据（或同名格式），用于种族、性别-年龄等子群的公平性审计。
- **模拟数据**：Table 1/2 与 Figs2–3 所用 CSV 由各目录下的 `Generate-*-Data.R` 在本地生成，需在对应目录下依次运行相应脚本。
- **未纳入仓库的大文件**：`Fig4-6/flag.pptx`、`Fig4-6/simhei.ttf` 已加入 `.gitignore`（避免推送超时），若本地需要可自行放入对应目录。

## 核心方法简述

- **EL Lagrange 计算**：在估计方程 \(g = (L - \epsilon) \cdot S\) 下求解 Lagrange 乘子 \(\lambda\)，为 EL 统计量、置信区间与 p 值的计算提供数值基础。本仓库采用 Chen et al. (2008)[^1] 第 432–433 页给出的修正 Newton–Raphson 算法实现 EL 的求解，代码见 `Glambda.R`；亦可使用 R 包 **emplik** 等现有实现。

[^1]: Chen, J., Variyath, A. M., and Abraham, B. (2008). Adjusted empirical likelihood and its properties. *Journal of Computational and Graphical Statistics*, 17(2): 426–443.

- **EEL**：基于经验协方差矩阵的二次型扩展，同样与卡方分布比较。
- **COMPAS 审计**：对“阳性预测”子群构造 PPV 差异 \(\varepsilon_G\) 的 EL 置信区间（见 `epsilonG_CI.R`），用于判断组间公平性。
- **FFR（Tab3）**：在分组损失下做 \(H_0: \varepsilon \leq \varepsilon_0\) 的 EL 检验，并结合 BH 程序控制 FFR。

## 提交代码

在本地修改后，可用以下命令提交并推送到 GitHub：

```bash
cd /path/to/EL-for-fairness-auditing
git add .
git commit -m "类型: 简短描述"   # 如：docs: 更新 README
git push origin main
```

建议 commit 信息前缀：`docs:`（文档）、`feat:`（新功能）、`fix:`（修复）、`chore:`（杂项）。

## 引用与参考

若使用本代码或 COMPAS 数据，请同时引用：

- 对应论文（EL for fairness auditing）
- ProPublica COMPAS 分析： [compas-analysis](https://github.com/propublica/compas-analysis)

## 许可证

见仓库根目录 LICENSE 文件（如有）。COMPAS 数据使用请遵循 ProPublica 相关说明。
