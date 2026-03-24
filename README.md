# spDNS

`spDNS` 是在 `scDNS` 基础上扩展的单细胞网络评分框架，用于在**单细胞分辨率**下分析**剪接因子（SF）对可变剪接事件的调控**。

> 核心升级：从传统基因-基因网络分析，扩展到**有向二部网络**（`SF -> splicing_event`）。

---

## 1. 方法概览

### 原始 scDNS 思路
- 输入：基因表达矩阵 + 基因网络
- 输出：基因/细胞层面的网络扰动评分

### spDNS 扩展思路
- 输入：
  1. 基因表达矩阵（cells × genes）
  2. PSI 矩阵（cells × splicing events）
  3. 有向调控网络（`SF -> target_gene`）
  4. 基因到剪接事件映射（`gene -> splicing_event`）
- 新增关键步骤：把 `SF -> gene` 扩展成 `SF -> splicing_event`
- 目标：
  - 计算每个细胞中每个 SF 的调控活性
  - 比较不同条件下 SF 活性差异
  - 识别调控重连（rewiring）

---

## 2. 安装

```r
# install.packages("devtools")
devtools::install_github("xiaolab-xjtu/scDNS")
```

> 当前仓库仍保持与 `scDNS` 的函数兼容，新增 `spDNS` 剪接分析函数已加入 `R/splicingDNS.R`。

---

## 3. 新增函数（spDNS）

### 3.1 网络扩展：`SF -> gene` 到 `SF -> event`

```r
sf_event_network <- expand_network_to_splicing(
  sf_gene_network = sf_gene_network,   # data.frame(SF, target_gene)
  gene_event_map  = gene_event_map     # data.frame(gene, splicing_event)
)
```

### 3.2 计算单细胞 SF 活性评分

```r
sf_activity <- compute_sf_activity_score(
  sf_expression    = sf_expr_mat,      # rows: SF, cols: cells
  psi              = psi_mat,          # rows: events, cols: cells
  sf_event_network = sf_event_network,
  center_psi       = TRUE,
  min_events       = 3,
  na_to_zero       = TRUE
)
```

默认评分公式：

\[
A_{s,c}=z(E_{s,c})\times
\frac{\sum_{e\in T(s)}w_{s,e}\tilde\psi_{e,c}}{\sum_{e\in T(s)}|w_{s,e}|}
\]

- \(E_{s,c}\)：SF 表达
- \(\tilde\psi_{e,c}\)：中心化后的 PSI
- \(w_{s,e}\)：边权重（默认 1）

### 3.3 条件差异分析（SF 活性）

```r
diff_sf <- differential_sf_activity(
  sf_activity = sf_activity,
  condition   = condition_vector,
  contrast    = c("Treatment", "Control"),
  method      = "wilcox"   # or "t.test"
)
```

输出包含：`mean_a`, `mean_b`, `delta`, `p_value`, `FDR`。

### 3.4 调控重连（rewiring）分析

```r
rewiring_res <- detect_sf_rewiring(
  sf_expression    = sf_expr_mat,
  psi              = psi_mat,
  sf_event_network = sf_event_network,
  condition        = setNames(condition_vector, colnames(sf_expr_mat)),
  contrast         = c("Treatment", "Control"),
  cor_method       = "spearman"
)
```

- 对每条 `SF -> event` 边，比较两个条件下 `cor(SF expression, event PSI)` 的差异。

---

## 4. 与原 scDNS 的兼容性

已保留原有流程；同时支持有向网络预处理：

- `preFilterNet(..., ignore_direction = TRUE)`
- `getDoubleDropout(..., ignore_direction = TRUE)`

当你做 `SF -> event` 有向分析时，建议设定：

```r
ignore_direction = FALSE
```

以避免把反向边合并。

---

## 5. 最小可运行示例

```r
# 1) 扩展网络
sf_event_network <- expand_network_to_splicing(sf_gene_network, gene_event_map)

# 2) 计算 SF 活性
sf_activity <- compute_sf_activity_score(
  sf_expression = sf_expr_mat,
  psi = psi_mat,
  sf_event_network = sf_event_network
)

# 3) 差异 SF 活性
diff_sf <- differential_sf_activity(
  sf_activity = sf_activity,
  condition = condition_vector,
  contrast = c("KO", "WT")
)

# 4) 网络重连
rewiring_res <- detect_sf_rewiring(
  sf_expression = sf_expr_mat,
  psi = psi_mat,
  sf_event_network = sf_event_network,
  condition = setNames(condition_vector, colnames(sf_expr_mat)),
  contrast = c("KO", "WT")
)
```

---

## 6. 注意事项

- 请确保表达矩阵与 PSI 矩阵的细胞列名可对齐。
- PSI 稀疏时建议保留 `na_to_zero = TRUE`。
- 若边数量很大，推荐使用稀疏矩阵输入以提升性能。

---

## 7. Citation

如在研究中使用，请引用 `scDNS/spDNS` 相关方法论文与仓库。

