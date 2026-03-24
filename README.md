# spDNS

`spDNS` 是从 `scDNS` 扩展而来的单细胞网络分析 R 包，面向**有向二部调控网络**场景，典型形式为：

- **Splicing Factor (SF) → Splicing Event**

当前版本重点支持：
1. 将 `SF -> gene` 网络扩展为 `SF -> event` 网络；
2. 在不改动原有 activity score / 下游统计逻辑的前提下，支持有向二部网络结构处理；
3. 提供**有向且度保持**（out-degree / in-degree）且**二部结构安全**的随机网络生成。

---

## 安装

> 仓库地址已变更为：`ztian608-oss/spDNS`

```r
# install.packages("devtools")
devtools::install_github("ztian608-oss/spDNS")
```

---

## 输入数据

### 1) 原始调控网络
- `SF -> target_gene`

### 2) 基因到事件映射
- `gene -> splicing_event`（一对多）

### 3) scDNS 主流程所需表达矩阵
- 保持与原始 scDNS 一致（不改变原有评分和统计接口）

---

## 新增：网络扩展预处理

### `expand_network_to_events()`

该函数用于在进入 scDNS 主流程前，先把网络从基因层扩展到事件层。

```r
sf_event_net <- expand_network_to_events(
  sf_gene_network = sf_gene_net,   # data.frame(SF, target_gene)
  gene_event_map  = gene_event_map # data.frame(gene, splicing_event)
)

head(sf_event_net)
#   source   target
# 1 SRSF1   SE:chr1:...
# 2 SRSF1   RI:chr3:...
```

函数会在网络上写入属性以标识该网络是 directed-bipartite，用于后续结构处理与随机网络生成。

---

## 有向二部网络支持（结构层）

### 1) 去重与过滤
`preFilterNet()` 在检测到 directed-bipartite 网络后，会按**有向边**去重（不再把 `A->B` 与 `B->A` 视为同一条边）。

### 2) 随机网络生成（关键）
新增函数：
- `randomize_directed_bipartite_network()`

该函数满足：
- 保持 SF 节点出度（out-degree）
- 保持 event 节点入度（in-degree）
- 仅生成 `SF -> event` 边（不产生 event->SF / SF->SF）

这会在 `NEAModel` 与 `NEAModel_v2` 路径中自动启用（当输入网络标记为 directed-bipartite 时）。

---

## 典型流程示例

```r
library(spDNS)

# Step 1: 网络扩展（SF->gene => SF->event）
sf_event_net <- expand_network_to_events(sf_gene_net, gene_event_map)

# Step 2: 构建 scDNS 对象（沿用原流程接口）
obj <- CreatScDNSobject(
  counts = counts_mat,
  data = expr_mat,
  Network = sf_event_net,
  GroupLabel = group_label
)

# Step 3: 计算网络散度（原函数）
obj <- scDNS_1_CalDivs(obj)

# Step 4: 构建 NEA 模型（原函数，内部已支持 directed-bipartite 随机网络）
obj <- scDNS_2_creatNEAModel_v2(obj)

# 后续打分、统计流程不变
# obj <- scDNS_3_GeneZscore_v2(obj)
# obj <- scDNS_4_scContribution(obj)
```

---

## 兼容性说明

- **未修改** activity score 计算逻辑；
- **未修改** 下游统计分析逻辑；
- 修改仅限：
  1) 网络结构处理；
  2) 随机网络生成。

因此原有 scDNS 使用方式基本兼容，新增能力可按需启用。

---

## License

沿用项目原有 License 约定。

