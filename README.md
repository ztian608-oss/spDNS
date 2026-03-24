# spDNS

`spDNS` 是在 `scDNS` 基础上做的结构层扩展版本，目标是支持 **有向二部调控网络** 分析（例如：**Splicing Factor, SF → Splicing Event**）。

> 关键原则：
> - ✅ 保留原有 activity score 与下游统计流程
> - ✅ 只扩展网络结构处理与随机网络生成

---

## 1. 安装

仓库地址：`ztian608-oss/spDNS`

```r
# install.packages("devtools")
devtools::install_github("ztian608-oss/spDNS")
```

---

## 2. 适用场景

原始网络通常是：
- `SF -> target_gene`

但剪接分析需要：
- `SF -> splicing_event`

因此在主流程前，需要先做网络展开（gene 映射到 event）。

---

## 3. 新增函数概览

### 3.1 `expand_network_to_events()`
将 `SF -> gene` 网络按 `gene -> splicing_event` 映射展开为 `SF -> event`。

```r
sf_event_net <- expand_network_to_events(
  sf_gene_network = sf_gene_net,   # columns: SF, target_gene
  gene_event_map  = gene_event_map # columns: gene, splicing_event
)
```

输出网络列为：
- `source`（SF）
- `target`（splicing event）

并自动写入 directed-bipartite 属性，供后续随机网络生成与去重使用。

### 3.2 `randomize_directed_bipartite_network()`
用于有向二部网络的随机化，满足：
- 保持 `source` 出度（out-degree）
- 保持 `target` 入度（in-degree）
- 保持二部方向（只允许 SF->event）

```r
rand_net <- randomize_directed_bipartite_network(sf_event_net, n.edge = nrow(sf_event_net))
```

---

## 4. 与原 scDNS 主流程的关系

`spDNS` 保持原 scDNS 主流程接口：

```r
# 1) 预处理：先扩展网络
sf_event_net <- expand_network_to_events(sf_gene_net, gene_event_map)

# 2) 构建对象（原函数）
obj <- CreatScDNSobject(
  counts = counts_mat,
  data = expr_mat,
  Network = sf_event_net,
  GroupLabel = group_label
)

# 3) 计算网络散度（原函数）
obj <- scDNS_1_CalDivs(obj)

# 4) 建模（原函数；内部已支持 directed-bipartite 随机网络）
obj <- scDNS_2_creatNEAModel_v2(obj)

# 5+) 下游步骤按原流程继续
# obj <- scDNS_3_GeneZscore_v2(obj)
# obj <- scDNS_4_scContribution(obj)
```

---

## 5. 本次扩展“改了什么 / 没改什么”

### 改了（仅结构层）
1. 网络展开（`SF->gene` 到 `SF->event`）
2. 有向二部图去重（避免把 `A->B` 与 `B->A` 合并）
3. 有向二部图度保持随机网络生成

### 没改（算法层）
1. activity score 计算公式
2. 下游统计与差异分析逻辑
3. 输出主结构格式

---

## 6. 输入数据建议

- `sf_gene_network`：至少包含 `SF`, `target_gene`
- `gene_event_map`：至少包含 `gene`, `splicing_event`
- `counts/data`：保持与原 scDNS 使用要求一致
- `Network` 列名建议使用 `source/target`（主流程内部默认前两列为边）

---

## 7. 常见问题

### Q1: 必须先调用 `expand_network_to_events()` 吗？
若你分析的是 SF→event 场景，建议必须先调用；否则主流程看到的仍是 gene-level 网络。

### Q2: 会不会影响原 scDNS 结果可比性？
不会。若仍使用原 gene-gene 网络，流程行为与原版保持一致（结构扩展逻辑仅在 directed-bipartite 场景触发）。

---

## 8. License

沿用项目原有 License 约定。

