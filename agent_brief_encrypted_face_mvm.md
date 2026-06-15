# Agent Brief: 全加密 1:N 人脸 Embedding Mat-Vec 路线与 Rhombus 改造要点

> 目的：供后续代码 Agent / 研究 Agent 快速查阅。  
> 背景：当前项目目标是 **加密人脸特征库 + 加密用户 Embedding 查询 + 第三方密文计算 + CKKS 阈值匹配 + indicator/tag 提取**。  
> 核心问题：1:N 场景下 `Enc(W) × Enc(q)` 的密文矩阵-向量乘法是主要瓶颈。

---

## 1. 当前业务目标

数据库：

```latex
W \in \mathbb{R}^{K \times d}
```

查询：

```latex
q \in \mathbb{R}^{d}
```

其中：

```text
K = 人脸库规模，例如 10^4 / 10^5 / 10^6
d = 人脸 embedding 维度，例如 128 或 512
```

目标：

```latex
s = Wq
```

其中：

```latex
s_i = \langle w_i, q \rangle
```

如果数据库模板和查询 embedding 都经过 L2 normalization，则：

```latex
s_i = \cos(w_i, q)
```

后续阈值匹配：

```latex
m_i = 1[s_i \ge 0.8]
```

再执行标签提取：

```latex
tag = \sum_i m_i \cdot tag_i
```

当前安全需求是：

```text
1. 用户 query embedding 加密；
2. 数据库人脸 embedding 也加密；
3. 第三方计算方不能看到 query 或 database；
4. 最终最好只释放 count/tag 或 yes/no，不释放完整 score 向量；
5. 后续希望接 CKKS threshold / BLEACH / indicator × tag。
```

---

## 2. 三种候选路线定位

### 2.1 HyDia-style diagonal ct-ct MVM

适用对象：

```text
Enc(DB) × Enc(query)
1:N face matching
K >> d
后续继续 CKKS threshold
```

核心结构：

```latex
C_g =
RelinRescale
\left(
\sum_{h=0}^{d-1}
DB_{g,h} \cdot Rot(ct_q,h)
\right)
```

其中：

```text
DB_{g,h} = 第 g 个 group 的第 h 条 diagonal database ciphertext
ct_q = 重复打包的 encrypted query
```

复杂度：

```latex
G = \lceil K / S \rceil
```

其中：

```text
S = CKKS slot 数
```

在线复杂度：

```text
ct-ct multiplication: G · d
relinearize/rescale: G
query rotations: d - 1
output score ciphertexts: G
```

关键优势：

```text
1. 数据库和 query 都可以加密；
2. score 直接落在 CKKS slots；
3. 后续可以直接执行 threshold polynomial；
4. d 个乘法先累加，再每个 group 只做一次 relin/rescale；
5. 最贴合当前项目主路线。
```

关键问题：

```text
1. ct-ct multiplication 数量仍然是 G · d；
2. threshold polynomial / cleaning 仍可能成为后续瓶颈；
3. encrypted DB 存储量大；
4. 如果输出 ID/tag，不能随意跨用户 folding。
```

---

### 2.2 IRIS/Park CCMM / large batched MatMul

适用对象：

```text
Enc(W) × Enc(Q)
large encrypted matrix-matrix multiplication
B = query batch size 足够大
```

Park CCMM 的核心思想：

```text
CC-MM → 4 × PP-MM + 3 × C-MT + relin/rescale
```

其中：

```text
PP-MM = plaintext/plain modular matrix multiplication
C-MT = ciphertext matrix transpose
```

原文 profile 结论：

```text
主要瓶颈是 PP-MM；
C-MT 是第二瓶颈；
Relin / Rescale 占比较小。
```

适合：

```latex
W_{K \times d} Q_{d \times B}
```

当：

```text
B 较大，或者像 IRIS 一样存在 batch eyes × rotations
```

不适合：

```text
单 query skinny MatVec: K × d by d × 1
```

原因：

```text
1. 单 query padding waste 非常大；
2. Park CCMM 更偏向大矩阵 block multiplication；
3. 如果只做 one-user authentication，batch 等待会增加首包延迟；
4. 如果需要 ID/tag，IRIS folding 不能直接使用，因为会丢位置。
```

---

### 2.3 Rhombus coefficient MVM

Rhombus 原始场景：

```text
Alice holds vector v
Bob holds plaintext matrix W
Bob computes W · Enc(v)
result converted to additive secret sharing
```

原始 Rhombus 不是 encrypted DB 方案。

原始复杂度：

```latex
m = rc/N
```

```text
HomMulPt: O(rc/N)
Rotation: O(sqrt(rc/N))
Communication: O(1)
```

其中：

```text
r = matrix rows
c = matrix cols
N = RLWE polynomial degree
```

Rhombus 关键优化：

```text
1. coefficient encoding；
2. input-output packing；
3. PackRLWEs；
4. split-point picking / SPP；
5. rotation reuse on common Enc(v)。
```

原始高效性的前提：

```latex
ct_i = \hat{y}_i \cdot Enc(v)
```

其中：

```text
\hat{y}_i 是明文矩阵预处理块。
```

如果数据库也加密，则变成：

```latex
Enc(y_i) \cdot Enc(v)
```

这会破坏 Rhombus 的 plaintext/ciphertext 分离优化。

---

## 3. 如果要把 Rhombus 改成全加密 MVM，需要什么？

### 3.1 最小可行改造思路

不要直接加密原始矩阵 `W`，而是：

```text
数据库方先运行 Rhombus 的矩阵预处理：
W
→ InputPacking
→ matrix_row_combine / SPP
→ 得到 z_{i1,tau}

然后加密：
z_{i1,tau} → Enc(z_{i1,tau})
```

原始在线计算：

```latex
ct_{i_1}
=
\sum_{\tau}
HomMulPt(ct_\tau, z_{i_1,\tau})
```

改造后：

```latex
ct_{i_1}
=
\sum_{\tau}
HomMult(ct_\tau, Enc(z_{i_1,\tau}))
```

其中：

```latex
ct_\tau = HomAut(Enc(v), \tau)
```

---

### 3.2 需要新增的数据结构

建议新增：

```cpp
struct EncryptedRhombusMatrix {
    size_t rows;
    size_t cols;
    size_t padded_rows;
    size_t padded_cols;

    size_t ell;
    size_t h;
    size_t u;

    size_t tau_count;  // 2^{u-h}
    size_t i1_count;   // usually r / 2^u or related packed count

    std::vector<Ciphertext> enc_z; // Enc(z_{i1,tau})
};
```

索引：

```cpp
enc_z[i1 * tau_count + tau_id]
```

对应：

```latex
Enc(z_{i_1,\tau})
```

---

### 3.3 需要新增的 API

数据库预处理/加密侧：

```cpp
EncodePreprocessedMatrix(...)
EncryptPreprocessedMatrix(...)
SaveEncryptedMatrix(...)
LoadEncryptedMatrix(...)
```

在线计算侧：

```cpp
EncryptedMatVecMul(...)
EncryptedHomInnerProd(...)
```

验证侧：

```cpp
DecryptCoeffMatVecResult(...)
CompareWithPlainMatVec(...)
```

如果要接生物识别后处理，还需要：

```cpp
CoeffToSlotResult(...)
HomomorphicThreshold(...)
CleanIndicator(...)
IndicatorTagSelect(...)
```

---

### 3.4 在线 hom_inner_prod 的关键改造

原始代码大致为：

```cpp
multiply_plain_ntt(cached_ct[j0], encoded_mat[idx], tmp);
add_inplace(acc, tmp);
```

全加密改成：

```cpp
multiply(cached_ct[j0], enc_z[idx], tmp);
add_inplace(acc_deg2, tmp);
```

但不要每次乘法后立刻 relinearize。

推荐模式：

```cpp
for (i1 = 0; i1 < i1_count; i1++) {
    Ciphertext acc_deg2;
    bool init = false;

    for (tau = 0; tau < tau_count; tau++) {
        idx = i1 * tau_count + tau;

        Ciphertext tmp;
        evaluator.multiply(cached_ct[tau], enc_matrix.enc_z[idx], tmp);

        if (!init) {
            acc_deg2 = tmp;
            init = true;
        } else {
            evaluator.add_inplace(acc_deg2, tmp);
        }
    }

    evaluator.relinearize_inplace(acc_deg2, relin_keys);
    evaluator.rescale_to_next_inplace(acc_deg2);

    matvec[i1] = acc_deg2;
}
```

这样 relinearization 次数从：

```latex
i1\_count \cdot tau\_count
```

降为：

```latex
i1\_count
```

这与 HyDia 的“先累加 degree-2 ciphertext，再统一 relin/rescale”思想一致。

---

## 4. 全加密 Rhombus 改造的关键新增要求

### 4.1 必须加入 relinearization key

因为从：

```text
ciphertext × plaintext
```

变成：

```text
ciphertext × ciphertext
```

乘法后 ciphertext degree 变成 2，需要：

```text
Relinearization key: sk^2 -> sk
```

---

### 4.2 必须管理 scale / level

CKKS ct-ct multiplication 后：

```latex
scale = \Delta^2
```

需要 rescale 回：

```latex
\Delta
```

所有待 PackRLWEs 的 `matvec[i1]` 必须满足：

```text
1. same modulus level
2. same scale
3. size-2 ciphertext
4. same NTT/INTT state
```

---

### 4.3 PackRLWEs 可复用，但调用前要规范化

保留：

```cpp
PackRLWEs(matvec, u, galois_keys, result)
```

调用前必须保证：

```cpp
mod_switch_to_same_level(matvec);
align_scale(matvec);
ensure_ciphertext_size_2(matvec);
ensure_expected_ntt_state(matvec);
```

---

### 4.4 ConvToSS 是否保留取决于目标

如果做 HE-MPC hybrid：

```text
保留 ConvToSS。
```

如果做生物识别 CKKS threshold：

```text
不要立刻 ConvToSS；
输出应保持 ciphertext，用于后续 threshold。
```

---

## 5. Rhombus 全加密版本的最大问题：输出 layout 不匹配

即使完成全加密 Rhombus，它的输出仍是 coefficient encoding。

也就是 score 落在 polynomial coefficients 上：

```text
coefficient domain: score_i 位于系数 index
```

但 CKKS threshold polynomial 是 slot-wise：

```latex
P(s_i - 0.8)
```

它要求 score 位于 slots：

```text
slot domain: score_i 位于 CKKS SIMD slot
```

因此如果继续用 Rhombus coefficient route，后面还要：

```text
CoeffToSlot / ring packing / bootstrapping layout conversion
```

或者改成：

```text
HE coefficient result → secret sharing → MPC comparison
```

这就是 Rhombus 不适合作为当前主线的核心原因。

---

## 6. 推荐实现路线

### 6.1 实验 baseline：Rhombus-coeff 全加密 MVM

目的：

```text
验证 coefficient encoding 下 ct-ct MVM 是否可跑；
作为对比 baseline；
分析 ct-ct mult + PackRLWEs 成本。
```

最小实现：

```text
1. EncryptPreprocessedMatrix
2. EncryptedHomInnerProd
3. delayed relinearize/rescale
4. PackRLWEs
5. decrypt coefficient result 验证 Wv
```

不强求：

```text
CoeffToSlot
threshold
indicator tag
```

---

### 6.2 主线：HyDia-style slot-domain diagonal ct-ct MVM

目的：

```text
服务于人脸 1:N 生物识别端到端路线。
```

核心：

```latex
C_g =
RelinRescale
\left(
\sum_{h=0}^{d-1}
DB_{g,h} \cdot Rot(ct_q,h)
\right)
```

优点：

```text
1. Enc(DB) and Enc(query) both supported;
2. output scores are directly in CKKS slots;
3. threshold polynomial can be applied directly;
4. indicator × tag route is natural;
5. compatible with ThFHE / threshold decryption architecture.
```

后续：

```text
1. threshold polynomial for s_i >= 0.8;
2. cleaning / BLEACH;
3. count = sum indicator;
4. tag_limb = sum indicator_i * tag_i_limb;
5. threshold decryption or authorized receiver decryption.
```

---

## 7. 三种方法的性能瓶颈总结

### 7.1 HyDia

原文 profile 结论：

```text
传统 bottleneck: rotate-sum + per-multiply relinearization
HyDia 后瓶颈: 大量 ct-ct multiplication + threshold polynomial
通信显著降低，但 encrypted DB 存储大
```

迁移判断：

```text
当前项目最贴近 HyDia；
应优先 profile:
G·d ct-ct mult
G relin/rescale
d-1 query rotations
threshold polynomial / cleaning
```

---

### 7.2 Park / IRIS CCMM

原文 profile 结论：

```text
CC-MM → 4 × PP-MM + 3 × C-MT + relin/rescale
主瓶颈是 PP-MM
C-MT 第二
Relin/rescale 占比较小
```

适合：

```text
large batched MatMat:
W_{K×d} Q_{d×B}
B sufficiently large
```

不适合：

```text
single-query skinny MatVec:
W_{K×d} q_{d×1}
```

---

### 7.3 Rhombus

原文 profile 结论：

```text
解决 PC-MVM 的 rotation / communication / Encrypt&Decrypt overhead
适合 plaintext matrix × encrypted vector
适合 HE-to-SS PPML linear layer
```

不适合直接作为：

```text
encrypted DB × encrypted query + CKKS threshold
```

原因：

```text
1. matrix is not plaintext;
2. SPP plaintext/ciphertext separation失效；
3. output is coefficient layout;
4. threshold route requires slot layout。
```

---

## 8. 安全模型提醒

全加密 DB 和 query 必须在同一密钥体系下。

### 模型 A：用户持 sk

```text
用户生成 pk/sk
数据库方用 pk 加密 DB
第三方计算
用户解密最终结果
```

问题：

```text
如果第三方把 Enc(DB) 给用户，用户可以解密 DB。
依赖 no-collusion。
```

HyDia 类似这种模型。

---

### 模型 B：数据库方持 sk

```text
数据库方生成 pk/sk
用户用 pk 加密 query
第三方计算
数据库方/授权方解密
```

问题：

```text
数据库方若拿到 Enc(query)，可解密用户 query。
依赖流程隔离。
```

---

### 模型 C：Threshold CKKS / ThFHE

```text
多方 DKG 生成 pk
数据库和 query 都用 pk 加密
没有任何单方持有完整 sk
最终由 threshold decryptors 解密
```

这是长期推荐的强安全模型。

---

## 9. Agent TODO 清单

### 9.1 若实现 Rhombus 全加密 baseline

1. 阅读并定位：
   ```text
   src/matvec.h
   src/matvec.cpp
   PackRLWEs
   hom_inner_prod
   matrix_row_combine64/128
   hom_aut_galois_group
   ```

2. 新增：
   ```text
   EncryptedRhombusMatrix
   EncryptPreprocessedMatrix
   EncryptedHomInnerProd
   EncryptedMatVecMul
   ```

3. 修改在线核：
   ```text
   multiply_plain_ntt → multiply
   accumulate degree-2 ciphertexts
   relinearize once per i1
   rescale once per i1
   level/scale alignment
   ```

4. 保留：
   ```text
   PackRLWEs
   decrypt coefficient result verification
   ```

5. 暂不实现：
   ```text
   CKKS threshold
   indicator tag
   unless CoeffToSlot is added
   ```

---

### 9.2 若实现主线 HyDia-style score engine

1. 实现 database diagonal packing：
   ```text
   split W into d×d blocks
   pack same diagonal from multiple blocks into one ciphertext
   ```

2. 实现 query repeated packing：
   ```text
   ct_q = Enc(q | q | ... | q)
   ```

3. 在线：
   ```text
   precompute/hoist Rot(ct_q,h), h=0..d-1
   for each group g:
       acc_deg2 = sum_h DB[g,h] * Rot(ct_q,h)
       relin/rescale once
       output Enc(scores_g)
   ```

4. 后处理：
   ```text
   threshold polynomial
   cleaning
   count
   tag_limb selection
   ```

---

## 10. 最终判断

如果问题是：

```text
“Rhombus 怎么改成全加密 MVM？”
```

答案是：

```text
加密 SPP 预处理后的矩阵块 z_{i1,tau}，
把 HomMulPt 改成 ct-ct HomMult，
并延迟 relinearize/rescale 到每个 i1 累加完成后。
```

如果问题是：

```text
“它是否适合作为当前加密人脸库 + CKKS threshold 的主线？”
```

答案是：

```text
不推荐作为主线。
```

原因：

```text
1. Rhombus 原始优势依赖明文矩阵；
2. 全加密后会引入大量 ct-ct multiplication；
3. 输出仍是 coefficient layout；
4. 后续 CKKS threshold 需要 slot layout；
5. HyDia-style diagonal ct-ct MVM 更贴合当前目标。
```

推荐：

```text
Rhombus-coeff 全加密版 = baseline / ablation
HyDia-style diagonal MVM = 主线 score engine
Park/IRIS CCMM = 大 batch / membership-only / GPU-BLAS 路线
```
