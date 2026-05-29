---
title: "CSAI_math_foundation : Generalized Inverse Matrices(广义逆矩阵)"
date: 2026-05-29
draft: false
tags: ["数学基础", "稀疏编码", "学习笔记"]
categories: ["学习笔记"]
math: true
---
# CSAI_math_foundation : Generalized Inverse Matrices(广义逆矩阵)
******************
## 广义逆矩阵
满足$LA=I$,但不满足$AL=I$的矩阵$L$称为矩阵$A$的左逆矩阵(n>=m possible)，反之有右逆矩阵(n<=m possible).

考查 $n>m$ 且 $A$ 具有满列秩 ($\text{rank } A = m$) 的情况。此时，$m \times m$ 矩阵 $A^T A$ 是可逆的，容易验证 $L = (A^T A)^{-1} A^T$ 满足 $LA = I$，这种左逆矩阵是唯一确定的，常称为**左伪逆矩阵**。

考查 $n<m$ 且 $A$ 具有满行秩 ($\text{rank } A = n$) 的情况。此时，$n \times n$ 矩阵 $AA^T$ 是可逆的，容易验证 $R = A^T (AA^T)^{-1}$ 满足 $AR = I$，这种右逆矩阵是唯一确定的，常称为**右伪逆矩阵**。


考虑一个 $n \times m$ 维的秩亏缺矩阵 $A$，满足 $\text{rank}(A) < \min\{m, n\}$。
$n \times m$ 维秩亏缺矩阵的逆矩阵称为**广义逆矩阵**，它是一个 $m \times n$ 矩阵。令 $A^-$ 表示 $A$ 的广义逆矩阵。

无论 $A^- A = I_{m \times m}$ 还是 $A A^- = I_{n \times n}$ 都不可能成立，因此需要用三个矩阵的乘积来定义秩亏矩阵的逆。

考虑线性矩阵方程 $Ax = y$ 的求解。两边左乘 $A A^-$，则有：
\[
A A^- A x = A A^- y
\]
若 $A^-$ 是 $A$ 的广义逆矩阵，则 $Ax = y \Rightarrow x = A^- y$，代入上式可得：
\[
A A^- A x = A x
\]
我们希望该式对任意非零向量 $x$ 均成立，因此必须满足约束条件：
\[
A A^- A = A
\]

满足 $A A^- A = A$ 的矩阵 $A^-$ 称为 $A$ 的广义逆矩阵。
但实际上$A^-$不是唯一的，有明确的缺陷。因此需多加限制条件得到Moore-Penrose逆矩阵

令 $A$ 是任意 $n \times m$ 矩阵，称 $A^+$ 是 $A$ 的广义逆矩阵，若 $A^+$ 满足以下四个条件（常称 **Moore-Penrose 条件**）：
1.  $AA^+A = A$
2.  $A^+AA^+ = A^+$
3.  $AA^+ = (AA^+)^T$
4.  $A^+A = (A^+A)^T$

与逆矩阵、左伪逆矩阵、右伪逆矩阵一样，Moore-Penrose 逆矩阵也是**唯一的**。

逆矩阵和之前介绍的各种广义逆矩阵都是 Moore-Penrose 逆矩阵的特例。
Moore-Penrose 逆矩阵可以用奇异值分解求出$A=U \Sigma V^T$
维数分别为$n \times n,n \times m,m \times m$,则$A^+=V \Sigma^+ U^T$

$Σ^+是对角阵Σ上元素取倒数然后行列数交换$
********************
## 矩阵方程
$ Ax=b $
在诸多领域中作为矩阵方程，根据数据向量 $b \in \mathbb{R}^n$ 和数据矩阵 $A \in \mathbb{R}^{m \times n}$ 的不同，矩阵主要有以下三种类型：方程数少于未知数的欠定系数矩阵方程，$A$未知的盲矩阵方程，方程数多于未知参数的超定矩阵方程。以下仅讨论超定矩阵和欠定矩阵，（盲矩阵似乎仅存在于复杂系统识别中，难以确定系数的情形中）
*****************

### 超定矩阵方程
考虑超定矩阵方程$Ax=b$,其中$b$为$n \times 1$数据向量，$A$为$n \times m$数据矩阵，且n>m.
由于在实际工业实践中对于$b$中的参数的测量总会受到一定不可抗因素的扰动，即存在观测误差或噪声：$b=b_{0} + e$
加入校正向量$\Delta b$,使得$Ax=b+\Delta b\\ b+\Delta b=b_{0}+e+\Delta b → b_{0}$
目标：$\Delta b→0\\
\min_{x} \|b\|^2 = \|Ax - b\|_2^2 = (Ax - b)^T(Ax - b)\\
\phi = x^T A^T A x - x^T A^T b - b^T A x + b^T b\\
\frac{d\phi}{dx} = 0 = 2A^T A x - 2A^T b$
$A^{T}Ax = A^{T}b$
情形一：超定方程 ($n>m$) 满列秩，即 $\text{rank}(A)=m$

由于 $A^T A$ 非奇异，所以方程有唯一解
\[
x_{\text{LS}} = (A^T A)^{-1} A^T b
\]

情形二：超定方程 ($n>m$) 秩亏缺，即 $\text{rank}(A)<m$
\[
x_{\text{LS}} = (A^T A)^\dagger A^T b
\]

其中 $B^\dagger$ 代表 $B$ 的广义逆矩阵（Moore-Penrose 逆）。
*****************
### 稀疏矩阵和稀疏编码
由线性代数的线性空间的知识可知，信号向量$y \in R^n$最多可分解为n个正交基，此时系数向量必为非稀疏的，若将$y \in R^n$分解为m个n维向量$d_{i} \in R^n,i=1,...,m$的线性组合，(m>n)则m个向量不可能是正交基的集合。这里的$d_{i}$称为原子，码字或基向量。这些原子的集合是过完备的。组成的矩阵$D=[d_{1},...,d_{m}] \in R^{n \times m}$称为字典或码本。
字典的特征：
（1） D的行数小于列数
（2） D具有满行秩
（3） D的列向量模长为1
当系数向量 $x$ 是稀疏向量时，信号分解 $y = Dx$ 称为（信号的）**稀疏分解**。其中，字典矩阵 $D$ 的列常称为解释变量或预测变量；向量 $y$ 称为响应变量或目标信号；而 $x$ 则可视为目标信号 $y$ 相对于字典 $D$ 的一种表示。

**稀疏编码问题**：给定一个 $n$ 维实值输入向量 $y \in \mathbb{R}^n$，确定 $m$ 个 $n$ 维基向量 $d_1, \dots, d_m \in \mathbb{R}^n$ 以及一个稀疏的 $m$ 维权向量（系数向量）$s \in \mathbb{R}^m$，使得部分基向量的加权线性组合可以充分逼近输入向量，即
\[
y \approx Ds, \quad \text{其中 } D = [d_1, \dots, d_m] \in \mathbb{R}^{n \times m}
\]

如果给定的是 $n$ 维实值输入向量的一组集合 $\{y_1, y_2, \dots, y_k\}$，则稀疏编码的目的就是确定基矩阵 $D$ 和系数矩阵 $S$，使得
\[
Y \approx DS
\]
其中系数矩阵的每一个列向量都是稀疏向量。

稀疏编码的主要特点是：系数向量只有少数元素不等于零，大多数元素为零。

稀疏编码可能是临界完备的或过完备的。**临界完备**是指基向量的个数 $m$ 等于输入向量的维数 $n$。
**稀疏矩阵方程的求解：**
贪婪求解算法
(1) 匹配追踪算法(MP算法)：迭代地构造一个稀疏解
从字典矩阵D（也称为过完备原子库中），选择一个与信号 y最匹配的原子(也就是某列)，构建一个稀疏逼近，并求出信号残
差，然后继续选择与信号残差最匹配的原子，反复迭代，信号y可以由这些原子来线性和，再加上最后的残差值来表示。
与信号y 最匹配的原子就是与信号y的相关系数的绝对值最大的原子。一般地，预先将信号y零均值化，将字典中的原子（即D中的
列向量）也零均值化，再将原子归一化后再进行求解。此时，与信号y 最匹配的原子就是与y的内积的绝对值最大的原子。
优点：计算速度快、收敛性差


### MP Python 实现
```python
import numpy as np

def mp(y, D, K):
    n, m = D.shape
    x = np.zeros(m)
    r = y.copy()
    for _ in range(K):
        corr = D.T @ r
        idx = np.argmax(np.abs(corr))
        x[idx] += corr[idx]
        r -= corr[idx] * D[:, idx]
    return x
```

(2) 正交匹配追踪算法(OMP算法)
在分解的每一步对残差与所选择的全部原子进行正交化处理,这使得在精度要求相同的情况下，OMP算法的收敛速度更快。
优点：计算速度快、收敛性好
一个原子不会被选择两次，因此结果会在有限的几步收敛。
```
def omp(y, D, K):
    n, m = D.shape
    x = np.zeros(m)
    r = y.copy()
    Lambda = []

    for _ in range(K):
        corr = D.T @ r
        idx = np.argmax(np.abs(corr))
        Lambda.append(idx)
        D_lam = D[:, Lambda]
        # 最小二乘更新
        x_lam = np.linalg.pinv(D_lam) @ y
        r = y - D_lam @ x_lam

    x[Lambda] = x_lam
    return x
```
字典训练算法K-SVD
K-SVD算法主要分为两个步骤，稀疏表示与字典更新：
稀疏表示：对于目标信号集$Y \in R^{n\times N}$，初始化字典$D\in R^{n\times m}$，采用各类稀疏编码的求解算法来求解出系数矩阵$X\in R^{m \times N}$；
输入：信号矩阵 Y，初始字典 D，稀疏度 K，迭代次数 J
输出：训练字典 D，稀疏系数 X

1. 初始化系数矩阵 X
2. for j = 1 to J:
    1. 稀疏编码：对每个信号 y_i，用 OMP 更新 x_i
    2. 逐列更新字典原子 D[:,k]
        - 找到使用该原子的样本集合 I_k
        - 计算去除当前原子的残差矩阵 E_k
        - 对 E_k 做 SVD 分解
        - 用主左奇异向量更新原子，主奇异值与右奇异向量更新系数
3. 返回 D, X
```
def omp_batch(Y, D, K):
    _, N = Y.shape
    m = D.shape[1]
    X = np.zeros((m, N))
    for i in range(N):
        X[:, i] = omp(Y[:, i], D, K)
    return X

def ksvd(Y, D, K, max_iter=10):
    n, N = Y.shape
    m = D.shape[1]
    X = np.zeros((m, N))

    for _ in range(max_iter):
        # 稀疏编码步骤
        X = omp_batch(Y, D, K)

        # 字典更新步骤（逐原子更新）
        for k in range(m):
            used_idx = np.nonzero(X[k, :])[0]
            if len(used_idx) == 0:
                continue
            # 去除第k个原子的贡献
            D_temp = D.copy()
            D_temp[:, k] = 0
            Ek = Y - D_temp @ X
            Ek_reduced = Ek[:, used_idx]

            # SVD 更新原子与系数
            U, s, Vt = np.linalg.svd(Ek_reduced, full_matrices=False)
            D[:, k] = U[:, 0]
            X[k, used_idx] = s[0] * Vt[0, :]
    return D, X
```

