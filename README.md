# 《Applied Stochastic Analysis》数学核心骨干

> Weinan E（鄂维南）、Tiejun Li（李铁军）、Eric Vanden-Eijnden  
> *Applied Stochastic Analysis*, Graduate Studies in Mathematics, Vol. 199, AMS
>
> 提炼口径：只保留**定义—定理—公式—方法**的主干链条，以及章与章之间的逻辑接口；习题、证明细节、文献注记、数值实验从略。目标是一万汉字以内能通读一遍、 thereafter 当手册查。

---

## 目录

- [0. 一句话与三条主线](#0-一句话与三条主线)
- [1. 概率空间、条件期望与收敛（第 1 章）](#1-概率空间条件期望与收敛第-1-章)
- [2. 极限定理与大偏差（第 2 章）](#2-极限定理与大偏差第-2-章)
- [3. Markov 链（第 3 章）](#3-markov-链第-3-章)
- [4. Monte Carlo 方法（第 4 章）](#4-monte-carlo-方法第-4-章)
- [5. 随机过程的一般框架（第 5 章）](#5-随机过程的一般框架第-5-章)
- [6. Wiener 过程（第 6 章）](#6-wiener-过程第-6-章)
- [7. Itô 随机微分方程（第 7 章）](#7-itô-随机微分方程第-7-章)
- [8. Fokker–Planck 方程（第 8 章）](#8-fokkerplanck-方程第-8-章)
- [9. 路径积分（第 9 章）](#9-路径积分第-9-章)
- [10. 随机场（第 10 章）](#10-随机场第-10-章)
- [11. 统计力学（第 11 章）](#11-统计力学第-11-章)
- [12. 稀有事件（第 12 章）](#12-稀有事件第-12-章)
- [13. 化学反应动力学（第 13 章）](#13-化学反应动力学第-13-章)
- [14. 附录：六个工具箱（A–F）](#14-附录六个工具箱af)
- [15. 收束：十个必背公式、依赖图与读法](#15-收束十个必背公式依赖图与读法)

---

## 0. 一句话与三条主线

这本书可以用一句话概括：**把随机现象翻译成方程，再用分析的工具把方程解出来。**

翻译的方向有两边。一边是"概率 → 方程"：用 Markov 性把过程的分布压进转移核，用 Itô 公式把轨道压进生成元，于是得到 Fokker–Planck 方程、后向方程、化学主方程。另一边是"方程 → 概率"：用 Feynman–Kac 公式把 PDE 的解写成路径平均，用大偏差把罕见事件的概率写成变分问题。全书就在这两边来回走。

三条主线把它们串起来：

1. **极限链**：大数律 → 中心极限定理 → 大偏差 → 极值理论。回答"微观随机性在宏观尺度上如何变成确定性，又如何偶尔失效"。
2. **动力学链**：Markov 链 → 扩散过程 → SDE → Fokker–Planck → 路径积分。回答"随机系统如何随时间演化，以及演化算子的谱是什么"。
3. **多尺度链**：平均/齐次化 → 统计力学 → 稀有事件 → 化学反应动力学。回答"快变量如何被消去，慢变量如何获得封闭的有效方程"。

贯穿全书的五个核心对象，值得先记住：

| 对象 | 符号 | 作用 |
|---|---|---|
| 生成元 | $L$ | 描述过程的一阶时间导数：$Lf=\lim_{t\downarrow0}\frac{\mathbb E_x f(X_t)-f(x)}{t}$ |
| Fokker–Planck 算子 | $L^*$ | $L$ 的对偶，驱动概率密度演化 |
| 半群 | $T_t=e^{tL}$ | 把初值演化到 $t$ 时刻 |
| 作用量 | $I[\phi]$ | 路径的"代价"，大偏差与稀有事件的货币 |
| 速率函数 | $I(x)$ | 大偏差指数，与熵、自由能同构 |

### 章间接口（谁给谁供料）

| 上游 | 下游 | 供的是什么 |
|---|---|---|
| 第 1–2 章 | 全书 | 收敛语言、特征函数、大偏差 |
| 第 3 章 | 第 4、5、13 章 | 转移矩阵、$Q$-过程、细致平衡 |
| 第 5–6 章 | 第 7–9 章 | 半群、生成元、Wiener 测度、KL 展开 |
| 第 7 章 | 第 8、11、12 章 | Itô 公式、SDE、Langevin |
| 第 8 章 | 第 11–13 章 | FP 方程、不变分布、多尺度平均 |
| 第 9 章 | 第 12 章 | 路径积分、Girsanov、作用量 |
| 第 10 章 | 第 11 章 | Gibbs 分布、Markov 随机场 |
| 第 11 章 | 第 12、13 章 | 涨落耗散、线性响应、Mori–Zwanzig |
| 第 12 章 | 第 13 章 | 大偏差、Kramers 速率、准势 |

这张表的用法：读到某一章卡住时，先回到它的上游章找工具；写论文时，先看下游章有没有现成的模板。

---

## 1. 概率空间、条件期望与收敛（第 1 章）

**概率空间**是三元组 $(\Omega,\mathcal F,\mathbb P)$：样本空间、$\sigma$-代数、概率测度。$\sigma$-代数的作用不是装饰，而是规定"哪些事件可以被赋值概率"；由此才有可测函数（随机变量）、可测过程、停时这些概念。

**条件期望**是全书的第一个技术支点。给定子 $\sigma$-代数 $\mathcal G\subset\mathcal F$，$X\in L^1$ 的条件期望 $\mathbb E[X\mid\mathcal G]$ 由两条性质唯一确定：

$$\mathbb E[X\mid\mathcal G]\ \text{是}\ \mathcal G\text{-可测},\qquad
\int_A \mathbb E[X\mid\mathcal G]\,d\mathbb P=\int_A X\,d\mathbb P,\quad \forall A\in\mathcal G.$$

存在性靠 Radon–Nikodym 定理：测度 $A\mapsto\int_A X\,d\mathbb P$ 关于 $\mathbb P|_{\mathcal G}$ 绝对连续，其 RN 导数就是条件期望。若 $X\in L^2$，它还有一个更直观的身份——**$L^2(\Omega,\mathcal G,\mathbb P)$ 上的正交投影**：

$$\mathbb E[X\mid\mathcal G]=\arg\min_{Y\in L^2(\mathcal G)}\mathbb E|X-Y|^2.$$

三条使用频率最高的性质：线性性、塔性质 $\mathbb E[\mathbb E[X\mid\mathcal G_2]\mid\mathcal G_1]=\mathbb E[X\mid\mathcal G_1]$（$\mathcal G_1\subset\mathcal G_2$）、"已知的可提出" $\mathbb E[YX\mid\mathcal G]=Y\mathbb E[X\mid\mathcal G]$（$Y$ 为 $\mathcal G$-可测）。这三条几乎支撑了全书所有鞅论证和 Markov 论证。

**四种收敛**必须能默画关系图：

$$X_n\xrightarrow{\text{a.s.}}X\ \Longrightarrow\ X_n\xrightarrow{\mathbb P}X\ \Longrightarrow\ X_n\xrightarrow{d}X,
\qquad X_n\xrightarrow{L^p}X\ \Longrightarrow\ X_n\xrightarrow{\mathbb P}X.$$

反向一般不成立；唯一的例外是"依分布收敛 + 极限为常数 ⇒ 依概率收敛"。依分布收敛等价于对一切有界连续函数 $f$ 有 $\mathbb E f(X_n)\to\mathbb Ef(X)$，也等价于特征函数逐点收敛（Lévy 连续性定理）。

**特征函数** $\varphi_X(t)=\mathbb E[e^{itX}]$ 是分布的傅里叶变换。它把独立和的卷积变成乘积，因此是证明 CLT 的标准工具。Lévy 连续性定理说：若 $\varphi_n(t)\to\varphi(t)$ 对每个 $t$ 成立，且 $\varphi$ 在 $0$ 连续，则 $\varphi$ 必是某个随机变量的特征函数，且 $X_n\Rightarrow X$。Bochner 定理给出反向判据：连续、正定、$\varphi(0)=1$ 的函数必是特征函数。

**两个组合工具**。Borel–Cantelli 引理：若 $\sum_n\mathbb P(A_n)<\infty$，则 $\mathbb P(A_n\ \text{i.o.})=0$；若 $A_n$ 独立且 $\sum_n\mathbb P(A_n)=\infty$，则 $\mathbb P(A_n\ \text{i.o.})=1$。Wick 定理：对中心高斯向量，

$$\mathbb E[X_{i_1}\cdots X_{i_{2n}}]=\sum_{\text{配对}}\prod_{(a,b)\in\text{配对}}\mathbb E[X_{i_a}X_{i_b}],$$

奇数阶矩为零。这条定理在 Wiener 混沌展开和量子场论的费曼图里会反复出现。

---

## 2. 极限定理与大偏差（第 2 章）

**大数律**回答"平均是否稳定"。独立同分布且 $\mathbb E|X|<\infty$ 时，弱大数律给 $\bar X_n\xrightarrow{\mathbb P}\mu$，强大数律给 $\bar X_n\xrightarrow{\text{a.s.}}\mu$。证明路线：弱律用 Chebyshev 或特征函数；强律用 Borel–Cantelli 加截断。

**中心极限定理**回答"涨落有多大"。独立同分布、均值 $\mu$、方差 $\sigma^2$ 时，

$$\sqrt n\,\frac{\bar X_n-\mu}{\sigma}\ \Longrightarrow\ N(0,1).$$

更一般的形式是 Lindeberg–Feller 定理：三角阵列独立、方差和满足 $\sum_k\sigma_{nk}^2\to\sigma^2$，且对每个 $\epsilon>0$ 有 Lindeberg 条件

$$\sum_k\mathbb E\big[|\xi_{nk}|^2\mathbf 1_{|\xi_{nk}|>\epsilon}\big]\to0,$$

则和依分布收敛到 $N(0,\sigma^2)$。Lindeberg 条件的含义是"没有单个项主导总和"。

**Edgeworth 展开**把 CLT 推进一阶：用累积量 $\kappa_3,\kappa_4,\dots$ 修正正态近似，得到 $F_n(x)=\Phi(x)+\frac{\kappa_3}{6\sqrt n}(1-x^2)\varphi(x)+O(1/n)$ 一类的展开。它是"高阶精度"的源头。

**稳定律**回答"不服从 CLT 时极限是什么"。若独立同分布的尾部分布满足 $\mathbb P(|X|>x)\sim x^{-\alpha}L(x)$，$\alpha\in(0,2)$，则适当中心化、用 $n^{1/\alpha}$ 归一化后，极限是 $\alpha$-稳定分布。$\alpha=2$ 时退回正态；$\alpha<2$ 时方差无穷，极限分布有厚尾。$\alpha$ 就是"吸引域"的标签。

CLT 失效有三种典型方式，对应三种不同的极限语言。第一种是**尾太厚**（$\alpha<2$），方差不存在，极限是稳定律；第二种是**相关性太强**（长记忆），方差和发散，极限可能是分数 Brownian 运动；第三种是**极大值主导**（少数超大项），极限是 Poisson 型或极值型。识别是哪一种，决定了你该用 CLT、稳定律还是大偏差。

Edgeworth 展开的实用价值在于"误差条"：当 $n$ 不够大时，正态近似的偏差由偏度 $\kappa_3/\sigma^3$ 主导；对称分布的三阶项消失，四阶项（峰度）接管。这解释了为什么重尾数据的小样本 $t$ 检验不可靠，也解释了为什么金融收益的"正态近似"总在尾部失效。

**Cramér 大偏差**回答"平均偏离均值有多难"。设矩母函数有限：

$$\Lambda(\lambda)=\log\mathbb E[e^{\lambda X}],\qquad
I(x)=\sup_{\lambda\in\mathbb R}\{\lambda x-\Lambda(\lambda)\}.$$

则对 $x>\mu$，

$$\lim_{n\to\infty}\frac1n\log\mathbb P(\bar X_n\ge x)=-I(x).$$

$I$ 称为**速率函数**，是 $\Lambda$ 的 Legendre–Fenchel 共轭：凸、非负、在 $x=\mu$ 处取最小值 $0$。这个"共轭"结构是全书反复出现的母题：自由能是内能的 Legendre 变换，Hamiltonian 是 Lagrangian 的 Legendre 变换，大偏差速率是累积量母函数的 Legendre 变换。

两个例子值得记住。Bernoulli$(p)$ 的速率函数是

$$I(x)=x\log\frac{x}{p}+(1-x)\log\frac{1-x}{1-p}=D(\text{Bern}(x)\,\|\,\text{Bern}(p)),$$

即 KL 散度。高斯 $N(0,1)$ 的速率函数是 $I(x)=x^2/2$。

**Sanov 定理**把大偏差从"均值"推广到"经验分布"：若 $X_i\sim p$，经验分布 $\hat\mu_n=\frac1n\sum\delta_{X_i}$，则对分布集合 $A$，

$$\mathbb P(\hat\mu_n\in A)\approx\exp\Big\{-n\inf_{q\in A}D(q\|p)\Big\},\qquad
D(q\|p)=\int q\log\frac{q}{p}.$$

这是"最大熵/最小相对熵"原则的概率论根基，也是统计力学里"熵 = 大偏差速率函数"的严格版本。

**极值理论**回答"最大值落在哪里"。独立同分布的最大值 $M_n=\max_{i\le n}X_i$，若适当归一化后收敛，极限只能是三型之一：

$$\text{Gumbel: } \Lambda(x)=\exp(-e^{-x}),\qquad
\text{Fréchet: } \Phi_\alpha(x)=\exp(-x^{-\alpha})\ (x>0),\qquad
\text{Weibull: } \Psi_\alpha(x)=\exp(-(-x)^\alpha)\ (x\le0).$$

选择哪一型由尾部的衰减速度决定：指数尾给 Gumbel，幂尾给 Fréchet，有界尾给 Weibull。

---

## 3. Markov 链（第 3 章）

**Markov 性**是"给定现在，未来与过去条件独立"。离散时间、有限状态的链由转移矩阵 $P=(p_{ij})$ 描述，Chapman–Kolmogorov 方程是矩阵乘法：

$$p_{ij}^{(n+m)}=\sum_k p_{ik}^{(n)}p_{kj}^{(m)},\qquad P^{(n)}=P^n.$$

**不变分布** $\pi$ 满足

$$\pi=\pi P,\qquad \pi_i\ge0,\quad \sum_i\pi_i=1.$$

存在唯一性由 **Perron–Frobenius 定理**保证：不可约非负矩阵有唯一正的主特征值（等于谱半径），对应正特征向量；若矩阵本原（非周期），其余特征值的模严格小于主特征值。于是 $P^n\to\mathbf 1\pi$，收敛速率由**谱隙** $1-\lambda_2$ 控制。

**遍历定理**把"时间平均"和"空间平均"划等号：

$$\frac1n\sum_{k=0}^{n-1}f(X_k)\ \longrightarrow\ \sum_i\pi_i f(i)\quad\text{a.s.}$$

这是 MCMC 能工作的理论依据：只要链不可约、正常返，样本的时间平均就是目标分布的期望。

收敛速率的定量工具是**全变差距离** $\|\mu P^n-\pi\|_{\mathrm{TV}}=\frac12\sum_i|\mu P^n(i)-\pi_i|$。证明套路是**耦合**：构造两个边缘分别为 $\mu P^n$ 和 $\pi$ 的联合过程，让它们在某个随机时刻之后重合；重合时刻的尾概率给收敛速率。Doeblin 条件（存在 $n$ 使 $P^n$ 的所有元素有正下界）给出几何收敛。

**Poisson 过程**是连续时间最简单的跳过程：$N_0=0$，独立增量，$N_t-N_s\sim\text{Poisson}(\lambda(t-s))$。等待时间独立且服从指数分布 $\text{Exp}(\lambda)$，这是"无记忆性"的直接后果，也是 Gillespie 算法和 KMC 的基石。

**Q-过程**（连续时间 Markov 链）由生成元 $Q$ 描述：$q_{ij}\ge0$（$i\ne j$）是 $i\to j$ 的跳速率，$q_{ii}=-\sum_{j\ne i}q_{ij}$。转移概率矩阵 $P(t)$ 满足

$$P(t)=e^{Qt},\qquad
\frac{dP}{dt}=P(t)Q\ \text{（前向/Forward）},\qquad
\frac{dP}{dt}=QP(t)\ \text{（后向/Backward）}.$$

前向方程演化"当前分布在何处"，后向方程演化"从某点出发的期望"。两者是同一枚硬币的两面，这个对偶结构在第 8 章会以 Fokker–Planck / 后向方程的形式重演。

**嵌入链**把连续时间链拆成"跳到哪里"和"跳得多快"：嵌入链的转移概率是 $q_{ij}/|q_{ii}|$，等待时间是 $\text{Exp}(|q_{ii}|)$。

**时间反演**：若链有不变分布 $\pi$，反演链的转移概率为

$$\hat p_{ij}=\frac{\pi_j}{\pi_i}p_{ji}.$$

若 $\hat p=p$，即

$$\pi_i p_{ij}=\pi_j p_{ji}\quad\forall i,j,$$

则称链满足**细致平衡**。细致平衡比"不变"更强：不变只要求每个节点的净流入为零，细致平衡要求每条边上的净流为零。这个区别在第 8 章的扩散过程里会变成"$\nabla\cdot j_s=0$"与"$j_s=0$"的区别。

**Dirichlet 型**把链的"能量"写成边的二次型：

$$\mathcal E(f,g)=-\langle f,Lg\rangle_\pi=\frac12\sum_{i,j}\pi_i p_{ij}(f_i-f_j)(g_i-g_j).$$

它是谱隙、收敛速率、以及第 8 章 Witten Laplacian 的离散原型。

**隐 Markov 模型（HMM）**是"观测不完整"的链：状态 $s_t$ 按 $P$ 演化，观测 $x_t$ 由发射概率 $B$ 生成。三个算法构成一条流水线：前向算法算 $\alpha_t(i)=\mathbb P(x_1,\dots,x_t,s_t=i)$；Viterbi 算法用 $\delta_t(i)=\max_{s_1,\dots,s_{t-1}}\mathbb P(s_1,\dots,s_{t-1},s_t=i,x_1,\dots,x_t)$ 回溯最可能路径；Baum–Welch 是 EM，用前向–后向的期望计数更新 $P,B$。

---

## 4. Monte Carlo 方法（第 4 章）

**误差的量级**是 Monte Carlo 的全部卖点。用 $N$ 个独立样本估计 $\mathbb Ef(X)$，估计量的标准差是 $\sigma/\sqrt N$，与维数无关；而确定性网格法的误差随维数指数恶化。这就是"高维积分只能靠抽样"的原因。

**均匀随机数**：线性同余发生器 $x_{n+1}=(ax_n+c)\bmod m$ 便宜但有结构。Marsaglia 超平面定理说，$d$ 维的 LCG 点列最多落在 $m^{1/d}$ 个超平面上；这是"伪随机数不能直接用在高维"的经典警告。

**从均匀到目标分布**的三条路：

- 逆变换：$X=F^{-1}(U)$，$U\sim\text{Unif}(0,1)$。要求 $F$ 可逆，一维首选。
- Box–Muller：$Z_1=\sqrt{-2\log U_1}\cos(2\pi U_2)$，$Z_2=\sqrt{-2\log U_1}\sin(2\pi U_2)$，得到两个独立标准正态。
- 接受–拒绝：从提议分布 $g$ 采样 $Y$，以概率 $f(Y)/(Mg(Y))$ 接受。接受率是 $1/M$，$M$ 越紧越好。

**方差缩减**是"用更少的样本换同样的精度"：

- 重要性采样：$\mathbb E_p[f(X)]=\mathbb E_q[f(X)p(X)/q(X)]$。最优提议 $q\propto|f|p$，但通常只能近似。
- 控制变量：用已知期望的 $g$ 构造 $f-c(g-\mathbb Eg)$，选 $c$ 使方差最小。
- Rao–Blackwellization：用 $\mathbb E[f(X)\mid Y]$ 代替 $f(X)$，条件期望的方差不超过原变量。
- 分层、对偶、公共随机数：针对具体问题的技巧。

**Markov 链 Monte Carlo（MCMC）**解决"目标分布只知道未归一化密度"的情形。Metropolis–Hastings 的提议 $q(y\mid x)$，接受概率

$$\alpha(x,y)=\min\left\{1,\ \frac{\pi(y)q(x\mid y)}{\pi(x)q(y\mid x)}\right\}.$$

它自动满足细致平衡 $\pi(x)q(x\mid y)\alpha(x,y)=\pi(y)q(y\mid x)\alpha(y,x)$，因此 $\pi$ 是不变分布。收敛速度由谱隙和"链是否混合"决定；实践中要看自相关时间、有效样本量（ESS）。ESS 的定义是 $N_{\rm eff}=N/(1+2\sum_{k\ge1}\rho_k)$，$\rho_k$ 是滞后 $k$ 的自相关；经验上 ESS 小于几百就意味着链没混合好，此时增加样本数没用，要改提议分布或重参数化。

**Ising 模型**是 MCMC 的试金石：$H=-J\sum_{\langle i,j\rangle}s_i s_j-h\sum_i s_i$，$s_i=\pm1$。Glauber 动力学以速率 $1/(1+e^{\beta\Delta H})$ 翻转自旋；Metropolis 以 $\min(1,e^{-\beta\Delta H})$ 接受翻转。二维 Ising 在 $T_c$ 附近出现临界 slowing down，单自旋更新效率骤降，这是 cluster 算法（Swendsen–Wang、Wolff）的动机。

**Kinetic Monte Carlo（KMC）/BKL** 处理"事件稀疏、速率已知"的跳过程。总速率 $R_{\rm tot}=\sum_i r_i$，时间增量

$$\Delta t=-\frac{\log u}{R_{\rm tot}},\qquad \mathbb P(\text{选事件 }i)=r_i/R_{\rm tot}.$$

它精确复现了指数等待时间，是 Gillespie 算法的前身。

**模拟退火与模拟回火**：把目标分布写成 $\pi_\beta\propto e^{-\beta H}$，让 $\beta=\beta(n)$ 随步数增大。若冷却足够慢，

$$\beta(n)\le\frac{C}{\log n}\ \Longrightarrow\ \text{依概率收敛到全局最优}.$$

太快会卡在局部极小；这是"探索–利用"权衡的数学版本。

---

## 5. 随机过程的一般框架（第 5 章）

**Kolmogorov 扩张定理**是"存在性"的总开关：若给出一族相容的有限维分布（任意有限个时刻的联合分布，且边缘化一致），则存在一个概率空间和过程，使这些分布就是它的有限维分布。它把"构造过程"降级为"写出一族分布"，代价是样本空间取为 $\mathbb R^{[0,\infty)}$、$\sigma$-代数取乘积 $\sigma$-代数，轨道性质（连续性等）需要另外证明。

**滤流与停时**：$\mathcal F_t$ 表示 $t$ 时刻已知的信息；停时 $\tau$ 满足 $\{\tau\le t\}\in\mathcal F_t$，即"是否已停"在 $t$ 时刻可知。首达时、首出时都是停时，这是 Dynkin 公式和强 Markov 性的前提。

**Markov 过程**的现代定义用转移核：$p(s,x;t,B)$ 是"$s$ 时刻在 $x$、$t$ 时刻落入 $B$"的条件概率，满足

$$\mathbb P(X_{t+s}\in B\mid\mathcal F_s)=\mathbb P(X_{t+s}\in B\mid X_s).$$

时齐情形下 $p(s,x;t,B)=p(t-s,x,B)$，于是有**半群** $T_t f(x)=\mathbb E_x[f(X_t)]$，满足 $T_{t+s}=T_tT_s$、$T_0=I$。**无穷小生成元**是

$$Af(x)=\lim_{t\downarrow0}\frac{T_tf(x)-f(x)}{t}.$$

对扩散过程，$A$ 就是第 8 章的 $L$；对跳过程，$A$ 是第 3 章的 $Q$ 或第 13 章的化学主方程生成元。整个"用 PDE 研究随机过程"的程序，就是把半群翻译成生成元。

**Gaussian 过程**由均值函数 $m(t)=\mathbb EX_t$ 和协方差函数 $K(s,t)=\mathrm{Cov}(X_s,X_t)$ 完全决定：任意有限维分布都是高斯。**Mercer 定理**说，正定核可以展开为

$$K(s,t)=\sum_{k=1}^\infty\lambda_k e_k(s)e_k(t),\qquad \lambda_k\ge0,$$

$\{e_k\}$ 是正交基。于是有 **Karhunen–Loève 展开**

$$X_t=\sum_{k=1}^\infty\sqrt{\lambda_k}\,\xi_k\,e_k(t),\qquad \xi_k\ \text{i.i.d.}\ N(0,1),$$

收敛在 $L^2$ 意义下。这是"把随机过程对角化"的标准手段，第 6 章的 Wiener 过程和第 10 章的高斯随机场都用它。

**Markov 性与强 Markov 性**要分清。Markov 性说"给定当前，未来与过去条件独立"；强 Markov 性说"把当前换成停时，结论仍成立"。后者不是前者的推论：它需要过程有"重新开始"的结构。扩散过程和 $Q$-过程都有强 Markov 性，这是 Dynkin 公式、首达时分解、以及"从边界重新注入"算法的基础。反例也存在：某些带边界的确定性流程在首达时处会"记住"进入方向，从而破坏强 Markov 性。

**生成元的定义域**是一个容易被忽略的技术点。$A$ 一般是无界算子，只在 $D(A)\subset B$ 上定义；$D(A)$ 在 $B$ 中稠密，$A$ 是闭算子。对扩散过程，$C_c^2\subset D(A)$，且在这些函数上 $A=L$；但 $D(A)$ 通常比 $C_c^2$ 大，包含一些二阶可微性更弱的函数。Hille–Yosida 定理（附录 F）给出"$A$ 生成压缩半群"的充要条件，这是所有"算子 ↔ 过程"构造的抽象底座。

---

## 6. Wiener 过程（第 6 章）

**定义**：$W_0=0$；独立增量；$W_t-W_s\sim N(0,t-s)$；轨道几乎必然连续。等价的高斯刻画是 $m(t)=0$、$K(s,t)=s\wedge t$。

**从随机游走到 Wiener**：设随机游走每步长 $\pm l$、时间步长 $\tau$，则 $n$ 步后的方差是 $nl^2$，对应扩散系数

$$D=\frac{l^2}{2\tau}.$$

取 $l=\sqrt{2D\tau}$、$\tau\to0$，则 $S_{\lfloor t/\tau\rfloor}\Rightarrow W_t$。这是**Donsker 不变原理**的雏形：不只是单点分布收敛，而是整条路径的分布收敛（在 $C[0,T]$ 上）。

**热方程**：$p(x,t)=\frac1{\sqrt{2\pi t}}e^{-x^2/(2t)}$ 满足 $\partial_t p=\frac12\partial_{xx}p$。更一般地，$u(x,t)=\mathbb E_x[f(W_t)]$ 满足后向热方程，这是第 8 章 Feynman–Kac 的最简情形。

**反正弦律**：$\mathbb P\big(\int_0^1\mathbf 1_{W_t>0}\,dt\le s\big)=\frac2\pi\arcsin\sqrt s$。它说明 Brownian 轨道"停留时间"的分布是 U 形的——大部分时间待在正侧或负侧，而不是均匀分配。这是"直觉失效"的经典例子。

**KL 展开**：Wiener 过程在 $[0,1]$ 上的协方差核 $s\wedge t$ 的特征值是

$$\lambda_k=\frac{1}{((k-\tfrac12)\pi)^2},\qquad
e_k(t)=\sqrt2\sin\big((k-\tfrac12)\pi t\big),\quad k=1,2,\dots$$

于是

$$W_t=\sum_{k=1}^\infty\frac{\sqrt2}{(k-\tfrac12)\pi}\sin\big((k-\tfrac12)\pi t\big)\,\xi_k.$$

这个展开把 Wiener 过程变成"可数多个独立高斯的叠加"，是路径积分有限维近似的严格版本。Haar 基给出另一组展开，每一项是区间上的分段常数，适合数值实现。

**路径性质**（必须记住的四条）：

- 二次变差：$\sum_i(W_{t_{i+1}}-W_{t_i})^2\to t$（依 $L^2$），记作 $[W,W]_t=t$。
- 总变差：几乎必然无穷。因此 $\int_0^t|\dot W_s|\,ds=\infty$，路径不可求导。
- Hölder 连续性：对任意 $\alpha<\frac12$，轨道几乎必然是 $\alpha$-Hölder 的；$\alpha=\frac12$ 不成立。
- 重对数律：$\limsup_{t\to\infty}\frac{W_t}{\sqrt{2t\log\log t}}=1$ a.s.

**反射原理**：$\mathbb P(\max_{0\le s\le t}W_s\ge a)=2\mathbb P(W_t\ge a)=2(1-\Phi(a/\sqrt t))$。它把"最大值分布"化为"端点分布"，是首达时、障碍期权、以及第 8 章边值问题的常用工具。

**Wiener 混沌展开**：$L^2(\Omega,\sigma(W))$ 可以正交分解为 $\bigoplus_n\mathcal H_n$，$\mathcal H_n$ 由 $n$ 阶 Hermite 多项式（Wick 乘积）张成。这就是多项式混沌（polynomial chaos）/gPC 的数学基础：把随机响应写成"高斯基 + 系数"的展开，用于不确定性量化。

---

## 7. Itô 随机微分方程（第 7 章）

**白噪声**是 Wiener 过程的形式导数：$\dot W_t$ 不存在为普通函数，但作为广义过程，它的协方差是 $\delta(t-s)$。SDE 的写法 $dX_t=b\,dt+\sigma\,dW_t$ 是"白噪声驱动"的紧凑记号。

**Itô 积分**的构造分两步。对简单过程 $H_t=\sum_i h_i\mathbf 1_{[t_i,t_{i+1})}$，定义

$$\int_0^T H_s\,dW_s=\sum_i h_i\,(W_{t_{i+1}}-W_{t_i}).$$

关键点是**被积函数取左端点**（非预见），这保证了积分的鞅性。然后由 **Itô 等距**

$$\mathbb E\Big[\Big(\int_0^T H_s\,dW_s\Big)^2\Big]=\mathbb E\int_0^T H_s^2\,ds$$

把定义扩展到 $H\in L^2(0,T)$。等距是"随机积分是 $L^2$ 等距映射"的精确表述，也是所有数值误差分析的起点。

为什么不能把 $\int H\,dW$ 当成普通的 Riemann–Stieltjes 积分？因为 Wiener 轨道的总变差无穷，而 Stieltjes 积分要求积分器有界变差。Itô 的做法是"绕开轨道、直接在 $L^2$ 里定义"：先用简单过程逼近，再用等距延拓。代价是链式法则要改（多出二阶项），收益是积分保持鞅性——这两件事是同一枚硬币的两面。

**Itô 公式**是随机分析的链式法则。对 $f(t,x)$，

$$df=\partial_t f\,dt+\nabla f\cdot dX+\frac12\sum_{i,j}\partial^2_{ij}f\,d\langle X_i,X_j\rangle.$$

对 $dX=b\,dt+\sigma\,dW$，$d\langle X_i,X_j\rangle=(\sigma\sigma^T)_{ij}dt$。形式乘法表是

$$dt\cdot dt=0,\qquad dt\cdot dW=0,\qquad dW_i\cdot dW_j=\delta_{ij}\,dt.$$

最重要的一行是 $(dW)^2=dt$：它不是"小量的高阶项"，而是一阶项，这正是随机分析与确定性微积分的分水岭。

基本例子：$\int_0^t W_s\,dW_s=\frac12W_t^2-\frac12t$。多出来的 $-\frac12t$ 就是 Itô 修正项；若用 Stratonovich 积分，结果是 $\frac12W_t^2$。

**BDG 不等式**（Burkholder–Davis–Gundy）把 $L^p$ 控制从端点搬到上确界：

$$\mathbb E\Big[\sup_{s\le t}\Big|\int_0^s H\,dW\Big|^p\Big]\le C_p\,\mathbb E\Big[\Big(\int_0^t H^2\,ds\Big)^{p/2}\Big].$$

它是证明 SDE 解的存在性、数值格式的强收敛、以及鞅的极限定理的核心工具。

**存在唯一性**：若 $b,\sigma$ 满足 Lipschitz 条件和线性增长条件，则 SDE 有唯一强解。证明用 Picard 迭代 $X^{(n+1)}_t=x+\int_0^t b(X^{(n)}_s)ds+\int_0^t\sigma(X^{(n)}_s)dW_s$，再用 Gronwall 不等式控制 $\mathbb E\sup|X^{(n+1)}-X^{(n)}|^2$。

**三个必须会的模型**：

- **Ornstein–Uhlenbeck**：$dX_t=-\gamma X_t\,dt+\sigma\,dW_t$，解为
  $$X_t=e^{-\gamma t}X_0+\sigma\int_0^te^{-\gamma(t-s)}dW_s,$$
  平稳分布 $N(0,\sigma^2/(2\gamma))$，自协方差 $\frac{\sigma^2}{2\gamma}e^{-\gamma|t-s|}$。
- **几何 Brownian 运动**：$dS_t=\mu S_t\,dt+\sigma S_t\,dW_t$，解为
  $$S_t=S_0\exp\Big((\mu-\tfrac12\sigma^2)t+\sigma W_t\Big).$$
  对应的 Black–Scholes PDE 是 $\partial_tV+\frac12\sigma^2S^2\partial_{SS}V+rS\partial_SV-rV=0$。
- **Langevin**：$m\,dV_t=-\gamma V_t\,dt-\nabla U(X_t)\,dt+\sqrt{2\gamma k_BT}\,dW_t$。它是第 11 章涨落耗散关系的原型。

**扩散过程的刻画**：$dX=b(X)dt+\sigma(X)dW$ 的生成元是

$$Lf=b\cdot\nabla f+\frac12A:\nabla^2f,\qquad A=\sigma\sigma^T.$$

反过来，给定一个二阶椭圆算子 $L$，就能构造对应的扩散过程。这是"算子 ↔ 过程"字典的核心条目。

**Stratonovich 积分**取中点，遵守普通链式法则；与 Itô 积分的换算关系是

$$b^{\text{It\^o}}_i=b^{\text{Strat}}_{i}+\frac12\sum_{j,k}\sigma_{jk}\partial_j\sigma_{ik}.$$

一维时 $b^{\text{It\^o}}=b^{\text{Strat}}+\frac12\sigma\sigma'$。物理建模常用 Stratonovich（因为白噪声是真实有色噪声的极限），金融与数值分析常用 Itô（因为鞅性质好）。

**数值格式的阶数**（强收敛：轨道逼近；弱收敛：期望逼近）：

| 格式 | 强阶 | 弱阶 |
|---|---|---|
| Euler–Maruyama | $1/2$ | $1$ |
| Milstein | $1$ | $1$ |
| Itô–Taylor（高阶） | $1.5$ 或 $2$ | $2$ |

Euler–Maruyama 是 $X_{n+1}=X_n+b(X_n)\Delta t+\sigma(X_n)\Delta W_n$；Milstein 多了一项 $\frac12\sigma\sigma'(\Delta W_n^2-\Delta t)$，用来消掉一阶强误差。

**多层 Monte Carlo（MLMC）**：把期望写成 telescoping sum

$$\mathbb Ef(X_T)=\mathbb Ef(X_T^{(0)})+\sum_{\ell=1}^L\mathbb E\big[f(X_T^{(\ell)})-f(X_T^{(\ell-1)})\big],$$

每层用不同步长、不同样本数，把方差大的层放在粗网格上。强阶 $1/2$ 时，总代价从标准 MC 的 $O(\varepsilon^{-3})$ 降到 $O(\varepsilon^{-2}(\log\varepsilon)^2)$。这是"用算法换精度"的典范。

---

## 8. Fokker–Planck 方程（第 8 章）

**Fokker–Planck（FP）方程**是概率密度的演化方程。对 $dX=b(X)dt+\sigma(X)dW$，

$$\partial_t p=L^*p=-\nabla\cdot(b\,p)+\frac12\nabla^2:(A\,p),\qquad A=\sigma\sigma^T.$$

它是前向 Kolmogorov 方程；对偶的**后向方程**是 $\partial_t u=Lu$，$u(x,t)=\mathbb E_x[f(X_t)]$。两者由分部积分（伴随）联系：$\langle Lf,g\rangle=\langle f,L^*g\rangle$。

**生成元** $Lf=b\cdot\nabla f+\frac12A:\nabla^2f$ 是"沿轨道的期望变化率"；FP 算子 $L^*$ 是"密度的输运 + 扩散"。$L$ 决定"从一点出发会怎样"，$L^*$ 决定"一团概率会怎样"。

**Smoluchowski（过阻尼 Brownian）方程**：$\gamma\dot X=-\nabla U+\sqrt{2\gamma k_BT}\,\dot W$，FP 为

$$\partial_t p=\nabla\cdot\Big(\frac1\gamma p\,\nabla U+\frac{k_BT}{\gamma}\nabla p\Big).$$

它可以写成**自由能梯度流**：

$$F[p]=\int U(x)p(x)\,dx+k_BT\int p\log p\,dx,\qquad
\mu=\frac{\delta F}{\delta p}=U+k_BT\log p,$$

$$\partial_t p=\nabla\cdot\Big(\frac{p}{\gamma}\nabla\mu\Big).$$

$F$ 是"内能 + 熵"的 Helmholtz 自由能；$\mu$ 是化学势。这个结构说明：扩散方程不是随便写的，它是自由能泛函在 Wasserstein 度量下的梯度流。平衡态 $p_s\propto e^{-U/k_BT}$ 使 $F$ 最小、$\mu$ 为常数。

**边界条件**：反射边界要求法向概率流为零，$\mathbf n\cdot j=0$，其中 $j=-\frac1\gamma p\nabla\mu$；吸收边界要求 $p=0$。前者保概率，后者"漏概率"，对应 Dirichlet 边值问题。

**不变分布与细致平衡**。不变分布 $p_s$ 满足 $L^*p_s=0$，即 $\nabla\cdot j_s=0$。若更强的 $j_s=0$ 成立，则称**细致平衡**：

$$b(y)p_s(y)=\frac12\nabla\cdot(A(y)p_s(y)).$$

对 Smoluchowski，这给出 Gibbs 分布 $p_s\propto e^{-U/k_BT}$。对 Langevin（$x$ 偶、$v$ 奇），需要**广义细致平衡**：存在对角符号矩阵 $s=\mathrm{diag}(I,-I)$，使

$$p_s(x)=p_s(sx),\qquad p(z,t\mid x,0)p_s(x)=p(sx,t\mid sz,0)p_s(z).$$

对 SDE 的等价条件是

$$s\,b(sy)p_s(y)=-b(y)p_s(y)+\nabla\cdot(A(y)p_s(y)),\qquad s\,A(sy)=A(y)\,s.$$

代入 Langevin 得到 $p_s(y)=\frac1Z\exp(-H(x,v)/k_BT)$，$H=\frac12m|v|^2+U(x)$。这是第 11 章"正则分布是不变分布"的严格版本。

**Feller 半群**：在 $b,\sigma$ 有界 Lipschitz 时，$T_tf(x)=\mathbb E_x[f(X_t)]$ 是 $C_0(\mathbb R^d)$ 上的强连续压缩半群，生成元 $A$ 在 $C_c^2$ 上与 $L$ 一致。**定理 8.5**：$u(x,t)=\mathbb E_x[f(X_t)]$ 满足 $\partial_tu=Au$，$u(0)=f$。

**Feynman–Kac 公式**：若 $v$ 解 $\partial_tv=Lv+qv$，$v(0)=f$，则

$$v(x,t)=\mathbb E_x\Big[\exp\Big(\int_0^tq(X_s)ds\Big)f(X_t)\Big].$$

它是"PDE 的解 = 路径上的指数权重平均"。$q$ 可以看作 killing rate、贴现率或势能；$q<0$ 时它给路径打折，$q>0$ 时它给路径加权。这个公式在第 9 章会用路径积分重新解释，在第 11 章会变成线性响应理论的出发点。

**Dynkin 公式**把 Feynman–Kac 从固定时刻推广到停时：

$$\mathbb E_x[f(X_\tau)]=f(x)+\mathbb E_x\Big[\int_0^\tau Lf(X_t)dt\Big],\qquad \mathbb E_x\tau<\infty.$$

**边值问题**：若 $Lu=-g$ 在 $U$ 内、$u=f$ 在 $\partial U$ 上，则

$$u(x)=\mathbb E_x\Big[f(X_{\tau_U})+\int_0^{\tau_U}g(X_t)dt\Big],\qquad
\tau_U=\inf\{t\ge0:X_t\notin U\}.$$

这是"解 PDE 不用网格，用随机游走"的理论依据；walk-on-sphere 算法就是它的蒙特卡洛实现。

**平均首达时间（MFPT）**：一维、$a$ 反射、$b$ 吸收，$\tau(x)=\mathbb E_x\tau_b$ 满足

$$L\tau=-1,\qquad \tau(b)=0,\qquad \tau'(a)=0.$$

对 $dX=-U'\,dt+\sqrt{2\varepsilon}\,dW$，解为

$$\tau(x)=\frac1\varepsilon\int_x^b e^{U(y)/\varepsilon}\int_a^y e^{-U(z)/\varepsilon}dz\,dy.$$

第 12 章的 Kramers 速率就是对这个积分做 Laplace 渐近得到的。

**谱理论**：当细致平衡成立时，$L$ 在加权空间 $L^2_\mu$（$\mu=p_s$）中自伴，且 $-L$ 正定：

$$(Lf,g)_\mu=(f,Lg)_\mu,\qquad (-Lf,f)_\mu=\|\nabla f\|_\mu^2\ge0.$$

因此特征值 $0=\lambda_0>\lambda_1\ge\lambda_2\ge\cdots$ 全为实数、非正；$L^*$ 的特征函数是 $\psi_\lambda=\mu\phi_\lambda$。谱隙 $\lambda_1$ 控制收敛到平衡的速率。

**Witten Laplacian**：相似变换

$$\Delta_U=e^{-U/2}Le^{U/2}=\Delta-\frac14|\nabla U|^2+\frac12\Delta U$$

把扩散算子变成"量子 Hamiltonian"，势能为 $V=\frac14|\nabla U|^2-\frac12\Delta U$。这把扩散的亚稳态问题与量子隧穿联系起来，是"WKB 分析"（第 12 章）的谱论视角。

**多尺度渐近**：设后向方程为

$$\partial_tu^\epsilon=\frac1{\epsilon^2}L_1u^\epsilon+\frac1\epsilon L_2u^\epsilon+L_3u^\epsilon.$$

若 $e^{L_1t}\to P$（投影到 $L_1$ 的零空间）、$PL_2P=0$（可解性条件）、$Pf=f$（初值一致性），则展开 $u^\epsilon=u_0+\epsilon u_1+\cdots$ 给出有效方程

$$\partial_tu_0=\bar Lu_0,\qquad
\bar L=PL_3P-PL_2L_1^{-1}L_2P.$$

这就是"消去快变量"的一般配方：投影掉快方向，用 $L_1^{-1}$ 解出快变量的响应，再代回慢方程。第 11 章的 Mori–Zwanzig、第 13 章的化学多尺度分析，都是这个配方的具体化。

**Euler–Maruyama 的弱收敛**：局部截断误差 $O(h^2)$，全局弱误差 $O(h)$。证明的关键是把误差写成"后向方程的解在两条轨道上的差"，再用 Itô 公式和 Gronwall 控制。这解释了为什么"弱阶 = 局部阶 − 1"。

---

## 9. 路径积分（第 9 章）

**形式 Wiener 测度**：把 Wiener 过程的有限维密度

$$p_n(w_1,\dots,w_n)=Z_n^{-1}\exp\Big(-\frac12\sum_{j=1}^n\frac{(w_j-w_{j-1})^2}{t_j-t_{j-1}}\Big)$$

取极限，形式上得到

$$d\mu_W=\delta(w_0)\,Z^{-1}\exp\Big(-\frac12\int_0^1\dot w_t^2\,dt\Big)\,\mathcal Dw.$$

这个表达式有三处"不成立"，必须知道：

1. 归一化常数 $Z_n=(2\pi/n)^{n/2}\to0$；
2. $\int_0^1\dot w_t^2dt=+\infty$ a.s.（Wiener 路径不可导）；
3. 路径空间上不存在平移不变、局部有限的 Lebesgue 测度 $\mathcal Dw$（用单位球内互相分离的小球可导出矛盾）。

所以路径积分是**形式工具**，不是严格测度。但它在物理直觉、渐近分析、以及"有限维近似再取极限"的计算中极其有效。

**有限维计算范式**：以 $\mathbb E\exp(-\frac12\int_0^1W_t^2dt)$ 为例。离散化后，被积函数是 $\exp(-\frac12\delta t\,x^TAx)$，Wiener 权重是 $\exp(-\frac12\delta tx^TBx)$，$B$ 是离散二阶差分矩阵。高斯积分给出

$$\mathbb E\exp\Big(-\frac12\int_0^1W_t^2dt\Big)=\Big(\frac{\det B}{\det(A+B)}\Big)^{1/2},$$

两个行列式都发散，但比值有定义；取极限后与 KL 展开的结果 $\sqrt{2e/(1+e^2)}$ 一致。这个"行列式比值"就是路径积分的通用算法：**先把无穷维行列式约成分式，再取极限**。

路径积分要严格化，通常走三条路。第一条是**有限维近似**：先固定分割，算有限维 Gaussian 积分，再证明极限与分割无关；KL 展开就是这条路的对角化版本。第二条是**测度论路线**：直接在路径空间上构造 Wiener 测度（Kolmogorov 扩张定理 + 连续性定理），把 $\mathcal Dw$ 解释为 Wiener 测度，把 $e^{-I}$ 解释为 Radon–Nikodym 导数；Girsanov 定理就是这条路的产物。第三条是**算子路线**：把路径积分看成半群 $e^{tL}$ 的核，用谱理论或热核估计来定义；Feynman–Kac 公式是这条路的代表。

三条路各有代价：有限维近似直观但收敛性琐碎；测度论严格但需要拓扑与 $\sigma$-代数；算子路线最优雅但把"路径"藏进了半群。做计算用第一条，做理论用第二条，做谱分析用第三条。

**Girsanov 变换**是路径空间的"换测度"公式。设

$$d\tilde W_t=\phi(t,\omega)dt+dW_t,$$

在**Novikov 条件** $\mathbb E\exp(\frac12\int_0^T|\phi|^2ds)<\infty$ 下，定义

$$Z_t=\exp\Big(-\int_0^t\phi\,dW_s-\frac12\int_0^t\phi^2ds\Big),\qquad
d\tilde{\mathbb P}=Z_T\,d\mathbb P.$$

则在 $\tilde{\mathbb P}$ 下 $\tilde W$ 是标准 Wiener 过程。等价的说法是：把漂移"吸进"测度，让过程在新测度下变成无漂移的 Brownian。$Z_T$ 就是 Radon–Nikodym 导数。

对一般扩散 $dX=b\,dt+\sigma\,dW$ 与 $dY=(b+\gamma)dt+\sigma\,dW$，取 $\phi=\sigma^{-1}\gamma$，则

$$\frac{d\mu_Y}{d\mu_X}=\exp\Big(-\frac12\int_0^1|\phi|^2dt+\int_0^1\phi\,dW_t\Big).$$

**Cameron–Martin 公式**是 $\phi=\dot\varphi$（确定性平移 $X_t=W_t+\varphi(t)$）的特例：

$$\frac{d\mu_X}{d\mu_W}=\exp\Big(-\frac12\int_0^1\dot\varphi^2dt+\int_0^1\dot\varphi\,dW_t\Big).$$

注意：只有当平移 $\varphi$ 属于 Cameron–Martin 空间（$H^1$）时，两个测度才互相绝对连续；一般的连续平移会得到互相奇异的测度。

**Feynman–Kac 的路径积分形式**：

$$v(x,t)=\int\delta(w_0-x)\,Z^{-1}\exp\Big(-\int_0^t\big(\tfrac12\dot w_s^2-q(w_s)\big)ds\Big)f(w_t)\,\mathcal Dw.$$

**Feynman 与 Schrödinger**：量子力学的形式解为

$$\psi(x,t)=\int\delta(w_0-x)Z^{-1}\exp\Big(\frac{i}{\hbar}I[w]\Big)\psi_0(w_t)\,\mathcal Dw,\qquad
I[w]=\int_0^t\Big(\frac m2\dot w_s^2-V(w_s)\Big)ds.$$

把 $m=1$、$\hbar=-i$ 代入，虚时间 Schrödinger 方程变成抛物方程，Feynman 的积分变成 Kac 的严格公式。这就是"Feynman–Kac"名字的来历：Kac 把 Feynman 的物理直觉翻译成了概率论。

---

## 10. 随机场（第 10 章）

**随机场**是"多元随机函数"，即定义在 $\mathbb R^d$ 或图上的随机函数/广义函数。四类典型来源：统计力学（Ising 自旋构型）、量子场论（函数积分）、图像处理（Markov 随机场）、结构可靠性与湍流。

**Gaussian 随机场**由三件东西完全决定：

- 均值泛函 $m(\psi)=\mathbb E(\psi,\omega)$（$\psi$ 是测试函数，$(\cdot,\cdot)$ 是对偶配对）；
- 协方差算子 $(\psi,K\varphi)=\mathbb E(\psi,\omega-m)(\varphi,\omega-m)$；
- 特征泛函 $S(\psi)=\mathbb E e^{i(\psi,\omega)}$。

严格地说，Gaussian 随机场通常是**广义函数空间上的测度**，样本是随机分布而不是随机函数；所以一切都要用测试函数来检验。

三个标准例子：

- **Gaussian 白噪声**：$m=0$，$K=I$，$S(\psi)=e^{-\frac12\|\psi\|^2}$。不同区域的限制独立；一维时它是 Wiener 过程的形式导数。
- **OU 测度**：$K=(-\Delta+I)^{-1}$，一维对应平稳 OU 过程，协方差 $\frac12e^{-|t-s|}$。
- **Brownian sheet**：二维版本，$\mathbb EB(s,t)=0$，$\mathbb EB(s,t)B(\tilde s,\tilde t)=(s\wedge\tilde s)(t\wedge\tilde t)$；由二维随机游走的部分和 $S_{m,n}=\sum_{i\le m}\sum_{j\le n}\xi_{ij}$ 归一化得到。

**Markov 随机场**（离散情形）：在图上，若把顶点分成 $V_1,V_2,V_3$，$V_1$ 与 $V_3$ 之间没有边，则

$$\mathbb P(X_1=x_1\mid X_2=x_2,X_3=x_3)=\mathbb P(X_1=x_1\mid X_2=x_2),$$

即"给定分隔集，两侧条件独立"。

**Gibbs 分布**：取团（clique）集合 $\mathcal C$，给每个团一个势函数 $U_C$，

$$U(x)=\sum_{C\in\mathcal C}U_C(x),\qquad \mathbb P(x)=\frac1Z e^{-U(x)}.$$

**Hammersley–Clifford 定理**：在 $\mathbb P(x)>0$ 的前提下，Gibbs 分布与 Markov 随机场**等价**。这是"局部交互 ⇒ 全局 Gibbs 形式"的严格表述，也是图像分割、语义标记里 MRF 建模的理论基础。

**协方差算子的谱**决定了随机场的"相关长度"与"光滑度"。若 $K$ 的特征值衰减慢，场在空间上粗糙；衰减快，场光滑。白噪声的 $K=I$ 所有特征值相等，因此完全不相关、样本几乎必然不连续；OU 测度的 $K=(-\Delta+I)^{-1}$ 特征值 $\sim1/(1+k^2)$，给出指数衰减的相关函数与连续样本。这个"谱 ↔ 光滑度"的对应，与第 6 章 Wiener 过程的 KL 展开是同一件事。

**Markov 随机场的推断**：给定观测 $y$，后验

$$p(x\mid y)\propto\exp\Big(-U(x)-\frac{1}{2\sigma^2}\|y-x\|^2\Big)$$

仍是 Gibbs 形式，只是势函数多了数据项。于是可以用 Gibbs 采样、ICM（iterated conditional modes）或图割（graph cut）求解。这是图像去噪、分割、立体匹配的标准套路。

**广义 Ising**：$U(x)=\sum_{\langle i,j\rangle}(x_i-x_j)^2+\lambda\sum_j(x_j-u_j)^2$，$\mathbb P\propto e^{-\beta U}$。$\beta\to0$ 趋于均匀分布，$\beta\to\infty$ 趋于 $U$ 的全局最小点；$\beta^{-1}$ 扮演温度。这个模型可以直接用于图像去噪/分割：第一项鼓励相邻像素平滑，第二项拉住观测数据。

---

## 11. 统计力学（第 11 章）

**三种系综**对应三种边界交换：

- 孤立系统：$(N,V,E)$，**微正则系综**，能量固定；
- 封闭系统：$(N,V,T)$，**正则系综**，与热库交换能量；
- 开放系统：$(\mu,V,T)$，**巨正则系综**，同时交换能量与粒子；
- 力学平衡：$(N,P,T)$，**等温等压系综**。

**熵**：Boltzmann 形式 $S=k_B\log g(E)$（$g$ 是能量 $E$ 的微观状态数）；Gibbs 形式 $S=-k_B\sum_xp(x)\log p(x)$。微正则时 $p(x)=1/g(E)$，两者一致。熵是广延量，因为取了对数。

**温度、化学势、压强**都是自由能对相应变量的偏导：

$$\frac1T=\Big(\frac{\partial S}{\partial E}\Big)^{-1},\qquad
\mu=\Big(\frac{\partial F}{\partial N}\Big)_{V,T},\qquad
P=-\Big(\frac{\partial F}{\partial V}\Big)_{N,T}.$$

热平衡要求 $T_1=T_2$；化学平衡要求 $\mu_1=\mu_2$；力学平衡要求 $P_1=P_2$。这三条都是"自由能最小"的推论。

**正则分布**：把系统 $A$ 与热库 $R$ 接触，总能量固定，对热库的微观状态数做 Taylor 展开：

$$p(x)\propto g_R(E-\epsilon(x))=\exp\Big(\sigma_R(E)-\frac{\epsilon(x)}{k_BT}+\cdots\Big)
\ \Longrightarrow\ p(x)=\frac1Z e^{-\beta\epsilon(x)}.$$

配分函数 $Z(T)=\sum_xe^{-\beta\epsilon(x)}=\int_0^\infty g(\epsilon)e^{-\beta\epsilon}d\epsilon$，即密度 of states 的 Laplace 变换。

**热力学量**：

$$U=\langle\epsilon\rangle=-\frac{\partial\log Z}{\partial\beta},\qquad
F=-k_BT\log Z,\qquad F=U-TS.$$

能量涨落：

$$\langle(\delta\epsilon)^2\rangle=\frac{\partial^2\log Z}{\partial\beta^2}=k_BT^2C_V,\qquad
\frac{\sqrt{\langle(\delta\epsilon)^2\rangle}}{\langle\epsilon\rangle}\sim O(N^{-1/2}).$$

最后这个 $N^{-1/2}$ 说明：大 $N$ 时微正则与正则等价，宏观量是确定的。

**巨正则**：$p(x)=\frac1{Z(\mu,T)}e^{\beta(\mu n(x)-\epsilon(x))}$，grand potential $\Phi=-k_BT\log Z=U-TS-\mu N$。**等温等压**：$G=F+PV=U-TS+PV$，且 $G(P,T)=\inf_V\{F(V,T)+PV\}$，这是 $F$ 关于 $V$ 的 Legendre–Fenchel 变换。

**理想气体**：单粒子配分函数 $Z_1=\frac V{h^3}\int e^{-\beta p^2/2m}dp=n_QV$，$n_Q=(2\pi mk_BT/h^2)^{3/2}$；$Z_N=Z_1^N/N!$（$N!$ 来自不可分辨性）。由此得

$$PV=Nk_BT,\qquad
S=k_BN\Big(\log\frac{n_Q}{n_V}+\frac52\Big),\qquad
P=C_0\rho^{5/3}\ \text{（绝热）}.$$

**Ising 模型与转移矩阵**：一维最近邻

$$H_N(x;\beta,h)=-\beta J_0\sum_{i=1}^{N-1}x_ix_{i+1}-\beta h\sum_{i=1}^Nx_i.$$

配分函数写成矩阵乘积 $Z_N=\frac12\sum_{x_1,x_N}a(x_1)B^{N-1}(x_1,x_N)a(x_N)$，$B$ 是 $2\times2$ 转移矩阵。由 Perron–Frobenius，

$$\frac1N\log Z_N\to\log\lambda_B,\qquad
\lambda_B=e^{\beta J_0}\Big(\cosh\beta h+\sqrt{\sinh^2\beta h+e^{-4\beta J_0}}\Big).$$

自由能 $F(y;\beta,h)$ 是序参量 $y=\frac1N\sum x_i$ 的速率函数，通过 Legendre–Fenchel 变换与累积量母函数 $H(\theta)$ 互换：

$$H(\theta)=\sup_y\{\theta y-F(y)\},\qquad F(y)=\sup_\theta\{\theta y-H(\theta)\}.$$

这是 Gärtner–Ellis 定理的实例。一维最近邻 Ising 的 $F$ 严格凸、唯一极小，因此**没有相变**；相变需要二维及以上。

**Langevin 与涨落耗散**：$m\dot V=-\gamma V+\sqrt{2\sigma}\dot W$ 的平稳动能是 $\langle\frac12mV^2\rangle=3\sigma/(2\gamma)$。热力学要求它等于 $\frac32k_BT$，于是

$$\sigma=k_BT\gamma.$$

这就是**涨落–耗散关系**：噪声强度 $\sigma$ 与摩擦系数 $\gamma$ 由同一微观起源决定，不能独立指定。由此还得 Einstein 关系

$$D=\lim_{t\to\infty}\frac{\langle(X_t-X_0)^2\rangle}{6t}=\frac{k_BT}{\gamma}.$$

过阻尼极限（$t\gg m/\gamma$）给出 $\gamma\dot X=-\nabla U+\sqrt{2k_BT\gamma}\dot W$。

**广义 Langevin**：带记忆核与有色噪声，

$$m\dot V_t=-\int_{-\infty}^t\gamma(t-s)V_sds-\nabla U+\xi_t,\qquad
\langle\xi_t\xi_s\rangle=k_BT\gamma(|t-s|).$$

这就是涨落–耗散关系的一般形式：噪声的协方差就是记忆核（乘 $k_BT$）。当 $\gamma(t)=2\gamma\delta(t)$ 时退回普通 Langevin。

**线性响应理论**：在平衡态上加小外场 $F(t)$，响应

$$\langle\delta A(t)\rangle=\int_0^t\chi_{AB}(t-s)F(s)ds,$$

响应函数与平衡态自相关函数的关系是

$$\chi_{AB}(t)=-\beta\frac{d}{dt}C_{AB}(t),\qquad C_{AB}(t)=\langle A(t)B(0)\rangle_{\rm eq}.$$

这就是**涨落–耗散定理**：左边是"被外力驱动时的响应"，右边是"没有外力时的自发涨落"。频域形式是 Green–Kubo 关系

$$\hat\mu(\omega)=\beta\int_0^\infty\langle A(t)\dot B(0)\rangle_{\rm eq}e^{-i\omega t}dt,\qquad
D=\hat\mu(0)k_BT.$$

**Mori–Zwanzig 约化**：只保留感兴趣的分量 $\hat x$，把其余变量投影掉。取 $Pf(\hat x)=\mathbb E(f(x)\mid\hat x)$，$Q=I-P$。由 Dyson 公式

$$e^{tL}=e^{tQL}+\int_0^t e^{(t-s)L}PL\,e^{sQL}ds,$$

得到广义 Langevin 方程

$$\frac{d\phi}{dt}=e^{tL}PL\phi_0+\int_0^t e^{(t-s)L}K(s)ds+R(t),$$

其中 $K(t)=PLR(t)$ 是记忆核，$R(t)=e^{tQL}QL\phi_0$ 是噪声项。三项分别对应：Markov 项、记忆项、噪声项。记忆与噪声的关联正是涨落–耗散关系的一般来源。

**Kac–Zwanzig 模型**：一个 distinguished 粒子（质量 1）通过 $N$ 个谐振子热浴耦合，

$$H=\frac12v^2+U(x)+\sum_{i=1}^N\frac12m_iv_i^2+\frac{k}{2N}\sum_{i=1}^N(x_i-x)^2.$$

解出热浴变量后代回，得到

$$\ddot X=-U'(X)-\int_0^t\gamma_N(t-s)\dot X(s)ds+\frac1{\sqrt\beta}\xi_N(t),$$

$$\gamma_N(t)=\frac kN\sum_{i=1}^N\cos\omega_it,\qquad
\mathbb E\xi_N(t)\xi_N(s)=\gamma_N(t-s).$$

当 $N\to\infty$、$\omega_i$ 独立同分布时，$\gamma_N\to\gamma(t)=k\int\cos(\omega t)\rho(\omega)d\omega$，噪声收敛为高斯过程。这是"从 Hamiltonian 微观模型严格导出广义 Langevin 方程"的完整示例。

---

## 12. 稀有事件（第 12 章）

**亚稳态**：考虑

$$dX_t^\epsilon=-\nabla U(X_t^\epsilon)dt+\sqrt\epsilon\,dW_t.$$

$\epsilon=0$ 时，ODE 收敛到 $U$ 的局部极小；每个极小有吸引域，边界是 separatrix。$\epsilon>0$ 但很小时，过程在 $O(1)$ 时间内待在某个吸引域内，但在 $O(e^{\delta U/\epsilon})$ 的时间尺度上会"跳"到另一个域。这就是**亚稳态**与**稀有事件**。

三个核心问题：最可能路径是什么？过渡态（saddle）在哪里？平均过渡时间多长？

**WKB 分析**：把转移密度写成 $p^\epsilon\propto\exp(-\epsilon^{-1}\varphi)$，代入 FP 方程，领头阶给出 Hamilton–Jacobi 方程

$$\partial_t\varphi=H(x,\nabla_x\varphi),\qquad
H(x,p)=b^Tp+\frac12p^TAp.$$

$p$ 是"动量"，$A=\sigma\sigma^T$。对应的 Lagrangian 是 $H$ 的 Legendre 变换：

$$L(x,z)=\frac12(z-b(x))^TA^{-1}(x)(z-b(x)),\qquad
L(x,z)=\sup_p\{p\cdot z-H(x,p)\}.$$

作用量 $I_t[\phi]=\int_0^tL(\phi,\dot\phi)ds$，且

$$\varphi(x,t\mid y)=\inf_{\phi(0)=y,\phi(t)=x}I_t[\phi].$$

**Kramers 速率**：一维双井 $U$，$x_-$ 与 $x_+$ 为极小、$x_s$ 为极大。MFPT 的 Laplace 渐近给出

$$\tau(x_-)\sim\frac{2\pi}{\sqrt{|U''(x_s)|U''(x_-)}}e^{\delta U/\epsilon},\qquad
k=\frac1\tau=\frac{\sqrt{|U''(x_s)|U''(x_-)}}{2\pi}e^{-\delta U/\epsilon}.$$

这就是 **Kramers/Arrhenius 定律**：速率指数由能垒 $\delta U=U(x_s)-U(x_-)$ 决定，前因子由曲率决定。多维、index-one saddle 的版本是

$$k=\frac{|\lambda_s|}{2\pi}\sqrt{\frac{\det H_-}{|\det H_s|}}e^{-\delta U/\epsilon},$$

$\lambda_s<0$ 是 saddle 处唯一的负特征值。

**大偏差原理（Freidlin–Wentzell）**：对 $dX^\epsilon=b\,dt+\sqrt\epsilon\sigma\,dW$，在 $b,\sigma$ 有界 Lipschitz、$A$ 一致椭圆时，对任意 $T>0$：

$$\limsup_{\epsilon\to0}\epsilon\log\mathbb P(X^\epsilon\in F)\le-\inf_{\phi\in F}I_T[\phi]\quad(F\ \text{闭}),$$

$$\liminf_{\epsilon\to0}\epsilon\log\mathbb P(X^\epsilon\in G)\ge-\inf_{\phi\in G}I_T[\phi]\quad(G\ \text{开}),$$

$I_T[\phi]=\int_0^T\frac12(\dot\phi-b)^TA^{-1}(\dot\phi-b)dt$。这就是"路径空间的 Laplace 原理"：概率的指数由作用量最小的路径决定。

**最小作用路径（MEP）**：对梯度系统 $b=-\nabla U$、$\sigma=I$，连接两个极小的最小作用路径由两段组成：

$$\dot\phi_1=\nabla U\ (\text{爬坡}),\qquad \dot\phi_2=-\nabla U\ (\text{下坡}),$$

最小作用 $I^*=2\delta U$。路径满足

$$\dot\phi=\pm\nabla U\iff(\nabla U(\phi))^\perp=0,$$

即"梯度方向与曲线切向平行"——这就是**最小能量路径**的定义。MEP 必过 index-one saddle。

**String method**：把 MEP 作为曲线演化，

$$\partial_t\phi(\alpha,t)=-(\nabla U(\phi))^\perp=-\nabla U(\phi)+r\hat\tau,\qquad
\phi(0)=x_-,\ \phi(1)=x_+.$$

数值上分两步：先让每个点沿梯度下降走一步，再按等弧长重新参数化。$r$ 是 Lagrange 乘子，用来固定参数化；等弧长条件 $(|\phi_\alpha|)_\alpha=0$。这是目前最常用的 MEP 算法之一。

**准势（quasipotential）**：

$$V(x;y)=\inf_{T>0}\inf_{\phi(0)=y,\phi(T)=x}I_T[\phi].$$

对梯度系统，$V(x;y)=2(U(x)-U(y))$（在 $y$ 的吸引域内）；对非梯度系统，它是"势能"的推广。它满足稳态 Hamilton–Jacobi 方程

$$H(x,\nabla_xV)=0,$$

且在极小路径上 $H=0$（"零能量路径"）。由 Maupertuis 原理，可以消去时间优化：

$$V(x;y)=\inf_{H(\phi,p)=0}\int_\gamma p\cdot d\phi
=\inf_{\psi(0)=y,\psi(1)=x}\int_0^1\big(|b(\psi)||\psi'|-b(\psi)\cdot\psi'\big)d\alpha.$$

这个形式在数值上更好用，因为它只在曲线空间里优化，不需要处理时间。

**Exit problem**：设 $0$ 是 $\dot x=b(x)$ 的稳定平衡点，$D$ 是其吸引域的一部分。则

$$\lim_{\epsilon\to0}\epsilon\log\mathbb E\tau^\epsilon=\min_{x\in\partial D}V(x;0),$$

且若最小值点 $x_0$ 唯一，则 $X^\epsilon_{\tau^\epsilon}\to x_0$（依概率）。这就是"从哪个门出去、要等多久"的严格答案。

**全局准势**：$U(x)=-\lim_{\epsilon\to0}\epsilon\log p^\epsilon_\infty(x)$，其中 $p^\epsilon_\infty$ 是不变测度。注意两个极限不可交换：

$$-\lim_{\epsilon\to0}\lim_{t\to\infty}\epsilon\log p^\epsilon(x,t\mid y)=U(x)
\ \ne\ -\lim_{t\to\infty}\lim_{\epsilon\to0}\epsilon\log p^\epsilon(x,t\mid y)=V(x;y).$$

局部准势 $V(x;y)$ 可以"剪接"成全局准势 $U(x)$。

---

## 13. 化学反应动力学（第 13 章）

**反应速率方程（RRE）**：质量作用定律说，反应速率正比于反应常数与参与物浓度的乘积。以 decaying-dimerization 为例：

$$S_1\xrightarrow{c_1}\varnothing,\qquad
2S_1\rightleftharpoons[c_2]{c_3}S_2,\qquad
S_2\xrightarrow{c_4}S_3,$$

$$\dot x_1=-c_1x_1-2c_2x_1^2+2c_3x_2,\qquad
\dot x_2=c_2x_1^2-c_3x_2-c_4x_2,\qquad
\dot x_3=c_4x_2.$$

这是"大量分子、忽略离散随机性"的宏观模型。

**化学主方程（CME）**：当分子数很少（单细胞、基因表达），必须逐分子建模。设 $N$ 个物种、$M$ 个反应通道，状态 $X_t\in\mathbb N^N$；通道 $j$ 有倾向函数 $a_j(x)$ 和状态改变向量 $\nu_j$。$a_j(x)dt$ 是"当前在 $x$、未来 $dt$ 内发生一次 $R_j$"的概率。CME 是

$$\partial_tP(x,t)=\sum_{j=1}^M\Big[a_j(x-\nu_j)P(x-\nu_j,t)-a_j(x)P(x,t)\Big].$$

它是 Kolmogorov 前向方程，但状态空间是 $\mathbb N^N$，所以是无穷维 ODE 而非 PDE。CME 的解析解几乎从不存在的：状态空间随物种数指数增长，即使 10 个物种、每个 100 个分子，也有 $100^{10}$ 个状态。这正是 SSA 与各种近似算法存在的理由。

倾向函数的形式来自分子碰撞理论：

$$a_j(x)=\frac{\kappa_j}{V^{|\nu_j^-|-1}}\binom{x}{\nu_j^-},\qquad
\binom{x}{\nu_j^-}=\prod_{m=1}^N\binom{x_m}{\nu_j^{-,m}}.$$

$V$ 是体积（乘 Avogadro 常数），$\kappa_j$ 是过渡态理论的反应速率。组合数来自"从 $x_m$ 个分子中选 $\nu_j^{-,m}$ 个"的概率。

**跳过程的 SDE 表示**：用独立单位速率 Poisson 过程 $P_j$，

$$X_t=X_0+\sum_{j=1}^M\nu_jP_j\Big(\int_0^ta_j(X_s)ds\Big).$$

这是"随机时间变换"：把真实时间 $t$ 换成"累积倾向" $\int a_jds$，每个通道按单位速率计数。更一般地，用 Poisson 随机测度 $\lambda(dt\times da)$ 和接受–拒绝函数 $c_j(a;X_t)=\mathbf 1_{a\in(h_{j-1}(X_t),h_j(X_t)]}$ 表示，$h_j=\sum_{i\le j}a_i$。

**Gillespie 算法（SSA）**：给定 $X_t=x$，下一次反应等待时间 $\tau\sim\text{Exp}(a_0(x))$，$a_0=\sum_ja_j$；反应通道 $j$ 被选中的概率是 $a_j(x)/a_0(x)$；更新 $X_{t+\tau}=X_t+\nu_j$。它是**精确**的，因为对应的前向方程正是 CME。

SSA 的两个瓶颈：速率悬殊（快反应控制步长，但我们关心慢反应）、分子数很大（离散性成为负担）。这催生了后面的近似与多尺度算法。

**大体积极限（Kurtz limit）**：设 $X_0,X_t\sim O(V)$，$a_j(Vx)=Va_j^V(x)$，$a_j^V\to\tilde a_j$。令 $X^V=X/V$（浓度），则

$$X^V_t=X^V_0+\int_0^tb(X^V_s)ds+M^V_t,\qquad
b(x)=\sum_j\nu_j\tilde a_j(x),$$

$M^V$ 是中心化 Poisson 的缩放。由 Poisson 的大数律与 Gronwall，

$$\sup_{s\le t}|X^V_s-x(s)|\to0\ \text{a.s.},\qquad
\dot x=b(x).$$

这就是"从 CME 回到 RRE"的严格桥梁。速率常数的换算关系是

$$c_j=\frac{\kappa_j}{\nu_j^-!}=\frac{k_j}{\nu_j^-!V^{|\nu_j^-|-1}}.$$

**扩散近似（化学 Langevin）**：CLT 给出涨落的线性扩散。令 $Y^V=\sqrt V(X^V-x)$，则

$$dY=\sum_j\nu_j\nabla\tilde a_j(x(t))\cdot Y\,dt+\sum_j\nu_j\sqrt{\tilde a_j(x(t))}\,dW_j.$$

直接把 Poisson 的 CLT 代入跳过程，得到

$$dZ_t=\sum_j\nu_ja_j(Z_t)dt+\sum_j\nu_j\sqrt{a_j(Z_t)}\,dW_j,$$

这就是**化学 Langevin 方程**。它在 $V$ 大时与 CME 同阶，但把离散跳变成了连续扩散。

**Tau-leaping**：SSA 逐事件模拟太慢，tau-leaping 在步长 $\tau$ 内"批量"计数：

$$X_{t+\tau}=X_t+\sum_{j=1}^M\nu_j\,\text{Poisson}(a_j(X_t)\tau).$$

步长选取的 **leap condition**：$\tau$ 要小到倾向函数在 $[t,t+\tau)$ 内没有可观变化，即 $|a_j(X+\tau\xi)-a_j(X)|\le\epsilon a_0(X)$，Taylor 展开给

$$\tau\approx\min_j\frac{\epsilon a_0(X)}{|\xi\cdot\nabla a_j|},\qquad \xi=\sum_j\nu_ja_j(X).$$

显式 tau-leaping 的强阶是 $1/2$、弱阶是 $1$（在 $\tau\to0$ 标度下），因为它本质上是跳过程 SDE 的 Euler–Maruyama 离散。

**平稳分布**：可逆反应 $S_1\rightleftharpoons[k_1]{k_2}S_2$，$x+y=n$。细致平衡 $a_1(x)p(x)=a_2(x-1)p(x-1)$ 给出

$$p(x)=\binom nx\Big(\frac{k_2}{k_1+k_2}\Big)^x\Big(\frac{k_1}{k_1+k_2}\Big)^{n-x},$$

即二项分布 $B(n,q)$，$q=k_2/(k_1+k_2)$。推广到一阶反应网络 $S_1\rightleftharpoons S_2\rightleftharpoons\cdots\rightleftharpoons S_N$，不变分布是多项分布

$$p(x)=\frac{n!}{x_1!\cdots x_N!}c_1^{x_1}\cdots c_N^{x_N},\qquad c\cdot Q=0,\ \sum_ic_i=1,$$

$Q$ 是单分子反应的生成元。这个"乘积形式"是化学网络理论（deficiency zero）的经典结果。

**多尺度分析**：把反应分成快、慢两组，$a=(a_s,\frac1\epsilon a_f)$。**慢变量**是快反应的守恒量：

$$v(x+\nu_j^f)=v(x)\iff b\cdot\nu_j^f=0.$$

取基 $\{b_1,\dots,b_J\}$，$z_j=b_j\cdot x$。在快时间尺度，$z$ 固定，快过程达到条件平衡 $\mu_z(x)$；在慢时间尺度，有效速率是

$$\bar a_i^s(z)=\mathbb E_{\mu_z}[a_i^s]=\sum_xa_i^s(x)\mu_z(x),\qquad
\bar\nu_i^s=(b_1\cdot\nu_i^s,\dots,b_J\cdot\nu_i^s).$$

有效方程是 $\partial_t\bar U=\sum_i\bar a_i^s(z)(\bar U(z+\bar\nu_i^s)-\bar U(z))$。这与第 8 章的多尺度渐近是同一个公式，只是把微分算子换成了跳算子。

**嵌套 SSA**：外层只跑慢反应，但慢反应的速率由内层 SSA 实时估计：内层跑 $K$ 个只含快反应的副本，计算

$$\tilde a_j^s=\frac1K\sum_{k=1}^K\frac1{T_f}\int_{T_0}^{T_0+T_f}a_j^s(X^k_\tau)d\tau.$$

这样不用预先算出有效速率，平均"在飞行中"完成。它可以递归扩展到多个时间尺度。

---

## 14. 附录：六个工具箱（A–F）

**A. Laplace 渐近与 Varadhan 引理**。若 $h$ 在 $0$ 取唯一全局最大、$h''(0)<0$，则

$$\int_{\mathbb R}e^{th(x)}dx\sim\sqrt{2\pi}\,(-th''(0))^{-1/2}e^{th(0)},\qquad t\to\infty.$$

等价地，$\lim_{t\to\infty}\frac1t\log\int e^{th(x)}dx=\sup_xh(x)$。**大偏差原理**：存在速率函数 $I$ 使闭集上有上界、开集上有下界：

$$\limsup_{\epsilon\to0}\epsilon\log\mathbb P_\epsilon(F)\le-\inf_FI,\qquad
\liminf_{\epsilon\to0}\epsilon\log\mathbb P_\epsilon(G)\ge-\inf_GI.$$

**Varadhan 引理**：若 $F$ 有界连续，则

$$\lim_{\epsilon\to0}\epsilon\log\int e^{F(x)/\epsilon}\mathbb P_\epsilon(dx)=\sup_x\{F(x)-I(x)\}.$$

它是"路径积分 Laplace 方法"的严格版本，也是第 12 章所有渐近的底层工具。

**B. Gronwall 不等式**。若 $f(t)\le a(t)+\int_0^tb(s)f(s)ds$，$a,b\ge0$，则

$$f(t)\le a(t)+\int_0^ta(s)b(s)\exp\Big(\int_s^tb(u)du\Big)ds.$$

常数情形 $f(t)\le ae^{bt}$。离散版本：若 $F_{n+1}\le(1+b_n\delta t)F_n+a_n$，则

$$F_n\le\exp\Big(\sum_{k=0}^{n-1}b_k\delta t\Big)F_0+\sum_{k=0}^{n-1}a_k\exp\Big(\sum_{l=k+1}^{n-1}b_l\delta t\Big).$$

这是所有 SDE 存在性、数值收敛性证明的"最后一步"。

**C. 测度与积分**。测度：$\sigma$-可加、$\mu(\varnothing)=0$。可测函数、简单函数、积分的三层构造（简单函数 → 非负可测 → 一般可测）。三条收敛定理：单调收敛、Fatou 引理、控制收敛。**Radon–Nikodym 定理**：若 $\eta\ll\mu$ 且 $\mu$ 是 $\sigma$-有限，则存在 $f=d\eta/d\mu$ 使 $\eta(A)=\int_Af\,d\mu$。这是条件期望存在的基石。

**D. 鞅**。$X_t$ 是 $\mathcal F_t$-鞅，若 $\mathbb E(X_t\mid\mathcal F_s)=X_s$（$s\le t$）。凸函数作用后得下鞅。**Doob 不等式**：

$$\mathbb P\Big(\sup_{u\in[s,t]}X_u\ge\lambda\Big)\le\frac{\mathbb EX_t^+}{\lambda},\qquad
\mathbb E\Big[\sup_{u\in[s,t]}X_u\Big]^p\le\Big(\frac p{p-1}\Big)^p\mathbb EX_t^p.$$

它是强收敛、首达时估计、以及 BDG 不等式的来源。

**E. 强 Markov 性**。Markov 性：给定 $X_m=i$，$\{X_{m+n}\}$ 是以 $\delta_i$ 为初分布、$P$ 为转移矩阵的链，且与过去独立。强 Markov 性把 $m$ 换成停时 $N$：给定 $N<\infty$、$X_N=i$，$\{X_{N+n}\}$ 仍是以 $\delta_i$ 为初分布的链。连续时间版本对 $Q$-过程同样成立。这是 Dynkin 公式和"首达时后重新开始"论证的前提。

**F. 算子半群**。强连续半群：$S(0)=I$、$S(t+s)=S(t)S(s)$、$\|S(t)f-f\|\to0$。压缩半群：$\|S(t)\|\le1$。生成元 $Af=\lim_{t\downarrow0}\frac1t(S(t)f-f)$，定义域 $D(A)$ 在 $B$ 中稠密，$A$ 是闭算子（可以无界）。**Hille–Yosida 定理**：$A$ 生成压缩半群，当且仅当 $\lambda\in\rho(A)$ 且 $\|R_\lambda\|=\|(\lambda I-A)^{-1}\|\le1/\lambda$ 对所有 $\lambda>0$ 成立。这是"从算子回到半群"的充要条件，也是第 8 章 Feller 半群存在的抽象保证。

这六个工具箱的分工可以这样记：Laplace/Varadhan 管"指数尺度的极限"；Gronwall 管"误差传播"；测度论管"存在性"；鞅管"上确界与停时"；强 Markov 管"重新开始"；半群管"算子与过程的对应"。证明一个 SDE 结果的标准流水线是：测度论构造概率空间 → 鞅/BDG 控制随机积分 → Gronwall 控制确定性误差 → 半群/Hille–Yosida 把生成元翻译回过程。少任何一环，证明都会漏气。

---

## 15. 收束：十个必背公式、依赖图与读法

### 10 个必背公式

1. **大偏差速率**：$I(x)=\sup_\lambda\{\lambda x-\Lambda(\lambda)\}$，$\Lambda=\log\mathbb E e^{\lambda X}$。
2. **不变分布**：$\pi=\pi P$；连续时间 $P(t)=e^{Qt}$。
3. **KL 展开**：$K(s,t)=\sum_k\lambda_ke_k(s)e_k(t)$，$X_t=\sum_k\sqrt{\lambda_k}\xi_ke_k(t)$。
4. **Itô 公式**：$df=f_tdt+f_xdX+\frac12f_{xx}(dX)^2$，$(dW)^2=dt$。
5. **Fokker–Planck**：$\partial_tp=-\nabla\cdot(bp)+\frac12\nabla^2:(Ap)$。
6. **Feynman–Kac**：$v(x,t)=\mathbb E_x[e^{\int_0^tq(X_s)ds}f(X_t)]$。
7. **Girsanov**：$Z_t=\exp(-\int_0^t\phi\,dW-\frac12\int_0^t\phi^2ds)$，$d\tilde{\mathbb P}=Z_Td\mathbb P$。
8. **Kramers 速率**：$k\simeq\frac{\sqrt{|U''(x_s)|U''(x_-)}}{2\pi}e^{-\delta U/\epsilon}$。
9. **化学主方程**：$\partial_tP=\sum_j[a_j(x-\nu_j)P(x-\nu_j)-a_j(x)P(x)]$。
10. **Mori–Zwanzig**：$\frac{d\phi}{dt}=e^{tL}PL\phi_0+\int_0^te^{(t-s)L}K(s)ds+R(t)$。

### 依赖图

```
概率基础(1) ──> 极限定理(2) ──> 大偏差(12)
     │                │
     v                v
 Markov链(3) ──> 随机过程(5) ──> Wiener(6) ──> SDE(7) ──> Fokker-Planck(8)
     │                              │              │            │
     v                              v              v            v
 Monte Carlo(4)              路径积分(9)      数值格式(7)   谱/多尺度(8)
     │                              │              │            │
     └────────────> 稀有事件(12) <──┘              v            v
                                     随机场(10) ──> 统计力学(11) ──> 化学反应(13)
```

三条最短阅读路径：

- **概率 → 分析路线**（想做理论）：1 → 2 → 3 → 5 → 6 → 7 → 8 → 9 → 12。
- **计算 → 应用路线**（想做算法）：1 → 2 → 3 → 4 → 7 → 8 → 12 → 13。
- **物理 → 多尺度路线**（想做建模）：3 → 5 → 7 → 8 → 11 → 12 → 13。

### 五个常见误区

1. **把"依分布收敛"当成"轨道收敛"**。CLT 只给分布，Donsker 才给路径；数值格式的强阶和弱阶是两回事。
2. **忘记 Itô 修正项**。$(dW)^2=dt$ 是一阶项，不是高阶小量；换到 Stratonovich 要把 $\frac12\sigma\sigma'$ 加回去。
3. **把"不变"当成"细致平衡"**。不变只要求 $\nabla\cdot j_s=0$，细致平衡要求 $j_s=0$；后者才有自伴谱理论。
4. **在稀有事件里用"典型路径"直觉**。最可能路径是"爬坡 + 下坡"两段，爬坡段由噪声驱动，下坡段由梯度驱动；$I^*=2\delta U$ 而不是 $\delta U$。
5. **把大偏差当成"小概率的普通渐近"**。它是指数尺度的律，$I(x)$ 是变分问题的解，不能由局部 Taylor 展开得到。

### 全书反复使用的六个证明技术

1. **截断 + Borel–Cantelli**（强收敛）：先把坏事件截掉，证明截断后的级数绝对收敛，再用 Borel–Cantelli 说明坏事件只发生有限次。
2. **耦合**（Markov 链收敛）：构造两个边缘正确的联合过程，让它们在随机时刻重合；重合时间的尾概率给全变差上界。
3. **Picard 迭代 + Gronwall**（SDE 存在唯一）：逐次逼近，用 Gronwall 控制相邻两项的差，得到 $L^2$ 收敛与唯一性。
4. **Itô 公式 + 鞅/BDG**（矩估计与强收敛）：把误差写成随机积分，用 BDG 把上确界的 $L^p$ 范数压回积分范数。
5. **Laplace 渐近/Varadhan**（大偏差与稀有事件）：把概率写成指数积分，用最速下降/变分求指数；Kramers 速率、准势、exit problem 都是这一招。
6. **投影 + Fredholm 可解性**（多尺度平均）：把快变量投影掉，用 $L_1^{-1}$ 解出快响应，可解性条件 $PL_2P=0$ 保证慢方程封闭。

把这六个技术认全，书中绝大多数证明的骨架就能一眼看出来。真正的难点通常不在技术本身，而在"该用哪一个"——选错工具，证明会绕远路；选对工具，三行就完。

### 这套数学在 AI 里的投影

- **扩散生成模型**：前向加噪是 OU 过程，逆向采样是时间反演的 SDE；FP 方程、Girsanov（ELBO 的路径测度比）、以及"score = ∇log p"三者是同一件事的三种写法。
- **Langevin 采样 / MCMC**：过阻尼 Langevin 的平稳分布是 Gibbs 分布；MCMC 的收敛速率由谱隙与 Poincaré 不等式控制。
- **稀有事件与重要性采样**：大偏差的作用量给出"最优提议分布"的理论极限；退火、tempering、以及强化学习里的探索噪声都可以在这个框架下读。
- **化学/生物动力学**：CME、Gillespie、tau-leaping 是单细胞噪声建模的标准工具；大体积极限给出确定性 ODE，扩散近似给出化学 Langevin。
- **神经网络的训练动力学**：SGD 的连续极限是 SDE，其 FP 方程与稳态涨落决定泛化；这属于本书第 7、8、11 章的直接延伸。

---

> **一句话收尾**：这本书的真正核心不是某个定理，而是**三组字典**——概率 ↔ PDE（生成元/半群）、路径 ↔ 作用量（大偏差/路径积分）、微观 ↔ 宏观（平均/齐次化/统计力学）。把这三组字典背熟，剩下的都是查表。
