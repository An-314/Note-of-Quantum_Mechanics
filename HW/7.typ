#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第7次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [7.1])[
  若 $A,B$ 为厄米算符，证明 $(A B+B A)/2$ 与 $(A B-B A)/(2i)$ 均为厄米算符，并把任意算符分解为两个厄米算符之和。
]
#proof[
  利用 $(A B)^dagger=B A$ 即得前两式。对任意 $F$，令
  $
    F_1=(F+F^dagger)/2, quad F_2=(F-F^dagger)/(2i),
  $
  则 $F_1,F_2$ 均为厄米算符，且 $F=F_1+i F_2$。
]

#exercise(subname: [7.2])[
  证明 $bold(p) times bold(L)+bold(L) times bold(p)=2 i hbar bold(p)$。
]
#proof[
  第 $i$ 个分量为
  $
    epsilon_(i j k)(p_j L_k+L_j p_k)
    =epsilon_(i j k)(p_j L_k-p_k L_j)
    =epsilon_(i j k)[p_j,L_k].
  $
  利用 $[L_k,p_j]=i hbar epsilon_(k j l)p_l$，并缩并两个 Levi-Civita 符号，得到 $2 i hbar p_i$。
]

#exercise(subname: [7.3])[
  证明处于离散能量本征态的粒子有 $〈 bold(p) 〉=0$。
]
#proof[
  对 $H=bold(p)^2/(2m)+V(bold(r))$，有
  $
    [H,bold(r)]=-i hbar bold(p)/m.
  $
  在能量本征态中 $〈 [H,bold(r)] 〉=E 〈 bold(r) 〉-E 〈 bold(r) 〉=0$，故结论成立。
]

#exercise(subname: [7.4])[
  判断并说明下列说法：非定态中所有平均值都随时间变化；定态中本征值的测量概率不变；哈密顿量中出现的物理量都是守恒量；中心势定态中角动量必有确定值；自由粒子定态中动量必有确定值；一维粒子能级从不简并；中心势束缚能级至少有 $2l+1$ 重简并。
]
#solution[
  依次为：错误、正确、错误、需区分、错误、需限定、正确。

  非定态也可能使某些算符的平均值保持不变；定态只积累整体相位，所以任一本征结果的概率不变。物理量守恒要求 $pdv(A,t)=0$ 且 $[A,H]=0$，仅在 $H$ 中出现并不充分。中心势可同时选取 $H,L^2,L_z$ 的共同本征态，此时 $L^2,L_z$ 确定，但任意简并叠加未必如此。自由粒子同一能量可由 $plus.minus bold(p)$ 等简并态叠加。非简并定理适用于一维束缚态，不适用于连续谱。中心势不含方向，能量与 $m$ 无关，故给定 $l$ 至少有 $2l+1$ 重磁量子数简并。
]

#exercise(subname: [7.5])[
  若算符 $A$ 不显含时间，证明系统处于能量本征态时 $〈 A 〉$ 不随时间变化。
]
#proof[
  定态为 $|psi(t)〉=exp(-i E t/hbar)|E 〉$，整体相位在矩阵元中相消，故 $〈 A 〉=〈 E|A|E 〉$ 与时间无关。等价地，Ehrenfest 公式中 $〈 [H,A] 〉=0$。
]

#exercise(subname: [7.6])[
  若 $U$ 为幺正矩阵，证明 $U^dagger,U^(-1),U^*,U^T$ 也都是幺正矩阵。
]
#proof[
  由 $U^dagger U=U U^dagger=1$ 可知 $U^(-1)=U^dagger$，故前两者幺正。对该式取复共轭得 $U^T U^*=1$，再取转置得 $U^*U^T=1$，所以 $U^*$ 与 $U^T$ 互为逆与伴随，二者也幺正。
]

#exercise(subname: [7.7])[
  将坐标表象的定态薛定谔方程变换到动量表象。
]
#solution[
  傅里叶变换下 $-i hbar dif_x$ 对应乘法 $p$，而 $x$ 对应 $i hbar dif_p$。因此
  $
    [p^2/(2m)+V(i hbar dif_p)] phi(p)=E phi(p).
  $
  若势能不能作幂级数理解，则更一般地写成积分核形式
  $
    p^2/(2m) phi(p)+integral 〈 p|V|p' 〉 phi(p') dif p'=E phi(p).
  $
]
