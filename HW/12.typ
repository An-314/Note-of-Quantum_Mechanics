#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第12次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [14.1])[
  把耦合态 $|l-1/2,m 〉$ 展开为非耦合态 $|l,m_l 〉|1/2,m_s 〉$ 的线性组合，并给出系数。
]
#solution[
  给定总磁量子数 $m$ 时只有两项：
  $
    |l-1/2,m 〉=c_1|l,m-1/2 〉|1/2,1/2 〉
    +c_2|l,m+1/2 〉|1/2,-1/2 〉.
  $
  令它与 $|l+1/2,m 〉$ 正交并归一化，可取
  $
    c_1=-sqrt((l-m+1/2)/(2l+1)), quad
    c_2=sqrt((l+m+1/2)/(2l+1)).
  $
  整体符号取决于 Clebsch—Gordan 系数的相位约定。
]

#exercise(subname: [14.2])[
  在 $L-S$ 反平行耦合态 $|l-1/2,m 〉$ 中，$S_z$、$L_z$ 以及 $J_z=L_z+S_z$ 是否具有确定值？若没有，求其平均值。
]
#solution[
  该态是两个不同 $(m_l,m_s)$ 非耦合态的叠加，所以 $S_z,L_z$ 一般没有确定值；$J_z$ 有确定值 $m hbar$。利用上一题系数得
  $
    〈 S_z 〉=-hbar m/(2l+1),
  $
  $
    〈 L_z 〉=m hbar-〈 S_z 〉
    =2(l+1)m hbar/(2l+1).
  $
]

#exercise(subname: [14.3])[
  在按 $|1,1 〉,|1,-1 〉,|1,0 〉,|0,0 〉$ 排列的双自旋耦合基底中，求 $S_(1+)+S_(2-)$ 的矩阵。
]
#solution[
  利用
  $
    |1,0 〉=1/sqrt(2)(|↑↓ 〉+|↓↑ 〉),
  $
  $
    |0,0 〉=1/sqrt(2)(|↑↓ 〉-|↓↑ 〉),
  $
  逐个作用升降算符，得到
  $
    S_(1+)+S_(2-)=hbar/sqrt(2)
    mat(0,0,1,-1;0,0,1,-1;1,1,0,0;1,1,0,0).
  $
]

#exercise(subname: [14.4])[
  证明 $exp(i lambda sigma_z)=cos lambda+i sigma_z sin lambda$，并求其迹。
]
#proof[
  因 $sigma_z^2=1$，指数级数中的偶次幂与奇次幂分别求和为余弦和正弦，故
  $
    exp(i lambda sigma_z)=cos lambda+i sigma_z sin lambda.
  $
  又 $"tr" sigma_z=0$，所以
  $
    "tr" exp(i lambda sigma_z)=2 cos lambda,
  $
  这也等于两个本征值 $exp(plus.minus i lambda)$ 之和。
]

#exercise(subname: [14.5])[
  电子的磁矩算符为 $bold(mu)=-e/(2m_e c)(bold(L)+2bold(S))$。求耦合态 $|l,s=1/2;j,m 〉$ 中 $mu_z$ 的平均值。
]
#solution[
  在定 $j,m$ 态中，矢量平均值沿 $bold(J)$。利用
  $
    〈 L_z 〉=m hbar
    (j(j+1)+l(l+1)-s(s+1))/(2j(j+1)),
  $
  以及对应的 $S_z$ 公式，得到
  $
    〈 mu_z 〉=-mu_B g_j m,
  $
  $
    g_j=1+(j(j+1)+s(s+1)-l(l+1))/(2j(j+1)), quad s=1/2.
  $
]

#exercise(subname: [14.6])[
  双电子自旋态为 $|↓↑ 〉$。在按 $|1,1 〉,|1,0 〉,|1,-1 〉,|0,0 〉$ 排列的耦合基底中，测得各耦合态的概率是多少？
]
#solution[
  由
  $
    |↓↑ 〉=1/sqrt(2)(|1,0 〉-|0,0 〉)
  $
  可知列矩阵为 $mat(0;1/sqrt(2);0;-1/sqrt(2))$，相应概率依次为
  $
    0,quad 1/2,quad 0,quad 1/2.
  $
]
