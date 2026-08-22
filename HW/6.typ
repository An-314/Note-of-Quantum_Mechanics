#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第6次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [6.1])[
  求 $0<x<a$ 的一维无限深势阱基态中的 $Delta x Delta p$。
]
#solution[
  对 $psi_1=sqrt(2/a)sin(pi x/a)$，直接积分得
  $
    〈 x 〉=a/2, quad
    〈 x^2 〉=a^2(1/3-1/(2 pi^2)),
  $
  $
    (Delta x)^2=a^2(1/12-1/(2 pi^2)).
  $
  又有 $〈 p 〉=0$、$〈 p^2 〉=pi^2 hbar^2/a^2$，故
  $
    Delta x Delta p=hbar sqrt((pi^2-6)/12)>hbar/2.
  $
]

#exercise(subname: [6.2])[
  粒子原来处于宽度为 $a$ 的无限深势阱基态。右壁突然移到 $2a$ 处，求随即测得新势阱第二能级的概率。
]
#solution[
  突变瞬间波函数来不及改变。新势阱第二本征态在 $0<x<2a$ 为
  $
    phi_2(x)=1/sqrt(a) sin(pi x/a).
  $
  它在原势阱内与旧基态 $psi=sqrt(2/a)sin(pi x/a)$ 成比例。因此
  $
    c_2=integral_0^a phi_2^* psi dif x=1/sqrt(2), quad P_2=1/2.
  $
]

#exercise(subname: [6.3])[
  利用谐振子升降算符推导相邻本征函数之间的导数关系。
]
#solution[
  令 $alpha=sqrt(m omega/hbar)$。由
  $
    a=sqrt(m omega/(2hbar))x+1/(sqrt(2)alpha) dif_x,
  $
  $
    a^dagger=sqrt(m omega/(2hbar))x-1/(sqrt(2)alpha) dif_x
  $
  以及 $a psi_n=sqrt(n)psi_(n-1)$、$a^dagger psi_n=sqrt(n+1)psi_(n+1)$，相减得
  $
    dif_x psi_n=alpha(sqrt(n/2)psi_(n-1)-sqrt((n+1)/2)psi_(n+1)).
  $
]

#exercise(subname: [6.4])[
  证明在 $L_z$ 的本征态 $|l m 〉$ 中 $〈 L_x 〉=〈 L_y 〉=0$。
]
#proof[
  由 $L_x=(L_++L_-)/2$、$L_y=(L_+-L_-)/(2i)$，而 $L_+|l m 〉$ 与 $L_-|l m 〉$ 分别正交于 $|l m 〉$，故两个平均值均为零。
]

#exercise(subname: [6.5])[
  系统处于两个同时为 $L^2,L_z$ 本征态的归一化叠加
  $
    |psi 〉=c_1|l_1 m_1 〉+c_2|l_2 m_2 〉.
  $
  求测量 $L^2$ 与 $L_z$ 的可能值、概率和平均值。
]
#solution[
  测得 $L^2$ 的可能值为 $hbar^2 l_1(l_1+1)$、$hbar^2 l_2(l_2+1)$，概率分别为 $abs(c_1)^2$、$abs(c_2)^2$；若两个本征值相同，应把对应概率相加。类似地，$L_z$ 的可能值为 $hbar m_1,hbar m_2$。平均值为
  $
    〈 L^2 〉=hbar^2 sum_(j=1)^2 abs(c_j)^2 l_j(l_j+1),
  $
  $
    〈 L_z 〉=hbar sum_(j=1)^2 abs(c_j)^2 m_j.
  $
]

#exercise(subname: [6.6])[
  若一维束缚势满足 $V(x)>=0$，证明任一束缚态能量 $E>0$。
]
#proof[
  对归一化本征态，
  $
    E=〈 H 〉
    =hbar^2/(2m) integral abs(psi')^2 dif x
    +integral V(x)abs(psi)^2 dif x>=0.
  $
  等号要求 $psi'=0$ 且在 $V>0$ 处 $psi=0$，这与非零、归一化的束缚态矛盾，故实际上 $E>0$。
]

#exercise(subname: [6.7])[
  求自由高斯波包随时间的演化，并说明波包的展宽。
]
#solution[
  设初态
  $
    psi(x,0)=1/(2 pi sigma_0^2)^(1/4)
    exp(-(x-x_0)^2/(4 sigma_0^2)+i p_0 x/hbar).
  $
  对各动量分量乘以 $exp(-i p^2 t/(2m hbar))$ 后反变换，得到
  $
    abs(psi(x,t))^2=1/(sqrt(2 pi) sigma_t)
    exp(-(x-x_0-p_0 t/m)^2/(2 sigma_t^2)),
  $
  其中
  $
    sigma_t=sigma_0 sqrt(1+(hbar t/(2m sigma_0^2))^2).
  $
  波包中心按经典速度 $p_0/m$ 运动，而宽度随时间增加。
]
