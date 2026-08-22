#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第3次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [3.1])[
  一维自由粒子的波函数为左右行平面波的叠加
  $
    psi(x,t)=c_1 exp(i(p x-E t)/hbar)+c_2 exp(-i(p x+E t)/hbar).
  $
  求概率流密度。
]
#solution[
  代入 $j=hbar/(2 m i)(psi^* dif_x psi-psi dif_x psi^*)$，交叉项相消，得到
  $
    j=p/m (abs(c_1)^2-abs(c_2)^2).
  $
]

#exercise(subname: [3.2])[
  从薛定谔方程导出概率守恒的连续性方程，并证明归一化随时间保持不变。
]
#proof[
  将薛定谔方程及其复共轭式分别乘以 $psi^*$、$psi$ 后相减，得
  $
    pdv(rho,t)+nabla dot bold(j)=0,
  $
  其中 $rho=abs(psi)^2$，$bold(j)=hbar/(2 m i)(psi^* nabla psi-psi nabla psi^*)$。对全空间积分并令无穷远处概率流为零，即有
  $
    dif_t integral rho dif V=-integral nabla dot bold(j) dif V=0.
  $
]

#exercise(subname: [3.3])[
  证明：一维对称势 $V(x)=V(-x)$ 的束缚态可选为具有确定宇称的态。
]
#proof[
  若 $psi(x)$ 是能量 $E$ 的本征函数，则 $psi(-x)$ 也是同一能量的本征函数。一维束缚能级非简并，因此 $psi(-x)=c psi(x)$。再次反演得 $c^2=1$，故 $c=plus.minus 1$，本征函数分别为偶函数或奇函数。
]

#exercise(subname: [3.4])[
  求三维无限深矩形势阱 $0<x<a,0<y<b,0<z<c$ 的能级和归一化本征函数，并讨论立方势阱的简并。
]
#solution[
  分离变量得到
  $
    psi_(n_x n_y n_z)=sqrt(8/(a b c))
    sin(n_x pi x/a) sin(n_y pi y/b) sin(n_z pi z/c),
  $
  $
    E_(n_x n_y n_z)=hbar^2 pi^2/(2 m)
    (n_x^2/a^2+n_y^2/b^2+n_z^2/c^2).
  $
  对立方势阱 $a=b=c=L$，能量只依赖 $n_x^2+n_y^2+n_z^2$。三个量子数互不相同时，排列通常给出六重简并；两个相同时给出三重简并，此外还可能出现偶然简并。
]

#exercise(subname: [3.5])[
  求中心位于原点、宽度为 $a$ 的一维无限深势阱基态的动量分布。
]
#solution[
  基态为 $psi_1(x)=sqrt(2/a) cos(pi x/a)$（$abs(x)<a/2$）。傅里叶变换给出
  $
    phi_1(p)=sqrt(pi a/hbar) dot
    cos(p a/(2 hbar))/(pi^2-(p a/hbar)^2).
  $
  因而动量概率密度为
  $
    w(p)=pi a/hbar dot cos^2(p a/(2 hbar))/(pi^2-(p a/hbar)^2)^2.
  $
  在 $p=plus.minus pi hbar/a$ 处取连续极限。
]

#exercise(subname: [3.6])[
  用索末菲量子化条件求宽度为 $a$ 的一维无限深势阱能级。
]
#solution[
  粒子在两壁之间往返一周的作用量为 $integral p dif x=2 p a=n h$，故 $p=n h/(2a)$，于是
  $
    E_n=p^2/(2m)=n^2 h^2/(8 m a^2)=n^2 pi^2 hbar^2/(2 m a^2).
  $
]

#exercise(subname: [3.7])[
  求吸引 $delta$ 势 $V(x)=-gamma delta(x)$（$gamma>0$）的束缚态。
]
#solution[
  令 $E=-hbar^2 kappa^2/(2m)$，则 $x != 0$ 时可归一化解为 $psi=A exp(-kappa abs(x))$。在原点积分薛定谔方程得跳跃条件
  $
    psi'(0^+)-psi'(0^-)=-2m gamma/hbar^2 psi(0),
  $
  所以 $kappa=m gamma/hbar^2$。归一化后
  $
    psi(x)=sqrt(m gamma/hbar^2) exp(-m gamma abs(x)/hbar^2),
  $
  $
    E=-m gamma^2/(2 hbar^2).
  $
  仅有这一个束缚态。
]
