#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第2次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [2.1])[
  利用玻尔模型的角动量量子化条件，求氢原子的轨道半径、能级和里德伯常数。
]
#solution[
  圆周运动满足
  $
    m v^2/r=e^2/(4 pi epsilon_0 r^2), quad m v r=n hbar.
  $
  联立得
  $
    r_n=n^2 a_0, quad a_0=(4 pi epsilon_0 hbar^2)/(m e^2),
  $
  $
    E_n=-m e^4/(2(4 pi epsilon_0)^2 hbar^2 n^2)=-13.6/n^2 "eV".
  $
  跃迁波数为 $1/lambda=R_H(1/n_1^2-1/n_2^2)$，其中
  $
    R_H=(m e^4)/(8 epsilon_0^2 h^3 c) approx 1.097 times 10^7 "m"^(-1).
  $
]

#exercise(subname: [2.2])[
  用狄拉克 $delta$ 函数写出坐标本征态和动量本征态的归一化、正交与完备性关系。
]
#solution[
  取
  $
    〈 x|x' 〉=delta(x-x'), quad
    〈 p|p' 〉=delta(p-p'),
  $
  $
    integral |x 〉 〈 x| dif x=1, quad
    integral |p 〉 〈 p| dif p=1.
  $
  在约定 $〈 x|p 〉=(2 pi hbar)^(-1/2) exp(i p x/hbar)$ 下，
  $
    phi(p)=1/sqrt(2 pi hbar) integral exp(-i p x/hbar) psi(x) dif x,
  $
  且帕塞瓦尔等式保证两种表象具有相同的归一化。
]

#exercise(subname: [2.3])[
  已知 $psi_1,psi_2$ 归一且相互正交，求 $psi=c_1 psi_1+c_2 psi_2$ 的归一化条件。
]
#solution[
  $
    〈 psi|psi 〉
    =abs(c_1)^2+abs(c_2)^2
  $
  因此归一化条件是 $abs(c_1)^2+abs(c_2)^2=1$。整体相位不影响物理状态。
]

#exercise(subname: [2.4])[
  粒子严格定域在 $x=x_0$，即 $psi(x)=delta(x-x_0)$，求其动量表象波函数。
]
#solution[
  作傅里叶变换得
  $
    phi(p)=1/sqrt(2 pi hbar) exp(-i p x_0/hbar).
  $
  因而 $abs(phi(p))^2$ 与 $p$ 无关：位置完全确定时动量完全不确定。严格的 $delta$ 态不是平方可积态，应理解为广义本征态。
]

#exercise(subname: [2.5])[
  已知归一化高斯波包
  $
    psi(x)=(alpha/pi)^(1/4) exp(-alpha x^2/2),
  $
  求动量表象波函数。
]
#solution[
  利用高斯积分得
  $
    phi(p)=1/(pi alpha hbar^2)^(1/4)
    exp(-p^2/(2 alpha hbar^2)).
  $
  它满足 $integral abs(phi(p))^2 dif p=1$，并有
  $
    Delta x=1/sqrt(2 alpha), quad Delta p=hbar sqrt(alpha/2),
  $
  故 $Delta x Delta p=hbar/2$。
]

#exercise(subname: [2.6])[
  用索末菲量子化条件 $integral p dif q=n h$ 求一维谐振子的能量。
]
#solution[
  谐振子的相轨道是椭圆
  $
    p^2/(2 m E)+m omega^2 q^2/(2 E)=1.
  $
  其面积即作用量
  $
    integral p dif q=2 pi E/omega=E/nu.
  $
  因而旧量子论给出 $E_n=n h nu=n hbar omega$。它没有包含现代量子力学中的零点能 $hbar omega/2$。
]
