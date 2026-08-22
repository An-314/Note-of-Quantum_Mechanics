#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第4次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)
#let hw-figure(path, width: 42%) = align(center, image(path, width: width))

#exercise(subname: [4.1])[
  说明无限深势阱的量子数为什么不能取零或负数。把势阱从 $[-a,a]$ 平移为 $[0,2a]$ 后，能级是否改变？本征函数是否仍有确定宇称？
]
#solution[
  边界条件要求非零解的波数为 $k_n=n pi/L$。$n=0$ 只给出恒等于零的波函数；$-n$ 与 $n$ 仅相差整体符号，不代表新态。因此取 $n=1,2,dots$。平移不改变宽度 $L=2a$，故
  $
    E_n=n^2 pi^2 hbar^2/(8 m a^2).
  $
  关于势阱中心 $x=a$，本征态仍有确定宇称；但关于坐标原点一般不再是奇函数或偶函数。
]

#exercise(subname: [4.2])[
  能量 $E>0$ 的粒子从左侧入射到势阶
  $
    V(x)=cases(-V_0 & x<0, 0 & x>0),
  $
  求反射系数和透射系数。
  #hw-figure("pic/4.2.png", width: 34%)
]
#solution[
  令
  $
    k_1=sqrt(2 m (E+V_0))/hbar, quad k_2=sqrt(2 m E)/hbar.
  $
  两侧波函数取为
  $
    psi_1=A exp(i k_1 x)+B exp(-i k_1 x), quad
    psi_2=C exp(i k_2 x).
  $
  由 $psi$ 与 $psi'$ 在原点连续，得
  $
    B/A=(k_1-k_2)/(k_1+k_2), quad C/A=2k_1/(k_1+k_2).
  $
  用概率流之比计算
  $
    R=((k_1-k_2)/(k_1+k_2))^2, quad
    T=4 k_1 k_2/(k_1+k_2)^2,
  $
  且 $R+T=1$。
]

#exercise(subname: [4.3])[
  电荷为 $q$ 的一维谐振子处在均匀电场 $cal(E)$ 中，势能为
  $
    V(x)=1/2 m omega^2 x^2-q cal(E) x.
  $
  求能量本征值和本征函数。
]
#solution[
  配方并令 $x_0=q cal(E)/(m omega^2)$，有
  $
    V(x)=1/2 m omega^2(x-x_0)^2-(q^2 cal(E)^2)/(2m omega^2).
  $
  因而
  $
    E_n=(n+1/2)hbar omega-(q^2 cal(E)^2)/(2m omega^2),
  $
  $
    psi_n(x)=sqrt(alpha/(sqrt(pi) 2^n n!))
    H_n(alpha(x-x_0)) exp(-alpha^2(x-x_0)^2/2),
  $
  其中 $alpha=sqrt(m omega/hbar)$。
]

#exercise(subname: [4.4])[
  求半谐振子势 $V(x)=infinity$（$x<0$）、$V(x)=m omega^2 x^2/2$（$x>0$）的能级与本征函数。
]
#solution[
  原点边界条件为 $psi(0)=0$，故只保留完整谐振子的奇宇称态。若以 $r=0,1,2,dots$ 编号，则
  $
    E_r=(2r+3/2)hbar omega,
  $
  且 $x>0$ 上的本征函数为完整谐振子 $n=2r+1$ 态的 $sqrt(2)$ 倍，$x<0$ 时为零。
]

#exercise(subname: [4.5])[
  在球坐标中写出轨道角动量算符，并导出 $L_z$。
]
#solution[
  由 $bold(L)=bold(r) times bold(p)=-i hbar bold(r) times nabla$ 及
  $
    nabla=bold(e)_r pdv(,r)+bold(e)_theta 1/r pdv(,theta)
    +bold(e)_phi 1/(r sin theta) pdv(,phi),
  $
  得
  $
    bold(L)=-i hbar (bold(e)_phi pdv(,theta)
    -bold(e)_theta 1/(sin theta) pdv(,phi)).
  $
  投影到 $z$ 轴即
  $
    L_z=-i hbar pdv(,phi).
  $
]

#exercise(subname: [4.6])[
  推导球坐标中的 $L^2$，并证明 $[L^2,L_z]=0$。
]
#solution[
  由 $nabla^2=1/r^2 pdv(,r)(r^2 pdv(,r))+1/r^2 nabla_Omega^2$，角向部分给出
  $
    L^2=-hbar^2 [1/(sin theta) pdv(,theta)(sin theta pdv(,theta))
    +1/(sin^2 theta) frac(partial^2,partial phi^2)].
  $
  该算符的系数与 $phi$ 无关，且各偏导可交换，因此它与 $L_z=-i hbar pdv(,phi)$ 对易。
]
