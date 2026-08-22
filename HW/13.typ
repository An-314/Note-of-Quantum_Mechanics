#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第13次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [15.1])[
  在双电子三重态基底 $|1,1 〉,|1,0 〉,|1,-1 〉$ 中写出总自旋分量 $S_z$。对单电子自旋基底作任意 $S U(2)$ 变换
  $
    U=mat(alpha,beta;-beta^*,alpha^*), quad abs(alpha)^2+abs(beta)^2=1,
  $
  求诱导到三重态空间的变换，并说明新表象中的 $S_z$。
]
#solution[
  旧表象中
  $
    S_z=hbar mat(1,0,0;0,0,0;0,0,-1).
  $
  把 $U|↑ 〉,U|↓ 〉$ 的张量积对称化，得到三重态中的自旋一表示
  $
    D^1(U)=mat(alpha^2,sqrt(2)alpha beta,beta^2;
    -sqrt(2)alpha beta^*,abs(alpha)^2-abs(beta)^2,sqrt(2)alpha^* beta;
    beta^(*)^2,-sqrt(2)alpha^* beta^*,alpha^(*)^2).
  $
  新表象中的矩阵由
  $
    S_z'=D^1(U)^dagger S_z D^1(U)
  $
  给出。这保留本征值 $hbar,0,-hbar$，但一般不再对角。
]

#exercise(subname: [15.2])[
  忽略自旋—轨道耦合以外的微扰。对 $l=1,s=1/2$ 的六维非耦合空间写出 $bold(L) dot bold(S)$ 的作用，并求其本征值与简并度。
]
#solution[
  利用
  $
    bold(L) dot bold(S)=L_z S_z+1/2(L_+S_-+L_-S_+)
  $
  可在任意指定的 $|m_l,m_s 〉$ 顺序中逐列写出矩阵。更直接地用
  $
    bold(L) dot bold(S)=1/2(J^2-L^2-S^2),
  $
  其中 $j=3/2,1/2$，得到本征值
  $
    hbar^2/2 quad (j=3/2,"四重"),
    quad -hbar^2 quad (j=1/2,"二重").
  $
  若微扰为 $H'=xi bold(L) dot bold(S)/hbar^2$，相应一级能移就是 $xi/2$ 与 $-xi$。
]

#exercise(subname: [15.3])[
  两个全同粒子处于一维无限深势阱 $0<x<a$ 中。分别对自旋为零的玻色子和自旋为 $1/2$ 的费米子，写出最低三条两粒子能量及空间本征函数。
]
#solution[
  记单粒子态为 $phi_n$，$epsilon_n=n^2 epsilon_1$，$epsilon_1=pi^2 hbar^2/(2m a^2)$。玻色子的最低三条能量为
  $
    2epsilon_1, quad epsilon_1+epsilon_2=5epsilon_1,
    quad 2epsilon_2=8epsilon_1,
  $
  空间波函数依次为 $phi_1 phi_1$、$(phi_1 phi_2+phi_2 phi_1)/sqrt(2)$、$phi_2 phi_2$。

  对自旋 $1/2$ 费米子，总波函数须反对称。最低能量 $2epsilon_1$ 只能配自旋单态；能量 $5epsilon_1$ 可由对称空间态配单态，或反对称空间态配三重态；下一条 $8epsilon_1$ 再由 $phi_2phi_2$ 配单态。各态的交换对称性由此完全确定。
]

#exercise(subname: [15.4])[
  证明任意纯的单电子自旋态都对应空间中的一个方向，使该方向上的自旋投影以概率一取得最大值 $hbar/2$。
]
#proof[
  任意归一化二分量旋量除去整体相位后都可写为
  $
    chi=mat(cos(theta/2);exp(i phi)sin(theta/2)).
  $
  令 $bold(n)=(sin theta cos phi,sin theta sin phi,cos theta)$，直接相乘得到
  $
    bold(n) dot bold(sigma) chi=chi.
  $
  所以 $bold(n) dot bold(S)$ 的测量值必为 $hbar/2$。这就是纯自旋态的 Bloch 球表示。
]

#exercise(subname: [15.5])[
  一维无限深势阱 $0<x<a$ 的粒子受到分段线性微扰
  $
    H'(x)=cases(2 h x/a & 0<x<a/2, 2 h(1-x/a) & a/2<x<a).
  $
  求基态能量的一级修正。
]
#solution[
  基态为 $psi_1=sqrt(2/a)sin(pi x/a)$。微扰关于 $a/2$ 对称，因此
  $
    E_1^((1))=2 integral_0^(a/2) 2/a sin^2(pi x/a) dot 2h x/a dif x.
  $
  积分得到
  $
    E_1^((1))=h(1/2+2/pi^2).
  $
]

#exercise(subname: [15.6])[
  在第一 Born 近似下，求球对称势的散射振幅，并用于有限球方势阱、指数势、屏蔽库仑势和三维 $delta$ 势。
]
#solution[
  令动量转移 $q=2k sin(theta/2)$，球对称势的一阶振幅为
  $
    f_B(theta)=-2m/(hbar^2 q) integral_0^infinity r V(r) sin(q r) dif r,
    quad dif sigma/dif Omega=abs(f_B)^2.
  $
  对 $V=-V_0$（$r<a$）有
  $
    f_B=2 m V_0/(hbar^2 q^3)[sin(q a)-q a cos(q a)].
  $
  对 $V=V_0 exp(-r/a)$ 有 $f_B=-4 m V_0/(hbar^2 a(q^2+a^(-2))^2)$；对 Yukawa 势 $V=g exp(-mu r)/r$ 有
  $
    f_B=-2m g/(hbar^2(q^2+mu^2)).
  $
  对 $V=g delta^3(bold(r))$ 应直接使用三维傅里叶公式，得 $f_B=-m g/(2 pi hbar^2)$（归一化约定改变时整体常数相应改变），散射各向同性。
]

#exercise(subname: [15.7])[
  从 Lippmann—Schwinger 方程推导三维散射波函数的远区形式，并写出精确到第一 Born 近似的散射振幅。
]
#solution[
  定态方程的出射格林函数为
  $
    G^+(bold(r),bold(r)')=-m/(2 pi hbar^2)
    exp(i k abs(bold(r)-bold(r)'))/abs(bold(r)-bold(r)').
  $
  当 $r$ 远大于势的作用范围时，
  $
    abs(bold(r)-bold(r)') approx r-hat(bold(r)) dot bold(r)',
  $
  因而
  $
    psi(bold(r)) tilde exp(i bold(k) dot bold(r))+f(theta,phi) exp(i k r)/r.
  $
  把积分方程中的精确波函数以入射平面波代替，得到第一 Born 振幅
  $
    f_B(theta,phi)=-m/(2 pi hbar^2)
    integral exp(-i bold(q) dot bold(r)')V(bold(r)') dif^3 r',
  $
  其中 $bold(q)=bold(k)'-bold(k)$。因此一阶 Born 振幅就是势的三维傅里叶变换。
]
