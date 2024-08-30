#import "@preview/physica:0.9.2": *
#import "@local/mytemplate:1.0.0": *

= 一维运动问题的一般分析

一维定态Schrödinger方程为：
$
- hbar^2/(2m) dd(""^2psi)/dd(x^2) + V(x) psi = E psi
$
写成二阶齐次线性微分方程的标准形式：
$
dd(""^2psi)/dd(x^2) + (2m)/hbar^2(E - V(x)) psi = 0
$
其中，$psi$是波函数，$V(x)$是势能函数，$E$是能量。

在经典力学的意义上$E = T + V$，因此$E - V >= 0$（*经典允许区*）。而在量子力学中，由于不确定关系，无法谈粒子“在某点处的动能”因此，$E - V < 0$（*经典禁戒区*）区的波函数仍然有非零解，粒子仍然会在这些区域出现。

方程重写为
$
1/psi dd(""^2psi)/dd(x^2) = -(2m)/hbar^2(E - V)
$
- 在经典允许区，$psi''$与$psi$符号相反，呈现振荡解。
- 在经典禁戒区，$psi''$与$psi$符号相同，呈现单调解。

== 一维定态Schrödinger方程解的性质

一维定态Schrödinger方程是：
$
dd(""^2psi)/dd(x^2) + (2m)/hbar^2(E - V(x)) psi = 0
$

=== Wronskian定理

设$psi_1$和$psi_2$是方程同一能量的两个解，且$Delta = W(psi_1, psi_2) = psi_1 psi_2' - psi_1' psi_2$是它们的Wronskian，则
$
W(psi_1, psi_2) = C 
$
其中$C$是常数。

_证明：_
$
&W(psi_1, psi_2)' \
=& psi_1 psi_2^'' - psi_1^'' psi_2 \
=& - psi_1 (((2m)/hbar^2)(E - V(x)) psi_2) + psi_2 (((2m)/hbar^2)(E - V(x)) psi_1) \
=& 0
$

当$Delta = 0$时, $phi_(1,2)(x)$ 是线性相关的，即它们只相差一个常数因子(同一个波函数)；当$Delta != 0$时，$phi_(1,2)(x)$ 是线性无关的。

=== 共轭定理

设$psi$是定态Schrödinger方程的一个解，$psi^*$是它的共轭，则$psi^*$也是也是该方程的解（且能量相同）。

=== 反射定理

设势能函数$U(x)$是关于原点对称的，即它满足
$
U(-x) = U(x)
$
那么若$psi(x)$是定态Schrödinger方程的解，则$psi(-x)$也是该方程的解（且能量相同）。

== 一维束缚态的一般性质

=== 一维定态的分类

*束缚态与非束缚态*

如果
$
psi(x) -> 0, x -> plus.minus oo
$
从而粒子在无穷远处出现的几率为零，那么这样的量子状态就称为*束缚态*，否则如
$
psi(x) -> plus.minus oo, x != 0
$
就称为*非束缚态*，或者称为*散射态*。

粒子处于束缚态还是非束缚态的另一判据：

假设 $U(x)$ 在 $x-> plus.minus oo$ 时有确定的极限，那么当
$
E < U(+oo, -oo)
$
时粒子处于束缚态，否则处于非束缚态。

=== 不简并定理
*简并度*：如果对一个给定的能量$E$，只有一个线性独立的波函数存在，则称该能级是非简并的，否则称它是简并的，其线性独立的波函数的个数称为它的*简并度*。

*不简并定理：一维束缚态必是非简并态*

_证明：_

设$psi_1$和$psi_2$是方程同一能量的两个解，且$Delta = W(psi_1, psi_2) = psi_1 psi_2' - psi_1' psi_2$是它们的Wronskian，则
$
W(psi_1, psi_2) = C = W(psi_1, psi_2) |_oo =  0
$
所以 
$
psi_1 = A psi_2
$

而这就表示 $psi_1,psi_2$ 代表相同的量子状态，所以它是非简并态。

注意:这个定理的两个前提*“一维”和“束缚态”*是缺一不可的。

波函数是复函数，可以写成下面的形式：
$
psi(x) = rho(x) e^(i theta(x))
$
其中，$rho(x)$为波函数的模，$theta(x)$为波函数的相位。

=== 一维束缚定态波函数的相位

*推论1：一维束缚定态波函数的相位必是常数*

_证明：显然。_

如果波函数$psi(x)$满足$psi(-x) = plus.minus psi(x)$，则称$psi(x)$有*正/负宇称*。
=== 宇称定理

*推论2（宇称定理）：如果 $U(-x) = U(x)$，则一维束缚定态波函数必有确定的宇称*

_证明：显然。_

束缚态（不只是一维束缚态）还有一个更重要的特征：它的能级是不连续地（离散地）变化的，即仅仅当取某些离散的数值时，定态Schrödinger方程才有符合单值、有限、连续条件的解。这就是通常意义的“量子化”。

#pagebreak(weak: true)

= ⼀维运动问题的解

== 一维自由粒子

一维自由粒子($U=0$)的定态方程：
$
- hbar^2/(2m)  dd(""^2psi)/dd(x^2)psi(x) = E psi(x)
$
解为：
$
psi(x) = c_1 e^(i/hbar p x) + c_2 e^(-i/hbar p x)
$
其中，$p = sqrt(2m E)$。

== 一维无限深方势阱

一维无限深方势阱的势能函数为：
$
U(x) = cases(
  0  &|x| < a,
  oo  &|x| > a
)
$
定态薛定谔方程的形式：
$
cases(
  - hbar^2/(2m) dd(""^2psi)/dd(x^2)psi(x) = E psi(x)  &|x| < a,
  - hbar^2/(2m) dd(""^2psi)/dd(x^2)psi(x) + U_0 psi(x) = E psi(x)  &|x| > a
)
$
#newpara()

阱外：$U=oo$要求$psi=0$粒子被束缚在势阱内——束缚态。

阱内：$U=0$，方程的形式类似于一维自由粒子。

解为：
$
psi(x) = cases(
  c_1 e^(i/hbar p x) + c_2 e^(-i/hbar p x)  &|x| < a,
  0  &|x| > a
)
$
其中，$p = sqrt(2m E)$。

$c_1, c_2$为待定常数，由波函数应满足的“单值、有限、连续”条件决定。“单值、有限”已经满足，下面看*连续条件*：
$
psi(-a) &= 0 &=> &c_1 e^(-i/hbar p a) + c_2 e^(i/hbar p a) = 0\
psi(a) &= 0 &=> &c_1 e^(i/hbar p a) + c_2 e^(-i/hbar p a) = 0
$

得到：
$
4i /hbar p a = 2 i pi n, n = 1, 2, 3, ...
$
得到能级：
$
E_n = p^2/(2m) = (n^2 pi^2 hbar^2) / (8 m a^2), n = 1, 2, 3, ...
$
称作*一维无限深势阱能量本征值*。其中$n$称为量子数，$n=1$代表基态，取其它值代表激发态。这表明，一维无限深方势阱中运动粒子的能量是量子化的。

得到波函数：
$
psi_n (x) = sqrt(1/a) sin(n pi (x+a)/(2a))
$
是*驻波形式*。如果加上时间项：
$
Psi_n (x)= sqrt(1/a) sin(n pi (x+a)/(2a)) e^(-i/hbar E_n t)
$
在势能无穷大的地方，一阶导数可以不连续，这个波函数就是例子。

*能级间隔*：
$
Delta E_n = E_(n+1) - E_n = (pi^2 hbar^2 (2n+1)) / (8 m a^2) prop 1/(m a^2)
$
$
Delta E_n/ E_n = (2n+1)/n^2 prop 1/n
$
宏观情况(能级变化$>>Delta E_n$)或量子数很大时$(n>>1)$，可认为能量连续。
*
最低能量（基态能量）*：$E_1 = pi^2 hbar^2 / (8 m a^2) >0$零点能——真空不空。

态的宇称是偶奇相间，基态为偶宇称。

波函数的节点数为 $n -1$。

*玻尔对应原理*：量子数$n$很大时，经典力学与量子力学的结果相近。

波函数的正交归一：
$
integral psi_n^* (x) psi_m (x) dd(x) = delta_(n,m)
$
其中，$delta_(n,m)$是Kronecker delta，当$n=m$时为1，否则为0。

波函数的完备性：
$
Psi(x) = sum_(n=1)^oo c_n psi_n (x)
$
其中，$c_n = integral psi_n^* (x) Psi(x) dd(x)$。

== 非对称无限深势阱

非对称无限深势阱的势能函数为：
$
U(x) = cases(
  0    & 0<x<a,
  oo  & x < 0 x> a
)
$
能级是：
$
E_n = (n^2 pi^2 hbar^2) / (2 m a^2), n = 1, 2, 3, ...
$

== 对称有限深方势阱

对称有限深方势阱的势能函数为：
$
U(x) = cases(
  0  & |x| < a,
  U_0    & |x| > a
)
$
定态薛定谔方程的形式：
$
cases(
   dd(""^2psi)/dd(x^2)psi(x) + k^2 psi(x) = 0  &|x| < a,\
    dd(""^2psi)/dd(x^2)psi(x) - alpha^2 psi(x) = 0  &|x| > a
)
$
其中，$k = sqrt(2m E/hbar^2)$，$alpha = sqrt(2m (U_0 - E)/hbar^2)$。

有限解为
$
psi(x) = cases(
  C e^(alpha x) &x<-a,\ 
  A cos(k x) + B sin(k x)  &|x| < a,\
  D e^(-alpha x)  &x > a
)
$

1. 偶宇称解

$
psi(x) = psi(-x) , B = 0, C = D
$
在$x=a$处$psi$和$psi'$连续，得到
$
k tan(k a) = alpha
$

2. 奇宇称解

$
psi(x) = -psi(-x) , A = D = 0
$
在$x=a$处$psi$和$psi'$连续，得到
$
k cot(k a) = -alpha
$

#newpara()

采用图解法，令
$
xi = k a, eta = alpha a (xi, eta > 0)
$
得到
$
cases(
  xi tan(xi) = plus.minus eta,
  xi^2+eta^2 = (2 mu U_0 a^2) / hbar^2
)

$

#figure(
  image("pic/2024-03-20-13-39-33.png", width: 40%),
  caption: [
    对称有限深方势阱的解
  ],
)

找出这两族曲线的交点，记交点的$xi$值为$xi_1,xi_2,...,$则能级就是
$
E_n = (hbar^2 xi_n^2) / (2 m a^2)
$

_讨论:_
- 能级的宇称是偶奇相间，最低的能级是偶宇称
- $0 < xi_1 < pi/2 < xi_2 < pi < xi_3 < ...$。每个能级都比无限深势阱的相应能级低一些。

    $U_0 -> oo $时，$xi_n -> n pi/2$，能级趋近于无限深势阱的能级。在$|x|>=a$处，波函数趋于0。
- 不论势阱多浅或多窄，至少存在一个束缚态，并且宇称为偶。
- 对于偶宇称，当
    $
    eta^2 + xi^2 = (2 mu U_0 a^2) / hbar^2 >= pi^2
    $
    时，才能出现第一个偶宇称的激发态。
- 对于奇宇称，当
    $
    eta^2 + xi^2 = (2 mu U_0 a^2) / hbar^2 >= (pi/2)^2
    $
    时，才能出现第一个奇宇称的激发态。
- 在给定的势阱中，能级的个数是
    $
    ceil((8 mu U_0 a^2) / (pi^2 hbar^2))
    $
== $delta$势阱

考虑$delta$势阱的势能：
$
V(x) = - gamma delta(x)
$
其中$gamma > 0$。

定态Schrödinger方程：
$
dd(""^2psi)/dd(x^2) + (2m)/hbar^2(E + gamma delta(x)) psi(x) = 0
$
在$x≠0$处：
$
dd(""^2psi)/dd(x^2) + (2m)/hbar^2 E psi(x) = 0
$
结合束缚态的要求，得到一般解为：
$
psi(x) = cases(
  C_1 e^(beta x) &x<0,
  C_2 e^(-beta x) &x>0
)
$
其中
$
beta = sqrt(2m |E|)/hbar
$
且在$x=0$处$psi$连续，得到：
$
C_1 = C_2
$
归一化条件：
$
integral |psi(x)|^2 dd(x) &= 1
$
令$beta = 1/L$，得到：
$
psi(x) = 1/sqrt(L) e^(-(|x|)/L)
$
$L$被称为特征长度。

进一步地，对Schrodinger方程两边积分，得到：
$
dd(psi)/dd(x) |_0^0 &= - (2m)/hbar^2 gamma psi(0)\
psi'(0^+) - psi'(0^-) &= - (2m)/hbar^2 gamma psi(0)\
-2 beta &= - (2m gamma)/hbar^2 
$
得到束缚态基态的能量：
$
E = - (gamma^2 m)/(2 hbar^2)
$
束缚态能级有且只有一个。


== 线性谐振子

线性谐振子的势能函数是：
$
U(x) = 1/2 mu omega^2 x^2
$
其中，$omega$是谐振子的角频率。

定态Schrödinger方程：
$
- hbar^2/(2m) dd(""^2psi)/dd(x^2) + 1/2 mu omega^2 x^2 psi = E psi
$
做如下的无量纲化变换：
$
cases(
  xi = sqrt((mu omega)/hbar) x = alpha x,\
  lambda = (2 E)/(hbar omega)
)
$
得到：
$
dd(""^2psi)/dd(xi^2) + (lambda - xi^2) psi = 0
$
束缚态的解的要求：
$
psi(xi) -> 0, xi -> plus.minus oo
$

在$xi -> plus.minus oo$时，方程近似为：
$
dd(""^2psi)/dd(xi^2) - xi^2 psi = 0
$
有近似解：
$
psi(xi) = e^(-xi^2/2)
$
进行变换：
$
psi(xi) = e^(-xi^2/2) H(xi)
$
得到*Hermite方程*：
$
dd(""^2H)/dd(xi^2) - 2 xi dd(H)/dd(xi) + (lambda - xi^2) H = 0
$

可以用级数解法求解Hermite方程，得到*Hermite多项式*。
$
H(xi) = sum_(n=0)^oo c_n xi^n
$
带入方程得到递推关系：
$
c_(k+2) = (2k + 1 - lambda)/((k+1)(k+2)) c_k
$
可以分析：
$
c_k tilde.op 1/(k/2)!, H(xi) tilde.op~ e^(xi^2)
$
这会使得
$
psi(xi) tilde.op~ e^(xi^2/2)
$
发散！这是不可能的，所以*级数必须截断*，即存在一个最大的$n$，使得$c_n = 0$。*这就对$lambda$有了限制：*
$
lambda = 2n + 1, n = 0, 1, 2, ...
$
首先得到了能量本征值：
$
E_n = (n + 1/2) hbar omega, n = 0, 1, 2, ...
$
这就是*线性谐振子的能级*。

*Hermite方程对于本征值$lambda$的多项式解就是Hermite多项式*。
$
H_0 (xi) &= 1,\
H_1 (xi) &= 2 xi,\
H_2 (xi) &= 4 xi^2 - 2,\
$
一般形式：
$
H_n (xi) &= (-1)^n e^(xi^2) dd(""^n )/dd(xi^n)e^(-xi^2)
$
正交归一：
$
integral e^(-xi^2) H_n (xi) H_m (xi) dd(xi) = sqrt(pi) 2^n n! delta_(n,m)
$
正交归一化常数为：
$
N_n = sqrt(alpha/(2^n n! sqrt(pi)))
$

#newpara()

得到的波函数为：
$
psi_n (x) = N_n e^(-xi^2/2) H_n (xi) = N_n H_n (alpha x) e^(-alpha^2 x^2/2)
$
$n$的奇偶性决定了谐振子波函数的奇偶性，即宇称的奇偶：
$
psi_n (-x) = (-1)^n psi_n (x)
$

#figure(
  image("pic/2024-03-20-14-34-15.png", width: 50%),
  caption: [
    粒子具有一定的概率处在经典允许区之外
  ],
)

经典情况下：
$
xi = a sin(omega t + delta)\
v = a omega cos(omega t + delta) = a omega sqrt(1 - xi^2/a^2)\
w(xi) prop 1/v prop 1/sqrt(1 - xi^2/a^2)
$

#figure(
  image("pic/2024-03-20-14-36-21.png", width: 40%),
  caption: [
    实线(虚线)表示谐振子在量子(经典)情况下的概率密度分布
  ],
)

_讨论：_
- $n$较大时能谱连续：
    $
    (Delta E_n)/E_n = 1/(n + 1/2)  ->0
    $
- $n$较小时，量子效应明显。比如在$n=0$时，量子力学告诉我们粒子在$x=0$处出现几率最大，但在经典力学中，粒子在此处速度极大，出现几率反而最小。
- $n$较大时，空间几率趋于经典分布（量子态$n$对应的典运动振幅为$a = sqrt(2n  +1)$。
- 在$xi>a$的经典力学粒子不能达到的区域，量子世界中仍有分布
- 同势阱问题一样，振子*基态能量(零点能)*不为0($E_0 = 1/2 hbar omega$)。零点能的实验现象：Casimir效应。
- 由于谐振子势有空间反射不变性，所以有确定的宇称。能级的宇称偶奇相间，基态是偶宇称。
- $psi_n (x)$有$n$个节点。

== 一维散射问题

为简单可以假设：
$
U(+ oo) = U(- oo) = 0, E>0
$
所以这时的量子状态是非束缚态，也就是散射态。

在$x -> plus.minus oo$时，$U = 0$
$
dd(""^2psi)/dd(x^2) + k^2 psi = 0
$
其中$k = sqrt((2m E)/hbar^2)$。

解为：
$
psi(x) = c_1 e^(i/hbar p x) + c_2 e^(-i/hbar p x)
$

在$psi = A e^(i/hbar p x)$时，几率流密度为：
$
j = (hbar p)/(2mu) (psi dd(psi^*)/dd(x) - psi^* dd(psi)/dd(x)) = (hbar p)/(mu) |psi|^2 = |A|^2 nu
$
其中，$nu = (hbar k)/mu = p/mu$是粒子的经典速度。

*散射情况*是：粒子从一边入射，被势场散射而分成了反射和透射两个部分。这给方程提出了一定的条件。以左方入射为例，边界条件是：
$
cases(
  psi(x) =A e^(i/hbar p x) + B e^(-i/hbar p x) & "入射加反射" &x -> -oo,\
  psi(x) = C e^(i/hbar p x) & "只有透射" &x -> +oo
)
$
由此定义*反射系数*和*透射系数*：
$
R = J_R/J_I = (|B|^2)/(|A|^2),\
D = J_T/J_I = (|C|^2)/(|A|^2)
$

#newpara()

*方势垒的穿透*

方势垒的势能函数为：
$
U(x) = cases(
  U_0  &0 < x < a,
  0  &x < 0 or x>a
)
$
定态Schrödinger方程：
$
cases(
  dd(""^2psi)/dd(x^2)psi(x) + (2m)/hbar^2(E - U_0) psi(x) = 0  &(0 < x < a),\
  dd(""^2psi)/dd(x^2)psi(x) + (2m)/hbar^2(E) psi(x) = 0  &(x < 0 or x>a)
)
$
解为：
$
psi_1 &= A e^(i k_1 x) + A' e^(-i k_1 x)  &(x < 0),\
psi_2 &= B e^(i k_2 x) + B' e^(-i k_2 x)  &(0 < x < a),\
psi_3 &= C e^(i k_1 x)  &(x > a)
$
其中，$k_1 = sqrt(2m E/hbar^2), k_2 = sqrt(2m (E - U_0)/hbar^2)$。表达式中第一项(第二项)代表从左向右(从右向左)传播的平面波。在$x>a$的区域只有向右的透射波，所以$C'=0$。

利用波函数及其导数在边界处的连续性，得出波函数中各系数的关系式：
$
psi_1 (0) = psi_2 (0), psi_1 ' (0) = psi_2 ' (0), psi_2 (a) = psi_3 (a), psi_2 ' (a) = psi_3 ' (a)
$
得到：
$
A' = (2i(k_1^2 - k_2^2) sin(k_2 a))/((k_1 - k_2)^2 e^(i k_2 a) - (k_1 + k_2)^2 e^(-i k_2 a)) A,\
C = (4k_1 k_2 e^(- i k_1 a))/((k_1 + k_2)^2 e^(- i k_2 a) - (k_1 - k_2)^2 e^(i k_2 a)) A
$
计算入射波($A$项)，反射波($A'$项)和透射波($C$项)的概率流密度得:
$
J_I = (hbar k_1)/mu |A|^2,\
J_R = (hbar k_1)/mu |A'|^2,\
J_D = (hbar k_1)/mu |C|^2
$
接着可以得到粒子的*反射系数和透射系数*：

$
R = abs(J_R)/abs(J_I) = ((k_1^2 - k_2^2)^2 sin^2(k_2 a))/((k_1^2-k_2^2)^2 sin^2(k_2 a) + 4k_1^2 k_2^2 ),\
D = abs(J_D)/abs(J_I) = (4k_1^2 k_2^2)/((k_1^2-k_2^2)^2 sin^2(k_2 a) + 4k_1^2 k_2^2 )
$

很容易看出$R+D=1$，这是粒子流守恒的必然结果。

$E<U_0$情形下,$k_2$是虚数,令$k_2 = i k_3 $,则
$
D = (4k_1^2 k_3^2)/((k_1^2+k_3^2)^2 sinh^2(k_3 a) + 4k_1^2 k_3^2 ) > 0
$
当$U_0$很大时时，$k_3a>>1$，透射系数可近似为:
$
D &= D_0 (k_1, k_3) e^(-2k_3 a)\
  &= D_0 (k_1, k_3) e^(-2 sqrt(2mu(U_0 - E)/hbar^2) a)
$
$E<U_0$时$D>0$，这就是*量子隧穿效应*。

#figure(
  image("pic/2024-03-22-13-32-46.png", width: 20%),
  caption: [
    量子隧穿效应
  ],
)

量子隧穿效应随着$a$的上升指数下降，在宏观尺度上是充分小的。理解量子隧道效应 – 在势垒内部，根据总能量守恒，粒子的动能将变为负数：
$
E_"kin" = E - U < 0
$
这在经典力学中是不可能的。

在量子力学里，粒子动能与势能不能同时具有确定值。而且，力学量的平均值是一个全域的积分平均，在某个局域内讨论是没有意义的。一旦讨论限制在某一局域(势垒内部)，粒子动量就在某一范围内不确定。

== 方势阱的共振隧穿

如果把方势垒改为方势阱，整个实轴都是经典允许区不存在势垒的隧穿，发生共振隧穿。

对于势阱情况下的透射系数：
$
T= (4k_1^2 k'^2)/((k_1^2-k'^2)^2 sin^2(k' a) + 4k_1^2 k'^2 )
$
其中$k' = sqrt(2m (V_0 + E)/hbar^2)$。

共振透射，共振能量的位置：
$
E = - V_0 + (n pi hbar)^2 / (2m a^2), n = 1, 2, 3, ...
$
以$-V_0$为势能起点宽度为$a$的无线深势阱中的能级表达式。

#figure(
  image("pic/2024-03-22-13-52-52.png", width: 80%),
  caption: [
    共振隧穿
  ],
)