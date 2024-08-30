#import "@preview/physica:0.9.2": *
#import "@local/mytemplate:1.0.0": *

= 波函数和薛定谔方程

== 波函数及其统计解释
=== 波函数和波粒二象性
==== 波函数

#figure(
  table(
  columns: (auto, auto, auto),
  inset: 10pt,
  align: horizon,
  [],[*保留经典概念的哪些特征*],[*不具有经典概念的哪些特征*],
  [粒子性],[有确定的质量、电荷、自旋等],[没有确定的轨道],
  [波动性],[有干涉、衍射等现象],[振幅不直接可测]
),
  caption: [
    波粒二象性
  ],
)


对于微观粒子或量子力学体系，可以用与时间和空间相关的复函数来描写：
$
psi = psi(arrow(r), t)
$

被称为波函数。

Schrodinger方程所描述的波函数，并不是象经典波那样的代表实在的物理量的波动，而是刻画粒子在空间的几率分布的几率波

$
|psi(arrow(r), t)|^2
$

波函数在某点的强度（绝对值的平方）与在该点找到粒子的*几率密度*成正比。波函数本身称为*几率振幅*。

例如电子的双缝干涉：
$
P_(12) &= |psi_1 + psi_2|^2 = |psi_1|^2 + |psi_2|^2 + (psi_1^*psi_2 + psi_2^*psi_1)\
&= P_1 + P_2 + 2sqrt(P_1P_2)cos(phi)
$
电子呈现的波动性反应了微观客体运动的一种统计特性，因此波函数也被称为概率波幅。
==== 波函数的归一

可以看到$psi(arrow(r))$和$c psi(arrow(r))$描述的相对几率分布是完全相同的。

对于空间中两点$arrow(r_1)$和$arrow(r_2)$，粒子出现在这两点的几率之比为：
$
(|psi(arrow(r_1))|^2)/(|psi(arrow(r_2))|^2) = (|c psi(arrow(r_1))|^2)/(|c psi(arrow(r_2))|^2)
$

对于$psi(arrow(r))$是某个波函数，按照几率解释，在点$(arrow(r),t)$附近的体积元 $d tau$ 中发现粒子的几率是：
$
dd(W(arrow(r),t)) = |psi(arrow(r),t)|^2 dd(tau)
$
粒子的空间几率密度是：
$
w(arrow(r),t) = dd(W(arrow(r),t)) / dd(tau) = |psi(arrow(r),t)|^2
$
因此在全空间发现粒子的几率是：
$
W(t) = integral w(arrow(r),t) dd(tau) = integral |psi(arrow(r),t)|^2 dd(tau)
$
根据波函数的统计诠释，很自然地要求粒子不产生，不湮灭，取归一化:
$
W(t) = 1
$
即：
$
integral |psi(arrow(r),t)|^2 dd(tau) = 1
$
这就是波函数的归一条件。

_注：_
- _即使要求波函数是归一的，它仍然有一个整体的相位因子$e^(i theta)$（$theta$为实常数）不能确定。_
- _如果积分_
  $
  integral |psi(arrow(r),t)|^2 dd(tau)
  $
  _是无穷大，相当于粒子的运动范围没有限制，粒子可以达到无穷远处。【正交归一化为delta函数】这样的波函数就是不能（有限地）归一的例如平面波（自由粒子的de Broglie波），此时_
  $
  |psi(arrow(r),t)|^2
  $
  _代表“相对几率密度”。_

==== 波的相干性

两个独立无关的粒子的波函数，由于各自存在相位不定性，所以两个粒子之间没有干涉（相干性coherence）。

例如，两束自然光之间是没有相干性的。为了具有相干性，可以采用同一个光源的波前分割，或者振幅分割的方法。或者利用量子力学原理，使得不同光源之间具有强相干性，例如激光。

两个不相干粒子的整体波函数不是写成相加的数学形式，而是相乘，即所谓的*直积态*。
$
Psi(x_1, x_2) = psi(x_1)Phi(x_2)
$
则
$
(P(x_1 = 1, x_2))/(P(x_1 = 2, x_2)) = (|psi(1)|^2)/(|psi(2)|^2) (|Phi(x_2)|^2)/(|Phi(x_2)|^2) = (|psi(1)|^2)/(|psi(2)|^2)
$
即两个粒子的几率分布是独立的。

N个粒子的系统，对于一个由N个粒子构成的系统，系统的波函数是N个粒子的坐标和时间的复函数：
$
Psi(arrow(r_1), arrow(r_2), ..., arrow(r_N), t)
$

==== 态叠加原理

如果$Psi_(1(2))$是体系的可能状态，那么
$
Psi = c_1 Psi_1 + c_2 Psi_2
$
也是体系的可能状态，其中$c_1$和$c_2$是复数。
$
|Psi|^2 = |c_1 Psi_1 + c_2 Psi_2|^2 = |c_1|^2 |Psi_1|^2 + |c_2|^2 |Psi_2|^2 + (c_1^*c_2 Psi_1 Psi_2 + c_2^*c_1 Psi_2 Psi_1)
$
其中$c_1^*c_2 Psi_1 Psi_2 + c_2^*c_1 Psi_2 Psi_1$是干涉项。

一般地说，叠加原理可以写成：
$
Psi = sum_n c_n Psi_n
$
对于一个指定的量子体系，如果我们找到了它的“完备的基本状态”，例如：
$
{Psi_n(n = 1, 2, 3, ...)}
$
那么任何一个体系的波函数都可以用这些基本状态的线性组合来表示，处于状态$n$的几率是$|c_n|^2$。

由于量子态满足叠加原理，一个量子系统的全部状态构成了线性空间，称为这个量子系统的Hilbert空间。

这里补充一些重要的公式：
- Fourier变换：
$
f(x) = 1/sqrt(2pi) integral_(-oo)^oo F(k) e^(i k x) dd(k)
$
- 逆变换：
$
F(k) = 1/sqrt(2pi) integral_(-oo)^oo f(x) e^(-i k x) dd(x)
$
- $delta$函数的Fourier变换：
$
delta(x) = 1/(2pi) integral_(-oo)^oo e^(i k x) dd(k)
$

_例：*动量几率分布*_

_一个自由粒子以动量 $arrow(p)$ 和能量_
$
E = p^2/(2mu)
$
_运动的状态是平面波：_
$
psi_arrow(p)(arrow(r), t) = e^(i/hbar (arrow(p)dot arrow(r) - E t))
$
_根据叠加原理，任何的波函数(不一定是自由粒子)都可以写成：_
$
psi(arrow(r), t) = sum_arrow(p) c_arrow(p) psi_arrow(p)(arrow(r), t)
$
_即是各种不同动量的平面波的叠加。这个例子在数学上就是函数的Fourier变换 – 任何函数（在分段光滑、绝对可积等条件下）都可以展开为正、余弦函数的线形叠加。_

_定义：_
$
psi_arrow(p) (arrow(r)) = sqrt(1/(2 pi hbar)^3) e^(i/hbar (arrow(p)dot arrow(r)))
$

_那么任何波函数（不一定是自由粒子的）都可以写成 :_

$
psi(arrow(r), t) = integral c(arrow(p), t) sqrt(1/(2 pi hbar)^3) psi_arrow(p) (arrow(r))  dd(""^3arrow(p))
$

_其中的展开系数由下式得出（Fourier变换）：_

$
c(arrow(p), t) = integral psi(arrow(r), t) sqrt(1/(2 pi hbar)^3) psi_arrow(p)^* (arrow(r)) dd(""^3arrow(r))
$

_$c(arrow(p), t)$的物理意义是“动量测量几率振幅”。_

_Fourier变换相当于作频谱分析，代表在$psi(arrow(r), t)$中出现动量为$arrow(p)$能量为$E$的平面单色波的概率为$|c(arrow(p), t)|^2$。_

_$|c(arrow(p), t)|^2$代表粒子在$t$时刻具有动量$arrow(p)$的概率密度，或者说$psi(arrow(p), t)$中包含平面波$e^(i/hbar arrow(p) dot arrow(r))$的成分有多少。_

_可以证明归一化：_

$
integral |c(arrow(p), t)|^2 dd(""^3arrow(p)) = integral |psi(arrow(r), t)|^2 dd(""^3arrow(r))
$

_对于一维情况，Fourier变换为：_

$
psi(x , t) = 1/sqrt(2 pi hbar) integral c(p, t) e^(i/hbar p x) dd(p)\
c(p, t) = 1/sqrt(2 pi hbar) integral psi(x, t) e^(-i/hbar p x) dd(x)
$

事实上$ψ(r,t)$是波函数在*坐标表象*的表示形式，$c(p,t)$是波函数在*动量表象*的表示形式。它们在描述粒子状态方面是等价的，只是具体描述方式不同。

== 薛定谔方程

薛定谔方程应满足:
1. 线性齐次 —— 波函数线性叠加原理。
2. 方程系数不应包含如动量,能量等系统参量 —— 方程的解应描述所有可能的粒子状态。

量子力学的基本定律是波函数所满足的偏微分方程。这个基本定律在本质上是一个假说。

=== 从De Broglie波到薛定谔方程

de Broglie波 $psi(arrow(r), t) = e^(i/hbar (arrow(p)dot arrow(r) - E t))$ 有等式：
$
diff/(diff t) psi = -i/hbar E psi
$
$
nabla psi = i arrow(p)/hbar  psi
$
$
nabla^2 psi = -p^2/hbar^2 psi
$

这可以看作是在经典能量-动量关系：
$
E = p^2/(2m)
$
中进行代换：
$
cases(
  E -> i hbar diff/(diff t),
  p -> -i hbar nabla
)
$
并且把它们作用于波函数得到的。

推广到粒子在外势场 $U(arrow(r))$ 中运动，其能量的表达式为：
$
E = p^2/(2m) + U(arrow(r))
$
得到单粒子运动的*薛定谔(Schrodinger)方程*：
$
i hbar diff/(diff t) psi = -hbar^2/(2m) nabla^2 psi + U(arrow(r)) psi
$

_例：相对论粒子所满足的Klein-Gorden方程_
静止质量为$m_0$的粒子自由运动时，能量-动量关系为：
$
E^2 = p^2c^2 + m_0^2c^4
$
因此，相对论波动方程为：
$
(c^2 hat(p)^2 + m_0^2c^4) psi = hat(E)^2 psi
$
可以描述相对论性玻色子的运动。但该方程可能解出负能量的解，有不可逾越的困难。

=== 几率守恒定律

在*非相对论*情况下，实物粒子没有产生和湮灭的现象，所以随时间演化的过程中，*粒子数目守恒*，即*几率守恒*。

粒子的空间几率密度是
$
w(arrow(r),t) = |psi(arrow(r),t)|^2 = psi^*psi
$
可以计算时间的变化率：
$
diff/(diff t) w = diff/(diff t) (psi^*psi) = (diff/(diff t) psi^*) psi + psi^* (diff/(diff t) psi)
$
利用Schrodinger方程：
$
diff/(diff t) psi &= (i hbar)/(2m) nabla^2 psi - i/hbar U psi\
diff/(diff t) psi^* &= (-i hbar)/(2m) nabla^2 psi^* + i/hbar U psi^*
$
认为势能是实函数，从而：
$
diff/(diff t) w = (i hbar)/(2m) (nabla^2 psi^* psi - psi^* nabla^2 psi)
$
设*几率流密度*为：
$
arrow(J) = (i hbar)/(2m) (psi^* nabla psi - psi nabla psi^*)
$
则：
$
diff/(diff t) w + nabla dot arrow(J) = 0
$
即*几率守恒定律*。

$arrow(J)$是“几率流密度”，而上式表现了几率守恒。将体积$V$扩展到全空间，在非相对论情况下，实物粒子没有产生和湮灭的现象，所以在随时间演化的过程中，粒子数目将保持不变。对于一个粒子而言，在全空间找到它的概率应不随时间改变，即：
$
dd("")/dd(t) integral_oo | psi(arrow(r),t) |^2 dd(tau) = - integral.cont_oo arrow(J) dot arrow(n) dd(S) = 0
$
*几率守恒也就是粒子数守恒。*

一个粒子既不可能凭空产生，也不可能凭空消失。当粒子在空间某处的概率减少时，必然在空间某个地方的概率增加了。

=== 波函数应满足的条件

从波函数的几率解释以及波函数满足二阶微分方程这一要求，一般地说，波函数应该满足以下三个条件：

+ *单值性*（波函数是处处单值的）
+ *有限性*（波函数是处处有限的）
+ *连续性*（波函数及其一阶导数是处处连续的）：由于有波函数对时间和空间的微商，波函数应在时空坐标$(arrow(r),t)$上连续，并且对空间的一阶微商$(nabla Y)$也连续，否则方程中空间二阶导数项发散

但是连续性允许有例外：在势能有无穷大跳跃的地方，波函数的一阶导数可以是不连续的

== 定态薛定谔方程

若势能$U$不显含时间，则薛定谔方程为：
$
i hbar diff/(diff t) psi = (-hbar^2/(2m) nabla^2 + U(arrow(r))) psi(arrow(r),t)
$
记Hamilton算符为：
$
hat(H) = -hbar^2/(2m) nabla^2 + U(arrow(r))
$
时间算符$i hbar diff/(diff t)$与Hamilton算符$hat(H)$都是*能量算符*。

分离变量法：
$
Psi(arrow(r),t) = psi(arrow(r)) f(t)
$
代入薛定谔方程得：
$
i hbar psi(arrow(r)) (dd(f(t)))/(dd(t)) &= f(t) hat(H) psi(arrow(r))\
(i hbar)/(f(t)) (dd(f(t)))/(dd(t)) &= (hat(H) psi(arrow(r)))/(psi(arrow(r))) = E
$
得到*本征方程*：
$
hat(H) psi(arrow(r)) = E psi(arrow(r))
$
和$f(t)$满足的方程：
$
(dd(f(t)))/(dd(t)) &= -i E f(t)/hbar\
f(t) &= C e^(-i/hbar E t)
$
为时间震动因子，$C$为常数。

因此，定态薛定谔方程的解为：
$
psi(arrow(r)) = sum C_n psi_n (arrow(r)) e^(-i/hbar E_n t)
$


#newpara()

*本征方程*就是*定态薛定谔方程*：
$
E psi(arrow(r)) = (-hbar^2/(2m) nabla^2 + U(arrow(r))) psi_n (arrow(r))
$

多粒子系统的Hamilton算符为：
$
hat(H) = sum_(i = 1)^N (-hbar^2/(2m_i) nabla_i^2 + U_i(arrow(r_i))) + V(arrow(r)_1, arrow(r)_2, ..., arrow(r)_N)
$
其中$V(arrow(r)_1, arrow(r)_2, ..., arrow(r)_N)$是粒子间相互作用势能。