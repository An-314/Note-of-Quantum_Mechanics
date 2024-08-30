#import "@preview/physica:0.9.2": *
#import "@local/mytemplate:1.0.0": *

= ⼒学量算符与波函数

== 量子力学的基本公设

+ 公设1：微观体系的状态由波函数描述，波函数满足单值、有限、连续条件
+ 公设2：波函数的动力学演化满足薛定鄂方程
+ 公设3：力学量用*厄密算符*表示，且有组成*完备集的本征函数系*
+ 公设4：任一波函数可以展开为力学量算符本征函数的线性叠加，测得力学量为本征值$lambda_n$的几率（密度）为展开式中对应本征函数系数的模方$|c_n|^2$#footnote[意思是测量的只能是本征值，而且是又概率的]

== 力学量的算符表示

在量子力学中力学量有完全不同于经典力学的表示方法，这就是用算符表示：
$
"基本的力学量算符" <=> "数学上的函数变换、算子"
$
算符就是可以作用于波函数把它变成另一个函数的运算。

代表力学量$F$的算符是$hat(F)$。

量子力学中基本的力学量算符是：

- *动量算符*
    $
    hat(arrow(p)) = -i hbar nabla
    $
    $
    hat(p)_x - i hbar (diff)/(diff x)
    $
- *位置算符*
    $
    hat(arrow(r)) = arrow(r)
    $
其它的力学量算符按下列规则来构成：若在经典力学中力学量$F$用坐标和动量表示出的关系式是
$
F = f(arrow(r), arrow(p))
$
则在量子力学中$F$的算符表示是
$
hat(F) = f(hat(arrow(r)), hat(arrow(p))) = f(arrow(r), -i hbar nabla)
$
$f$代表相同的关系函数。

总能量（动能加势能）在分析力学中称为Hamiltonian（哈密顿量），记为H。对于单粒子，
$
H = T + U = (arrow(p)^2)/(2m) + U(arrow(r))
$
对应的*Hamilton算符* ：
$
hat(H) = hat(p)^2/(2m) + V(arrow(r)) = -hbar^2/(2m) nabla^2 + V(arrow(r))
$

轨道角动量的经典表达式是
$
arrow(L) = arrow(r) crossproduct arrow(p)
$
对应的$L$*角动量算符*是
$
hat(L) = hat(arrow(r)) crossproduct hat(arrow(p)) = - i hbar arrow(r) crossproduct nabla
$

更准确地说，上面所定义的算符应该称作是“坐标表象”中的算符。用算符来代替经典力学中的力学量，是把经典力学模型“量子化”的步骤的重要部分。

在量子力学中有一些量是没有经典力学的对应物的，比如宇称和自旋角动量。那时我们就要直接从量子力学的分析出发来引进它们的算符。

== 不同坐标系下的微分算符表示

=== 直角坐标系

在量子力学中，我们经常需要在不同坐标系下表示微分算符。

在直角坐标系下，微分（梯度）算符是#footnote[有$partial f = f' + f partial$]
$
nabla = (partial)/(partial x) arrow(i) + (partial)/(partial y) arrow(j) + (partial)/(partial z) arrow(k) = partial_x arrow(i) + partial_y arrow(j) + partial_z arrow(k)
$
Laplacian算符是
$
nabla^2 = partial^2/(partial x^2) + partial^2/(partial y^2) + partial^2/(partial z^2)
$

=== 柱坐标系

在柱坐标系下，微分算符是
$
nabla = (partial)/(partial r) arrow(r) + 1/r (partial)/(partial theta) arrow(theta) + (partial)/(partial z) arrow(z)
$
Laplacian算符是
$
nabla^2 &= partial^2/(partial r^2) + 1/r (partial)/(partial r) + 1/r^2 partial^2/(partial theta^2) + partial^2/(partial z^2)\
&= 1/r partial/(partial r) (r partial/(partial r)) + 1/r^2 partial^2/(partial theta^2) + partial^2/(partial z^2)
$
并且有Jacobi行列式
$
dd(x)dd(y)dd(z) = r dd(r) dd(theta) dd(z)
$

=== 球坐标系

在球坐标系下，微分算符是
$
nabla = (partial)/(partial r) arrow(r) + 1/r (partial)/(partial theta) arrow(theta) + 1/(r sin(theta)) (partial)/(partial phi) arrow(phi)
$
Laplacian算符是
$
nabla^2 &= partial^2/(partial r^2) + 1/r (partial)/(partial r) + 1/r^2 partial^2/(partial theta^2) + 1/(r^2 sin(theta)^2) partial^2/(partial phi^2)\
&= 1/r^2 partial/(partial r) (r^2 partial/(partial r)) + 1/(r^2 sin(theta)) partial/(partial theta) (sin(theta) partial/(partial theta)) + 1/(r^2 sin(theta)^2) partial^2/(partial phi^2)
$
并且有Jacobi行列式
$
dd(x)dd(y)dd(z) = r^2 sin(theta) dd(r) dd(theta) dd(phi)
$

*定义*
$
hat(arrow(Y)) = hat(r) crossproduct nabla
$
其在球坐标系的表示是
$
hat(arrow(Y)) = partial/(partial theta) arrow(phi) - 1/sin(theta) partial/(partial phi) arrow(theta)
$
且有
$
hat(Y)^2 &= partial^2/(partial theta^2) + 1/sin(theta) partial/(partial theta) + 1/(sin(theta)^2) partial^2/(partial phi^2)\
&= 1/sin(theta) partial/(partial theta) (sin(theta) partial/(partial theta)) + 1/(sin(theta)^2) partial^2/(partial phi^2)
$
于是有
$
nabla^2 = 1/r^2 (hat(Y)^2 + partial/(partial r)(r^2 partial/(partial r)))
$

*角动量算符：*
$
hat(arrow(L)) = hat(arrow(r)) crossproduct hat(arrow(p)) = -i hbar hat(Y)
$
$
hat(L)^2 = - hbar^2 hat(arrow(Y))^2
$

== 算符的一般性质和运算规则

量子力学中的算符，代表着对波函数(量子态)的一种运算(或操作)。

=== 线性算符

算符$hat(A)$是线性的，就是说对于任意两个波函数$psi_1$和$psi_2$和任意两个复数$a$和$b$，有
$
hat(A)(a psi_1 + b psi_2) = a hat(A) psi_1 + b hat(A) psi_2
$
例如：$hat(arrow(p)) = -i hbar nabla$是线性算符。

描述可观测量的算符都是线性算符，这是态叠加原理的体现。

=== 单位算符

单位算符$hat(I)$是一个恒等算符，对任意波函数$psi$有
$
hat(I) psi = psi
$

=== 算符的相等

若对于体系的任何波函数，都有
$
hat(A) psi = hat(B) psi
$
则称算符$hat(A)$和$hat(B)$相等，记作$hat(A) = hat(B)$。

=== 算符之和

若$hat(A)$和$hat(B)$是两个算符，定义它们的和为
$
(hat(A) + hat(B)) psi = hat(A) psi + hat(B) psi
$
例如：Hamilton算符$hat(H) = hat(T) + hat(U)$。

显然算符求和*满足交换率和结合率*：
$
hat(A) + hat(B) = hat(B) + hat(A)\
hat(A) + (hat(B) + hat(C)) = (hat(A) + hat(B)) + hat(C)
$

可以证明：两个线性算符之和仍为线性算符。

=== 算符的乘积

若$hat(A)$和$hat(B)$是两个算符，定义它们的乘积为
$
(hat(A) hat(B)) psi = hat(A) (hat(B) psi)
$
一般说来，算符之积不满足交换率：
$
hat(A) hat(B) != hat(B) hat(A)
$
或者说对易关系
$
[hat(A), hat(B)] = hat(A) hat(B) - hat(B) hat(A) != 0
$
这里$[hat(A), hat(B)]$称为算符$hat(A)$和$hat(B)$的*对易子*。

=== 算符的复共轭算符、转置算符

算符$hat(arrow(P))$的复共轭算符是
$
hat(arrow(P))^* = (- i hbar nabla)^* = i hbar nabla = - hat(arrow(P))
$
写成内积的形式：
$
(psi, hat(A) phi) = (hat(A)^* phi^*, psi^*)
$

#newpara()

算符$hat(A)$的转置算符$hat(A)^T$：
$
integral psi^* hat(A) phi dd(V) = integral phi hat(A)^T psi^* dd(V)
$
#newpara()
定义两个波函数的*内积*或标积为
$
(psi_1, psi_2) = integral psi_1^* psi_2 dd(tau)
$
有：
$
(psi_1, psi_2) = (psi_2^*, psi_1^*)
$
那么转置算符的定义又可写为
$
(psi, hat(A)^T phi) = (phi^*, hat(A) psi^*)
$
#newpara()

转置有性质：
$
(hat(A) hat(B))^T = hat(B)^T hat(A)^T\
(hat(A) + hat(B))^T = hat(A)^T + hat(B)^T
$
给出第一个性质的证明：
$
(psi, (hat(A) hat(B))^T phi) = (phi^*, hat(A) hat(B) psi^*) = (hat(B)^* psi, hat(A)^* phi) = (hat(A)^(T*) phi^*, hat(B)^T psi^*) = (phi, hat(B)^T hat(A)^T psi)
$



=== 算符的逆算符

若算符$hat(A)$满足
$
hat(A) psi = phi\
hat(A)^(-1) phi = psi
$
也就是说根据$phi$可以唯一确定$psi$，则称$hat(A)$是可逆的，$hat(A)^(-1)$称为$hat(A)$的逆算符。但不是所有的算符都有逆算符，如投影算符。

有性质：
$
hat(A) hat(A)^(-1) = hat(A)^(-1) hat(A) = hat(I)\
(hat(A)hat(B))^(-1) = hat(B)^(-1) hat(A)^(-1)\
((hat(A))^(-1))^(-1) = hat(A)
$

在散射微扰问题中，算符$E - hat(H)_0$的逆算符定义为$(E - hat(H)_0)^(-1)$，也就是与传播子相关的格林函数算符。

=== 算符的厄密共轭算符

算符$hat(A)$的厄密共轭算符$hat(A)^dagger$定义为
$
integral psi^* hat(A) phi dd(V) = integral (hat(A)^dagger psi)^* phi dd(V)
$
写成内积的形式：
$
(psi, hat(A) phi) = (hat(A)^dagger psi, phi)
$
#newpara()
有性质：
$
(hat(A))^dagger = hat(A)^(*T) = hat(A)^(T*)\
(hat(A) + hat(B))^dagger = hat(A)^dagger + hat(B)^dagger\
(hat(A) hat(B))^dagger = hat(B)^dagger hat(A)^dagger
$
#newpara()
算符$hat(A)$为*厄密算符*的条件：
$
hat(A)^dagger = hat(A)
$

== 厄密算符的本征值和本征函数

可以定义*厄密算符*
$
integral psi^* hat(F) psi dd(V) = integral (hat(F) psi)^* psi dd(V)
$

=== 厄密算符的本征值

*Hermitian算符的本征值都是实数*
$
hat(F) psi_lambda &= lambda psi_lambda\
(hat(F) psi_lambda)^* &= lambda^* psi_lambda^*\
$
代入定义式
$
lambda integral psi_lambda^* psi_lambda dd(V) &= lambda^* integral psi_lambda psi_lambda^* dd(V)\
lambda &= lambda^*
$
从而得到$lambda$是实数。由于这个定理，我们*要求所有的物理量（或者称为“可测量量”）的算符都是Hermitian算符*（但是反过来不一定）。

不难证明坐标算符和动量算符都是Hermitian算符。在一定条件下，它们的函数也是Hermitian算符。

假设已经证明了$hat(arrow(P))^T = - hat(arrow(P))$，那么有
$
hat(arrow(P))^dagger = (hat(arrow(P))^T)^* = - hat(arrow(P))^* = hat(arrow(P))
$
而对于$hat(arrow(L)) = arrow(r) crossproduct hat(arrow(P))$，有
$
hat(arrow(L))_x = y hat(p)_z - z hat(p)_y \
hat(arrow(L))_x^dagger = hat(p)_z^dagger y - hat(p)_y^dagger z = hat(p)_z y - hat(p)_y z = hat(arrow(L))_x
$

径向动量算符$hat(p)_r = hat(r) dot hat(arrow(p))$*不是Hermitian算符*，因为
$
(hat(r) dot hat(arrow(p)))^dagger =  hat(arrow(p))^dagger dot hat(r)^dagger = hat(arrow(p)) dot hat(r) = - i hbar (partial/(partial r) hat(r)+ 1/r hat(theta) + 1/(r sin(theta)) hat(phi)) dot hat(r) = - i hbar (partial/(partial r) + 2/r)
$
于是有
$
(hat(r) dot hat(arrow(p)))^dagger = hat(r) dot hat(arrow(p)) - 2 i hbar/r
$
所以径向动量算符不是Hermitian算符。但我们可以构造一个厄密算符——*坐标表象*表达式：
$
hat(p)_r = 1/2 ((hat(r) dot hat(arrow(p)))^dagger + hat(r) dot hat(arrow(p))) = - i hbar (partial/(partial r) + 1/(r))
$

=== 厄密算符的本征函数

正交：若两个函数$psi_1$和$psi_2$满足
$
integral psi_1^*psi_2 dd(tau) = 0
$
则称$psi_1$和$psi_2$是正交的。 

*Hermitian算符的本征函数对应于不同本征值的本征函数是正交的*：
$
hat(F) psi_lambda &= lambda psi_lambda\
hat(F) psi_mu &= mu psi_mu\
lambda integral psi_lambda^* psi_mu dd(V) &= mu integral psi_lambda psi_mu^* dd(V)\
(lambda - mu) integral psi_lambda^* psi_mu dd(V) &= 0\
integral psi_lambda^* psi_mu dd(V) &= 0
$
说明了Hermitian算符的本征函数是正交的。

- 若本征值谱是非简并的和离散的，本征值为${lambda_i}$，本征函数为${phi_i}$, 那么波函数是平方可积的，因而可以归一化，所以正交和归一可统一写为
$
integral phi_i^* phi_j dd(tau) = delta_(i j)
$
- 若$F$的本征值谱是非简并的和连续的，本征函数可以按$delta$函数正交“归一”化，即
$
integral phi^*_(lambda) phi_(mu) dd(tau) = delta(lambda - mu)
$
或者是箱归一化。

=== 平面波的箱归一化

对于平面波：
$
psi(arrow(r)) = A e^(i/hbar arrow(p) dot arrow(r))
$
箱归一化要求：粒子波函数在任意边长为$L$的正方体内正交归一化。
- 粒子在三维空间自由运动
- 周期性边界条件
- 箱边长 $L -> oo$

$
integral_V psi^*_(arrow(p)_1)(arrow(r),t) psi_(arrow(p)_2)(arrow(r),t) dd(V) = |C|^2 integral_V e^(i/hbar (arrow(p)_1 - arrow(p)_2) dot arrow(r)) dd(V)
$

在$arrow(p)_1=arrow(p)_2$时，归一化为($a,b,c$为任意常数)：
$
|C|^2 integral_a^(a+L) integral_b^(b+L) integral_c^(c+L) e^(i/hbar (arrow(p)_1 - arrow(p)_2) dot arrow(r)) dd(r) = 1
$
即
$
|C|^2 L^3 = 1 => C = 1/sqrt(L^3)
$

在$arrow(p)_1!=arrow(p)_2$时，
$
integral_0^L e^(i/hbar (p_1 - p_2) dot x) dd(x) = (i/hbar Delta p) ^(-1) (e^(i/hbar Delta p L) - 1)
$
正交要求为：
$
(-i hbar)^3/(L^3 Delta p_x Delta p_y Delta p_z) e^(i/hbar (Delta p_x a + Delta p_y b + Delta p_z c)) (e^(i/hbar Delta p_x L) - 1) (e^(i/hbar Delta p_y L) - 1) (e^(i/hbar Delta p_z L) - 1) = 0
$
于是有：
$
e^(i/hbar Delta p_x L) =  e^(i/hbar Delta p_y L) =  e^(i/hbar Delta p_z L) = 1
$
即
$
Delta p_i = (2 pi hbar n_i) / L\
$
从而
$
p_i = (2 pi hbar n_i) / L + delta_i
$
其中$n_i$是整数，$delta_i$是初始相位，一般不妨设为0。

于是得到$p_i$
$
p_i = (2 pi hbar n_i) / L
$
系统的动量是分立的，但当$L -> oo$时，又过渡到连续的动量谱。

最终的波函数为
$
psi_(arrow(p))(arrow(r),t) = 1/sqrt(L^3) e^((2 pi i)/ L arrow(n) dot arrow(r))
$

== 算符的本征方程

对于Hermitian算符$hat(F)$，有*本征方程*
$
hat(F) psi_lambda = lambda psi_lambda
$
其中$lambda$是*本征值*，$psi_lambda$是$hat(F)$属于本征值$lambda$的*本征函数*。

量子力学关于测量问题的基本假设是：

算符$hat(F)$的本征值集${lambda}$就是力学量$hat(F)$的*测量值集*。

$hat(F)$的本征函数$psi_lambda$代表力学量$F$有确定值$lambda$的量子状态$psi = sum_lambda c_lambda psi_lambda$中的一个分量，概率为$|c_lambda|^2$。

_例：动量本征函数_

动量算符$hat(arrow(p))$的本征方程是
$
hat(arrow(p)) psi_(arrow(p)) = - i hbar nabla psi_(arrow(p)) = arrow(p) psi_(arrow(p))
$
于是
$
psi_arrow(p) = A e^(i arrow(p) dot arrow(r)/hbar)
$
其中$A = 1/(2 pi hbar)^(3/2)$是归一化系数。但是，在无穷空间中它们是平方不可积的，这时它们正交归一化为$delta$函数。

_例：位置算符_

位置算符$hat(arrow(r))$的本征方程是
$
hat(arrow(r)) psi_(arrow(r)) = arrow(r) psi_(arrow(r))
$
解出的本征函数是$psi_(arrow(r)_0) = delta(arrow(r) - arrow(r_0))$，即位置算符的本征函数是位置本身。

对波函数，可以按照位置算符的本征函数展开：
$
psi(arrow(r)) &= sum_arrow(r_0) c_(arrow(r_0)) psi_(arrow(r_0)) psi_(arrow(r)) \
&= integral c_(arrow(r_0)) delta(arrow(r) - arrow(r_0)) dd(r_0) = c_(arrow(r))
$
这就是波函数的位置表象。


== 简并波函数的正交化

如果出现简并（即一个本征值有若干个线性独立的本征函数）的情形，则*正交性定理不能保证同一本征值的不同本征函数是彼此正交的*。

经过对本征函数进行适当的重新组合，可以使它们仍然是正交的。这个过程称为*正交化*。

=== Schmidt正交化

*Gram-Schmidt正交化*方法是一种常用的正交化方法。

设函数$rho_1, rho_2, ..., rho_n$是线性独立的，但不正交的函数。我们要把它们正交化。

1. 先把第一个函数归一化：
2. 第二个函数减去第一个函数在第二个函数方向上的投影，得到新的函数，再归一化：
3. 第三个函数减去前两个函数在第三个函数方向上的投影，得到新的函数，再归一化：
4. 以此类推，直到最后一个函数。
$
rho_1 = rho_1\
rho_2 = rho_2 - ((rho_2, rho_1))/((rho_1, rho_1)) rho_1\
rho_3 = rho_3 - ((rho_3, rho_1))/((rho_1, rho_1)) rho_1 - ((rho_3, rho_2))/((rho_2, rho_2)) rho_2\
...
$
这样得到的函数就是正交的。

=== 共同本征函数

在量子力学中，一个更为物理的解决简并本征函数的办法是考虑*两个算符的共同本征函数*。

*对易：*若$hat(F)$和$hat(G)$是两个算符，若它们的*对易子*是
$   
[hat(F), hat(G)] = hat(F) hat(G) - hat(G) hat(F)
$
若$[hat(F), hat(G)] = 0$，则称$hat(F)$和$hat(G)$是*对易*的。

*共同本征函数：*若$[hat(F), hat(G)] = 0$则$hat(F)$和$hat(G)$有共同的本征函数。即存在$phi$使得
$
hat(F) phi = lambda phi\
hat(G) phi = mu phi
$

该定理也很容易推广到多个算符的情形。

共同本征函数描写的就是几个力学量同时有确定值的状态。

这样，如果算符$hat(F)$的本征值$l$有简并，我们就再引进另一算符$hat(G)$，使得$hat(F)$和$hat(G)$有共同的本征函数，这样就可以把简并的本征函数正交化。

如果对于$F$简并的本征函数对于$G$不是简并的，那么*正交性定理就保证了它们是正交的*。但也可能$F$和$G$的共同本征函数仍然有简并，我们就再引进第三个算符，如此等等，直到所有的简并完全去除为止。这时，*一组量子数*就完全确定了一个量子态。

这种情形多半出现在多自由度体系中。对这种体系，一组两两对易的、完全去除简并的算符集称为它的*对易可观测量完全集(CSCO)*。完备算符集中算符的数目就是体系的*自由度数*。

如果这些量子数都是分立量子数，共同本征函数的正交归一关系就是：
$
(phi_(n l m), phi_(n' l' m')) = delta_(n n') delta_(l l') delta_(m m')
$

#newpara()

*例如：动量算符*

动量算符的各个分量是彼此对易的：
$
[p_x, p_y] = [p_x, p_z] = [p_y, p_z] = 0
$
所以动量算符的三个分量有共同的本征函数：
$
phi_arrow(p) (x,y,z) = (1/sqrt(2 pi hbar))^3 e^(i/hbar p_x x) e^(i/hbar p_y y) e^(i/hbar p_z z)
$
即是三维平面波，任何波函数都可以用它们来展开（函数的Fourier变换）。

*例如：动量算符和哈密顿算符*

对一维自由粒子来说：
$
hat(p)_x = - i hbar partial/(partial x)\
hat(H) = hat(p)^2/(2m) = - hbar^2/(2m) partial^2/(partial x^2)
$
这两个算符是对易的：
$
[hat(p)_x, hat(H)] = 0
$
所以它们有共同的本征函数：
$
phi_(p_x) (x) = 1/sqrt(2 pi hbar) e^(i/hbar p_x x)
$
其中参数$p_x$为$hat(p)_x$算符的某个本征值。

*例：对于氢原子，考察下面三个算符*
$
hat(L) &= - i hbar hat(Y)_z = - i hbar partial/(partial phi)\
hat(L)^2 &= - hbar^2 hat(Y)^2 = - hbar^2 (1/sin(theta) partial/(partial theta) (sin(theta) partial/(partial theta)) + 1/(sin(theta)^2) partial^2/(partial phi^2))\
hat(H) &= - hbar^2/(2m) nabla^2 + U(r) = - hbar^2/(2m r^2) partial/(partial r) (r^2 partial/(partial r)) + hat(L)^2/(2 m r^2) + U(r)
$
有：
$
[hat(L)^2, hat(H)] = [hat(L), hat(H)] = [hat(L), hat(L)^2] = 0
$
所以它们有共同的本征函数——氢原子能量本征函数：
$
psi_(n l m) (r, theta, phi) = R_(n l) (r) Y_(l m) (theta, phi)
$
满足正交归一：
$
(psi_(n l m), psi_(n' l' m')) = delta_(n n') delta_(l l') delta_(m m')
$

#newpara()

已知算符$hat(A),hat(B)$满足$[hat(A),hat(B)] = 0$ ，$phi_B$是$hat(B)$的一个本征态（对应本征值为$B$），则$hat(A) phi_B$也是$hat(B)$的一个本征态（对应本征值为$B$）。从而$hat(A) phi_B = A phi_B$。其中$A$也是该函数的本征值。这就证明了$hat(A)$和$hat(B)$有共同的本征函数。

=== 对易括号的运算

1. 对易括号是交换反对称的，即
$
[hat(A), hat(B)] = - [hat(B), hat(A)]
$
2. 对易括号的运算满足线性性，即
$
[hat(A), hat(B) + hat(C)] = [hat(A), hat(B)] + [hat(A), hat(C)]\
[hat(A) + hat(B), hat(C)] = [hat(A), hat(C)] + [hat(B), hat(C)]\
[c hat(A), hat(B)] = [hat(A), c hat(B)] = c [hat(A), hat(B)]
$
3. 对易括号的运算满足Leibniz法则，即
$
[hat(A), hat(B) hat(C)] = [hat(A), hat(B)] hat(C) + hat(B) [hat(A), hat(C)]\
[hat(A) hat(B), hat(C)] = hat(A) [hat(B), hat(C)] + [hat(A), hat(C)] hat(B)
$
4. 量子力学的基本对易括号是
$
[hat(x)_i, hat(p)_j] = i hbar delta_(i j)
$
其中$hat(p)_i = - i hbar partial/(partial x_i)$（坐标表象）。

_证明：_
$
[hat(x)_i, hat(p)_j] psi &=  hat(x)_i hat(p)_j psi - hat(p)_j hat(x)_i psi\
&= - i hbar (x_i partial/(partial x_j) psi -  partial/(partial x_j) (x_i psi))\
&= i hbar delta_(i j) psi
$

#newpara()

利用上面给出的对易括号的性质和运算法则：
$
[hat(x), hat(F)] = i hbar hat((partial F)/(partial p_x))\
[hat(p), hat(F)] = - i hbar hat((partial F)/(partial x))
$
其中$hat(F) = hat(F)(hat(x), hat(p)) = sum_(m,n = 0)^oo c_(m n) hat(x)^m hat(p)^n$，是算符$hat(x)$和$hat(p)$的多项式。

*例：角动量算符的对易括号*
$
[hat(L)_i, hat(L)_j] = i hbar epsilon_(i j k) hat(L)_k
$
其中：
$
hat(L)_i = hat(r) crossproduct hat(p)_i = - i hbar epsilon_(i j k) hat(r)_j hat(p)_k
$
而角动量平方算符的对易括号：
$
[hat(L)^2, hat(L)_i] = 0
$
其中：
$
hat(L)^2 = hat(L)_x^2 + hat(L)_y^2 + hat(L)_z^2
$
角动量各分量之间互相不对易有深刻的物理结果。

== 波函数按本征函数系展开

一维情形。假设力学量算符$hat(F)$的本征值集是${λ_n, n=1,2,...}$,(离散的、非简并的)，本征函数系是${phi_n(x), n = 1, 2,...}$按叠加原理，
$
psi(x) = sum_n c_n phi_n (x)
$
注意到${phi_n(x)}$是正交归一的，
$
(phi_n, phi_m) = delta_(n m)
$
所以
$
c_n = (phi_n, psi) = integral phi_n^* psi dd(x)
$

只有当${phi_n (x)}$是完备的函数系时，才能用它来展开任意的连续函数：
$
psi(x) = sum_n integral phi_n (x') psi(x') dd(x') phi_n (x) = integral (sum_n phi_n (x) phi_n (x')) psi(x') dd(x')\
psi(x) = integral delta(x - x') psi(x') dd(x')
$
从而得到
$
sum_n phi_n (x) phi_n (x') = delta(x - x')
$
这个条件就称为函数系${phi_n (x)}$的完备性条件。

_注：_
- 本征值是连续谱，本征函数系是$phi_lambda (x)$（$lambda$连续变化）
    $
    psi (x) = integral c_lambda phi_lambda (x) dd(lambda)\
    integral phi_lambda^* phi_lambda' dd(x) = delta(lambda - lambda')\
    c_lambda = integral phi_lambda^*(x) psi(x) dd(x)\
    integral phi_lambda^* psi (x) phi_lambda (x') dd(lambda) = delta(x - x')
    $

- 多自由度体系（例如三维运动）。这时要按CSCO算符集的共同本征函数系展开。系数的计算方法是类似的。
- 与时间有关的波函数$c_n -> c_n (t)$，展开系数也是时间的函数。
- 根据完备力学量集的定义和态的叠加原理，完备力学量集的全体算符的共同本征函数构成了表示该系统量子状态的正交归一的完备基底，即系统的任何状态都可以展开为这些基底的线性组合。

== 量子力学量的测量-波包坍缩

量子力学的测量结果是几率性的，比如我们测一个非定态系统的能量，其波函数为：
$
psi(x, t) = sum_n c_n(t) phi_n(x) e^(-i E_n t/hbar)
$
在测量以前，系统的状态是许许多多本征态的叠加。测量之后，系统坍缩为某一个本征态：
$
sum_n c_n(t) phi_n(x) e^(-i E_n t/hbar) ->^"测量并读数" phi_n(x) e^(-i E_n t/hbar)
$
这一过程称为*“波包坍缩”*（von Neumann，1932年）。

波包坍缩的动力学过程至今仍在研究（不服从薛定谔方程）。量子力学关于测量的假定是理论的基本假定之一，是量子力学目前无法解释的。比如，在对粒子做空间位置测量后的一刻，其波函数坍缩为
$
psi(x) = delta(x - x_0)
$

=== 力学量的测量几率

一维离散情形：假设力学量算符$hat(F)$的本征值集是${lambda_n, n=1,2,...}$，本征函数系是${phi_n (x), n = 1, 2,...}$，波函数为
$
psi(x) = sum_n c_n phi_n (x)
$
则测量力学量$F$的本征值为$lambda_n$的几率是
$
w(lambda_n) = |c_n|^2
$

#newpara()

总几率不变的验证: 测量$hat(F)$得到各种可能测量值的总几率为
$
sum_n w(lambda_n) &= sum_n |integral phi_n^* psi dd(x)|^2 \
&= sum_n integral phi_n (x) psi^*_n (x) dd(x) integral psi (x') phi_n^*(x') dd(x') \
&= integral.double (sum_n phi_n (x)psi_n^* (x')) phi_n^*(x) psi (x') dd(x)dd(x') \
&= integral.double delta(x - x') psi^* (x) psi (x') dd(x)dd(x') \
&= integral delta(x - x) (integral psi^* (x) psi (x') dd(x)) dd(x') \
&= integral psi^* (x) psi (x) dd(x) = 1
$

推广：
- 本征值是连续的。此时要引入几率密度：记测量值在$lambda->lambda+dd(lambda)$之间的几率为$dd(W(lambda))$，则
    $
    w(lambda) = dd(W(lambda))/dd(lambda) = |c_lambda|^2
    $
    是$lambda$的测量几率密度，它的计算公式是
    $
    c_lambda = integral phi_lambda^*(x) psi(x) dd(x)
    $
- 对多自由度体系，只问某一个力学量的测量几率，经常会有简并。这时要找到一个包含$hat(F)$的CSCO完全集，并求出它们的共同本征函数。设力学量$F$的离散的本征值$lambda_n$的简并度为$k$，简并的本征态为$phi_(n 1), phi_(n 2), ..., phi_(n k)$，则
  $
  psi(x) = sum_n sum_(i = 1)^k c_(n i) phi_(n i) (x)
  $
    测量$F$得到$lambda_n$的几率是
    $
    w(lambda_n) = sum_(i = 1)^k |c_(n i)|^2
    $
    对连续谱的情况也做类似的推广。

=== 力学量的平均值

可以定义力学量$F$的平均值为

$
macron(F) &=  sum_n lambda_n w(lambda_n) \
&= sum_n lambda_n |c_n|^2 \
&= sum_n integral phi_n (x) psi^*_n (x) dd(x) integral psi (x') (hat(F)phi_n (x'))^* dd(x') \
&= integral.double (sum_n phi_n (x)psi_n^* (x')) (hat(F)phi_n (x'))^* dd(x)dd(x') \
&= integral.double delta(x - x') psi^* (x) hat(F) psi (x') dd(x)dd(x') \
&= integral psi^* (x) hat(F) psi (x) dd(x)
$
这个计算式的条件是$psi(x)$已经归一，即
$
integral psi^* (x) psi (x) dd(x) = 1
$
如果没有归一：
$
macron(F) =( integral psi^* (x) hat(F) psi (x) dd(x) )/( integral psi^* (x) psi (x) dd(x))
$

#pagebreak(weak: true)

= 算符对易和不确定关系

== 对易算符与本征函数

定理：两个力学量算符$hat(F)$和$hat(G)$有一组完备的共同本征函数的充要条件是它们彼此对易
$
[hat(F), hat(G)] = 0
$

*必要性：*设$hat(F)$和$hat(G)$有一组完备的共同本征函数${phi_(f g)}$，则
$
hat(F) phi_(f g) = f phi_(f g)\
hat(G) phi_(f g) = g phi_(f g)
$
对于任意波函数$psi$，都可以展开为这组共同本征函数的线性组合：
$
[hat(F), hat(G) psi] = hat(F) hat(G) psi - hat(G) hat(F) psi = (f g - g f)psi = 0
$
根据$psi$的任意性
$
[hat(F), hat(G)] = 0
$

*充分性：*如果$hat(F)$和$hat(G)$中有一个是非简并的（比如$hat(F)$），也就是说对应每个本征值只有一个线性无关本征态
$
hat(F) phi_f = f phi_f
$
那么
$
0 = [hat(F), hat(G)] phi_f = hat(F) hat(G) phi_f - hat(G) hat(F) phi_f = hat(F) hat(G) phi_f - f hat(G) phi_f\
hat(F)(hat(G) phi_f) = f (hat(G) phi_f)
$
也就是说，$phi_f$和$hat(G)phi_f$都属于$hat(F)$的对应本征值$f$的本征函数，非简并有：
$
hat(G) phi_(f g) = g phi_(f g)
$
这样就证明了$hat(F)$和$hat(G)$有一组完备的共同本征函数。

对于二者都简并的情况，设
$
hat(F) phi_(f k) = f phi_(f k) (k = 1, 2, ..., n)
$
有
$
hat(F)(hat(G) phi_(f k)) = f (hat(G) phi_(f k))
$
那么必有
$
hat(G) phi_(f k) = sum_(m = 1)^n c_(k m) phi_(f m)
$
写成矩阵形式，同时略去角标$f$：
$
hat(G) 
mat(
    phi_1;
    phi_2;
    dots.v;
    phi_n
)
= 
mat(
    c_11 , c_12, ..., c_(1n);
    c_21 , c_22, ..., c_(2n);
    dots.v, dots.v, dots.down, dots.v;
    c_(n 1) , c_(n 2), ..., c_(n n)
)
mat(
    phi_1;
    phi_2;
    dots.v;
    phi_n
)
$
*把系数矩阵对角化(厄密矩阵可以对对角化)，就得到了$hat(G)$的本征值和本征函数。*

存在$n times n$的变换矩阵$U$，使得
$
U mat(
    c_11 , c_12, ..., c_(1n);
    c_21 , c_22, ..., c_(2n);
    dots.v, dots.v, dots.down, dots.v;
    c_(n 1) , c_(n 2), ..., c_(n n)
) U^(-1) = diag(g_1, g_2, ..., g_n)
$
那么：
$
U hat(G) mat(
    phi_1;
    phi_2;
    dots.v;
    phi_n
) = U 
mat(
    c_11 , c_12, ..., c_(1n);
    c_21 , c_22, ..., c_(2n);
    dots.v, dots.v, dots.down, dots.v;
    c_(n 1) , c_(n 2), ..., c_(n n)
)
mat(
    phi_1;
    phi_2;
    dots.v;
    phi_n
)
$
即
$
hat(G) 
mat(
    phi_1^';
    phi_2^';
    dots.v;
    phi_n^'
) =
mat(
    g_1 , 0, ..., 0;
    0 , g_2, ..., 0;
    dots.v, dots.v, dots.down, dots.v;
    0 , 0, ..., g_n
)
mat(
    phi_1^';
    phi_2^';
    dots.v;
    phi_n^'
)
$
其中
$
mat(
    phi_1^';
    phi_2^';
    dots.v;
    phi_n^'
) = U
mat(
    phi_1;
    phi_2;
    dots.v;
    phi_n
)
$
这样就得到了(再把角标$f$放回)：
$
hat(F) phi_(f k)^' = f phi_(f k)^'\
hat(G) phi_(f k)^' = g_k phi_(f k)^'
$
这样就证明了$hat(F)$和$hat(G)$有一组完备的共同本征函数。

这个证明过程有两点启示：
- 提供了一种消除简并的方法（找对易算符或CSCO算符集）；
- 算符的作用可以用矩阵来表示－海森堡矩阵力学。

两个力学量算符不对易，则它们不能有共同的完备的本征函数集，但不排除它们碰巧有个别共同本征函数。例：
$
[hat(L)_x, hat(L)_y] = i hbar hat(L)_z\
psi(arrow(r)) = 1/(2 pi hbar)^(3/2)\
$
有共同本征函数：
$
hat(L)_z psi(arrow(r)) = 0 psi(arrow(r))\
hat(L)_y psi(arrow(r)) = 0 psi(arrow(r))
$

== 对易算符与不确定性原理

若
$
[hat(F), hat(G)] != 0
$
则一般来说$F$和$G$不能同时有确定值。例如:
$
hat(p)_x = - i hbar partial/(partial x)\
hat(x) = x
$
它在本质上是波粒二象性的反映。例如在粒子的单缝衍射实验中，$Delta x$越小,$Delta p_x$越大,
$
Delta x Delta p_x approx hbar
$
二者不能同时有确定值。所以，运动轨道的概念对微观粒子是不适用的。



对这种不确定性的定量描写如下。定义偏差算符为：
$
Delta hat(F) = hat(F) - macron(F)
$
有性质：
$
[Delta hat(F), Delta hat(G)] = [hat(F) - macron(F), hat(G) - macron(G)] = [hat(F), hat(G)]
$
$
macron(Delta hat(F)) = macron(hat(F) - macron(F)) = macron(F) - macron(F) = 0
$
$
macron(Delta hat(F)^2) = macron((hat(F) - macron(F))^2) = macron(hat(F)^2 - 2 macron(F) hat(F) + macron(F)^2) = macron(hat(F)^2) - 2 macron(F) macron(F) + macron(F)^2 = macron(hat(F)^2) - macron(F)^2
$

#newpara()

如果
$
[hat(F), hat(G)] = i hat(C) != 0
$
其中$hat(C)$是厄密算符，考虑积分不等式：
$
I(xi) = integral |(xi Delta hat(F) - i Delta hat(G))psi|^2 dd(tau) >= 0
$
其中$psi$为体系的任一态，$xi$为任意实数。化简得到：
$
I(xi) &= integral (xi(Delta hat(F) psi)^* + i (Delta hat(G) psi)^*) (xi Delta hat(F) psi - i Delta hat(G) psi) dd(tau) \
&= xi^2 integral (Delta hat(F) psi)^* (Delta hat(F)) psi dd(tau) + i xi integral (Delta hat(G) psi)^* (Delta hat(F) psi) dd(tau) - i xi integral (Delta hat(F) psi)^* (Delta hat(G) psi) dd(tau) + integral (Delta hat(G) psi)^* (Delta hat(G) psi) dd(tau) \
&=^"Hermit性质" xi^2 integral psi^* (Delta hat(F))^2 psi dd(tau) +  xi integral psi^* (Delta hat(G) Delta hat(F) + Delta hat(F) Delta hat(G)) psi dd(tau) + integral psi^* (Delta hat(G))^2 psi dd(tau) \
&= xi^2 macron(Delta hat(F)^2) + xi macron(hat(C)) + macron(Delta hat(G)^2) >= 0
$
这是一个关于$xi$的二次函数，所以它的判别式小于等于0，得到*Heisenberg不确定关系*：
$
macron(Delta hat(F)^2) macron(Delta hat(G)^2) >= 1/4 |macron(hat(C))|^2
$
也是Schwartz不等式。

若$[hat(F),hat(G)] != 0$则一般说来$Delta F$和$Delta G$不能同时为零，即$F$和$G$不能同时有确定值（但是注意$macron([hat(F),hat(G)]) = 0$的特殊态例外），或者说它们不能有共同的本征态。

反之，若$[hat(F),hat(G)] = 0$，则可以找到这样的态，使得$Delta F$和$Delta G$同时为零，即$F$和$G$可以同时有确定值，或者说它们有共同的本征态。


#newpara()

对于坐标和动量算符：
$
macron(Delta x)^2 macron(Delta p)^2 >= 1/4 hbar^2
$
即：
$
Delta x Delta p_x >= hbar/2
$
上面的不等式中取“=”的量子态，被称为“最小不确定态”。

对于谐振子的基态：
$
hat(H) = 1/(2m) hat(p)^2 + 1/2 m omega^2 hat(x)^2\
macron(E) = 1/2 macron(p)^2/(2m) + 1/2 m omega^2 macron(x)^2\
$
对于谐振子任意本征态，$macron(p) = 0$，$macron(x) = 0$，所以
$
macron(hat(p)^2) = macron(Delta hat(p)^2) \
macron(hat(x)^2) = macron(Delta hat(x)^2) \
$
从而得到
$
macron(E) = 1/(2m) macron(Delta p)^2 + 1/2 m omega^2 macron(Delta x)^2
$
在极限情况下，即最小不确定态，有
$
macron(E) = hbar omega/2
$
谐振子的基态是最小不确定态。（不是所有量子系统的基态都是最小不确定态）


== 联级Stern-Gerlach实验

$hat(L)_x$和$hat(L)_z$不对易，所以不能同时有确定值。

#figure(
  image("pic/2024-04-25-17-01-26.png", width: 80%),
  numbering: none
)


#figure(
  image("pic/2024-04-25-17-02-28.png", width: 80%),
  numbering: none
)

#figure(
  image("pic/2024-04-25-17-03-13.png", width: 80%),
  numbering: none
)


#figure(
  image("pic/2024-04-25-17-04-31.png", width: 80%),
  numbering: none
)

#pagebreak(weak: true)

= 角动量算符

== 角动量算符的本征值和本征态

角动量算符（轨道角动量）的定义是：
$
hat(L) = hat(r) crossproduct hat(p) = - i hbar hat(r) crossproduct nabla\
hat(L)_z = - i hbar (x partial/(partial y) - y partial/(partial x))\
hat(L)^2 = hat(L)_x^2 + hat(L)_y^2 + hat(L)_z^2
$
角动量算符的球坐标表示为：
$
hat(L)_x &= i hbar (sin(phi) partial/(partial theta) + cot(theta) cos(phi) partial/(partial phi))\
hat(L)_y &= - i hbar (cos(phi) partial/(partial theta) - cot(theta) sin(phi) partial/(partial phi))\
hat(L)_z &= - i hbar partial/(partial phi)\
hat(L)^2 &= - hbar^2 (1/sin(theta) partial/(partial theta) (sin(theta) partial/(partial theta)) + 1/(sin(theta)^2) partial^2/(partial phi^2))
$

== $hat(L)_z$的本征值和本征函数

$hat(L)_z$的本征值方程是：
$
hat(L)_z psi_(m) = m hbar psi_(m)
$
解出的本征函数是：
$
psi_(m) = C e^(i m phi)
$
由波函数的连续性，必须有：
$
psi_m (phi) = psi_m (phi + 2 pi)
$
周期性边界条件要求$m$是整数。归一化条件是：
$
1 = integral_0^(2 pi) |C|^2 d phi = 2 pi |C|^2\
C = 1/sqrt(2 pi)
$

== $hat(L)^2$的本征值和本征函数

$hat(L)^2$的本征值方程是：
$
hat(L)^2 Y = lambda hbar^2 Y , lambda = l (l + 1)\
1/(sin theta) partial/(partial theta) (sin theta partial/(partial theta) Y) + 1/(sin^2 theta) partial^2/(partial phi^2) Y = - lambda Y
$
分离变量法，设$Y(theta, phi) = Theta(theta) Phi(phi)$，得到两个方程：
$
1/(sin theta) partial/(partial theta) (sin theta partial/(partial theta) Theta) + lambda Theta = 0\
partial^2/(partial phi^2) Phi = - m^2 Phi
$
第一个方程是Legendre方程，解为连带Legendre多项式$P_l^m (cos theta)$，第二个方程是周期函数，解为$e^(i m phi)$。

详细地，缔合Legendre方程是：
$
dd("")/dd(w) (1 - w^2) dd(P)/dd(w) + (lambda - m^2/(1 - w^2)) P = 0
$
$w=plus.minus 1$是这个方程的“奇点”，除非$l$取某些特定值，方程的解将在$w=plus.minus 1$处变成无穷大。这些使得级数截断的$l$满足：
$
l = |m|, |m| + 1, |m| + 2, ...
$
它的解是连带Legendre函数：
$
P_l^m (cos theta) = 1/(2^l l!)(1 - cos^2 theta)^(m/2) dd("")^(l+m)/(d cos^(l+m) theta) (w^2 - 1)^l , |m| <= l
$
其正交归一性是：
$
integral_(-1)^1 P_l^m (cos theta) P_(l')^m (cos theta) dd(cos theta) = 2/(2 l + 1) (l + m)!/(l - m)! delta_(l ,l')
$
根据球谐函数的正交归一性：
$
integral_0^(2 pi) integral_0^pi Y_(l m)^* Y_(l' m') sin theta dd(theta) dd(phi) = delta_(l l') delta_(m m')
$
得到$Y(theta, phi)$球谐函数：
$
Y_(l m) (theta, phi) = sqrt((2 l + 1)/(4 pi) (l - m)!/(l + m)!) P_l^m (cos theta) e^(i m phi)
$
其中$l$是角量子数，$m$是磁量子数。

$l = 0,1,2,...$，对应SPDF态。对于给定的$l$，$m = -l, -l + 1, ..., l$，这证明$hat(L)^2$是$2l + 1$度简并的。

== 球谐函数的基本性质

1. $Y_(l m) (theta, phi)$是角动量算符$hat(L)^2$和$hat(L)_z$的共同本征函数，对应本征值分别是$l (l + 1) hbar^2$和$m hbar$。
   $
   hat(L)^2 Y_(l m) = l (l + 1) hbar^2 Y_(l m)\
    hat(L)_z Y_(l m) = m hbar Y_(l m)
   $
2. $Y_(l m) (theta, phi)$是正交归一的：
   $
   integral_0^(2 pi) integral_0^pi Y_(l m)^* Y_(l' m') sin theta dd(theta) dd(phi) = delta_(l l') delta_(m m')
   $
3. 关于宇称的定义也可以推广到三维空间。变换
   $
   cal(P) : arrow(r) -> - arrow(r)
   $
    称为宇称变换。对于球谐函数，有
    $
    cal(P) Y_(l m) (theta, phi) = (-1)^l Y_(l m) (theta, phi)
    $
    即$Y_(l m) (theta, phi)$的宇称是$(-1)^l$。
4. $Y_(l m) (theta, phi)$是单位球面$r=1$，上的完备正交函数系。任意函数$psi(theta, phi)$都可以展开为球谐函数的级数：
   $
   psi(theta, phi) = sum_(l = 0)^oo sum_(m = -l)^l c_(l m) Y_(l m) (theta, phi)
   $
   其中系数$c_(l m)$是：
   $
   c_(l m) = integral_0^(2 pi) integral_0^pi Y_(l m)^* psi(theta, phi) sin theta dd(theta) dd(phi)
   $
5. $Y_(l m) (theta, phi)$满足微分方程：
   $
   - r^2 nabla^2 Y_(l m) = l (l + 1) hbar^2 Y_(l m)
   $
   如果一个函数$f(theta, phi)$满足
    $
    - r^2 nabla^2 f(theta, phi) = lambda f(theta, phi)
    $
    那么$f(theta, phi)$就位于角动量量子数$l$的子空间里，相应的其展开式为
    $
    f(theta, phi) = sum_(m = -l)^l c_m Y_(l m) (theta, phi)
    $
6. 加法定理：设球坐标系中有两个矢量
   $
   arrow(r) = (r, theta, phi)\
    arrow(r') = (r', theta', phi')
   $
    则
    $
    cos gamma = cos theta cos theta' + sin theta sin theta' (cos(phi - phi'))
    $
    其中$gamma$是两个矢量之间的夹角。这个公式可以用来证明球谐函数的加法定理：
    $
    P_l (cos gamma) = (4 pi)/(2 l + 1) sum_(m = -l)^l Y_(l m) (theta, phi) Y_(l m)^* (theta', phi')
    $
7. 两点距离倒数的展开：
   $
   1/abs(arrow(r) - arrow(r')) = 1/r (1 + r'^2/r^2 - 2 r'/r cos gamma)^(-1/2) = cases(
    1/r sum_(l = 0)^oo (r'/r)^l P_l (cos gamma)  &"if" r > r'\
    1/r' sum_(l = 0)^oo (r/r')^l P_l (cos gamma) &"if" r < r'
   )
   $
   以及相应的积分：
    $
    integral_0^(2 pi) integral_0^pi 1/abs(arrow(r) - arrow(r')) sin theta dd(theta) dd(phi) = cases(
        4pi/r "if" r > r'\
        4pi/r' "if" r < r'
    )
    $

#pagebreak(weak: true)

= 量⼦系统的时间演化

== 不含时$hat(H)$的本征函数

定态薛定谔方程：
$
hat(H) psi = E psi
$
其中$hat(H)$是哈密顿算符，$E$是能量本征值，$phi$是能量本征函数。薛定谔方程的解是：
$
psi(arrow(r), t) = sum_n a_n (t) phi_n (arrow(r)) 
$
带入薛定谔方程，得到：
$
i hbar partial/(partial t) psi(arrow(r), t) = hat(H) psi(arrow(r), t)\
sum_n i hbar dd(a_n (t))/dd(t) phi_n (arrow(r)) = sum_n a_n (t) hat(H) phi_n (arrow(r))\
$
即：
$
i hbar dd(a_n (t))/dd(t) = E_n a_n (t)
$
解得：
$
a_n (t) = a_n(0) e^(-i/hbar E_n t)
$
所以：
$
psi(arrow(r), t) = sum_n a_n e^(-i/hbar E_n t) phi_n (arrow(r))
$
这就是量子系统的时间演化。其中
$
a_n (0) = a_n
$
与时间无关。并且可以得到：
$
|a_n (t)|^2 = |a_n|^2
$
即*系统在任意时刻的能量几率分布都和初始时刻的能量几率分布相同*。

在$hat(H)$和时间无关的情况下，只要我们完全地解决了定态薛定谔方程的问题，那么一旦知道了波函数的初始值，与时间有关的Schrödinger方程的解就可以很方便地得出。形式地说，这个解是
$
Psi(arrow(r), t) = e^(-i/hbar hat(H) t) Psi(arrow(r), 0)
$
其中$Psi(arrow(r), 0)$是初始波函数，而$e^(-i/hbar hat(H) t)$是*时间演化算符*。

我们可以用时间演化算符作用在初始波函数上来得到此后任一时刻系统的波函数：
$
psi(arrow(r), t) &= sum_n a_n (0) e^(-i E_n t/hbar) phi_n (arrow(r))\
&= sum_n a_n (0) e^(-i/hbar hat(H) t) phi_n (arrow(r))\
&= e^(-i/hbar hat(H) t) sum_n a_n (0) phi_n (arrow(r))\
&= e^(-i/hbar hat(H) t) psi(arrow(r), 0)\
&= hat(U)(t) psi(arrow(r), 0)
$
其中$hat(U)(t) = e^(-i/hbar hat(H) t)$是时间演化算符。

_例：一维自由粒子高斯波包的时间演化。_

$t=0$时刻波函数为
$
psi(x, 0) = (1/(2 pi sigma^2))^(1/4) e^(-x^2/(4sigma^2)) e^(i k_0 x)\
|psi(x, 0)|^2 = (1/(2 pi sigma^2))^(1/2) e^(-x^2/(2sigma^2))"是Gauss函数"
$
初始时刻的波包是一个高斯函数，它的宽度是$sigma$，时间$t$之后，波函数演化为
$
psi(x, t) &= e^(-i/hbar hat(H) t) psi(x, 0)\
$
其中*自由粒子的哈密顿算符的本征函数是平面波*：
$
psi(x, t) &= e^(-i/hbar hat(H) t) psi(x, 0)\
&= e^(-i/hbar hat(H) t) 1/sqrt(2pi) integral e^(-i/hbar hat(H) t) psi(k) dd(k)\
&= 1/sqrt(2pi) integral e^(-i/hbar hat(H) t) phi(k) e^(i k x) dd(k)\
&= 1/sqrt(2pi) integral phi(k) e^(i k x - i (hbar k^2)/(2m) t) dd(k)\
&= (sqrt(2 pi)sigma(1 + (i hbar t)/(2 m sigma^2)))^(-1/2) e^(i (k_0 x - (hbar k_0^2)/(2 m) t) - (x - (hbar k_0 t)/m)^2/(2(2 sigma^2 + i hbar t/m)))
$
其中
$
phi(k) = ((2 sigma^2)/pi)^(1/4) e^(-sigma^2 (k - k_0)^2)
$
于是得到模方：
$
|psi(x, t)|^2 prop e^(-(x - (hbar k_0 t)/m)^2/(2(2 sigma^2 + (hbar^2 t^2)/(4 m^2 sigma^2))))
$
对照高斯分布表达式，可见时间$t$之后，波包中心运动速度为
$
x/t = (hbar k_0)/m = p_0/m = v_0
$
波包宽度变为
$
sigma_t = sigma sqrt(1 + (hbar^2 t^2)/(4 m^2 sigma^2))
$
波包宽度随时间增大，即波包发散。这反驳了薛定谔对其波函数即是粒子的夸大解释。

对于位置算符的本征态$delta(x)$，对其做时间演化：
$
e^( - i/hbar hat(H) t) delta(t) &= e^(- i/hbar hat(H) t) (1/sqrt(2 pi))^2 integral e^(i k x) dd(k)\
&= 1/(2 pi) integral e^(i k x - i (hbar k^2)/(2m) t) dd(k)\
&= sqrt(m/(2 pi i hbar t)) e^(-m x^2/(2 i hbar t)) e^(-i pi/4)
$
似乎在无穷远处也能有概率，这是超光速的。这就证明Schrodinger方程是非相对论的，不能描述高速粒子，在此时不再适用。




== 力学量平均值随时间的变化

对于力学量$hat(F)$，其平均值为：
$
macron(F) = integral psi^* hat(F) psi dd(x)
$
在推导过程中，并未要求$psi$的展开系数与时间无关。实际上这个公式对任意时刻$t$的波函数都适用
$
macron(F) = integral psi^* (x,t) hat(F) psi(x,t) dd(x)
$
如此便产生了平均值随时间变化的结果。平均值随时间变化的原因一般有：
- 力学量算符本身显含时间
- 力学量算符与哈密顿算符不对易
我们假设系统的哈密顿量是与时间无关的
$
dd(macron(F))/dd(t) = integral dd(psi^* (arrow(r),t))/dd(t) hat(F) psi(arrow(r),t) + psi^* (arrow(r),t) hat(F) dd(psi(arrow(r),t))/dd(t) dd(tau)
$
将薛定谔方程代入，得到
$
dd(macron(F))/dd(t) &= (1/i hbar) integral psi^* (arrow(r),t) hat(H) hat(F) psi(arrow(r),t) - psi^* (arrow(r),t) hat(F) hat(H) psi(arrow(r),t) dd(tau)\
&= (1/i hbar) integral psi^* (arrow(r),t) [hat(H), hat(F)] psi(arrow(r),t) dd(tau)
$
如果$hat(H)$和$hat(F)$对易，则平均值不随时间变化。这是量子力学中的一个重要结论，称为*Ehrenfest定理*：
$
dd(macron(F))/dd(t) = 1/(i hbar) macron([hat(H), hat(F)]) + macron( dd(hat(F))/dd(t))
$
从上面公式容易发现，力学量F平均值不变的充分条件是力学算符$hat(F)$本身不显含时间，同时$hat(F)$与哈密顿算符$hat(H)$对易。

这个条件也保证了力学量$F$可以和$H$（也就是能量）同时有确定值，或者说算符$hat(F)$和$hat(H)$有共同本征函数系。这时，描写算符$hat(F)$的本征值的那个量子数被称为*“好量子数”*。

*守恒量*的严格定义：*在任意态下（不管是否是定态）平均值和所有可能测值的概率都不随时间改变的物理量称为守恒量*。

== 对易与守恒量

*定理：如果力学量算符$hat(A)$不显含时间，且与也不显含时间的哈密顿算符$hat(H)$对易，则$A$为守恒量（$hat(A)$在任意态下的平均值和所有可能测值的几率都守恒，不随时间变化）*

我们可以选取其共同本征函数系$(psi_n)$来展开任一波函数
$
psi(arrow(r), t) = sum_n a_n (t) psi_n (arrow(r))\
hat(H) psi_n = E_n psi_n\
hat(A) psi_n = a_n psi_n
$

于是在$n$态的概率随时间的变化为：
$
dd("")/dd(t) (a_n^* (t) a_n (t))
$
而$a_n$可以由本征函数正交归一化性质得出：
$
a_n (t) = integral psi_n^* psi dd(tau)
$
所以
$
dd("")/dd(t) a_n (t) &= integral psi^*_n dd("")/dd(t) psi dd(tau)\
&= integral psi^*_n 1/(i hbar) hat(H) psi dd(tau)\
&= 1/(i hbar) integral (hat(H) psi_n)^* psi_n dd(tau)\
&= 1/(i hbar) E_n integral psi_n^* psi_n dd(tau)\
&= E_n / (i hbar) a_n (t)
$
同理有：
$
dd("")/dd(t) a_n^* (t) = - E_n / (i hbar) a_n^* (t)
$
所以
$
dd("")/dd(t) (a_n^* (t) a_n (t)) = 0
$
也就是说，不论体系处于什么量子态下，如果$hat(A)$与$hat(H)$对易，则不仅$hat(A)$的力学量平均值为常数，其可能测值的概率分布也不随时间变化。

这一点与$hat(H)$本身的概率分布守恒是一致的。

设$t=0$时刻波函数满足$hat(A) psi(x,0) = A psi(x, 0)$，即在$psi(x,0)$态下的测量值为$A$的概率为$100%$。如果是$hat(A)$守恒量，那么在后续$t$时刻也还是这样，即（其中$hat(U)_t = e^(-i/hbar hat(H) t)$）：
$
hat(A) hat(U)_t psi(x, 0) = A hat(U)_t psi(x, 0) =  hat(U)_t A psi(x, 0) = hat(U)_t hat(A) psi(x, 0)
$
从而有：
$
hat(A) hat(U)_t = hat(U)_t hat(A)
$
对$t->0$可以推出：
$
[hat(A), hat(U)_t] = [hat(A), e^(-i/hbar hat(H) t)] = [hat(A), hat(H)] = 0
$

== 守恒量与能级简并

1.*如果系统有两个彼此不对易的守恒量，则系统能级一般简并*

设$hat(A)$和$hat(B)$是守恒量，则它们分别都和$hat(H)$对易。如果$psi$是$hat(H)$的对应于能量$E$的本征态，则$hat(A) psi$和$hat(B) psi$也都是属于$E$的本征态。如果$psi$是$hat(A)$和$hat(H)$的共同本征态：
$
hat(A) psi = A psi\
$
则一般$psi$不会是$hat(B)$的本征态（不对易力学量算符不能拥有共同完备本征函数集），也就是说
$
hat(B) psi != B psi = B/A hat(A) psi
$
也就是说$hat(A) psi$和$hat(B) psi$线性无关，即$E$能级是简并的。

例如$hat(L)_x$和$hat(L)_y$是守恒量，但它们不对易，所以氢原子的能级是简并的。

2.*如果体系有一个守恒量$A$和一个非简并的能级$E$，则此能级对应的本征态也是$hat(A)$的本征态*

这一能级$E$对应的本征态为$psi$，则因为守恒量$hat(A)$和哈密顿算符$hat(H)$对易，所以$hat(A) psi$也是$E$的本征态。又因为能级$E$无简并，所以$hat(A) psi$和$psi$线性相关，即$psi$是$hat(A)$的本征态。

例如一维线性谐振子能级无简并，宇称算符与$hat(H)$对易，所以谐振子的能量本征态必有确定的宇称

$
[hat(P), hat(H)] = 0\
hat(P) hat(H)(x) psi(x) =  hat(H)(-x) psi(-x) = hat(H)(x) psi(-x) = hat(H)(x) hat(P) psi(x)
$

== 守恒量与定态

*守恒量是指任意态下平均值及其概率分布不随时间变化的力学量。*

*定态是指系统处于某一特定的能量本征态。*

- 如果力学量$A$是守恒量，那么不管系统处于定态与否，$A$都守恒。守恒量不一定取确定值。
- 如果$[hat(A),hat(H)]!=0$，则$A$的平均值一般会随时间变化，但在某些特殊态下也可能不变，例如一维谐振子基态的动量平均值。
- 如果系统处于定态，一切不显含时间的力学量都是守恒的。
- 如果$[hat(A),hat(H)]!=0$，而系统又处于非定态，则 A 的平均值一般随时间变化。

== 定态下力学量平均值——Virial定理

如果系统处于定态，则不显含时间的力学量都守恒。由于能量本征值单一，于是
$
macron(T) + macron(V) = macron(H) = E
$
其中$T$和$V$分别是哈密顿量中的动能和势能项，假设$V$仅仅是位置的函数。通过Virial定理，我们可以进一步定出$T$和$V$之间的关系。设$psi$为能级$E$的本征函数，则
$
(hat(T) + hat(V)) psi &= E psi\
sum_i hat(x)_i hat(p)_i (hat(T) + hat(V) - E) psi &= 0\
sum_i (hat(x)_i hat(T) hat(p)_i - hat(T) hat(x)_i hat(p)_i +hat(T) hat(x)_i hat(p)_i +  hat(x)_i hat(p)_i hat(V) - hat(x)_i hat(V) hat(p)_i + hat(V) hat(x)_i hat(p)_i - E hat(x)_i hat(p)_i) psi &= 0\
sum_i ([hat(x)_i, hat(T)] hat(p)_i + hat(T) hat(x)_i hat(p)_i + hat(x)_i [hat(p)_i, hat(V)] + hat(V) hat(x)_i hat(p)_i - E hat(x)_i hat(p)_i) psi &= 0\
sum_i ([hat(x)_i, hat(T)] hat(p)_i + hat(x)_i [hat(p)_i, hat(V)] + (hat(T) + hat(V) - E) hat(x)_i hat(p)_i) psi &= 0\
$
第三步假设了$hat(V)$和$hat(p)_i$是对易的、$hat(T)$和$hat(x)_i$是对易的，利用：
$
[hat(x)_i, hat(T)] = i hbar hat(p)_i/m \
[hat(p)_i, hat(V)] = - i hbar (partial hat(V))/(partial x_i)
$
得到：
$
sum_i (i hbar hat(p)^2_i/m - i hbar hat(x)_i (partial hat(V))/(partial x_i) + (hat(T) + hat(V) - E) hat(x)_i hat(p)_i) psi &= 0\
(hat(p)^2/m - sum_i hat(x)_i (partial hat(V))/(partial x_i)) psi = (E - hat(T) - hat(V)) i/hbar sum_i hat(x)_i hat(p)_i psi\
integral psi^* (hat(p)^2/m - sum_i hat(x)_i (partial hat(V))/(partial x_i)) psi dd(tau) = integral psi^*(E - hat(T) - hat(V)) i/hbar  sum_i hat(x)_i hat(p)_i psi dd(tau)\
macron(hat(p)^2/m) - macron(sum_i hat(x)_i (partial hat(V))/(partial x_i)) = (E - macron(T) - macron(V)) i/hbar macron(sum_i hat(x)_i hat(p)_i) = 0
$

最终得到*Virial定理*：
$
macron(T) = (1/2 macron(sum_i hat(x)_i partial/(partial x_i) hat(V)))
$

#newpara()

*一维谐振子利用Virial定理：*

$
V(x) = 1/2 m omega^2 x^2\
macron(T) = (1/2 macron(sum_i hat(x)_i partial/(partial x_i) hat(V))) = macron(V)
$

从而
$
macron(T) = macron(V) = 1/2 macron(H) = 1/2 E\
1/m macron(hat(p)^2) = m omega^2 macron(hat(x)^2) = E
$

对一维简谐振子定态
$
macron(hat(x)) = macron(hat(p)) = 0
$
有
$
1/m macron(Delta hat(p)^2) = m omega^2 macron(Delta hat(x)^2) = E\
1/m (Delta p)^2 = m omega^2 (Delta x)^2 = E\
Delta p = sqrt(m E), Delta x = sqrt(E/(m omega^2))\
Delta p Delta x = sqrt(m E) sqrt(E/(m omega^2)) = (n + 1/2) hbar
$
再一次印证了坐标-动量不确定关系，以及最小不确定态。


*氢原子利用Virial定理：*

对氢原子来说
$
V(arrow(r)) = - 1/(4 pi epsilon_0) e^2/arrow(r)
$
用Virial定理得到
$
macron(T) = (1/2 macron(sum_i hat(x)_i partial/(partial x_i) hat(V))) = - macron(V)
$
结合
$
macron(T) + macron(V) = macron(H) = E
$
得到
$
macron(V) = 2E , macron(T) = -E
$
这和Bohr轨道量子化计算结果一致。


== 波包的时间演化、Ehrenfest定理

知道了量子力学中力学量平均值随时间的变化规律，我们会问，这与经典力学中物体的运动规律有何不同？

在经典力学中，物体的运动规律体现于牛顿第二运动定律：
$
arrow(F) = m arrow(a) = dd(arrow(p))/dd(t)\
arrow(p) = m arrow(v) = m dd(arrow(r))/dd(t)
$
量子力学中哈密顿算符：
$
hat(H) = (hat(p)^2)/(2m) + V(hat(x))
$
于是有：
$
dd(arrow(r))/dd(t) = 1/(i hbar) macron([ hat(arrow(r)) , hat(H) ]) = 1/(i hbar) macron([hat(arrow(r)), hat(p)^2/(2m)]) = macron(arrow(p)/m)\
m dd(""^2 macron(arrow(r)))/dd(t)^2 = dd(macron(arrow(p)))/dd(t)\
dd(macron(arrow(p)))/dd(t) = 1/(i hbar) macron([hat(H), hat(arrow(p))]) = - macron(nabla V(arrow(r))) = macron(arrow(F) (arrow(r)))
$
于是得到Ehrenfest定理：
$
m dd(""^2 macron(arrow(r)))/dd(t)^2 = macron(arrow(F) (arrow(r)))
$
Ehrenfest定理与经典的牛顿方程极为相似。考虑一个波包的运动，如果其空间范围很窄，则方程右边可近似为
$
m dd(""^2 macron(arrow(r)))/dd(t)^2 = arrow(F) (macron(arrow(r)))
$
这是经典的结果。如果波包范围不是很窄，则与经典物理偏离较大。

=== 原子结构的探测

考虑$α$粒子被原子散射来探测原子结构，就需要：

- 在整个散射过程中$α$粒子的波包的大小远小于原子的尺度（1埃），粒子的德布罗意波长也要远小于原子的尺度
- 原子势场在波包大小的范围内变化不大，这样原子核对波包的平均作用力可以用原子核对波包中心的作用力代替
- 由于波包会随着时间扩散变大，所以由要求散射过程所需时间极短，使得在散射过程中波包本身大小变化不大
上述条件都满足的情况下，$α$粒子的散射就可以用经典力学的方法来处理（卢瑟福$α$粒子散射实验），或者说Ehrenfest定理适用，微观粒子波粒二象性中的粒子性性质占主导地位。
$
lambda = h/p = h/sqrt(2 m E)
$

=== 近似条件

在粒子波包足够窄的情况下，如果使定理的近似条件成立，还必须满足势场变化缓慢的条件。
$
F(x) = - (partial V(x))/(partial x) = - (partial V(macron(x)))/(partial macron(x)) - (partial^2 V(macron(x)))/(partial macron(x)^2) (x - macron(x)) + 1/2 (partial^3 V(macron(x)))/(partial macron(x)^3) (x - macron(x))^2 + ...
$
如果要满足
$
macron(F(x)) = F(macron(x))
$
则要求
$
1/2 abs((partial^3 V(macron(x)))/(partial macron(x)^3) (x - macron(x))^2) << abs((partial V(macron(x)))/(partial macron(x)))
$
如果势能函数最多包含坐标的二次幂（线性、谐振子势场），则这个条件是满足的。

#pagebreak(weak: true)

= ⼳正算符和体系对称性

== 幺正算符（酉算子）

如果线性算符$hat(U)$的逆算符$hat(U)^(-1)$存在，且对于任意两个波函数$psi$和$phi$，有
$
integral psi^* phi dd(tau) = integral (hat(U) psi)^* (hat(U) phi) dd(tau)
$
则称算符$hat(U)$为幺正算符（Unitary）。

幺正算符相当于对波函数做幺正变换，而不改变波函数的内积，保持了波函数的正交归一性(经典物理中坐标旋转变换)。

有性质：
$
hat(U)^dagger hat(U) = hat(U) hat(U)^dagger = I\
hat(U)^(-1) = hat(U)^dagger
$
- 若$hat(U)$幺正，则$hat(U)^dagger$也是幺正的。
- 若$hat(U)$和$hat(V)$都是幺正的，则$hat(U) hat(V)$也是幺正的。

== 生成元

幺正算符包括单位算符$I$，如果幺正算符依赖于一个连续变化的参量$epsilon$（如空间旋转、时间平移等），即$hat(U) = hat(U)(epsilon)$，有如下性质
$
hat(U)(0) &= I\
hat(U)(epsilon_1) hat(U)(epsilon_2) &= hat(U)(epsilon_1 + epsilon_2)\
$
则在$epsilon→0$时，$hat(U)$能展开为
$
hat(U)(epsilon) = I + i epsilon hat(F) + O(epsilon^2)
$
利用$hat(U)$的幺正性，可以得到$hat(F)$的厄密性：
$
hat(U)^dagger hat(U) = (I - i epsilon hat(F)^dagger) (I + i epsilon hat(F)) = I - i epsilon (hat(F) - hat(F)^dagger)= I\
hat(F) = hat(F)^dagger
$
称$hat(F)$为$hat(U)$的*生成元*。

如果$hat(U)$不是无限接近$I$的，我们可以通过$n$次无限小的幺正操作实现任意有限大小参量$a$的幺正变换：
$
hat(U)(a) = lim_(n→oo) (hat(U)(a/n))^n = lim_(n→oo) (I + i a/n hat(F))^n = e^(i a hat(F))
$
算符出现在指数上也可以通过泰勒展开式来理解：
$
hat(U) = e^(i a hat(F)) = sum_(n = 0)^oo (i a hat(F))^n/n!
$

#newpara()

例：时间演化算符就是一个幺正算符
$
hat(U)(t) = e^(-i/hbar hat(H) t), psi(arrow(r), t) = hat(U)(t) psi(arrow(r), 0)
$
通过幺正算符的性质，可以得到
$
hat(U)^dagger (t) = hat(U)^(-1) (t) = hat(U)(-t)
$
从而得到
$
integral psi^* (arrow(r), t) psi(arrow(r), t) dd(tau) &= integral (hat(U)(t) psi)^* (hat(U)(t) psi) dd(tau) \
&= integral psi^*(arrow(r), 0) hat(U)^dagger (t) hat(U)(t) psi(arrow(r), 0) dd(tau) \
&= integral psi^*(arrow(r), 0) psi(arrow(r), 0) dd(tau) \
&= 1
$

== 量子不可克隆原理

*不可能有这样一台设备，能够完美复制任意量子比特，而不对此量子比特产生干扰。*

证：设任意量子比特$A$的波函数为$psi(x)$，复制前某量子比特$B$的波函数为$phi(y)$，复制之后变为和$A$相同，即$psi(y)$。要改变$B$的量子状态，可以通过测量，或者调整系统哈密顿算符，使得经过时间演化后$B$和$A$的状态相等。测量可能会改变$A$的状态，所以我们使用时间演化算符$hat(U)(t)$ ，对系统$( A + B )$做一个幺正变换。

#grid(
    columns: 3,
    inset: 8pt,
    [初态],[$U(U) (t)$],[末态],
    [$psi(x) phi(y)$],[$->$],[$psi(x) psi(y)$],
    [$psi'(x) phi(y)$],[$->$],[$psi'(x) psi'(y)$]
)

其中$psi(x)$和$psi'(x)$是$A$的两个任意量子态。两个初态波函数内积
$
integral integral psi^*(x) psi'(x) phi^*(y) phi(y) dd(x) dd(y) = integral psi^*(x) psi'(x) dd(x) = a
$
而末态内积为
$
integral integral psi^*(x) psi'(x) psi^*(y) psi'(y) dd(x) dd(y) = a^2
$
一般来说$a^2 != a$ ，除非$a = 0,1$，与$psi$和$psi'$的任意性相悖，于是证明了量子不可克隆。

== 幺正算符和幺正变换

用幺正算符实现的波函数和算符的变换称为幺正变换：
$
psi &->  &psi' = hat(U) psi\
hat(A) &->  &hat(A') = hat(U) hat(A) hat(U)^dagger
$
与经典物理中的坐标变换相似，幺正变换不改变系统的物理规律（算符方程、对易关系、平均值及概率）：
$
hat(A) psi = phi => hat(A') psi' = phi'\
hat(A)' psi' = hat(U) hat(A) hat(U)^dagger hat(U) psi = hat(U) hat(A) psi = hat(U) phi =  phi'
$
强调：这里的变换是*波函数和算符同时变换*，如果只变换其中一个，则量子系统的物理就有可能改变。

== Fourier变换

傅里叶变换也可以看作是一种幺正变换：
$
hat(U) psi(x) = 1/sqrt(2 pi hbar) integral dd(x) e^(-i/hbar p x) psi(x) = phi(p)\
hat(U)^dagger phi(p) = 1/sqrt(2 pi hbar) integral dd(p) e^(i/hbar p x) phi(p) = psi(x)
$
可以证明：
$
integral psi^*(x) phi(x) dd(x) = integral phi^*(p) phi(p) dd(p) = 1
$
_证明：_
$
hat(U) hat(U)^(-1) &= 1/sqrt(2 pi hbar) integral dd(x) e^(-i/hbar p x) 1/sqrt(2 pi hbar) integral dd(p') e^(i/hbar p' x)\
&= integral dd(p') 1/(2 pi hbar) integral dd(x) e^(i/hbar (p' - p) x)\
&= integral dd(p') delta(p' - p)\
&= 1_(p' -> p)
$
注意这事实上是两个算符分别操作，后面的$hat(U)^(-1)$中的变量是$p'$而不是$p$， $1_(p' -> p)$表示单位算符，但是要把自变量$p'$换为$p$。

傅里叶幺正变换对动量算符的变换：
$
hat(U) hat(p^2)/(2m) hat(U)^dagger &= 1/sqrt(2 pi hbar) integral dd(x) e^(-i/hbar p x) hat(p^2)/(2m) 1/sqrt(2 pi hbar) integral dd(p') e^(i/hbar p' x)\
&= 1/sqrt(2 pi hbar) integral dd(x) e^(-i/hbar p x) (-hbar^2)/(2m) dd("")^2/dd(x)^2 1/sqrt(2 pi hbar) integral dd(p') e^(i/hbar p' x)\
&= 1/sqrt(2 pi hbar) integral dd(x) e^(-i/hbar p x) 1/sqrt(2 pi hbar) integral dd(p') p'^2/(2m) e^(i/hbar p' x)\
&= integral dd(p') p'^2/(2m) 1/(2 pi hbar) integral dd(x) e^(i/hbar (p' - p) x)\
&= integral dd(p') p'^2/(2m) delta(p' - p)\
&= p^2/(2m) 1_(p' -> p)
$

#newpara()

傅里叶幺正变换对坐标算符的变换：
$
hat(U) hat(x) hat(U)^dagger &= 1/sqrt(2 pi hbar) integral dd(x) e^(-i/hbar p x) x 1/sqrt(2 pi hbar) integral dd(p') e^(i/hbar p' x)\
&= 1/sqrt(2 pi hbar) i hbar dd("")/dd(p) integral dd(x) e^(-i/hbar p x) 1/sqrt(2 pi hbar) integral dd(p') e^(i/hbar p' x)\
&= i hbar dd("")/dd(p) integral dd(p') delta(p' - p)\
&= i hbar dd("")/dd(p) 1_(p' -> p)
$
推广得到
$
hat(U) hat(x)^n hat(U)^dagger = (i hbar dd("")/dd(p))^n 1_(p' -> p)
$
#newpara()
傅里叶幺正变换对哈密顿算符的变换：
$
hat(U) hat(H) hat(U)^dagger &= hat(U) (hat(p)^2/(2m) + V(hat(x))) hat(U)^dagger\
&= p^2/(2m) + V(i hbar dd("")/dd(p))
$
也就是说，在*坐标表象*中，哈密顿算符形式为
$
hat(H) = - hbar^2/(2m) dd("")^2/dd(x)^2 + V(x)
$
幺正变换到*动量表象*中，哈密顿算符形式为
$
hat(H) = p^2/(2m) + V(i hbar dd("")/dd(p))
$
#newpara()

在坐标表象下算符替换为：
$
cases(
    hat(p) -> - i hbar dd("")/dd(x),
    hat(x) -> x
)
$
而在动量表象下算符替换为：
$
cases(
    hat(p) -> p,
    hat(x) -> i hbar dd("")/dd(p)
)
$

== 态和力学量的表象

在量子力学中，描写量子态和力学量算符的方式不是唯一的。一种具体的方式称为一种表象。

一维空间态的表象：

用$psi(x,t)$来描写量子态是坐标表象。按动量本征函数展开
$
psi(x,t) = integral c(p,t) phi_p (x) dd(p) , phi_p (x) = 1/sqrt(2 pi hbar) e^(i/hbar p x)
$
就变换到了动量表象，$c(p, t)$称为动量表象中的波函数
$
c(p,t) = integral psi^*_p (x) psi(x,t) dd(x), psi^*_p (x) = 1/sqrt(2 pi hbar) e^(-i/hbar p x)
$

#newpara()
坐标表象的优点：
- 容易根据具体的物理问题的要求写出波函数满足的边界条件，分束缚态和散射态；根据粒子的入射方向写出入射波、透射波和反射波
- 一些常见的势在坐标表象下是定域的
- 容易讨论量子力学和经典力学的关系

坐标表象中的定态薛定鄂方程：
$
(- hbar^2/(2m) dd("")^2/dd(x)^2 + V(x)) psi (x) = E psi (x)
$
动量表象中的定态薛定鄂方程：
$
(p^2/(2m) + V(i hbar dd("")/dd(p))) phi (p) = E phi (p)
$
*表象之间的变换是一种幺正变换*。

=== 简谐振子的傅里叶变换

一维简谐振子的哈密顿算符为
$
hat(H) = hat(p)^2/(2m) + 1/2 m omega^2 hat(x)^2
$
在坐标表象中，哈密顿算符形式为
$
hat(H) = - hbar^2/(2m) dd("")^2/dd(x)^2 + 1/2 m omega^2 x^2
$
幺正变换到动量表象中后，其形式为
$
hat(H) = p^2/(2m) + 1/2 m omega^2 (i hbar dd("")/dd(p))^2 = p^2/(2m) - 1/2 m omega^2 hbar^2 dd("")^2/dd(p)^2
$

=== 坐标表象和动量表象

#figure(
  three-line-table[
    | 符号| 坐标表象| 动量表象|
    | --| --| --|
    |$hat(x)$ | $x$ | $i hbar dd("")/dd(p)$|
    |$hat(p)$ | $- i hbar dd("")/dd(x)$ | $p$|
    | $hat(x)$本征态 | $delta(x - x')$| $1/sqrt(2 pi hbar) e^(i/hbar p x)$|
    | $hat(p)$本征态 | $1/sqrt(2 pi hbar) e^(-i/hbar p x)$ | $delta(p - p')$| 
  ],
  caption: [
    坐标表象和动量表象
  ],
  kind: table
)

== 幺正变换与系统对称性

前面说过，对波函数和算符同时进行幺正变换，量子力学规律不变。但如果只变换波函数，则量子力学规律可能改变。

_把假设条件加强，如果只对波函数或算符二者其一进行幺正变换，而量子力学规律不变，会有什么物理结果?_

首先证明二者是等价的。薛定鄂方程：
$
i hbar partial/(partial t) psi = hat(H) psi
$
对$psi$进行幺正变换$psi -> psi' = hat(U) psi$，得到
$
i hbar partial/(partial t) psi' &= hat(H) psi'\
i hbar partial/(partial t) (hat(U) psi) &= hat(H) (hat(U) psi)\
$
用算符$hat(U)^(-1)$从左边作用于方程两边。因为我们一般考虑的幺正算符都是与时间无关的，所以$hat(U)^(-1)$可以越过时间偏导算符作用于右方
$
i hbar hat(U)^(-1) partial/(partial t) psi' &= hat(U)^(-1) hat(H) hat(U) psi\
$
与原薛定鄂方程作对比，同时注意到$psi$是薛定鄂方程的任意解，所以有
$
hat(H) = hat(U)^(-1) hat(H) hat(U)
$ 
也就是说，*只对波函数进行幺正变换而量子力学规律不变，可以等效为只对系统算符进行幺正变换而量子力学规律不变*。

哈密顿算符$hat(H)$幺正变换不变的意义：
$
hat(U) hat(H) hat(U)^(-1)=  hat(H) \
[hat(U), hat(H)] = 0\
[1 + i epsilon hat(F), hat(H)] = 0\
[hat(F), hat(H)] = 0
$
也就是说，如果*哈密顿算符幺正变换不变，那么此幺正变换对应的生成元是守恒量*。

*Noether定理的量子版本*：每当量子系统存在一种对称性（$hat(H)$幺正不变性），就相应的存在一个守恒律和守恒量。

=== 守恒和$hat(H)$幺正不变性

前面讲过，如果$hat(A)$是守恒量，则$[hat(A), hat(H)] = 0$。

#figure(
  image("pic/2024-04-12-14-52-21.png", width: 80%),
  caption: [
    量子系统的对称性与守恒量
  ],
)

在以为生成元的幺正算符 Ü =的幺正变换下不变。和以为生成元的幺正算符 Ü = e “对易百和对易。时间演化算符 Ü =声每 IIÄ 对易。在时间演化算符 Ü = e [的幺正变换下不变。

=== 时间均匀性和能量守恒

设时间幺正算符把波函数时间平移（主动）$tau$：
$
hat(U) (tau) psi(x, t) = psi(x, t - tau)
$
对变化后的波函数做泰勒展开：
$
psi(x, t - tau) &= sum^oo_(n=0) 1/n! (- tau dd("")/dd(t))^n psi(x, t)\
&= sum^oo_(n=0) 1/n! ((i tau)/hbar hat(H))^n psi(x, t)\
&= e^(i/hbar tau hat(H)) psi(x, t)
$
用到了薛定鄂方程($hat(H)$不含时)：
$
dd(psi)/dd(t) = hat(H)/(i hbar) psi\
(dd("")/dd(t))^n psi = (hat(H)/(i hbar))^n psi\
$

从而得到
$
psi(t - tau) = sum_(n=0)^oo (i/hbar tau hat(H))^n/n! psi(t) = e^(i/hbar tau hat(H)) psi(t)
$

得到时间平移算符：
$
hat(U) (tau) = e^(i/hbar tau hat(H))
$

另外时间演化算符和时间平移算符没有必然的联系，但是有时二者可以达到相同的效果。

设$t_1$时刻的波函数是能量本征态$psi(x,t_1) = phi_n(x) e^(-i E_n t_1/hbar)$，在这个态下测得$E_n$的概率为100%，在$t_2$时刻下：
$
psi(x, t_2) &=^"时间演化算符" e^(-i/hbar hat(H) (t_2 - t_1)) psi(x,t_1)\ 
&=^"时间平移算符" e^(i E_n (t_1 - t_2)/hbar) psi(x,t_1)\
&= psi(x, t_1 - (t_1 - t_2))\
&= psi(x,t_2)
$
根据$hat(A) psi = A psi$时，$f(hat(A)) psi = f(A) psi$，有
$
psi(x , t_2) = e^(i/hbar E_n (t_2 - t_1)) psi(x, t_1)
$
相当于乘了一个常数相位因子，所以在$psi(x,t_2)$下测得的能量$E_n$概率仍然是100%，系统能量守恒。


时间平移算符的生成元为$hat(H)$，它当然是与自身对易的，也就是说（$hat(H)$不显含时间时）系统的能量是个守恒量。

*时间平移不变性等价于系统能量守恒。*


=== 空间均匀性和动量守恒

设空间幺正算符把波函数坐标平移（主动）$arrow(a)$：
$
hat(U) (arrow(a)) psi(arrow(r)) = psi(arrow(r) - arrow(a))
$
对变化后的波函数做泰勒展开：
$
psi(arrow(r) - arrow(a)) &= sum^oo_(n=0) 1/n! (- arrow(a) dot nabla)^n psi(arrow(r))\
&= e^(- i arrow(a) dot nabla) psi(arrow(r))\
&= e^(- i/hbar arrow(a) dot hat(arrow(p))) psi(arrow(r))
$
这里$hat(p)$是动量算符。得到空间平移算符：
$
hat(U) (arrow(a)) = e^(- i/hbar arrow(a) dot hat(arrow(p)))
$
空间平移算符的生成元为动量算符。

*空间平移不变性等价于系统动量守恒。*

平移不变不是数学上波函数不变，而是带入Schrodinger方程后，波函数的形式不变。

#figure(
  image("pic/2024-04-17-13-46-55.png", width: 80%),
  caption: [
    空间均匀性和动量守恒
  ],
)

显然，一般情况下$hat(p)$和$hat(H)$算符并不对易（氢原子中的电子、一维谐振子），但如果考虑的是孤立系统，则系统总的$hat(p)$和$hat(H)$算符一定对易（自由粒子、电子+氢原子核总系统）。

自由粒子波包：平均动量守恒，动量的分布概率也守恒（也就是说$Δ p$不变），但不同动量平面波的传播速度不同，导致波包的$Δ x$随时间增大，形成波包的弥散（色散）。

相比之下，光在真空中传播则没有色散现象。

=== 空间各向同性和角动量守恒

设空间幺正算符把波函数绕$arrow(e)_n$旋转（主动）小角度$alpha$：
$
hat(U) (arrow(alpha)) psi(arrow(r)) = psi(arrow(r) - Delta arrow(r)), arrow(alpha) = alpha arrow(e)_n
$
对变化后的波函数做泰勒展开：
$
psi(arrow(r) - Delta arrow(r)) &= sum^oo_(n=0) 1/n! (- Delta arrow(r) dot nabla)^n psi(arrow(r))\
&= sum_(i=0)^oo 1/n! ( - (arrow(alpha) crossproduct arrow(r) ) dot nabla)^n psi(arrow(r))\
&= sum_(i=0)^oo 1/n! ( - arrow(alpha) dot (arrow(r) crossproduct nabla))^n psi(arrow(r))\
&= e^(- i/hbar arrow(alpha) dot  hat(arrow(L)))psi(arrow(r))
$
这里$hat(L)$是角动量算符。得到空间旋转算符：
$
hat(U) (arrow(alpha)) = e^(- i/hbar arrow(alpha) dot  hat(arrow(L)))
$

空间旋转算符的生成元为角动量算符。

*空间旋转不变性等价于系统角动量守恒。*

例：
$
e^(-i/hbar arrow(alpha) dot  hat(arrow(L))) Y_(l m) (theta, phi) = sum_(m = -l)^l c_m Y_(l m) (theta, phi)
$

对于中心力场问题（氢原子），哈密顿算符在空间转动变换下不变，因而角动量的三个分量都是守恒量。

能量守恒、动量守恒、角动量守恒都是时空对称性的体现，这在经典物理学中都有。但是，量子物理学还有经典中没有的更丰富的对称性，如空间反射和全同粒子交换对称性等——系统内禀对称性。

=== 空间反射对称性和宇称守恒

宇称（parity）算符$hat(P)$
$
hat(P) psi(arrow(r)) = psi(- arrow(r))
$
则有
$
hat(P)^2 = I\
hat(P)^(-1) = hat(P)\
hat(P)^dagger = hat(P)
$
宇称算符是一个厄米算符也是幺正算符。

由于$hat(P)^2 = 1$，$hat(P)$只有两个本征值$plus.minus 1$，对应的本征函数分别是对称波函数和反对称波函数。

根据本征函数完备性，任一波函数都可以展开为$hat(P)$的本征函数的叠加（对称和反对称部分）：
$
psi(arrow(r)) = psi_s(arrow(r)) + psi_a(arrow(r))
$
由于宇称是内禀的，所以没有经典对应力学量。

如果系统宇称守恒，且系统能级是非简并的，则系统能量本征态必有确定的宇称。从而一维束缚态势能对称情况下系统本征态必有确定的宇称宇称守恒的系统并不一定处于宇称的本征态。

对于多粒子系统，系统的总宇称是各部分相乘的，而能量等力学量的本征值是相加的。

对于核或粒子物理反应：
$
a + b -> c + d
$
系统初末态的总宇称为$P_a P_b P_(a b)$和$P_c P_d P_(c d)$，其中$P_a$为内禀宇称，而$P_(a b)$为轨道宇称。如果系统处于$Y_(l m)$态，则
$
P_(a b) = (-1)^l\
$
若反应过程宇称守恒（哈密顿量中相关势能项与宇称算符对易），则
$
P_a P_b (-1)^l = P_c P_d (-1)^l'
$
人们一般期待所有自界的基本相互作用力都是宇称不变的。但是在弱相互作用中，宇称守恒恰恰被彻底打破了。

*弱相互作用的宇称破坏*
$
arrow(r) ->^hat(P) - arrow(r)\
arrow(p) ->^hat(P) - arrow(p)\
arrow(L) ->^hat(P) arrow(L)\
arrow(mu) ->^hat(P) arrow(mu)
$

我们注意下面的几个命题。

系统处于$psi = c_1 phi_1 + c_2 phi_2 + c_3 phi_3 + ...$态，其中$phi_i$是$hat(P)$的本征态。那么系统处于$psi' = c'_1 phi_1 + c'_2 phi_2$态（$|c'_1|^2 + |c'_2|^2 = 1$）的概率为
$
|c'_1^* c_1 + c'_2^* c_2|^2
$
证明可以用基变换，把$psi' = c'_1 phi_1 + c'_2 phi_2$变成第一维度的基矢量。

另外若系统处于$psi$态，则$t$时间后，系统处于$psi'$态的概率为
$
abs(integral psi'^* e^(- i/hbar hat(H) t) psi dd(tau))^2
$

#newpara()

在$beta$衰变过程发生的概率幅就是
$
A = integral psi_f^* e^(i/hbar hat(H) t) psi_i dd(tau)
$
如果宇称守恒则有$[hat(P), hat(H)] = 0$，也就是$hat(P) hat(H)^n hat(P) = hat(H)^n$，所以
$
hat(P) e^(i/hbar hat(H) t) hat(P) = e^(i/hbar hat(H) t)
$
从而
$
A' = integral (hat(P) psi_f)^*  e^(i/hbar hat(H) t) (hat(P) psi_i )dd(tau) = integral psi_f^* e^(i/hbar hat(H) t) psi_i dd(tau) = A
$
其中$psi_i,psi_f$不必是本征态。从而$|A|^2 = |A'|^2$，是宇称守恒的必要条件。

#pagebreak(weak: true)

= 全同粒子体系

== 多粒子体系的描写

由$N$个离子组成的体系。体系的波函数应该和所有粒子的坐标以及时间有关：
$
psi(arrow(r)_1, arrow(r)_2, ..., arrow(r)_N, t)
$
“坐标”$q$包括粒子的空间坐标和自旋量子数（也许还有其它的“内部”量子数）。体系的Hamiltonian算符是：
$
hat(H) = sum_(i=1)^N (- hbar^2/(2m_i) nabla_i^2 + U_i(q_i)) + V(q_1, q_2, ..., q_N)
$
其中包括了各个粒子的动能之和，在外场中的势能$U$，以及粒子间的相互作用势能$V$，由此即可写出体系的Schrödinger方程。

== 全同粒子的不可区分性

假设多粒子体系中的 N 个粒子是全同粒子，即质量、电荷、总自旋等内在性质完全相同的粒子。全同粒子体系例如多电子原子中的电子、固体中的“公用"电子、原子核中的核子等在量子力学中，全同粒子体系与非全同粒子体系有更多的区别。

在经典力学中，即使两个粒子是全同的，它们也仍然是可区别的，因为它们各自有自己的轨道。

但是在量子力学中，粒子的状态用波函数描写，当两个粒子的波函数在空间中发生重叠的时候，我们无法区分哪个是“第一个”粒子，哪个是“第二个”粒子。所以在量子理论中有*全同粒子不可区别性原理*：当一个全同粒子体系中各粒子的波函数有重叠的时候，这些全同粒子是不可区别的。

== 波函数的交换对称性和粒子的统计性

粒子交换算符$hat(P)_(i j)$，交换粒子$i$和粒子$j$的坐标：
$
hat(P)_(i j) psi(arrow(r)_1, arrow(r)_2, ..., arrow(r)_N) = psi(arrow(r)_1, arrow(r)_2, ..., arrow(r)_j, ..., arrow(r)_i, ..., arrow(r)_N)
$
全同粒子的不可区别性告诉我们：这样交换以后的状态与原来的状态是不可区别的，所以，按照量子力学的基本原理，这两个状态应该是相同的，即
$
hat(P)_(i j) psi = C psi
$
得到$C = ±1$，这个$C$称为粒子的统计性。

- 如果$C = 1$，则称为玻色子，玻色子的波函数是对称的，满足波函数交换对称性。
- 如果$C = -1$，则称为费米子，费米子的波函数是反对称的，满足波函数交换对称性。

交换对称性或反对称性是全同粒子体系波函数的特殊的、固有的性质，因此也是（微观）粒子的特殊的、固有的性质。它决定了粒子所服从的统计规律。

- 自旋为整数的粒子，波函数是交换对称的，服从Bose-Einstein统计，称为玻色子。例如光子（自旋为1）、介子（自旋为0）
- 自旋为半整数的粒子，波函数是交换反对称的，服从Fermi-Dirac统计，称为费米子。例如电子、质子、中子（自旋都是ℏ/2）

原子核、原子、分子这样的粒子是由质子、中子、电子这些更“基本的”粒子组成的，我们把它们称为“复合粒子”。如果复合粒子的内部自由度是“冻结”的，我们也可以把它们看做是“基本”粒子。如果一个复合粒子包含偶数个费米子，那么它是玻色子；如果它包含奇数个费米子，那么它还是费米子。它所包含的玻色子的数目对此毫无影响。

事实上，这正是因为偶数个费米子的总自旋一定是整数，而奇数个费米子的总自旋一定是半整数，这一点可以由角动量的合成规则得到说明。

== 交换对称或反对称波函数的构成

一般地说，一个全同粒子体系的波函数是解 schrödinger 方程得到的，未必有确定的交换对称性。所以我们要对它进行“对称化”或“反对称化”。这里只考虑比较简单的情形：*无耦合体系*，即体系的总波函数是单个粒子波函数的乘积：
$
psi(q_1, ..., q_N) = psi_1 (q_1) psi_2 (q_2) ... psi_N (q_N)
$
这称为单粒子近似。

以二粒子体系为例，单粒子近似的波函数是：
$
psi(q_1, q_2) = psi_1(q_1) psi_2(q_2)
$
对称化的波函数是：
$
psi_s(q_1, q_2) = 1/sqrt(2) (psi_1(q_1) psi_2(q_2) + psi_2(q_1) psi_1(q_2))
$
反对称化的波函数是：
$
psi_a(q_1, q_2) = 1/sqrt(2) (psi_1(q_1) psi_2(q_2) - psi_2(q_1) psi_1(q_2))
$
对于可区别粒子（波函数$psi_1$或$psi_2$），我们可以说系统状态是“第一个粒子处于状态$psi_1$，第二个粒子处于状态$psi_2$”。但对于不可区别粒子（波函数$psi_S$或$psi_A$），我们只能说“有一个粒子处于状态$psi_1$，一个粒子处于状态$psi_2$”。

类似的做法可以推广到N个粒子的体系。特别是，一般的反对称化波函数是*Slater行列式：*
$
psi_a (q_1, q_2, ..., q_N) = 1/sqrt(N!) det mat(
    psi_1 (q_1), psi_2 (q_1), ..., psi_N (q_1);
    psi_1 (q_2), psi_2 (q_2), ..., psi_N (q_2);
    dots.v, dots.v, dots.down, dots.v;
    psi_1 (q_N), psi_2 (q_N), ..., psi_N (q_N)
)
$
在$psi_i$中有两个是相同的函数时：
$
psi_a (q_1, q_2, ..., q_N) = 0
$
*Pauli不相容原理*：不可能有两个或更多的费米子处于完全相同的量子状态中。这是量子力学基本公理之一，它在统计物理中起重要的作用。

例如，对于两个粒子经典中：
$
ket(<->)ket(arrow.t.b) , ket(arrow.t.b)ket(<->), ket(arrow.t.b)ket(arrow.t.b), ket(<->)ket(<->)
$
而在量子力学中的两个光子：
$
ket(<->)ket(arrow.t.b) + ket(arrow.t.b)ket(<->), ket(arrow.t.b)ket(arrow.t.b), ket(<->)ket(<->)
$
只有上面三种等概率的状态出现。

== 微观粒子波动性的表现

到目前为止，我们已经了解了量子力学的一系列与经典物理不同的表现：

- 粒子的运动由波函数决定，是几率性的，其动力学演化由薛定鄂方程决定
- 力学量测量值由波函数本征值决定，其平均值是相应力学量算符在波函数中的积分平均
- 力学量之间能否同时取确定值由力学量算符之间的对易关系决定；不能同时取确定值的情况导致不确定关系

这些特性都是粒子的波动性在*单粒子*身上的表现。*粒子的波动性反映在多粒子系统中就是全同性原理。*

正是由于全同性原理植根于波动性原理，它比其它原理（如系统反射对称性）显得更为基本。例：系统能量本征态不一定就是宇称本征态，但是任何多粒子系统不管处于什么态，在波函数叠加区域它一定处于全同粒子交换算符的本征态。

这个原理也可以带来超距作用的量子纠缠。

== 多粒子系统的算符

我们认为不同粒子的算符是对易的。

多粒子系统本身就带来了某些单粒子系统中没有的复杂性。例如对单粒子来说：
$
[hat(x), hat(p)] = i hbar
$
而对多粒子系统来说，由于不同粒子的坐标和动量是对易的，所以
$
[hat(x)_1 - hat(x)_2, hat(p)_1 + hat(p)_2] = 0
$
所以可以构造一个波函数，使得它是$hat(x)_1 - hat(x)_2$的本征函数，同时也是$hat(p)_1 + hat(p)_2$的本征函数。

在多粒子系统中，虽然属于不同粒子的力学量算符都相互对易，也就是说可以同时取确定的值，但是这些可能的测值之间却可以有某种关联（量子纠缠、非定域性、爱因斯坦之问），对这些问题的研究一直处于量子力学的前沿领域。

== Pauli不相容原理

最初不相容原理是Pauli综合反常塞曼效应、原子不同壳层电子数为偶数等现象归纳得到的，Pauli引入了电子自旋的概念来解释壳外电子的填充规律。

粒子物理后来发展中遇到了$Delta^(++)$粒子：
$
Delta^(++) = u u u
$
这个粒子由三个同样的顶夸克组成，电荷一样，自旋相同（同向），同时局限于一个狭小的空间之中（波函数重叠），似乎违背了Pauli不相容原理。

后来发现它们还有一个量子数不同-色电荷（color）。色是量子色动力学的基础，它有三种（不妨设为红、绿、蓝）。三个夸克拥有不同的色荷量子数，所以Pauli不相容原理没有被打破。

== 全同粒子的干涉效应   

两个自由粒子的空间波函数：

1. 不考虑全同性（非全同粒子波函数）
$
psi(arrow(r)_1, arrow(r)_2) = 1/(2 pi)^3 e^(i arrow(k)_1 dot arrow(r)_1) e^(i arrow(k)_2 dot arrow(r)_2)
$

引入质心系坐标$arrow(R) = (arrow(r)_1 + arrow(r)_2)/2$和相对坐标$arrow(r) = arrow(r)_1 - arrow(r)_2$，总波矢$arrow(K) = arrow(k)_1 + arrow(k)_2$，相对波矢$arrow(k) = (arrow(k)_1 - arrow(k)_2)/2$，则
$
psi(arrow(R), arrow(r)) = 1/(2 pi)^3 e^(i arrow(K) dot arrow(R)) e^(i arrow(k) dot arrow(r))
$
在以一个粒子为中心，半径$r→r+dd(r)$的球壳内找到另一个粒子的几率密度为：
$
P(r) = integral |psi(arrow(R), arrow(r))|^2 dd(""^3arrow(R)) r^2 dd(omega) = A/(4 pi) r^2 integral dd(omega) = A r^2
$
2. 两个全同玻色子
$
psi_+ (arrow(r)_1, arrow(r)_2) = 1/sqrt(2) (psi(arrow(r)_1) psi(arrow(r)_2) + psi(arrow(r)_2) psi(arrow(r)_1))
$
从而
$
psi_+ (arrow(R), arrow(r)) = 1/(2 pi)^3 e^(i arrow(K) dot arrow(R))1/sqrt(2) (e^(i arrow(k) dot arrow(r)) + e^(i arrow(k) dot arrow(r))) = 1/(2 pi)^3 e^(i arrow(K) dot arrow(R)) sqrt(2) cos(arrow(k) dot arrow(r))
$
在以一个粒子为中心，半径$r→r+dd(r)$的球壳内找到另一个粒子的几率密度为：
$
P(r) = integral |psi_+ (arrow(R), arrow(r))|^2 dd(""^3arrow(R)) r^2 dd(omega) = A/(4 pi) r^2 integral 2 cos^2(k r cos theta)dd(omega) = A r^2 (1 + (sin 2 k r)/(2 k r))
$
3. 两个全同费米子
$
psi_- (arrow(r)_1, arrow(r)_2) = 1/sqrt(2) (psi(arrow(r)_1) psi(arrow(r)_2) - psi(arrow(r)_2) psi(arrow(r)_1))
$
从而
$
psi_- (arrow(R), arrow(r)) = 1/(2 pi)^3 e^(i arrow(K) dot arrow(R))1/sqrt(2) (e^(i arrow(k) dot arrow(r)) - e^(i arrow(k) dot arrow(r))) = 1/(2 pi)^3 e^(i arrow(K) dot arrow(R)) sqrt(2) sin(arrow(k) dot arrow(r))
$
在以一个粒子为中心，半径$r→r+dd(r)$的球壳内找到另一个粒子的几率密度为：
$
P(r) = integral |psi_- (arrow(R), arrow(r))|^2 dd(""^3arrow(R)) r^2 dd(omega) = A/(4 pi) r^2 integral 2 sin^2(k r cos theta)dd(omega) = A r^2 (1 - (sin 2 k r)/(2 k r))
$

#figure(
  image("pic/2024-04-19-14-29-51.png", width: 80%),
  caption: [
    全同粒子的干涉效应
  ],
)
- 对称空间波函数 → 两粒子相互靠近的几率增大
- 反对称空间波函数 → 两粒子相互排斥的几率增大
似乎在全同粒子间存在一种作用力，对玻色子来说是吸引力，对费米子来说是排斥力。这种力称为交换力，它不是一种真正意义上的力，无施力者。在$r→∞$时，这种交换力消失。

== 全同粒子系统的量子特性

- 全同*玻色子*系统在低温下呈现*超流*效应——具有量子特性的宏观物体（玻色-爱因斯坦凝聚）
- 全同*费米子*系统在低温下呈现*超导*效应——电子之间两两结成*库派对*（复合玻色子）
- Pauli不相容原理使全同费米子体系无法聚集——导致日常物体占有的空间尺度
- 电子在白矮星内部提供简并压力抵抗重力崩塌，但当其质量大于1.4倍太阳质量时电子被压入质子内部形成中子星，中子星内部压强改由中子的简并提供

=== BSC理论

核心：电子与晶格振动的相互作用（电声子耦合）

- 两个动量相等、方向和自旋相反的电子，通过晶格振动的相互作用产生吸引，形成电子对的束缚态
- 电子对也受散射，但是成对出现的散射不改变总动量
- 温度$T$升高，或电流$I$增大，超过电子对的束缚能：变为正常态

#pagebreak(weak: true)


= 量子力学的矩阵形式与狄拉克(Dirac)符号

== 波函数的矩阵表示

力学量$hat(Q)$，设它的本征值是离散的，本征值集为${q_n}$本征函数系（不含时）为$u_n (x)$。

假设所有本征值都非简并，这个本征函数系的正交归一性就是：
$
(u_m, u_n) = delta_(m n)
$
如果是连续本征值系统，就是
$
(u_q ,u_q') = delta(q - q')
$
在$hat(Q)$表象中，态函数$ψ$可表示为态展开系数的列矩阵形式
$
psi = mat(
    a_1; a_2; dots.v; a_n; dots.v
)
$
它称为$hat(Q)$表象中的*态矢量*或*表示*。这就是系统在$hat(Q)$表象中的波函数。

也可以加入时间因子：
$
psi(t) = mat(
    a_1 (t) ; a_2 (t) ; dots.v ; a_n (t) ; dots.v
)
$
厄密共轭态矢量排成行矩阵的形式：
$
psi^dagger = mat(
    a_1^* , a_2^* , dots , a_n^* , dots
)
$
内积则可定义为：
$
(psi, phi) = psi^dagger phi = sum_n a_n^* b_n
$
我们称：
- $psi$为态矢量
- $u$为表象的基底或基矢
- $a$是态矢量的分量或投影

== 算符矩阵表示

坐标表象中算符表示为
$
hat(F) (x, - i hbar partial/(partial x))
$
算符作用式$phi(x) = hat(F) psi(x)$变换到$hat(Q)$表象中为
$
psi(x ,t )= sum_n a_n (t) u_n (x)\
phi(x ,t )= sum_n b_n (t) u_n (x)
$
带入上面的方程有：
$
sum_n b_n (t) u_n (x) = sum_n a_n (t) hat(F) u_n (x)
$
左乘$u_m^* (x)$再做内积
$
b_m (t) = sum_n a_n (t) (u_m^* , hat(F) u_n)
$
记
$
F_(m n) = (u_m^* , hat(F) u_n)
$
则有
$
b_m (t) = sum_n F_(m n) a_n (t)
$
写成矩阵形式有
$
mat(
    b_1; b_2; dots.v
)
= mat(
    F_(1 1), F_(1 2), dots.v;
    F_(2 1), F_(2 2), dots.v;
    dots.v, dots.v, dots.down;
)
mat(
    a_1; a_2; dots.v
)
$
即
$
phi = F psi
$
这就是算符$hat(F)$在$hat(Q)$表象中的矩阵表示。

算符的Hermitian性质要求
$
F_(m n)^* = (u_n , hat(F) u_m) = (hat(F) u_m , u_n) = F_(n m) = F_(m n)^TT
$
即：
$
F^* = F^TT , F^dagger = F
$
这就是说，*算符$hat(F)$在$hat(Q)$表象中的矩阵是厄密的*。

恒等算符在$hat(Q)$表象中的矩阵表示是单位矩阵。

_$hat(L)_x$在动量表象中的矩阵元_

在坐标表象下，考虑$hat(L)_x$的矩阵元，将其作用于$hat(p)$的本征函数$psi_arrow(p) = e^(i arrow(p) dot arrow(r))$上：
$
(L_x)_(arrow(p') arrow(p)) &= integral psi_(arrow(p'))^* hat(L)_x psi_(arrow(p)) dd(arrow(r))\
&= integral psi_(arrow(p'))^* (y hat(p)_z - z hat(p)_y) psi_(arrow(p)) dd(arrow(r))\
&= 1/(2 pi hbar)^3 integral e^(-i/hbar arrow(p') dot arrow(r)) (y (- i hbar partial/(partial z)) - z (- i hbar partial/(partial y))) e^(i/hbar arrow(p) dot arrow(r)) dd(arrow(r))\
&= p_z (-i hbar partial/(partial p_y)) integral e^(-i/hbar arrow(p') dot arrow(r)) e^(i/hbar arrow(p) dot arrow(r)) dd(arrow(r)) - p_y (-i hbar partial/(partial p_z)) integral e^(-i/hbar arrow(p') dot arrow(r)) e^(i/hbar arrow(p) dot arrow(r)) dd(arrow(r))\
&= - i hbar (p_z partial/(partial p_y) - p_y partial/(partial p_z)) delta(arrow(p') - arrow(p))\
&= i hbar (p'_z partial/(partial p'_y) - p'_y partial/(partial p'_z)) delta(arrow(p') - arrow(p))
$
当然，我们可以考虑直接在动量表象下计算$hat(L)_x$的矩阵元，动量表现下的$hat(arrow(p))$的本征函数是$psi_arrow(p) = delta(arrow(p) - arrow(p'))$，所以
$
(L_x)_(arrow(p') arrow(p)) &= integral delta(arrow(tilde(p)) - arrow(p')) i hbar (tilde(p)_z partial/(partial tilde(p)_y) - tilde(p)_y partial/(partial tilde(p)_z)) delta(arrow(tilde(p)) - arrow(p)) dd(arrow(tilde(p))) \
&= i hbar (p'_z partial/(partial p'_y) - p'_y partial/(partial p'_z)) delta(arrow(p') - arrow(p))
$
虽然计算的是动量表象的矩阵元，计算公式中算符和波函数的表象可以任意选择，最终结果也不依赖于这些选择。


== 表象变换

仍以一维情形为例。

设我们再取另一个与算符$hat(Q)$函数独立的算符$hat(R)$，求出它的本征值集${r_n}$和本征函数系${u'_n (x)}$，我们就构造了$hat(R)$表象。

原来的基底${u_n (x)}$也可以用新的基底${u'_n (x)}$来展开：

$
mat(
    u_1, u_2, ...
)
=
mat(
    u'_1, u'_2, ...
)
mat(
    S_(1 1), S_(1 2), ...;
    S_(2 1), S_(2 2), ...;
    dots.v, dots.v, dots.down;
)
$
其中
$
S_(m n) = (u'_m , u_n)
$
如果一个态矢量$psi$在$hat(Q)$表象中的分量为${a_n}$，在$hat(R)$表象中的分量为${a'_n}$ ，则有
$
psi = sum_n a_n u_n = sum_n sum_m a_m u'_m S_(m n) = sum_m (sum_n a_m S_(m n)) u'_m = sum_m a'_m u'_m
$
所以有
$
a'_m = sum_n S_(m n) a_n
$
写成矩阵有
$
mat(
    a'_1; a'_2; dots.v;
)
=
mat(
    S_(1 1), S_(1 2), ...;
    S_(2 1), S_(2 2), ...;
    dots.v, dots.v, dots.down;
)
mat(
    a_1; a_2; dots.v;
)
$
即
$
psi' = S psi
$
注意*基底变换是行矩阵的形式，而态矢量是列矩阵的形式*。


下面讨论$S$应该满足的条件。

考虑到*态矢量*的模方为可观测量，应该要求其*内积*在表象变换下保持不变。

矢量$psi$和$phi$在$hat(Q)$表象中的分量分别是${a_n}$和${b_n}$，在$hat(R)$表象中的分量分别是${a'_n}$和${b'_n}$，则有
$
(psi, phi ) = sum_n a_n^* b_n = sum_m sum_j sum_k S^*_(m j) a_j^* S_(m k) b_k = sum_j sum_k (sum_m S^*_(m j) S_(m k)) a_j^* b_k
$
又因为
$
(psi, phi ) =  sum_k a_k^* b_k = sum_j sum_k delta_(j k) a_j^* b_k
$
对比得到
$
sum_m S^*_(m j) S_(m k) = delta_(j k)
$
或者
$
sum_m S_(j m)^dagger S_(m k) = (S^dagger S)_(j k) = delta_(j k)
$
即
$
S^dagger S = I
$
这就是说，*表象变换矩阵$S$是幺正的*。

也可以得到
$
mat(
    u'_1; u'_2; dots.v;
)
= 
mat(
    S_(1 1), S_(1 2), ...;
    S_(2 1), S_(2 2), ...;
    dots.v, dots.v, dots.down;
)^*
mat(
    u_1; u_2; dots.v;
)
$
即
$
u'_m = sum_n S_(m n)^* u_n\
u_m = sum_n S_(n m)^TT u'_n\
a'_m = sum_n S_(m n) a_n\
a_m = sum_n S_(n m)^dagger a'_n
$
其中
$
S_(m n) = (u'_m , u_n)
$
注意：*基底的变换矩阵和态矢的变换矩阵互为复共轭或转置。*

在表象变换下，一个算符所对应的矩阵的变换是
$
F' = S F S^dagger = S F S^(-1)
$
幺正变换不改变任何量了力学方程。即，如果$phi = F psi$，则$phi' =S phi = S F psi = F' psi'$。


== 量子力学的矩阵形式

坐标表象与离散表象的关系和对比如下表

#figure(
  three-line-table[
    || 坐标表象 | 离散表象|
    |--|---|---|
    |态| 波函数$psi(x,t)$ 复共轭波函数$psi^*(x,t)$ | 行矢量$ket(psi)$ 列矢量$bra(psi)$ |
    |算符| $hat(F)(x, - i hbar partial_x)$ | 矩阵$F = (F_(m n))$ |
    |算符作用到态| $hat(F) psi(x,t)$ | $F ket(psi)$ |
    |态的内积| $(psi, phi) = integral psi^* phi dd(x)$ | $bra(psi) ket(phi)$ |
  ],
  caption: [
    坐标表象与离散表象的关系和对比
  ],
  kind: table
)

1. 态的归一：$psi^dagger psi=1$, 两态正交：$phi^dagger psi=0$
2. 力学量的平均值（若 $psi$ 已归一）：$ macron(F) = psi^dagger F psi$
3. 本征方程：$hat(F) psi = lambda psi$
4. 含时间的薛定鄂方程：$i hbar partial/(partial t) psi = H psi$

=== 离散表象中的本征方程的解法

设
$
psi = mat(
    a_1; a_2; dots.v;
)
$
$
F = mat(
    F_(1 1), F_(1 2), dots.v;
    F_(2 1), F_(2 2), dots.v;
    dots.v, dots.v, dots.down;
)
$
本征方程：
$
F psi = lambda psi\
(F - lambda I) psi = 0
$
这是一个齐次线性方程组，它有非零解的充要条件是：
$
det(F - lambda I) = 0
$
即*久期方程*。

如果$F$是$n×n$矩阵，则是关于$lambda$的$n$次多项式方程。根据“代数基本定理”，在复数域内，$n$次代数方程一定有$n$个根，这些根就是本征值。另外，矩阵$F$的的厄密性保证了这些根都是实数。

把这些本征值记为${lambda_i}$, 再代回方程，假设没有重根
$
mat(
    F_(1 1) - lambda_1, F_(1 2), dots.v;
    F_(2 1), F_(2 2) - lambda_1, dots.v;
    dots.v, dots.v, dots.down;
)
mat(
    a_1; a_2; dots.v;
)
= 0
$
就可以对各个本征值求出${a_i}$，但有一个整体的常数因子未定，再利用归一化条件把它定出，就得到了完全归一化的本征态矢量。

【补充，对角化矩阵的一些性质】

我们已经求出了属于本征值${lambda_i}$的本征态矢量${psi_i}$，把他们排成一个矩阵
$
S = mat(
    psi_1, psi_2, dots;
)
$
它的Hermitian共轭矩阵$S^dagger$是
$
S^dagger = mat(
    psi_1^dagger;
    psi_2^dagger;
    dots.v;
)
$
则有
$
S^dagger S = mat(
    psi_1^dagger psi_1, psi_1^dagger psi_2, dots;
    psi_2^dagger psi_1, psi_2^dagger psi_2, dots;
    dots.v, dots.v, dots.down;
)
= mat(
    1, 0, dots;
    0, 1, dots;
    dots.v, dots.v, dots.down;
) = I
$
即$S$是幺正矩阵。

对$F$做幺正变换：
$
S^dagger F S = mat(
    psi_1^dagger;
    psi_2^dagger;
    dots.v;
)
F mat(
    psi_1, psi_2, dots;
)
=  mat(
    psi_1^dagger;
    psi_2^dagger;
    dots.v;
)
mat(
    lambda_1 psi_1, lambda_2 psi_2, dots;
)
\
= mat(
    lambda_1 psi_1^dagger psi_1, lambda_2 psi_1^dagger psi_2, dots;
    lambda_1 psi_2^dagger psi_1, lambda_2 psi_2^dagger psi_2, dots;
    dots.v, dots.v, dots.down;
)
= mat(
    lambda_1, 0, dots;
    0, lambda_2, dots;
    dots.v, dots.v, dots.down;
)
$
求$F$本征值问题归结为寻找一个幺正变换把$F$从某个表象变换到其自身表象，使$F$的矩阵表示对角化。

幺正变换不改变矩阵F的秩、迹和行列式。

== Dirac符号

不同的量子力学表象所表达的物理内容是完全相同的，但是从表面上看来，不同表象中的算符和量子态具体表达式却可能很不一样。为了避免不同表象带来的形式上的差异，Dirac引入了一种与表象无关的符号体系，被称为Dirac符号。

=== 态

量子体系的状态用态矢量表示。*态矢量*有
- 左矢 bra $bra(psi)$
- 右矢 ket $ket(psi)$
有关系
$
bra(psi) = ket(psi)^dagger, ket(psi) = bra(psi)^dagger
$
这里可把$dagger$看成一种“形式运算符号”，即矩阵力学中的转置加复共轭。

两个态的*内积*（即过去定义的$(psi, phi)$）用
$
braket(psi,phi)
$
是一个数，满足关系
$
braket(psi,phi)^* = braket(phi,psi)^dagger = braket(psi,phi)
$
互为共轭复数。所以内积的性质可以写成：
$
braket(psi,psi) >= 0
$
等号成立当且仅当$ket(psi) = 0$。并且，态的归一是
$
braket(psi,psi) = 1
$
态的正交是
$
braket(psi,phi) = 0
$

=== 算符

算符（例如$hat(F)$）对右矢的作用直接写为
$
hat(F) ket(psi) = ket(phi)
$
结果还是一个右矢。对左矢的作用写为
$
bra(psi) hat(F) = ((bra(psi) hat(F))^dagger)^dagger = (hat(F)^dagger ket(psi))^dagger = bra(phi)
$
结果还是一个左矢。其中$hat(F)^dagger$是$hat(F)$的Hermitian共轭，满足
$
(bra(psi) hat(F) ket(phi))^dagger = bra(phi) hat(F)^dagger ket(psi) 
$
对于任意的态矢量$ket(psi)$和$ket(phi)$成立。

而
$
ket(phi) bra(psi)
$
是一个算符。

如果有
$
hat(F) ket(psi) = ket(phi)
$
则有
$
bra(psi) hat(F)^dagger = bra(phi)
$
#newpara()

算符乘积的Hermitian共轭是
$
(hat(F) hat(G))^dagger = hat(G)^dagger hat(F)^dagger
$
如果算符$hat(F)$具有性质
$
hat(F) = hat(F)^dagger
$
则称$hat(F)$是Hermitian算符。

对于Hermitian算符，有
$
(bra(phi) hat(F) ket(psi) )^dagger = (bra(phi) hat(F) ket(psi) )^* = bra(psi) hat(F) ket(phi)
$
从而力学量的平均值
$
macron(F) = (bra(psi) hat(F) ket(psi)) / (braket(psi,psi))
$
是实数。

基底集合$\{ket(n)\}$是正交归一的，即
$
braket(m,n) = delta_(m n)
$
完备性可以写为
$
sum_n ket(n) bra(n) = I
$
上面中的一项被称作*投影算符*：
$
P_n = ket(n) bra(n)
$
称为处于态$ket(n)$的投影算符。有性质
$
P_n^2 = P_n
$
对于连续谱，狄拉克态矢的正交归一表示为
$
braket(lambda_1, lambda_2) = delta(lambda_1 - lambda_2)
$
比如坐标算符$x$的本征方程为：
$
x delta(x - x_0) = x_0 delta(x - x_0)
$
狄拉克符号表示
$
x ket(x_0) = x_0 ket(x_0)
$
有：
$
braket(x_0, x_1) = delta(x_0 - x_1)
$

#newpara()

一个抽象的态$ket(n)$在坐标表象中的函数表示
$
braket(x, n) = psi_n (x)
$

_例：动量为$p'$的平面波在坐标和动量表象的函数表示为_
$
braket(x, p') = 1/(2 pi hbar) e^(i p' x / hbar)\
braket(p, p') = delta(p - p')
$
#newpara()
基底完备性条件用狄拉克符号的表达：
$
sum_n ket(n) bra(n) = I
$

== 态矢量在具体表象中的表示

在$F$表象中（基矢量$ket(k)$），任何一个态矢量$ket(psi)$都可以用基矢量展开：
$
ket(psi) = sum_k ket(k) braket(k, psi) = sum_k a_k ket(k)
$
其中$a_k = braket(k, psi)$是态矢量在$ket(k)$表象中的分量。

${a_k} = {braket(k, psi)}$是态矢量$ket(psi)$在$ket(k)$表象中的表示
$
mat(
    a_1; a_2; dots.v;
)=
mat(
    braket(1, psi); braket(2, psi); dots.v;
)
$
在基底的量子数为连续谱时，完备性关系表示为
$
integral dd(x) ket(x) bra(x) = I\
integral dd(p) ket(p) bra(p) = I
$

在具体的$F$表象下，态矢量展开为：
$
ket(psi) = sum_k a_k ket(k) = sum_k braket(k, psi) ket(k)\
ket(phi) = sum_j b_j ket(j) = sum_j braket(j, phi) ket(j)
$
态矢量的内积为
$
braket(phi, psi) = sum_k braket(phi, k) braket(k, psi) = sum_k b_k^* a_k
$

== 算符在具体表象下的表示

算符代表着对态的一种运算：
$
hat(L) ket(psi) = ket(phi)
$
在$F$表象中，
$
bra(j) hat(L) ket(psi) = sum_k bra(j) hat(L) ket(k) braket(k, psi) = braket(j, phi)
$
即：
$
sum_k L_(j k) a_k = b_j
$
其中
$
L_(j k) = bra(j) hat(L) ket(k)
$
就是算符$hat(L)$在$F$表象中的矩阵元。

算符$hat(L)$的狄拉克符号表示为：
$
hat(L) = sum_(j k) L_(j k) ket(j) bra(k) = sum_(j k) ket(j) bra(j) hat(L) ket(k) bra(k)
$
算符$hat(F)$在其自身$F$表象中的矩阵元和狄拉克符号表示为：
$
F_(m n) = bra(m) hat(F) ket(n) =  bra(m) f_n ket(n) = f_n delta_(m n)\
hat(F) = sum_n f_n ket(n) bra(n)
$
其中$f_n$是$hat(F)$在$ket(n)$表象中的本征值。*任何算符在其自身表象中自然就是对角化的*。

= 中心力场中的运动和氢原子

== 牛顿力学中心力场问题

中心力场中运动的粒子角动量守恒：
$
dd("")/dd(t) arrow(L) &= dd("")/dd(t) (arrow(r) crossproduct arrow(p))\
&= (dd("")/dd(t) arrow(r)) crossproduct arrow(p) + arrow(r) crossproduct (dd("")/dd(t) arrow(p))\
&= arrow(v) crossproduct arrow(p) + arrow(r) crossproduct arrow(F)\
&= arrow(r) crossproduct arrow(F)\
$
在中心力场中，力$arrow(F)$与$arrow(r)$同向（排斥力）或反向（吸引力），所以叉乘结果为0，即
$
U(arrow(r)) = U(r) => arrow(F) = -nabla U(r) = F hat(r) => dd("")/dd(t) arrow(L) = 0
$
又由于$r × L = 0$，所以粒子运动平面恒垂直于$arrow(L)$，即粒子在中心力场中的运动是平面运动。

== 中心力场两体问题化为单体问题

中心力场中两体问题的定态薛定谔方程：
$
(-hbar^2/(2 m_1) nabla^2_1 -hbar^2/(2 m_2) nabla^2_2+ U(abs(arrow(r)_1 - arrow(r)_2))) Psi(arrow(r)_1, arrow(r)_2) = E_"tot" Psi(arrow(r)_1, arrow(r)_2)
$
两体问题系统总能量$E_"tot"$可分为整体*平动*动能和*相对运动*能量（相对动能+势能）两部分。为了进行这种分解，把两粒子坐标进行转化：
$
arrow(R) = (m_1 arrow(r)_1 + m_2 arrow(r)_2)/(m_1 + m_2) "（质心系坐标）"\ 
arrow(r) = arrow(r)_1 - arrow(r)_2 "（相对坐标）"
$
质心系坐标满足：
$
M = m_1 + m_2\
M dd(""^2)/dd(t^2) arrow(R) = m_1 dd(""^2)/dd(t^2) arrow(r)_1 + m_2 dd(""^2)/dd(t^2) arrow(r)_2\
$
进行坐标转化：
$
arrow(r)_1 , arrow(r)_2 => arrow(R), arrow(r)\
x_1 , x_2 => X, x\
y_1, y_2 => Y, y\
z_1, z_2 => Z, z
$
就有
$
partial/(partial x_1) = (partial X)/(partial x_1) partial/(partial X) + (partial x)/(partial x_1) partial/(partial x) = m_1 /M partial/(partial X) + partial/(partial x)\
partial^2/(partial x_1^2) = (m_1 /M partial/(partial X) + partial/(partial x))(m_1 /M partial/(partial X) + partial/(partial x)) = m_1^2 /M^2 partial^2/(partial X^2) + 2 m_1 /M partial^2/(partial X partial x) + partial^2/(partial x^2)
$
同理有：
$
partial/(partial x_2) = (partial X)/(partial x_2) partial/(partial X) + (partial x)/(partial x_2) partial/(partial x) = m_2 /M partial/(partial X) - partial/(partial x)\
partial^2/(partial x_2^2) = (m_2 /M partial/(partial X) - partial/(partial x))(m_2 /M partial/(partial X) - partial/(partial x)) = m_2^2 /M^2 partial^2/(partial X^2) - 2 m_2 /M partial^2/(partial X partial x) + partial^2/(partial x^2)
$
于是：
$
1/m_1 partial^2/(partial x_1^2) + 1/m_2 partial^2/(partial x_2^2) &= 1/m_1 m_1^2 /M^2 partial^2/(partial X^2) + 1/m_2 m_2^2 /M^2 partial^2/(partial X^2) + 1/m_1 partial^2/(partial x^2) + 1/m_2 partial^2/(partial x^2) \
&= 1/M partial^2/(partial X^2) + 1/((m_1 m_2)/M) partial^2/(partial x^2) \
&= 1/M partial^2/(partial X^2) + 1/mu partial^2/(partial x^2)
$
其中$mu = (m_1 m_2) /M$是约化质量。三维的情况则是：
$
1/m_1 nabla^2_1 + 1/m_2 nabla^2_2 = 1/M nabla^2_R + 1/mu nabla^2_r
$
于是定态薛定谔方程转换为：
$
(-hbar^2/(2 M) nabla^2_R - hbar^2/(2 mu) nabla^2_r + U(r)) psi(arrow(R), arrow(r)) = E_"tot" psi(arrow(R), arrow(r))
$
分离变量求特解：
$
psi(arrow(R), arrow(r)) = phi(arrow(R)) psi(arrow(r))
$
代入原方程得：
$
- hbar^2/(2 M) nabla^2_R phi(arrow(R)) = E_"cm" phi(arrow(R))\
(- hbar^2/(2 mu) nabla^2_r + U(r)) psi(arrow(r)) = (E_"tot" - E_"cm") psi(arrow(r))
$
其中$E_"cm"$是质心运动的能量，$E = E_"tot" - E_"cm"$是相对运动的能量。

第一个方程的解即自由粒子平面波：
$
phi(arrow(R)) = c_1 e^(i/hbar arrow(P) dot arrow(R)) + c_2 e^(-i/hbar arrow(P) dot arrow(R))
$
其中$P = sqrt(2 M E_"cm")$是质心动量。第二个方程的解才是我们关心的中心势场问题的解，其中能量$E$表示相对运动的总能量（相对运动动能+势能）：
$
(- hbar^2/(2 mu) nabla^2_r + U(r)) psi(arrow(r)) = E psi(arrow(r))
$
如果$U(r -> oo) > E$，那么这个方程的能量本征值是分立的。

- 对于氢原子来说$U(r) tilde - 1/r$，也就说如果电子处于原子核库伦势吸引下的束缚态，则相对运动总能量$E<0$
- 对于库伦散射来说$U(r) tilde ± 1/r$，也就说如果带电粒子处于散射态，则相对运动总能量$E>0$
- 对于三维无限深球势阱：
  $
  U(r) = cases(
    0 (r < R),
    oo (r > R)
  )
  $
  其能级为正，并且分立
- 对于三维各向同性谐振子$U(r) tilde r^2$其能级为正，且将一定处于束缚态（能级分立）

== 中心力场中的运动

=== 中心力场中Schrödinger方程的约化

$
(-hbar^2/(2 mu) nabla^2_r + U(r)) psi(arrow(r)) = E psi(arrow(r))
$
中心力场的势能函数与方向无关：
$
U(arrow(r)) = U(r)
$
在球坐标系中，
$
nabla^2_r = 1/r^2 partial/(partial r) (r^2 partial/(partial r)) + 1/(r^2 sin(theta)) partial/(partial theta) (sin(theta) partial/(partial theta)) + 1/(r^2 sin^2(theta)) partial^2/(partial phi^2)\
- hbar^2/(2 mu) nabla^2 = - hbar^2/(2 mu r^2) partial/(partial r) (r^2 partial/(partial r)) + 1/(2 mu  r^2) hat(L)^2
$
其中$hat(L)^2$是角动量平方算符。方程化为
$
(- hbar^2/(2 mu r^2) partial/(partial r) (r^2 partial/(partial r)) + 1/(2 mu  r^2) hat(L)^2 + U(r)) psi(arrow(r)) = E psi(arrow(r))
$
第一项是径向动能$hat(p)^2_r/(2 mu)$，$hat(p)_r = - i hbar (partial/(partial r) + 1/r)$，第二项是离心势能，角动量越大，离心势能就越大。

在球坐标系中分离变量：
$
psi(arrow(r)) = R(r) Y(theta, phi)\
hat(L)^2 Y_(l m) (theta, phi) = l (l + 1) hbar^2 Y_(l m) (theta, phi)
$
径向方程：
$
(- hbar^2/(2 mu r^2) dd("")/(dd(r)) (r^2 dd("")/(dd(r))) + l (l + 1) hbar^2/(2 mu r^2) + U(r)) R_l (r) = E R_l (r)
$
设$R(r) = u(r)/r$，则有
$
- hbar^2 / (2 mu) dd(""^2u) /dd(r^2)  + ((l (l + 1) hbar^2) / (2 mu r^2) + U(r)) u = E u
$
它很像是一维的方程，但是有$2l+1$重简并。有两点讨论：
1. “势能”是
   $
   U_"eff" (r) = U(r) + (l (l + 1) hbar^2) / (2 mu r^2)
   $
   它称为有效势能，包括*库伦势能*和*离心势能*。
2. 自变量区间是$0 <= r < +oo$而且波函数满足边界条件：
   $
   lim_(r -> 0) r R_l (r) = lim_(r ->0) u(r) = 0
   $
   下面解释这个边界条件的物理意义。

首先，在原点附近，$U(r) → ∞$的速度不能比$1/r^2$更快，即当$r→0$时，$r^2 U(r) -> 0$。通常碰到的中心力场都满足这个条件。例如
- 谐振子势：$U(r) prop r^2$
- 线性中心势：$U(r) prop r$
- 无限深球势阱
- coulomb势：$U(r) prop 1/r$
- 汤川势：$U(r) prop e^(-alpha r)/r$
径向方程：
$
dd(""^2R)/dd(r^2) + 2/r dd(R)/dd(r) + ((2 mu) / hbar^2 (E - U(r)) - (l(l+1))/r^2 )R= 0
$
当$r→0$时，上式的渐进式是：
$
dd(""^2R)/dd(r^2) + 2/r dd(R)/dd(r) - (l(l+1))/r^2 R= 0
$
在正则奇点$r=0$的邻域内，设$R∝ r^s$，代入后得
$
s(s-1) r^(s-2) + 2 s r^(s-1) - l(l+1) r^(s-2) = 0
$
称为指标方程（characteristic equation），根为
$
s = l, -l-1
$
后者发散，从而而当$l = 0$时
$
R prop r^l
$
因此要求在求解径向方程时，解满足边界条件：
$
u(r) = r R(r) prop r^(l + 1)
$
在散射问题中，入射$lambda$相对于势场$U(r)$的作用范围很大时，$l=0$的贡献最大，只考虑$l=0$的S波（即正碰）就行了。

=== 氢原子和类氢离子的能级和波函数

氢原子或类氢离子的核电荷是$Z e$($Z$是原子序数)，核外有一个电子，所以势能是：
$
U(r) = - 1/(4 pi epsilon_0) (Z e^2) / r
$
约化质量是
$
mu = (m_e m_N) / (m_e + m_N)
$
分离变量后径向方程是：
$
- hbar^2 / (2 mu) dd(""^2u) /dd(r^2)  + ((l (l + 1) hbar^2) / (2 mu r^2) - k_1 (Z e^2) / r) u = E u\
dd(""^2 u)/dd(r^2) + ((2 mu) / hbar^2( E + k_1 (Z e^2) / r )- (l (l + 1)) /  (r^2)) u = 0
$
对于束缚态，$E<0$。定义一个无量纲的新变量
$
rho = alpha r , alpha = sqrt(8 mu abs(E)) / (hbar)
$
以及一个无量纲的新参数
$
beta = (2 mu k_1 Z e^2)/(alpha hbar^2) = (k_1 Z e^2) / hbar sqrt(mu/(2 abs(E)))
$
于是方程化为：
$
dd(""^2 u)/dd(rho^2) + (- 1/4 + beta/rho - (l (l + 1)) / rho^2 ) u = 0
$
讨论两种极限情形：
1. $rho -> oo$
   
   $
   dd(""^2 u)/dd(rho^2) - 1/4 u = 0
   $
    有解
    $
    u(rho) -> e^(-rho/2)
    $
    正指数舍去。
2. $rho -> 0$

    $
    dd(""^2 u)/dd(rho^2) - (l (l + 1)) / rho^2 u = 0
    $
    有解
    $
    u(rho) -> rho^(l+1)
    $

#figure(
  image("pic/2024-05-20-20-56-40.png", width: 80%),
  numbering: none
)
得到量子化的能级是：
$
E_n = (mu k_1^2 Z^2 e^4)/(2 hbar^2 n^2)
$
能量本征态由量子数组合$(n, l, m)$表征，它们的意义是：
- 主量子数$n = 1, 2, 3, ...$，能级数$E = E_n$
- 角量子数$l = 0, 1, 2, ..., n-1$，角动量量子数$L^2 = l (l + 1) hbar^2$
- 磁量子数$m = -l, -l+1, ..., l-1, l$，$L_z = m hbar$
能级$E_n$只和$n$有关，所以对$l$和$m$是简并的，简并度是：
$
g_n = sum_(l=0)^(n-1) (2l + 1) = n^2
$
对应的波函数是：
$
psi_(n l m) (r, theta, phi) = R_(n l) (r) Y_(l m) (theta, phi)\
R_(n l) (r) = (u_(n l) (r))/r
$
#figure(
  image("pic/2024-05-20-21-00-03.png", width: 80%),
  numbering: none
)
#figure(
  image("pic/2024-05-20-21-00-26.png", width: 80%),
    numbering: none
)
#figure(
  image("pic/2024-05-20-21-03-06.png", width: 80%),
  numbering: none
)
#figure(
  image("pic/2024-05-20-21-03-24.png", width: 80%),
  numbering: none
)

=== 碱金属（正1价金属）原子的能谱

氢原子的能级只和量子数$n$有关而和$l$无关，这是Coulomb势场所特有的结果。碱金属中的价电子是在原子实（即原子核加上内壳层电子）的作用下运动，它受到的势场就不再是Coulomb势场，所以碱金属原子的能级实际上与$n, l$都有关：
$
E = E_(n l)
$
举例来说，在钠（Na）原子中，$n=1,2$的轨道都被内壳层电子填满，所以价电子的最小主量子数是$n=3$。对于这些状态，能量的高低顺序是$E_(3S)<E_(3P)<E_(3D)$。只要原子没有受到外磁场的作用，它的能量总是和量子数$m$无关的。

== 径向位置概率分布

利用径向波函数，可以求得电子的径向位置分布概率，即不管方向如何，找到电子在$(r, r+dd(r) )$球壳中的概率为：
$
r^2 dd(r) integral sin theta dd(theta) dd(phi) |psi(r, theta, phi)|^2\
=  |R_(n l) (r)|^2 r^2 dd(r)\
= |u_(n l) (r)|^2 dd(r) 
$
在量子力学中，电子没有严格的轨道概念，只能研究其位置分布概率。$|u_(n l) (r)|^2$的极大值位置称为*最可几半径*。

== 概率密度随角度的变化

电子在$(theta , phi)$方向的立体角$dd(Omega) = sin theta dd(theta) dd(phi)$中的概率密度是：
$
|Y_(l m) (theta, phi)|^2 dd(Omega) prop | P_l^m (cos theta) |^2 dd(Omega)
$
和$phi$角无关。

== 粒子流密度和磁矩

设氢原子中电子处于$psi_(n l m)$态，则其粒子概率流密度为：
$
(i hbar)/(2 mu) (psi grad psi^* - psi^* grad psi)
$
乘以电子电荷（$-e$）得到电流密度：
$
arrow(j) = - e (i hbar)/(2 mu) (psi grad psi^* - psi^* grad psi)
$
而
$
psi_(n l m) tilde R_(n l) (r) P^m_l (cos theta) e^(i m phi)
$
其中$R,P$是实函数，从而$j_r = j_theta =0$。
$
j_phi &= (i e hbar)/(2 mu) 1/(r sin theta) (psi^* partial/(partial phi) psi - psi partial/(partial phi) psi^*)\
&= (i e hbar)/(2 mu) (2 i m)/(r sin theta) abs(psi)^2\
&= - (e hbar m)/(mu r sin theta) |R_(n l) (r) P^m_l (cos theta)|^2
$
电流密度$j_phi$对应的磁矩为：
$
arrow(M) = 1/2 integral arrow(r) crossproduct arrow(j) dd(V)
$
$
M_z = 1/2 integral r sin theta j_phi dd(V) = - (e hbar m)/(2 mu) |R_(n l) (r) P^m_l (cos theta)|^2 dd(r) dd(theta) dd(phi) = - (e hbar m)/(2 mu)
$
定义Bohr磁子：
$
mu_B = (e hbar) / (2 mu)
$
则氢原子磁矩为：
$
M_z = - mu_B m
$
所以也把$m$称为磁量子数。定义回转磁比值：
$
g = M_z / L_z = (- m mu_B)/(m hbar) = - e/(2 mu)
$

= 中心力场问题——三维各向同性谐振子

三维各向同性谐振子的势能函数
$
V(r) = 1/2 mu omega^2 r^2 = 1/2 mu omega^2 (x^2 + y^2 + z^2)
$
哈密顿量可写为：
$
H = sum_i H_i , H_i = - hbar^2/(2 mu) nabla_i^2 + 1/2 mu omega^2 r_i^2
$
类似多粒子系统，系统波函数分离变量：
$
psi(x, y, z) = psi_(n_x) (x) psi_(n_y) (y) psi_(n_z) (z)
$
其中$psi_n$为*一维谐振子*与量子数$n$的本征函数。

系统能级为：
$
E = (n_x + n_y + n_z + 3/2) hbar omega = (N + 3/2) hbar omega
$
其中$N = n_x + n_y + n_z$是量子数。

Virial定理仍然成立：
$
macron(T) = macron(V) = 1/2 macron(E)
$
对于给定的能级量子力学数$N$，其简并度为
$
((N + 1)(N + 2))/2
$

#newpara()

在球坐标系中，定态薛定鄂方程的径向部分：
$
(1/r^2 dd("")/dd(r) (r^2 dd("")/dd(r)) + (2 mu) / hbar^2 (E - 1/2 mu omega^2 r^2) -( l (l + 1)) / r^2) R_l (r) = 0
$
若令$R(r) = u(r)/r$，则有
$
(dd("")/dd(r) + (2 mu)/hbar^2 (E - 1/2 mu omega^2 r^2) -( l (l + 1) hbar )/(2 mu r^2)) u = 0 
$
在$l=0$时，可以看到$u$的方程和一维谐振子方程非常相像，但是又有本质不同：$u(r)$的自变量定义在$0≤r<∞$范围内，而一维谐振子范围是$-∞<x<∞$，这直接导致了它们的基态能量也不相同。

引入无量纲常量$ρ、λ$：
$
rho = alpha r, alpha = sqrt(mu omega / hbar)\
lambda = 2 E / (hbar omega)
$
则径向方程化为：
$
dd("")/dd(rho) R + 2/rho dd(R)/dd(rho) + (lambda - rho^2 - (l (l + 1)) / rho^2) R = 0
$
在$rho  -> 0$时，$R(rho) -> rho^l$；在$rho -> ∞$时，$R(rho) -> e^(-rho^2/2)$。得到渐近形式后可设：
$
R(rho) = e^(-rho^2/2) rho^l v(rho)
$
代入原方程后，再做变量代换$ξ = ρ^2$得到*合流超几何方程*：
$
dd(""^2v)/dd(ξ^2)  + ((2l + 3)/(2 ξ) - 1) dd(v)/dd(ξ) + (lambda - 2l -3)/(4 eta) v = 0
$
和氢原子的合流超几何方程做对比
$
dd(""^2v)/dd(rho^2)  + ((2(l + 1))/(rho) - 1) dd(v)/dd(ξ) + (beta - l -1)/(rho) v = 0
$
类似于氢原子方程有解条件
$
n_r = beta - l - 1 = 0, 1, 2, ...
$
当前方程的有解条件为
$
n_r = (lambda - 2l - 3)/4 = 0, 1, 2, ...
$
解得
$
lambda = 2N + 3 (N = 2n_r + l)
$
代入λ表达式
$
E_N = (N + 3/2) hbar omega
$
在$N$给定以后，$l$可以取值
$
l = N , N - 2, ..., 0"或"1
$
最后系统径向波函数为
$
R(r) = C L_(n_r)^(l+1/2) (rho^2) rho^l e^(-rho^2/2) (rho = sqrt((mu omega) / hbar) r)
$
其中$L_(n_r)^(l+1/2)$是缔合Laguerre多项式。

#figure(
  image("pic/2024-05-24-00-57-38.png", width: 80%),
    numbering: none
)

实质上说，对于同一个$N$，直角坐标系中的的波函数$psi_(n_x) (x) psi_(n_y) (y) psi_(n_z) (z)$和球坐标系中的波函数$R(r) Y_(l m) (theta, phi)$（不同对易力学量完全集的基底）可以通过幺正变换互相联系，这是*表象变换*的一个实际例子。

#figure(
  image("pic/2024-05-24-00-59-09.png", width: 80%),
    numbering: none
)

#figure(
  image("pic/2024-05-24-01-08-16.png", width: 80%),
    numbering: none
)

#figure(
  image("pic/2024-05-24-01-14-46.png", width: 80%),
    numbering: none
)

#figure(
  image("pic/2024-05-24-01-16-13.png", width: 80%),
    numbering: none
)

下面求$hat(L)_z$在$(hat(H), hat(L)^2, hat(L)_z)$和$(hat(H)_x, hat(H)_y, hat(H)_z)$表象下的矩阵表示：
#figure(
  image("pic/2024-05-24-01-22-14.png", width: 80%),
    numbering: none
)

#figure(
  image("pic/2024-05-24-01-25-11.png", width: 80%),
    numbering: none
)

#figure(
  image("pic/2024-05-24-01-25-35.png", width: 80%),
    numbering: none
)
#figure(
  image("pic/2024-05-24-01-29-04.png", width: 80%),
    numbering: none
)

= 带电粒子在电磁场中的运动

== 薛定鄂方程的幺正变换

一般量子力学问题的薛定鄂方程：
$
(- hbar^2/(2 mu) nabla^2 + U(r)) psi(arrow(r)) = i hbar partial/(partial t) psi(arrow(r))
$
对于库仑势能来说：
$
U(r) = q Phi(arrow(r))
$
其中$Phi(arrow(r))$为外界电场力所对应的静电势（注意不一定是中心势，如恒定均匀外场势），q为粒子电荷。对氢原子（中心势场）中的电子来说
$
q = - e\
Phi(arrow(r)) = k_1 e / r
$
加入磁场后，经典哈密顿量变为（参考分析力学）:
$
H = 1/(2 m) (arrow(p) - q arrow(A))^2 + q Phi(arrow(r))
$
磁力不是保守力，不像库仑力那样有一个标量势能项。但是我们知道*电磁场是包含动量的*，电荷$q$产生的$arrow(E)$与外磁场$arrow(B)$结合产生动量密度$epsilon_0 arrow(E) crossproduct arrow(B)$，这反映在动量的改变量中。$arrow(p)$是*正则动量*，而$arrow(pi) = arrow(p) - q arrow(A)$是*机械动量*。

相应的，带电粒子在外电磁场作用下的哈密顿算符：
$
hat(H) = 1/(2 mu) (hat(p) - q hat(A))^2 + q Phi(arrow(r))
$
同氢原子问题（只有$q Phi$项）一样，这个算符对应的薛定鄂方程的适用范围是低速运动的粒子。对于高能问题，需对波函数和电磁场进行量子化（所谓二次量子化）。

对任意势场，对方程做幺正变换：
$
psi -> psi' = e^(i theta) psi\
hat(H) -> hat(H)' = e^(i theta) hat(H) e^(- i theta)
$
其中$theta = theta(arrow(r))$为不显含时间的任意实函数，显然这一变换是幺正变换，也不改变薛定鄂方程：
$
hat(H)' psi' &= e^(i theta) hat(H) e^(- i theta) e^(i theta) psi = e^(i theta) hat(H) psi\
&= i hbar partial/(partial t) psi'\
$
现在
$
hat(H) = 1/(2 mu) (- i hbar grad - q arrow(A))^2 + q Phi(arrow(r))
$
幺正变换的结果是：
#figure(
  image("pic/2024-05-24-01-48-06.png", width: 80%),
    numbering: none
)

#figure(
  image("pic/2024-05-24-01-49-56.png", width: 80%),
    numbering: none
)

$
arrow(A) -> arrow(A)' = arrow(A) + hbar/q grad theta\
$
根据幺正变换的性质，在量子力学中，这一代换不会引起任何物理上的变化。

在经典力学中：
$
arrow(B) = nabla crossproduct arrow(A) = nabla crossproduct arrow(A)' = nabla crossproduct (arrow(A) + hbar/q grad theta)
$
也就是说，这一代换在经典电磁学中同样不会产生任何物理上的不同。

以上考虑的是静电场和静磁场的情况，在变化的电磁场中，哈密顿量显含时间，相应的幺正变换则为
$
psi -> psi' = e^(i theta) psi\
hat(A) -> hat(A)' = hat(A) + hbar/q grad theta\
Phi -> Phi' = Phi - hbar/q (partial theta) / (partial t)
$
这一套变换又称为规范（gauge）变换。规范变换不改变系统物理学性质——系统具有规范不变性。

== 规范不变性与Yang-Mills理论

在量子场论中，初等量子力学中的波函数演变为经典场进入到哈密顿量中，与经典电磁场一起进行二次量子化。

在这种情况下，规范变换就不像初等量子力学那样对波函数和算符同时变换，而是仅对哈密顿量（或拉格郎日量）进行变换-这就体现为系统的一种对称不变性。

根据Noether定理，每一种对称性的背后都有一个守恒量。在量子场论中，如果系统规范不变，将带来深刻的物理结果，比如：
- 系统电荷守恒（诸如$e -> nu gamma$不可能发生）
- 光子质量为0
- 光子自旋投影只能是$±hbar$，没有0分量

#figure(
  image("pic/2024-05-24-01-56-05.png", width: 80%),
    numbering: none
)

规范不变对系统的拉格郎日量（或哈密度量）的形式做出了强烈的限制，基于规范不变思想的量子场论最早由Pauli提出（1953），可惜没有正式发表。

Yang、Mills二人发展了这一思想，把规范不变的公设从电磁U(1)相互作用延伸到了SU(N)相互作用，对SU(N)相互作用的形式作出了限定。美中不足的是，规范不变性成立的前提是所有基本粒子质量为0，这显然与实验不符，这也是Pauli没有发表他早期研究的原因。

后来希格斯等提出了基于自发性对称破缺的机制来解释为什么粒子质量不为0。2012年希格斯粒子的发现，使得基于规范不变和希格斯质量机制这两大支柱的粒子物理“标准模型”得以最终确立。

= 塞曼效应和郎道能级

== 在外场中的原子

带电粒子在外场中的定态薛定鄂方程：
$
(1/ (2m) (- i hbar grad - q arrow(A))^2 + q Phi(arrow(r))) psi(arrow(r)) =E psi(arrow(r))\
1/(2m) (-hbar^2 nabla^2 psi + i hbar q (grad dot (arrow(A) psi) + arrow(A) dot grad psi) + q^2 arrow(A)^2 psi) = (E - q Phi) psi\
1/(2m) (-hbar^2 nabla^2 psi + i hbar q ((grad dot arrow(A)) psi +2 arrow(A) dot grad psi) + q^2 arrow(A)^2 psi) = (E - q Phi) psi\
$
取*库仑规范*：
$
div arrow(A) = 0
$
则有
$
1/(2m) (-hbar^2 nabla^2 psi + 2i hbar q (arrow(A) dot grad psi) + q^2 arrow(A)^2 psi) = (E - q Phi) psi
$
利用Stokes定理：
$
integral.double_(S) (curl arrow(A)) dd(arrow(S)) = integral.cont arrow(A) dot dd(arrow(l))
$
对恒定静磁场$arrow(B)$来说：
$
B pi r^2 = A 2 pi r\
arrow(A) = 1/2 arrow(B) crossproduct arrow(r)
$
*库伦规范*就是：
$
div arrow(A) = 1/2 div (arrow(B) crossproduct arrow(r)) = 0
$
代入原式
$
1/(2m) (-hbar^2 nabla^2 psi + i hbar q  (arrow(B) crossproduct arrow(r)) dot grad psi + 1/4q^2  (arrow(B) crossproduct arrow(r))^2 psi) = (E - q Phi) psi\
1/(2m) (-hbar^2 nabla^2 psi + i hbar q  arrow(B) dot (arrow(r) crossproduct grad psi) + 1/4q^2  (r^2B^2 - (arrow(r)dot arrow(B))^2) psi) = (E - q Phi) psi
$

#figure(
  image("pic/2024-05-24-02-08-07.png", width: 80%),
    numbering: none
)

#figure(
  image("pic/2024-05-24-02-08-25.png", width: 80%),
    numbering: none
)

所以相对第二项（的变化量）来说，第三项可以忽略不计。于是：
$
(- hbar^2/(2m) nabla^2 - q arrow(B) dot hat(arrow(L)) + q Phi) psi = E psi
$
如果选择$z$轴的方向为$arrow(B)$的方向，则
$
(- hbar^2/(2m) nabla^2 - (q B)/(2 m) hat(L)_z + q Phi) psi = E psi
$
电子电荷$q = - e$，所以
$
(- hbar^2/(2 mu) nabla^2 + (e B)/(2 mu) hat(L)_z - e Phi) psi = E psi
$
设已求得未加磁场$(B=0)$时碱金属原子的能级与波函数
$
E_(n l) , psi_(n l) (r, theta, phi) = R_(n l) (r) Y_(l m) (theta, phi)
$
每一个能级是$(2l+1)$度简并的。那么加上外磁场后(相当于$hat(H)' = hat(H) + (e B)/(2 mu)hat(L)_z$)，本征波函数不变，本征值发生改变，简并将被打破。
$
E_(n l m) = E_(n l) + (e B)/(2 mu)hbar m 
$
其中$m$为磁量子数。而波函数的形式仍旧不变

#figure(
  image("pic/2024-05-24-13-41-40.png", width: 80%),
  caption: [
    塞曼效应
  ],
)

*碱金属原子的能级在强磁场中分裂的现象称为正常塞曼(Zeeman)效应。*

== 自由粒子在磁场中运动

在*对称规范*中，我们约定均匀磁场$arrow(B)$的矢势为
$
arrow(A) = 1/2 arrow(B) crossproduct arrow(r)
$
如果取$B$沿z轴方向，则：
$
arrow(A) = vec(-1/2 B y, 1/2 B x, 0)
$
进行规范变换(规范变换不改变$arrow(B)$)，新的磁矢势仍满足库伦规范
$
arrow(A) -> arrow(A)' = arrow(A) + grad f , f = 1/2 B x y
$
这时磁矢势变为（此即*朗道规范*）
$
arrow(A)' = vec(0, B x, 0)
$
取电子电荷$-e$，则在均匀磁场中运动电子的定态薛定鄂方程为
$
1/(2m) (hat(arrow(p)) + e arrow(A))^2 psi = E psi
$
设$arrow(B)$沿$z$轴方向，电子运动限制在$x-y$平面内（二维电子气模型），则在朗道规范下方程为
$
1/(2 m) (hat(p)_x^2 + (hat(p)_y + e B x)^2) psi = E psi
$
易见
$
[hat(p)_y, hat(H)] = 0
$
分离变量求特解：
$
psi(x, y) = e^(i k_y y) phi(x)
$
代入原方程得：
$
(- hbar^2/(2 m) dd(""^2)/dd(x^2) + 1/2 m omega_c^2 (x + x_0)^2 ) phi(x) = E phi(x)
$
其中
$
omega_c = (e B) / m, x_0 = k_y l_c^2 , l_c = sqrt(hbar /(m omega_c)) = 1/alpha , alpha^2 = (e B) / hbar
$
$ω_c$是回旋角频率，$l_c$是最小回旋半径
$
(m v)/R = e B\
T = (2 pi R)/v = (2 pi m) / (e B)\
omega_c = (2 pi) / T =  (e B)/m
$
$
2 pi l_c = lambda = h/p\
e B = p/l_c\
=> l_c = sqrt(hbar / (e B)) = sqrt(hbar / (m omega_c))
$
这个方程的解即是一维谐振子方程的解，只是坐标平移了$x_0$:
$
phi(x) = phi_n (x+x_0), psi(x, y) = e^(i k_y y) phi_n (x+x_0)\
E_n = (n + 1/2) hbar omega_c "朗道能级"
$
量子观点：粒子在$x-y$平面内绕$z$轴转动。*粒子能量就是这种转动产生的磁矩*与磁场的相互作用能
$
E = (n + 1/2) hbar omega_c = (n + 1/2) (e hbar) / m B = -mu_z B
$
则
$
mu_z = - (e hbar) / (2 m) = - mu_B
$
即磁矩方向与磁场方向相反——*朗道抗磁性*。朗道抗磁性与电荷正负无关，是自由粒子在磁场中运动的量子效应。

朗道能级$n$对应的波函数$e^(i k_y y) phi_n (x+x_0)$是一种平面波和谐振子波函数的乘积，简并度是无穷大的：对于每个能级$E_n$，对应波函数中的$k_y$可以任意取值。

考虑电子气局限于$L_x$宽的长条中，则必须有
$
0 < x_0 < L_x => 0 < k_y < L_x alpha^2
$
考虑$y$轴方向周期性边界条件：长条内每$L_y$长度内有一个电子（即一维箱归一化），得
$
k_y = (2 pi N)/L_y , N = 0, 1, 2, ...
$
所以
$
0 < N < (L_x L_y alpha^2)/(2 pi) = (e B L_x L_y) / (h)
$
于是*单位面积内的能级简并度为*
$
g = (e B) / (h)
$
这是一个重要的结果，对于理解量子霍尔效应很有用。

如果使用对称规范，则电子绕$z$轴转动的物理图像更加一目了然。但是物理结论不依赖于规范选择（如同三维谐振子在直角坐标和球坐标系表象中的解一样），这个问题中两个不同规范对应的波函数解可以通过幺正变换联系起来。

== 量子霍尔效应

= 电子自旋及其描述

== 角动量

轨道角动量算符$hat(arrow(L))$的各个分量满足对易关系
$
[hat(L)_i, hat(L)_j] = i hbar epsilon_(i j k) hat(L)_k
$
其中$epsilon_(i j k)$是三维Levi-Civita符号，即三阶完全反对称张量。只有
角标$i,j,k$各不相同时，$epsilon_(i j k)$才不为0，否则为0。约定$epsilon_(1 2 3) = 1$，再任意排列的情况下，$epsilon_(i j k)$的值由排列的奇偶性决定。

这些对易关系是*角动量算符的定义以及量子力学基本对易关系*（$[hat(x) , hat(p)] = i hbar$）所导致的结果。

可以把它们推广为量子力学中的一般角动量应该满足的对易关系，也就是说，我们假设若 $hat(arrow(J))$是一个角动量算符，那么它的各个分量算符要满足
$
[hat(J)_i, hat(J)_j] = i hbar epsilon_(i j k) hat(J)_k
$
在量子力学里，上式可以看作是*角动量算符的一般定义*。
$
[hat(J)^2 , hat(J)_i] = 0, hat(J)^2 = hat(J)_x^2 + hat(J)_y^2 + hat(J)_z^2
$
*角动量本征态*是$hat(J)^2$和$hat(J)_z$的共同本征态，本征值分别是$eta hbar^2$和$m hbar$，其中$j$是角动量量子数，$m$是磁量子数。
$
hat(J)^2 ket(eta"," m) = eta hbar^2 ket(eta","  m)\
hat(J)_z ket(eta","  m) = m hbar ket(eta","  m)
$
注意，在Dirac符号的形式下，我们只是说存在满足角动量对易关系的力学量算符和它们的本征态，但是并不需要把它们写成任何具体的函数形式。

== 阶梯算符

引进*阶梯算符*：
$
hat(J)_(plus.minus) = hat(J)_x ± i hat(J)_y
$
不难证明
$
[hat(J)_z , hat(J)_(plus.minus)] = [hat(J)_z, hat(J)_x] ± i [hat(J)_z, hat(J)_y] = ± hbar hat(J)_(plus.minus)
$
从而
$
hat(J)_z hat(J)_(plus.minus) = hat(J)_(plus.minus) hat(J)_z ± hbar hat(J)_(plus.minus)\
hat(J)_z hat(J)_(plus.minus) ket(eta"," m) =  hat(J)_(plus.minus) hat(J)_z ket(eta"," m) ± hbar hat(J)_(plus.minus) ket(eta"," m)\
hat(J)_z hat(J)_(plus.minus) ket(eta"," m) = (m ± 1) hbar hat(J)_(plus.minus) ket(eta"," m)
$
从而
$
hat(J)_(plus.minus) ket(eta"," m) prop ket(eta'"," m ± 1)
$
设
$
hat(J)_(plus.minus) ket(eta"," m) = c ket(eta'"," m ± 1)
$
再利用
$
[hat(J)^2 , hat(J)_(plus.minus)] = [hat(J)^2 , hat(J)_x] ± i [hat(J)^2 , hat(J)_y] = 0
$
令$hat(J)^2$作用于等式的左端：
$
hat(J)^2 hat(J)_(plus.minus) ket(eta"," m) = hat(J)_(plus.minus) hat(J)^2 ket(eta"," m) = eta hbar^2 hat(J)_(plus.minus) ket(eta"," m)
$
作用于右端：
$
hat(J)^2 c ket(eta'"," m ± 1) = c eta' hbar^2 ket(eta'"," m ± 1)
$
就有
$
eta = eta', hat(J)_(plus.minus) ket(eta"," m) = c ket(eta"," m ± 1)
$
这就是*阶梯算符*的含义（只改变m，不改变$eta$）。

可以证明，若$m$的极大值为$j$，则$eta=j(j+1)$。于是我们用$ket(j "," m)$来表示一个角动量本征态，而不再用$ket(eta "," m)$。有
$
hat(J)_(plus.minus) ket(j "," m) = hbar sqrt(j(j+1) - m(m ± 1)) ket(j "," m ± 1)
$
同时可证，$j$的可能取值为非负半整数或整数，即
$
j = 1/2 , 3/2 , 5/2 ... 或 0 , 1 , 2 ... 
$
这是从角动量算符对易关系得出的一般结果，与中心力场问题求解过程中得到的$(hat(L)^2, hat(J)_z)$的本征值和本征函数（即球谐函数）间的关系是完全一致的。但是，这里用的代数解法更加普适，把$j$为半整数的情形也推导出来了，这恰恰就是电子自旋的情况，球谐函数是无能为力的。


#figure(
  image("pic/2024-05-27-00-48-10.png", width: 80%),
    numbering: none
)

#figure(
  image("pic/2024-05-27-19-08-50.png", width: 80%),
numbering: none
)

#figure(
  image("pic/2024-05-27-19-10-08.png", width: 80%),
numbering: none
)

非相对论量子力学在解释许多实验现象上都获得了成功，例如氢原子的能谱结构，但是更进一步的实验发现，还有许多实验现象，例如光谱线在磁场下的分裂、光谱线的精细结构，用前面讲述的理论无法解释，原因在于，以前的理论只涉及到轨道角动量。而新的实验表明，电子还具有自旋角动量。

在非相对论量子力学中，自旋是作为一个新的附加的量子数引入的，只是在薛定鄂方程中加入自旋。

在相对论量子力学中，电子的自旋将自然地包含在相对论的波动方程Dirac方程中。

== 电子自旋的发现

碱金属元素有特征光谱的双线结构。

*Stern-Gerlach*实验（1922）：测量银原子的磁矩。

#figure(
  image("pic/2024-05-27-19-26-50.png", width: 60%),
  caption: [
    Stern-Gerlach实验
  ],
)

让银原子通过不均匀的磁场
$
arrow(B) = B(z) arrow(e)_z
$
根据银原子的磁矩$arrow(mu)$在这个磁场中的势能：
$
U = - arrow(mu) dot arrow(B) = - mu_z B(z)
$
受力
$
arrow(F) = - grad U = - mu_z (partial B(z))/(partial z) arrow(e)_z
$
这个力使它的飞行轨迹发生偏转，*偏转大小和原子磁矩在磁场方向上的投影$mu_z$成正比*。实验结果是
$
mu_z = ± mu_B, mu_B = (e hbar) / (2 m) "Bohr磁子"
$
电子有磁矩，其投影是量子化的。推论：*电子有自旋（内禀角动量），其投影也是量子化的*。

Uhlenbeck-Goudsmit假设(1925)：*电子有自旋角动量*，其投影只能取两个值：
$
S_z = ± hbar/2
$
不能将电子的自旋按照经典的图像看作是带电的小球绕自身轴的自转。要达到实验观测到的磁矩，小球表面的线速度将超过光速。

== 电子自旋的描述

自旋这个新的自由度的特点：
- 是个内禀的物理量，不能用坐标、动量、时间等变量表示
- 完全是一种量子效应，没有经典的对应：$hbar -> 0$时，自旋角动量消失
- 是角动量，满足角动量算符的最一般的对易关系
- 电子自旋在空间任何方向上的投影只取
  $
  plus.minus hbar/2
  $
  $hat(S)_x, hat(S)_y, hat(S)_z$的本征值是$± hbar/2$，
  $
  hat(S)_x^2 = hat(S)_y^2 = hat(S)_z^2 = (hbar/2)^2\
  hat(S)^2 = hat(S)_x^2 + hat(S)_y^2 + hat(S)_z^2 = (3/4) hbar^2
  $
- 一般总自旋为$j$的粒子，自旋算符要用$(2j + 1)times(2j + 1)$维矩阵表示。对于电子，$j = 1/2$，自旋算符是$2times 2$矩阵。

自旋角动量又导致电子有*自旋磁矩*，其$z$轴投影为
$
mu_z / S_z = - e/m_e , mu_z = minus.plus mu_B
$
*自旋角动量算符*记为$hat(arrow(S))$，*自旋磁矩算符*记为$hat(arrow(mu))$，它们之间的关系是
$
hat(arrow(mu)) = - e /m_e hat(arrow(S))
$
$hat(S)_x, hat(S)_y, hat(S)_z$都是2×2的矩阵，通常选择$hat(S)_z$为对角阵，即在$(hat(S)^2, hat(S)_z)$表象下，$hat(S)_z$的矩阵形式为
$
hat(S)_z = hbar/2 mat(
  1,0;
  0,-1
)
$
根据递推关系
$
hat(S)_plus.minus ket(m) = hbar sqrt(3/4 - m(m ± 1)) ket(m ± 1)
$
可以求其他算符的矩阵形式
$
hat(S)_plus mat(0;1) = hbar mat(1;0), hat(S)_plus mat(1;0) = 0 => hat(S)_plus = hbar mat(0,1;0,0)
$
同理有：
$
hat(S)_minus = hbar mat(0,0;1,0)
$
从而$hat(S)_x, hat(S)_y$的矩阵形式为
$
hat(S)_x = hbar/2 mat(
  0,1;
  1,0
), hat(S)_y = hbar/2 mat(
  0,-i;
  i,0
)
$
我们记Pauli矩阵
$
sigma_x = mat(
  0,1;
  1,0
), sigma_y = mat(
  0,-i;
  i,0
), sigma_z = mat(
  1,0;
  0,-1
)
$
则
$
hat(S)_x = hbar/2 sigma_x, hat(S)_y = hbar/2 sigma_y, hat(S)_z = hbar/2 sigma_z\
hat(arrow(S)) = hbar/2 hat(arrow(sigma))
$
这就是电子自旋的描述。

$hat(S)_z$对应的本征值$± hbar/2$，对应的本征态是
$
ket(hbar/2) = mat(1;0), ket(-hbar/2) = mat(0;1)
$
电子的任何自旋态$v$可以表示成
$
ket(v) = c_1 ket(hbar/2) + c_2 ket(-hbar/2) = mat(c_1;c_2)
$
其中$c_1$是某自旋态下$S_z$取$+hbar/2$的概率，$c_2$是某自旋态下$S_z$取$-hbar/2$的概率。自旋波函数归一化的条件：
$
v^dagger v = 1 => |c_1|^2 + |c_2|^2 = 1
$

== Pauli矩阵的主要性质

Pauli矩阵的主要性质：
- Pauli矩阵是厄密矩阵
- $i != j $时满足：
  $
  sigma_i sigma_j = -   sigma_j sigma_i = i epsilon_(i j k) sigma_k
  $
  Pauli矩阵是彼此反对易的，满足对易关系
  $
  [sigma_i, sigma_j] = i epsilon_(i j k) sigma_k
  $
- 满足
  $
  sigma_i^2 = I
  $
  其中$I$是单位矩阵，Pauli矩阵是*幺正矩阵*。

  上面两个式子可以合并为
  $
  sigma_i sigma_j = delta_(i j) I + i epsilon_(i j k) sigma_k
  $
  具有这类关系的矩阵称为*Clifford代数*。

- 满足：
  $
  (arrow(a) dot arrow(sigma)) (arrow(b) dot arrow(sigma)) = arrow(a) dot arrow(b) I + i arrow(a) crossproduct arrow(b) dot arrow(sigma)
  $

== 带有自旋的电子波函数和算符

现在电子的波函数应该同时描写它的自旋状态。由叠加原理
$
Psi(arrow(r), t) = Psi_1 (arrow(r), t) v_+ + Psi_2 (arrow(r), t) v_- = mat(Psi_1 (arrow(r), t); Psi_2 (arrow(r), t))
$
这称为电子的*二分量波函数*，又称为*二分量旋量*。其中
$
Psi_1 (arrow(r), t) = Psi(arrow(r) , t, S_z = hbar/2), Psi_2 (arrow(r), t) = Psi(arrow(r) , t, S_z = -hbar/2)
$
如果$Ψ_1$和$Ψ_2$不具有固定比例（即轨道和自旋波函数不能分离变量），则粒子处于轨道和自旋的*耦合态*。

粒子轨道和自旋无耦合时，粒子的总波函数可以写为空间部分$psi(x)$和自旋部分$ket(chi)$的直积态：
$
Psi(arrow(r), t) = psi(arrow(r), t) ket(chi)
$
设$ket(chi) = c_1 ket(arrow.t) + c_2 ket(arrow.b) = mat(c_1; c_2)$，其中$ket(arrow.t)$和$ket(arrow.b)$是自旋向上和向下的本征态，$c_1$和$c_2$是自旋向上和向下的概率振幅。则在空间任意位置$x$找到粒子自旋向上和向下概率之比不依赖于$x$：
$
abs(psi braket(arrow.t , chi))^2/abs(psi braket(arrow.b , chi))^2 = abs(c_1)^2/abs(c_2)^2
$
相反的，如果这个比值依赖于$x$，则说明轨道和自旋耦合：
$
psi_1 ket(arrow.t) + psi_2 ket(arrow.b) = psi_1 mat(1;0) + psi_2 mat(0;1) = mat(psi_1; psi_2)
$
对于这样的波函数和算符，原先的公式需要稍加修正
- 波函数的归一化是：
  $
  integral Psi^dagger Psi dd(arrow(r)) = integral (abs(Psi_1)^2 + abs(Psi_2)^2) dd(arrow(r)) = 1
  $
- 电子的空间几率密度是：
  $
  W(arrow(r)) = Psi^dagger(arrow(r)) Psi(arrow(r)) = abs(Psi_1)^2 + abs(Psi_2)^2
  $
- 电子的两种自旋状态的几率是：
  $
  W_arrow.t (arrow(r)) = abs(Psi_1)^2, W_arrow.b (arrow(r)) = abs(Psi_2)^2
  $
  如果自旋和轨道*非耦合*（即没有自旋-轨道相互作用）的状态，此时$Psi_1$和$Psi_2$函数形式呈固定比例：
  $
  Psi(arrow(r), t) = Psi_0 (arrow(r), t) mat(a; b)
  $
  其中$Psi_0$是轨道波函数，$a$和$b$是常数，$|a|^2 + |b|^2 = 1$。自然界的电子当然是带有自旋的，但是我们以前不考虑电子自旋也做过许多计算，在实质上，那等于是假设了电子是处在上述的自旋和轨道非耦合的状态下，所以电子的自旋自由度不带来可观察到的影响。
- 算符的平均值是：
  $
  macron(A) = integral Psi^dagger hat(A) Psi dd(arrow(r)) = integral (mat(Psi_1^*, Psi_2^*) hat(A) mat(Psi_1; Psi_2)) dd(arrow(r))
  $
  $Psi^dagger hat(A) Psi$在一般情况下既包括坐标函数的运算又包括矩阵运算。

  如果算符$hat(A)$和自旋无关（$hat(A)$不作用在自旋态$ket(chi)$上，如动量算符等），则$hat(A)$在自旋表象中的矩阵表示是对角化的：
  $
  A_(1 1) = braket(arrow.t , hat(A) , arrow.t) = hat(A) braket(arrow.t) = hat(A)\
  hat(A) = mat(
    hat(A), 0;
    0, hat(A)
  )
  $
== 静磁场中电子的自旋

自旋角动量与磁场耦合能对应的算符：
$
- hat(arrow(mu)) dot arrow(B) = - g hat(arrow(S)) dot arrow(B)
$
如果不考虑粒子空间运动，只考虑内禀自旋运动，那么：
$
hat(H) = - g hat(arrow(S)) dot arrow(B) = - (g hbar B)/2 arrow(sigma) dot arrow(e)_B
$
$arrow(sigma) dot arrow(e)_B$的本征值是$±1$。

$
omega_L = - g B
$
为*拉莫频率*。

== 自旋量子态的时间演化与量子跃迁

设粒子在$t=0$时处于自旋量子态$ket(chi)$，则其在后续任意时刻$t$的自旋波函数可表示为
$
ket(chi(t)) &= e^(- i/hbar t hat(H)) ket(chi(0))\
&= e^(- i (w_L t)/2 arrow(sigma) dot arrow(e)_B) mat(a_0; b_0) sigma_z"表象"\
&= (cos(w_L/2 t) - i sin(w_L/2 t) arrow(sigma) dot arrow(e)_B ) mat(a_0; b_0)\
$
如果时间演化算符具有非0非对角矩阵元，则有可能出现自旋向上和向下的部分相互“*跃迁*”。

例：取$arrow(B)$沿$x$轴方向，$ket(chi(0))$为自旋向上的$sigma_z$本征态，则
$
ket(chi(t)) &= (cos(w_L/2 t) - i sin(w_L/2 t) sigma_x) mat(1;0)\
&= mat(
  cos(w_L/2 t), - i sin(w_L/2 t);
  - i sin(w_L/2 t), cos(w_L/2 t)
  ) mat(1;0)\
&= mat(
  cos(w_L/2 t);
  - i sin(w_L/2 t)
  )
$
这种系统周期性地在两种不同量子态间来回跃迁又称为*振荡*(oscillation)。粒子自旋出现振荡现象的原因是：$sigma_z$和哈密顿算符不对易，自旋态$ket(chi)$不是定态。

相反，如果$arrow(e)_B$沿$z$轴方向，则$hat(H) = 1/2 hbar omega_L sigma_z$和$sigma_z$对易，那么就不存在两个自旋量子态之间的跃迁，这时$ket(chi)$可为任意自旋态
$
ket(chi(t)) = ( cos(omega_L t/2) - i sin(omega_L t/2) sigma_z ) ket(chi(0)) = mat(e^(- i omega_L t/2) a_0 ; e^(i omega_L t/2) b_0)
$
这时候自旋态的概率就不发生震荡了。

= 角动量的合成、角动量耦合表象、反常塞曼效应(Bell基)

== 角动量的合成

=== 角动量合成的一般规则

设$hat(arrow(J))_1$和$hat(arrow(J))_2$是两个互相独立的角动量，它们的分量分别满足角动量算符的对易关系（$[hat(J)_i, hat(J)_j] = i hbar epsilon_(i j k) hat(J)_k$），而它们互相之间是对易的：$hat(arrow(J))_1$和$hat(arrow(J))_2$是对易的，即$[hat(J)_(1 i), hat(J)_(2 j)] = 0$。

矢量和$arrow(J) = arrow(J)_1 + arrow(J)_2$也就是各个分量对应相加：
$
hat(J)_i = hat(J)_(1 i) + hat(J)_(2 i)
$
得到的$hat(arrow(J))$也是角动量算符，满足角动量算符的对易关系。

从*未耦合（也就是和还未相加）*的角度看来，这个体系的对易可观测量完全集CSCO是
$
hat(J)_1^2, hat(J)_2^2, hat(J)_1^2, hat(J)_2^2
$
共同本征态是
$
ket(j_1","m_1","j_2","m_2) = ket(j_1","m_1) ket(j_2","m_2)
$
是态的直积，并矢。Hilbert空间的总维数是$(2j_1 + 1)(2j_2 + 1)$。

而*角动量耦合*以后，体系的CSCO变成了
$
hat(J)^2, hat(J)_z, hat(J)_1^2, hat(J)_2^2
$
它们是两两对易的，把耦合以后的共同本征态记做
$
ket(j","m";"j_1","j_2)
$
根据态叠加原理和本征函数完备性：
$
ket(j","m";"j_1","j_2) = sum_(m_1,m_2) C(j ,m; j_1, m_1, j_2, m_2) ket(j_1","m_1) ket(j_2","m_2)
$
用$hat(J)_z = hat(J)_(1z) + hat(J)_(2z)$的本征态展开：
$
m ket(j m j_1 j_2) = sum_(m_1,m_2) (m_1 + m_2) C(j ,m; j_1, m_1, j_2, m_2) m_1 ket(j_1 m_1) ket(j_2 m_2)
$
得到
$
sum_(m_1,m_2) (m - m_1 -m_2) C(j ,m; j_1, m_1, j_2, m_2) ket(j_1 m_1) ket(j_2 m_2) = 0
$
其中$ket(j_1 m_1) ket(j_2 m_2)$两两正交，所以有
$
(m - m_1 - m_2) C(j ,m; j_1, m_1, j_2, m_2) = 0
$
只有在
$
m = m_1 + m_2
$
才可能有
$
C(j ,m; j_1, m_1, j_2, m_2) != 0
$
于是
$
ket(j","m";"j_1","j_2) = sum_(m = m_1 + m_2) C(j ,m; j_1, m_1, j_2, m_2) ket(j_1","m_1) ket(j_2","m_2)
$
#newpara()

在一个特殊的状态下，直积的本征态和耦合的本征态是相同的，那就是“最大投影态”$ket(j_1 j_1) ket(j_2 j_2)$。$m_1 = j_1, m_2 = j_2$，$m = j_1 + j_2$为其最大值。$m$的最大值显然也是$j$的最大值。注意到：
$
hat(J)^2 = (hat(arrow(J))_1 + hat(arrow(J))_2)^2 = hat(arrow(J))_1^2 + hat(arrow(J))_2^2 + 2 hat(arrow(J))_1 dot hat(arrow(J))_2 = hat(arrow(J))_1^2 + hat(arrow(J))_2^2 + hat(J)_(1 +) hat(J)_(2 -) + hat(J)_(1 -) hat(J)_(2 +) + 2 hat(J)_(1 z) hat(J)_(2 z)
$
作用到最大投影态上：
$
hat(J)^2 ket(j_1 j_1) ket(j_2 j_2) &= hat(arrow(J))_1^2 + hat(arrow(J))_2^2 + hat(J)_(1 +) hat(J)_(2 -) + hat(J)_(1 -) hat(J)_(2 +) + 2 hat(J)_(1 z) hat(J)_(2 z) ket(j_1 j_1) ket(j_2 j_2)\
&= (j_1(j_1 + 1) + j_2(j_2 + 1) + 2 j_1 j_2) ket(j_1 j_1) ket(j_2 j_2)\
&= (j_1 + j_2)(j_1 + j_2 + 1) ket(j_1 j_1) ket(j_2 j_2)
$
即
$
ket(j_1 j_1) ket(j_2 j_2) = ket(j = j_1 + j_2","m = j_1 + j_2","j_1","j_2)
$
对于态$j = j_1 +j_2$，它是$j$的最大可能值
$
j_max = j_1 + j_2
$
让$m_1$或$m_2$减小1，未耦合的本征态就是$ket(j_1 j_1 -1)ket(j_2 j_2)$或$ket(j_1 j_1)ket(j_2 j_2 -1)$，耦合以后的本征态就是$ket(j = j_1 + j_2","m = j_1 + j_2 - 1","j_1","j_2)$或$ket(j = j_1 + j_2-1","m = j_1 + j_2 - 1","j_1","j_2)$。所以$j$的下一个可能值是
$
j = j_1 + j_2 - 1
$

至于说$j$的最小可能值$j_min$，我们可以通过总自由度的分析得出。从未耦合的角度来看，体系的总自由度（即CSCO表象基底的个数）是$(2j_1 + 1)(2j_2 + 1)$，而从耦合表象的角度来看，总自由度是
$
sum_(j = j_min)^(j_max) (2j + 1) = (2j_1 + 1)(2j_2 + 1)
$
得到$j_min = |j_1 - j_2|$。

所以，最后结论是
$
j = |j_1 - j_2|, |j_1 - j_2| + 1, ..., j_1 + j_2
$
或者记为
$
abs(j_1 - j_2) <= j <= j_1 + j_2
$
以及
$
m = m_1 + m_2
$
这个合成法则是量子力学中的重要法则。在直观上，这是矢量相加的三角形法则的结果，上述关系又称为*三角形关系*。

== 合成角动量的本征态

我们还需要知道角动量合成以后的本征态是什么，也就是说，要找出在式子
$
ket(  j","m";"j_1","j_2) = sum_(m = m_1 + m_2) C(j ,m; j_1, m_1, j_2, m_2) ket(j_1","m_1) ket(j_2","m_2)
$
中的系数$C(j ,m; j_1, m_1, j_2, m_2)$被称为*Clebsch-Gordan*（克莱布希-戈尔丹）系数（CG系数）

首先，从前面的分析中知道：只有在
$
j = j_1 + j_2, j_1 + j_2 - 1, ..., |j_1 - j_2|  
$
和
$
m = m_1 + m_2
$
才有$C(j ,m; j_1, m_1, j_2, m_2) != 0$：
$
ket(j","m";"j_1","j_2) = sum_(m = m_1 + m_2) C(j ,m; j_1, m_1, j_2, m_2) ket(j_1","m_1) ket(j_2","m_2)
$
#newpara()

其次，求CG系数的基本出发点就是让各量子态满足
$
hat(J)^2 ket(j","m";"j_1","j_2) = j(j + 1) hbar^2 ket(j","m";"j_1","j_2)\
hat(J)_z ket(j","m";"j_1","j_2) = m hbar ket(j","m";"j_1","j_2)\
hat(J)_1^2 ket(j_1","m_1) = j_1(j_1 + 1) hbar^2 ket(j_1","m_1)\
hat(J)_(1 z) ket(j_1","m_1) = m_1 hbar ket(j_1","m_1)\
hat(J)_2^2 ket(j_2","m_2) = j_2(j_2 + 1) hbar^2 ket(j_2","m_2)\
hat(J)_(2 z) ket(j_2","m_2) = m_2 hbar ket(j_2","m_2)
$
求出CG系数的一般表达式是一件相当困难的工作，其结果也相当复杂。我们只在下面给出两个具体的（也是很重要的）例子。

=== 电子的旋-轨耦合总角动量

子的总角动量就是它的轨道角动量和自旋角动量之和：
$
hat(arrow(J)) = hat(arrow(L)) + hat(arrow(S))
$
其中
$
j_1 = l = 0, 1, 2, ...\
j_2 = s = 1/2
$
所以
$
j = l+1/2 或 l-1/2
$
从而$j$一定是半整数。

轨道角动量的本征态是球谐函数$Y_(l m)$，而电子自旋本征态是$v_+$和$v_-$，，所以电子总角动量的本征态应该写成如下的形式：
$
ket(j","m) = c_1 ket(l","m-1/2) ket(1/2","1/2) + c_2 ket(l","m+1/2) ket(1/2","-1/2)
$
或者等价地写为二分量旋量的形式：
$
psi_(j m) = c_1 Y_(l,m-1/2) v_+ + c_2 Y_(l,m+1/2) v_- = mat(
  c_1 Y_(l,m-1/2);
  c_2 Y_(l,m+1/2)
)
$

#figure(
  image("pic/2024-05-28-14-07-29.png", width: 80%),
  numbering: none
)

#figure(
  image("pic/2024-05-28-14-07-55.png", width: 80%),
  numbering: none
)

从而对于$j = l + 1/2$来说：
$
ket(l + 1/2","m) = sqrt((j + m)/(2l+1)) ket(l","m-1/2) ket(1/2","1/2) + sqrt((j - m)/(2l+1)) ket(l","m+1/2) ket(1/2","-1/2)
$

*CG系数的特点及符号约定*：CG系数都为*实数*，同时在$m=j,m_1=j_1$时，系数为*非负实数*。

把这两组本征态用幺正变换联系起来就是：
$
mat(
  ket(l + 1/2","m);
  ket(l - 1/2","m)
) = 1/sqrt(2l + 1) 
mat(
  sqrt(l + 1/2 + m), sqrt(l + 1/2 - m);
  - sqrt(l + 1/2 + m), sqrt(l + 1/2 - m)
)
mat(
  ket(l","m - 1/2) ket(1/2","1/2);
  ket(l","m + 1/2) ket(1/2","-1/2)
)
$

其逆变换：
$
mat(
  ket(l","m - 1/2) ket(1/2","1/2);
  ket(l","m + 1/2) ket(1/2","-1/2)
) = 1/sqrt(2l + 1)
mat(
  sqrt(l + 1/2 + m), - sqrt(l + 1/2 + m);
  sqrt(l + 1/2 - m), sqrt(l + 1/2 - m)
)
mat(
  ket(l + 1/2","m);
  ket(l - 1/2","m)
)
$

下面通过几个例题来研究耦合表象$(j m j_1 j_2)$的一些性质。

_L-S耦合表象本征态是否也是算符$hat(S)_z$的本征态？如果不是求其在L-S平行耦合态下的平均值。_

考虑到$hat(S)_z$与$hat(J)_z, hat(L)^2, hat(S)^2$对易，但对于$hat(J)^2$：
$
[hat(J)^2, hat(S)_z] &= [hat(L)^2 + hat(S)^2 + hat(L)_+ hat(S)_- + hat(L)_- hat(S)_+ + 2 hat(L)_z hat(S)_z, hat(S)_z] \
& = hat(L)_+ [ hat(S)_-, hat(S)_z] + hat(L)_- [ hat(S)_+, hat(S)_z]\
& =^([hat(J)_plus.minus , hat(J)_z] = minus.plus hbar hat(J)_plus.minus) hbar( hat(L)_+ hat(S)_- - hat(L)_- hat(S)_+)\
&!= 0
$
从而L-S耦合表象本征态不是$hat(S)_z$的本征态。我们可以求出其平均值：
$
macron(S)_z = braket(j m , hat(S)_z , j m )
$
方法1：对平行耦合态（j=l+1/2），量子态表示为二分量旋量形式
$
ket(j m) = 1/sqrt(2l + 1) mat(
  sqrt(j + m) Y_(l,m-1/2);
  sqrt(j - m) Y_(l,m+1/2)
)
$
从而
$
macron(S)_z &= hbar/(2(2l+1)) integral mat(
  sqrt(j + m) Y_(l,m-1/2);
  sqrt(j - m) Y_(l,m+1/2)
)^dagger sigma_z mat(
  sqrt(j + m) Y_(l,m-1/2);
  sqrt(j - m) Y_(l,m+1/2)
) dd(tau)\
&= hbar/(2(2l+1)) integral mat(
  sqrt(j + m) Y_(l,m-1/2);
  sqrt(j - m) Y_(l,m+1/2)
)^dagger mat(
  sqrt(j + m) Y_(l,m-1/2);
  - sqrt(j - m) Y_(l,m+1/2)
) dd(tau)\
&= (m hbar)/(2l+1)
$
方法2：直接分解到本征态的形式
$
ket(j m) = C_1 ket(l m - 1/2) ket(1/2 1/2) + C_2 ket(l m + 1/2) ket(1/2 -1/2)
$
从而
$
macron(S)_z &= abs(C_1)^2 hbar/2 - abs(C_2)^2 hbar/2\
&= ((j+m)/(2l+1) - (j-m)/(2l+1)) hbar/2\
&= (m hbar)/(2l+1)
$
#newpara()

_例(2)：求下面算符$hat(S)_z$作用在耦合表象基底上的展开式_
$
hat(S)_z ket(j m l 1/2)
$
首先由$hat(S)_z$与$hat(J)_z, hat(L)^2, hat(S)^2$对易，展开式中$m,l,1/2$量子数固定不变。而$hat(S)_z$与$hat(J)^2$不对易，用$hat(S)_z$的本征态展开：
$
hat(S)_z ket(j m l 1/2) &= C_(l + 1/2) ket(l+1/2","m","l","1/2) + C_(l - 1/2) ket(l-1/2","m","l","1/2)
$
其实这是一个表象变换的问题，在旧基底（非耦合表象）中，$hat(S)_z$是对角化的。即在基底
$
mat(
  ket(l","m-1/2) ket(1/2","1/2);
  ket(l","m+1/2) ket(1/2","-1/2)
)
$
下，$hat(S)_z$的矩阵形式是
$
hbar/2 mat(
  1, 0;
  0, -1
)
$
而新基底（耦合表象）用旧基底表示的展开式是
$
mat(
  ket(l + 1/2","m);
  ket(l - 1/2","m)
) &= 1/sqrt(2l + 1) 
mat(
  sqrt(l + 1/2 + m), sqrt(l + 1/2 - m);
  - sqrt(l+ 1/2 + m), sqrt(l + 1/2 - m)
)
mat(
  ket(l","m - 1/2) ket(1/2","1/2);
  ket(l","m + 1/2) ket(1/2","-1/2)
)\
& = U^* mat(
  ket(l","m - 1/2) ket(1/2","1/2);
  ket(l","m + 1/2) ket(1/2","-1/2)
)
$
于是$hat(S)_z$在新基底的矩阵形式是
$
U hbar/2 mat(
  1, 0;
  0, -1
) U^dagger &= hbar/(2(2l + 1)) mat(
  sqrt(l + 1/2 + m), sqrt(l + 1/2 - m);
  - sqrt(l + 1/2 + m), sqrt(l + 1/2 - m)
) mat(
  1, 0;
  0, -1
) mat(
  sqrt(l + 1/2 + m), - sqrt(l + 1/2 + m);
  sqrt(l + 1/2 - m), sqrt(l + 1/2 - m)
)\
&= hbar/(2(2l+1)) mat(
  2m,-sqrt((2l+1)^2 - 4m^2);
  -sqrt((2l+1)^2 - 4m^2), -2m
)
$
也就是说
$
hat(S)_z mat(
  ket(l + 1/2","m);
  ket(l - 1/2","m)
) &= hbar/(2(2l+1)) mat(
  2m,-sqrt((2l+1)^2 - 4m^2);
  -sqrt((2l+1)^2 - 4m^2), -2m
) mat(
  ket(l + 1/2","m);
  ket(l - 1/2","m)
)\
$
即有
$
braket(l plus.minus 1/2 "," m "," l "," 1/2 , hat(S)_z , l plus.minus 1/2 "," m "," l "," 1/2) = plus.minus (hbar m)/(2l+1)
$
这就是例1的结果。

方法二：先让$hat(S)_z$作用于直积态基底，然后把直积态基底转换为耦合态基底

#figure(
  image("pic/2024-05-30-16-46-15.png", width: 80%),
  numbering: none
)

#figure(
  image("pic/2024-05-30-16-47-05.png", width: 80%),
  numbering: none
)

_例(3)：算符$hat(arrow(L)) dot hat(arrow(S))$的本征态是非耦合表象的基底，还是L-S耦合表象的基底？_

$
hat(arrow(L)) dot hat(arrow(S)) = 1/2(hat(J)^2 - hat(L)^2 - hat(S)^2)
$
所以$hat(arrow(L)) dot hat(arrow(S))$和$hat(J)^2, hat(L)^2, hat(S)^2$对易，同时
$
[hat(J)^2, hat(S)_z] != 0, [hat(J)^2, hat(L)_z] != 0
$
所以$hat(arrow(L)) dot hat(arrow(S))$的本征态*是L-S耦合表象的基底*。

=== 碱金属原子光谱双线结构

电子轨道角动量产生磁场必定与电子本身自旋产生的磁距发生相互作用，从而改变原子能级，使光谱线产生分裂。

电子绕原子核旋转，在电子静止坐标系中看，等效于原子核$(+Z e)$绕电子旋转。把原子核绕电子转动想象成一个半径为$a$的电流环，则电流环圆心处磁感应强度$B$及环电流$I$为
$
B = (mu_0 I)/(2 a), I = (Z e v)/(2 pi a)\
B = (mu_0 Z e v)/(4 pi a^2), arrow(B) = (mu_0 Z e)/(4 pi a^3) arrow(a) crossproduct arrow(v)
$
于是电子旋-轨耦合能量为
$
- arrow(mu) dot arrow(B) = - ((-e)/m) arrow(S) dot (mu_0 Z e)/(4 pi a^3) arrow(a) crossproduct arrow(v) = 1/(m^2 c^2) 1/(4 pi epsilon_0) (Z e^2)/a^3 arrow(S) dot (arrow(a) crossproduct arrow(p)) = 1/(m^2 c^2) 1/a dd(V)/dd(a) arrow(S) dot arrow(L)
$
其中
$
c= 1/sqrt(epsilon_0 mu_0),V= - (Z e^2)/(4 pi epsilon_0 a), dd(V)/dd(a) = (Z e^2)/(4 pi epsilon_0 a^2) >0
$
正确的表达式还应加入Thomas进动修正（相对论修正，1926），所以最后（用$r$代替$a$）
$
- arrow(mu) dot arrow(B) = 1/(2 m^2 c^2) 1/r dd(V)/dd(r) arrow(S) dot arrow(L)
$
这一结果也可由狄拉克方程在非相对论极限下给出。由此看出，当旋-轨角动量平行时，耦合能量为正，反之为负。考虑到*旋-轨耦合*后的哈密顿算符为：
$
hat(H) = hat(p)^2/(2 mu) + V(r) + xi(r) arrow(L) dot arrow(S)
$
其中
$
xi(r) = 1/(2 m^2 c^2) dd(V)/dd(r)
$
这里$V(r)$应理解为库仑屏蔽势（碱金属原子内层电子对核有屏蔽作用）。由于有$hat(arrow(L))dot hat(arrow(S))$项，所以能级一般与量子数$n,l, j$都有关系：
$
hat(arrow(L)) dot hat(arrow(S)) = 1/2(hat(J)^2 - hat(L)^2 - hat(S)^2)
$
不考虑微扰项$xi(r) arrow(L) dot arrow(S)$，系统本征量子态的角度部分为（取L-S耦合表象）$ket(j","m_j","l","1/2)$，空间径向部分为$ket(n","l)$，整体就是
$
ket(n","j","m_j","l","1/2) = ket(n","l) ket(j","m_j","l","1/2)
$
如果把耦合项看作微扰，则耦合项引起的附加能量近似为
$
Delta E &= braket(n","j","m_j","l","1/2 , xi(r) arrow(L) dot arrow(S) , n","j","m_j","l","1/2) \
&= braket(n l, xi(r), n l) braket(j","m_j","l","1/2 , arrow(L) dot arrow(S) , j","m_j","l","1/2)\
&= xi_(n l) hbar^2/2 (j(j+1) - l(l+1) - 3/4) = cases(
  1/2 hbar^2 xi_(n l) 当(j= l+1/2) ,  - (l + 1)/2 hbar^2 xi_(n l)当 (j= l-1/2)
)
$
其中$xi_(n l) = braket(n l, xi(r), n l)$，所以由于旋-轨耦合作用使原来的每条能级分裂成了两条。

钠黄线的双线分裂：
#figure(
  image("pic/2024-05-30-17-17-05.png", width: 80%),
  caption: [
    钠黄线的双线分裂
  ],
)
在考虑旋-轨耦合作用后，钠原子$3P$能级分裂为$3P_(3/2)$和$3P_(1/2)$。其中前者的简并度为4，后者的简并度为2。

=== 反常塞曼效应

前面我们讲到了由于*旋-轨耦合*$xi(r) arrow(L) dot arrow(S)$产生的*碱金属原子的双线结构*。

由于*磁-轨耦合*$(e B)/(2 mu)hat(L)_z$产生的*正常塞曼效应*，以及*电子在外磁场中的能量*$(e B)/(mu)hat(S)_z$。

现在考虑*旋-轨耦合*和*磁-轨耦合*的共同作用，即*反常塞曼效应*。同时外加磁场$B$较弱，后两项与旋-轨耦合能量相当的情况。这时哈密顿算符的形式为
$
hat(H) &= hat(p)^2/(2 mu) + V(r) + xi(r) arrow(L) dot arrow(S) + (e B)/(2 mu) (hat(L)_z + 2 hat(S)_z)\
&= hat(p)^2/(2 mu) + V(r) + xi(r)/2 (hat(J)^2 - hat(L)^2 - hat(S)^2) + (e B)/(2 mu) hat(J)_z + (e B)/(2 mu) hat(S)_z
$
如果没有最后一项$(e B)/(2 mu) hat(S)_z$，根据前面的讨论，可以使用旋-轨耦合表象来表示系统的本征态。

设无外磁场时系统本征能量和本征态为
$
E_(n l j) , ket(n "," j "," m_j "," l "," 1/2)
$
每条能级是$(2j+1)$重简并。现在考虑加入$(e B)/(2 mu) hat(J)_z$项，则因为这一项与原哈密顿算符对易，系统量子态不变，但能级会多出一项变为：
$
E_(n l j m_j) = E_(n l j) + (e B)/(2 mu) m_j hbar
$
这样$(2j+1)$重简并就被完全消除了。

现在再考虑加入最后一项$(e B)/(2 mu) hat(S)_z$，由于这一项与原哈密顿算符不对易，所以新的本征态函数很难求出。但是如果仍沿用原有的波函数态$ket(n "," j "," m_j "," l "," 1/2)$，同时把最后一项看作*微扰*，则其对原能级的微扰修正为
$
Delta E = braket(n "," j "," m_j "," l "," 1/2 , (e B)/(2 mu) hat(S)_z , n "," j "," m_j "," l "," 1/2)
$
根据前面例2的计算，其结果为
$
Delta E = plus.minus (e B hbar)/(2 mu (2l + 1) ) m_j
$
于是最后修正后的能级为
$
E = cases(
  E_(n l j) + (e B)/(2 mu) (1  + 1/(2l + 1)) hbar m_j 当(j = l + 1/2) , E_(n l j) + (e B)/(2 mu) (1  - 1/(2l + 1)) hbar m_j 当(j = l - 1/2)
)
$

钠双黄线在弱磁场下的分裂：
#figure(
  image("pic/2024-06-07-12-47-42.png", width: 80%),
  caption: [
    钠双黄线在弱磁场下的分裂
  ],
)

与正常塞曼效应相比，反常塞曼效应是光谱线分裂为*偶数条*。

== 两个电子自旋的合成

设$hat(arrow(S))_1$和$hat(arrow(S))_2$是两个电子自旋，它们的和是
$
hat(arrow(S)) = hat(arrow(S))_1 + hat(arrow(S))_2
$
现在$j_1 = j_2 = 1/2$，所以
$
S = 1 , 0
$
当两个电子的自旋互相平行的时候$S=1$，而当它们是反平行的时候$S=0$。

$S=1$是三重态因为$m = 1,0,-1$，$S=0$是单态。

设
$
nu_+ = ket(arrow.t) , nu_- = ket(arrow.b)
$
于是
$
S_z ket(arrow.t) = hbar/2 ket(arrow.t) , S_z ket(arrow.b) = - hbar/2 ket(arrow.b)
$
还有
$
S_+ ket(arrow.t) = 0 , S_+ ket(arrow.b) = hbar ket(arrow.t)\
S_- ket(arrow.t) = hbar ket(arrow.b) , S_- ket(arrow.b) = 0
$
和
$
S^2 ket(arrow.t) = 3/4 hbar^2 ket(arrow.t) , S^2 ket(arrow.b) = 3/4 hbar^2 ket(arrow.b)
$
对于两个自旋$hat(arrow(S))_1$和$hat(arrow(S))_2$未耦合的本征态，在$(sigma_(z 1), sigma_(z 2))$表象记为
$
nu_(1 +) nu_(2 +) = ket(arrow.t "," arrow.t), nu_(1 +) nu_(2 -) = ket(arrow.t "," arrow.b), nu_(1 -) nu_(2 +) = ket(arrow.b "," arrow.t), nu_(1 -) nu_(2 -) = ket(arrow.b "," arrow.b)
$
而耦合后的本征态记为$ket(S","m)$，取
$
ket(1","1) , ket(1","0) , ket(1","-1) , ket(0","0)
$
对于$S=1$的情形，当$m=1,-1$时，耦合的本征态也就是未耦合的本征态$m = m_1 + m_2$，即
$
ket(1","1) = ket(arrow.t "," arrow.t) , ket(1","-1) = ket(arrow.b "," arrow.b)
$
但是$m=0$应该是线性组合，即
$
ket(1","0) = C_1 ket(arrow.t "," arrow.b) + C_2 ket(arrow.b "," arrow.t)
$
从最大投影态出发$ket(1","1) = ket(arrow.t "," arrow.t)$，两边作用$S_-$：
$
S_- ket(1","1) &= S_(1-) ket(arrow.t "," arrow.t) + S_(2-) ket(arrow.t "," arrow.t)\
sqrt(2) ket(1","0) &= ket(arrow.b "," arrow.t) + ket(arrow.t "," arrow.b)\
ket(1","0) &= 1/sqrt(2) (ket(arrow.b "," arrow.t) + ket(arrow.t "," arrow.b))
$
两边再作用一次降算符：
$
S_- ket(1","0) &= S_(1-) ket(arrow.b "," arrow.t) + S_(2-) ket(arrow.t "," arrow.b)\
sqrt(2) ket(1","-1) &= 1/sqrt(2)(ket(arrow.b "," arrow.b) + ket(arrow.b "," arrow.b))\
ket(1","-1) &= ket(arrow.b "," arrow.b)
$
同理可以得到$S=0$的情形：
$
ket(0","0) = 1/sqrt(2) (ket(arrow.b "," arrow.t) - ket(arrow.t "," arrow.b))
$
交换两电⼦的自旋后（设自旋交换算符为$hat(P)_(12)$）：
$
hat(P)_(12) ket(1","0) = ket(1","0) , hat(P)_(12) ket(0","0) = - ket(0","-0)
$

#newpara()

得到基底转换矩阵：
$
mat(
  ket(1","1);
  ket(1","-1);
  ket(1","0);
  ket(0","0)
) = mat(
  1,0 ,0 ,0 ;
  0,1 ,0 ,0 ;
  0 ,0 , 1/sqrt(2), 1/sqrt(2);
  0  ,0 , 1/sqrt(2), -1/sqrt(2)
)
mat(
  ket(arrow.t "," arrow.t);
  ket(arrow.b "," arrow.b);
  ket(arrow.t "," arrow.b);
  ket(arrow.b "," arrow.t);
)
$
即
$
U^* = mat(
  1,0 ,0 ,0 ;
  0,1 ,0 ,0 ;
  0 ,0 , 1/sqrt(2), 1/sqrt(2);
  0  ,0 , 1/sqrt(2), -1/sqrt(2)
)
$
算符$hat(S)_+ = hat(S)_(1+) + hat(S)_(2+)$在旧表象$mat(
  ket(arrow.t "," arrow.t);
  ket(arrow.b "," arrow.b);
  ket(arrow.t "," arrow.b);
  ket(arrow.b "," arrow.t);
)$下的矩阵形式是
$
mat(
  0, 0, hbar, hbar;
  0, 0, 0, 0;
  0, hbar, 0, 0;
  0, hbar, 0, 0
)
$
在新表象$mat(
  ket(1","1);
  ket(1","-1);
  ket(1","0);
  ket(0","0)
)$下的矩阵形式是
$
S_+ = U mat(
  0, 0, hbar, hbar;
  0, 0, 0, 0;
  0, hbar, 0, 0;
  0, hbar, 0, 0
) U^dagger = hbar mat(
  0, 0, sqrt(2), 0;
  0, 0, 0, 0;
  0, sqrt(2), 0, 0;
  0, 0, 0, 0
)
$
算符$hat(S)_x,hat(S)_y,hat(S)_z$在新表象中的表示是模块对角化的，也就是说具有如下形式
$
mat(
  F_11,F_12,F_13,0;
  F_21,F_22,F_23,0;
  F_31,F_32,F_33,0;
  0,0,0,F_44
)
$
从非耦合表象转换到耦合表象，就是把角动量算符的矩阵表示转变为*模块对角化（block diagonal）*的形式。

算符模块对角化的意义：每一个模块形成一个独立的子空间，不同子空间内的态矢量在转动变换下不会相互转化，只在各自所属子空间内部相互转化。如：双电子自旋三重态矢量通过角动量算符作用只在三重态内相互转换，而不会转化为单态。

#figure(
  image("pic/2024-06-14-13-11-38.png", width: 80%),
  numbering: none,
)

#figure(
  image("pic/2024-06-14-13-12-19.png", width: 80%),
  numbering: none,
)

#figure(
  image("pic/2024-06-14-13-16-27.png", width: 80%),
  numbering: none,
)

#figure(
  image("pic/2024-06-14-13-16-41.png", width: 80%),
  numbering: none,
)

= 定态微扰论

== 非简并情形

可以精确求解的量子力学问题是不多的，所以近似方法有重要的作用。微扰论是主要的近似方法之一（其它还有变分法、WKB法等）。

零级定态薛定鄂方程：
$
hat(H)^((0)) psi_n^((0)) = E_n^((0)) psi_n^((0))
$

其中$hat(H)^((0))$是容易解出的哈密顿算符，如氢原子系统、自由粒子

假设加入微扰能$hat(H)^'$，则薛定鄂方程形式应该为：
$
hat(H) psi_n = E_n psi_n, hat(H) = hat(H)^((0)) + hat(H)^'
$
$hat(H)^'$是$hat(H)^((0))$的小修正，方程的解可以用级数形式逐级展开：
$
hat(H)^' << hat(H)^((0))\
E_n = E_n^((0)) + E_n^((1)) + E_n^((2)) + ...\
psi_n = psi_n^((0)) + psi_n^((1)) + psi_n^((2)) + ...
$
其中$E_n^((0) )$和$psi_n^((0))$与$hat(H)^'$无关，而$E_n^((1))$和$psi_n^((1))$是一级微扰，和$hat(H)^'$的一次方成正比，$E_n^((2))$和$psi_n^((2))$是二级微扰......一般情况下，越高次的项越小，所以可以只保留最低的几阶，便有足够的精度。

把上述展开式代入原方程，得：
$
(hat(H)^((0)) + hat(H)^') (psi_n^((0)) + psi_n^((1)) + psi_n^((2)) + ...) = (E_n^((0)) + E_n^((1)) + E_n^((2)) + ...) (psi_n^((0)) + psi_n^((1)) + psi_n^((2)) + ...)\
(hat(H)^((0)) - E_n^((0)) + H' - E_n^((1)) - E_n^((2)) + ...) (psi_n^((0)) + psi_n^((1)) + psi_n^((2)) + ...) = 0
$
逐级比较方程两端就得到一系列方程，解出能量修正及本征函数修正。

零级方程就是无微扰时$hat(H)^((0))$的本征方程
$
(hat(H)^((0)) - E_n^((0))) psi_n^((0)) = 0
$
一级方程是
$
(hat(H)^((0)) - E_n^((0))) psi_n^((1)) = (E_n^((1)) - hat(H)^') psi_n^((0))
$
先处理非简并情形，即$hat(H)^((0))$的属于$E_n^((0))$的本征态只有一个。把$psi_n^((1))$按$hat(H)^((0))$表象的基底$psi_n^((0))$展开：
$
psi_n^((1)) = sum_m a_(n m)^((1)) psi_m^((0)), "其中"  sum_m abs(a_(n m)^((1)))^2 != 1
$
再代入方程中得：
$
sum_m a_(n m)^((1)) (hat(H)^((0)) - E_n^((0))) psi_m^((0)) = - (hat(H)^' - E_n^((1))) psi_n^((0))\
sum_m a_(n m)^((1)) (E_m^((0)) - E_n^((0))) psi_m^((0)) = - (hat(H)^' - E_n^((1))) psi_n^((0))
$
利用$psi_m^((0))$的正交性，得到
$
a_(n m )^((1)) (E_m^((0)) - E_n^((0)) )= - braket(psi_m^((0)) , hat(H)^' , psi_n^((0))) + E_n^((1)) delta_(m n)
$
由于$m$的任意性，取$m=n$，得到
$
E_n^((1)) = braket(psi_n^((0)) , hat(H)^' , psi_n^((0))) = H_(n n)^'
$
其中$H_(n n)^'$就是$hat(H)'$在$hat(H)^((0))$表象中的对角矩阵元，或者说在$psi_n^((0))$中的平均值。

考虑$m != n$的情况，可以求出系数$a_(n m)^((1))$：
$
a_(n m)^((1)) = - braket(psi_m^((0)) , hat(H)^' , psi_n^((0))) / (E_m^((0)) - E_n^((0))) = H_(m n)^' / (E_n^((0)) - E_m^((0)))
$
一级微扰波函数：
$
psi_n^((1)) = sum_(m != n) H_(m n)^'  / (E_n^((0)) - E_m^((0))) psi_m^((0))
$
由此，我们发现微扰论适用的条件是：
$
abs(H_(m n)^' / (E_n^((0)) - E_m^((0)))) << 1
$
利用波函数归一化：
$
braket(psi_n) = braket(psi_n^((0)) + psi_n^((1)) + psi_n^((2)) + ...) = 1
$
已知
$
braket(psi_n^((0))) = 1
$
所以逐级比较可得
$
braket(psi_n^((0)), psi_n^((1))) + braket(psi_n^((1)), psi_n^((0))) = 0\
braket(psi_n^((0)), psi_n^((2))) + braket(psi_n^((1)), psi_n^((1))) + braket(psi_n^((2)), psi_n^((0))) = 0
$
其中
$
psi_n^((1)) = sum_m a_(n m)^((1)) psi_m^((0))\
$
对$psi_n^((1))$来说，这就要求
$
a_(n n)^((1)) + a_(n n)^((1)*) = 0
$
有
$
a_(n n)^((1)) = i a_n^((1))
$
其中$a_n^((1))$是实数，且数量级与其它$a_(m n)^((1))$一致。

也能得到
$
a_(m n)^((1)) + a_(n m)^((1)*) = 0
$
而
$
a_(m n)^((1)) = H_(m n)^' / (E_n^((0)) - E_m^((0)))
$
所以
$
H' = H'^dagger
$
于是，展开到一级修正的波函数为（其中$O(a^((2)))$代表所有$a$的二次项及更高次项的集合，在每行等式里此项不一定相同）
$
psi_n &= psi_n^((0)) + i a_n^((1)) psi_n^((0)) + sum_(m != n) H_(m n)^' / (E_n^((0)) - E_m^((0))) psi_m^((0)) + O(a^((2)))\
&= e^(i a_n^((1))) psi_n^((0)) + e^(i a_n^((1))) sum_(m != n) H_(m n)^' / (E_n^((0)) - E_m^((0))) psi_m^((0)) + O(a^((2)))\
&= e^(i a_n^((1))) (psi_n^((0)) + sum_(m != n) H_(m n)^' / (E_n^((0)) - E_m^((0))) psi_m^((0)) ) + O(a^((2)))
$
可见$a_n^((1))$的效果就是给波函数乘上一个相位因子，波函数所包含的物理信息不变，不妨设$a_n^((1)) = 0$。

二级微扰方程是：
$
(hat(H)^((0)) - E_n^((0))) psi_n^((2)) = - (hat(H)^' - E_n^((1))) psi_n^((1)) + E_n^((2)) psi_n^((0))
$
将$psi_n^((1))$代入，得
$
(hat(H)^((0)) - E_n^((0))) psi_n^((2)) = - (hat(H)^' - E_n^((1))) sum_(m != n) H_(m n)^' / (E_n^((0)) - E_m^((0))) psi_m^((0)) + E_n^((2)) psi_n^((0))
$
然后方程两边左乘以$psi_n^((0))$并积分，左边为0，右边第二项为0，可以得到二阶修正能量$E_n^((2))$：
$
E_n^((2)) = sum_(m != n) H_(m n)^' / (E_n^((0)) - E_m^((0))) braket(psi_n^((0)) , hat(H)' ,psi_m^((0))) = sum_(m != n) (H_(m n)^'  H_(m n)^')/ (E_n^((0)) - E_m^((0))) = sum_(m != n) abs(H_(m n)^')^2 / (E_n^((0)) - E_m^((0)))
$
最终可以得到修正的波函数和能量：
$
psi_n = psi_n^((0)) + sum_(m != n) H_(m n)^' / (E_n^((0)) - E_m^((0))) psi_m^((0)) + ...\
E_n = E_n^((0)) + H_(n n)^' + sum_(m != n) abs(H_(m n)^')^2 / (E_n^((0)) - E_m^((0))) + ...
$

=== 在静电场中的一维谐振子

假设一维谐振子还带有电荷$q$，并处在外加恒定电场$E$（沿$x$轴正向）中，那么哈密顿量是
$
hat(H) = hat(H)^((0)) + hat(H)^'\
hat(H)^((0)) = hat(p)^2/(2 mu) + 1/2 mu omega^2 hat(x)^2\
hat(H)^' = - q E hat(x)
$
一级微扰能是
$
E_n^((1)) = braket(psi_n^((0)) , hat(H)^' , psi_n^((0)) )= - q E braket(psi_n^((0)) , hat(x) , psi_n^((0))) = 0
$
继续检查能级的二阶修正：
$
H'_(m n) = - q E braket(psi_m^((0)) , hat(x) , psi_n^((0)))
$
可利用递推关系：
$
hat(x) psi_n^((0)) = sqrt(hbar/(2 mu omega)) (sqrt(n+1) psi_(n+1)^((0)) + sqrt(n) psi_(n-1)^((0)))
$
于是
$
H'_(m n) &= - q E sqrt(hbar/(2 mu omega)) (sqrt(n+1) braket(psi_m^((0)) , psi_(n+1)^((0)) ) + sqrt(n) braket(psi_m^((0)) , psi_(n-1)^((0)) ))\
& = - q E sqrt(hbar/(2 mu omega)) (sqrt(n+1) delta_(m , n+1) + sqrt(n) delta_(m , n-1))
$
所以二级微扰能是
$
E_n^((2)) &= sum_(m != n) abs(H_(m n)^')^2 / (E_n^((0)) - E_m^((0)) ) = abs(H'_(n-1 ,n))^2 / (E_n^((0)) - E_(n-1)^((0)) ) + abs(H'_(n+1 ,n))^2 / (E_n^((0)) - E_(n+1)^((0)) )\
&= q^2 E^2 hbar / (2 mu omega) (n + 1) / (hbar omega n - hbar omega (n+1)) + q^2 E^2 hbar / (2 mu omega) n / (hbar omega n - hbar omega (n-1))\
&= - (q^2 E^2) / (2 mu omega^2)
$
所以微扰以后的能级是（准确到二级微扰）
$
E_n = (n + 1/2) hbar omega - (q^2 E^2) / (2 mu omega^2)
$
这个微扰能与$n$无关。实际上，这个问题是有精确解的
$
V(x) &= 1/2 mu omega^2 x^2 - q E x\
&= 1/2 mu omega^2 (x - (q E) / (mu omega^2))^2 -( q^2 E^2) / (2 mu omega^2)
$
它的第一项只不过是把原来的谐振子势能平移了一段距离，这个移动不会影响谐振子的能级，而它的第二项正是前面求出的与$n$无关的能级修正。

== 简并情形

一级微扰能和零级波函数：
$
E_n = E_n^((0)) + E_n^((1)) + ...\
psi_n = psi_n^((0)) + psi_n^((1)) + ...\
hat(H)^((0)) psi_n^((0)) = E_n^((0)) psi_n^((0))
$
$E_n^((0))$有$k$度简并，本征波函数为
$
psi_(n 1)^((0)) , psi_(n 2)^((0)) , ... , psi_(n k)^((0))\
hat(H)^((0)) psi_(n i)^((0)) = E_n^((0)) psi_(n i)^((0)) , i = 1,2,...,k
$
其中不同的本征态之间正交。在引入微扰后应设：
$
psi_n^((0)) = sum_i c_(n i) psi_(n i)^((0))
$
代入一级微扰方程
$
(hat(H)^((0)) - E_n^((0))) psi_n^((1)) = (E_n^((1)) - hat(H)^') psi_n^((0))
$
得：
$
(hat(H)^((0)) - E_n^((0))) psi_(n)^((1)) = (E_n^((1)) - hat(H)^') sum_i c_(n i)^((0)) psi_(n i)^((0))
$
两端左乘以$psi_(n j)^((0)*)$并积分，得：
$
sum_i c_(n i)^((0)) ( braket(psi_(n j)^((0)) , hat(H)^' , psi_(n i)^((0))) - E_n^((1)) braket(psi_(n j)^((0)) , psi_(n i)^((0))) ) = 0\
sum_i c_(n i)^((0)) ( H_(j i)^' - E_n^((1)) delta_(j i) ) = 0
$
这里的$H_(j i)^' = braket(psi_(n j)^((0)) , hat(H)^' , psi_(n i)^((0)) )$注意这里$n$是固定的。这和矩阵形式的本征方程完全一样，$c_(n i)^((0))$有非0解的条件是
$
det(H' - E_n^((1))I) = 0
$
久期方程
$
det
  mat(
  H_(1 1)^' - E_n^((1)) , H_(1 2)^' , ... , H_(1 k)^' ;
  H_(2 1)^' , H_(2 2)^' - E_n^((1)) , ... , H_(2 k)^' ;
  dots.v, dots.v, dots.down, dots.v ;
  H_(k 1)^' , H_(k 2)^' , ... , H_(k k)^' - E_n^((1))
) = 0
$
从中可以解出$E_n^((1))$和对应的展开系数$c_(n i)^((0))$。这就决定了一级微扰能和零级波函数。

一般来说$E_n^((1))$不仅和对角线元素$H_(i i)^'$有关，也和非对角线元素$H_(i j)^'$有关，但总有$k$个解。

假如$E_n^((1))$的$k$个解各不相同（方程没有重根），则$E_n^((0))$的简并度被完全消除，否则只可能是部分被消除。

== 原子能级在静电场中的分裂

原子能级在静电场中的分裂称为*Stark效应*。作为例子，让我们考虑氢原子。

设均匀外静电场$E$沿着正$z$轴方向，那么氢原子就受到了如下的附加势能：
$
hat(H)^' = e E z = e E r cos theta
$
在未加微扰时，氢原子的能级是
$
E_n^((0)) = - (m k_1^2 e^4)/(2 hbar^2 n^2) =  -(k_1 e^2)/(2 a_0 n^2)
$
对应波函数是
$
psi_n^((0)) = R_(n l) Y_(l m)
$

能级$n$的简并度为$n^2$，对$n=2$来说是4度简并（不考虑自旋）

能级$n=2$的简并态量子数$l, m$可以取值00、10、11、1-1，并依次简记为第1、2、3、4个态。于是要计算
$
H' = braket(psi_(2 l' m')^((0)), hat(H)^', psi_(2 l m)^((0))) = e E integral r^2 R_(2 l') R_(2 l) Y_(l' m') Y_(l m) r cos theta r^2 sin theta dd(theta) dd(phi) dd(r)
$

#figure(
  image("pic/2024-06-19-12-15-30.png", width: 80%),
  numbering: none,
)

#figure(
  image("pic/2024-06-15-01-22-21.png", width: 80%),
  numbering: none,
)

我们要求出下面这个久期方程从而得到$E_2^((1))$
$
det
  mat(
     - E_2^((1)) , -3e E a_0, 0 , 0 ;
    -3e E a_0 ,  - E_2^((1)) , 0 ,0;
    0 , 0 ,  - E_2^((1)) , 0 ;
    0 , 0, 0 ,  - E_2^((1))
  )
= 0
$
得到
$
E_2^((1)) = 0, plus.minus 3 e E a_0
$
这就是说，原来简并在$n=2$上的$4$个能级，现在有一个向上移动了$3 e E a_0$，一个向下移动了$3 e E a_0$，还有两个没有移动，简并是部分地消除了。

#figure(
  image("pic/2024-06-19-13-39-48.png", width: 80%),
  caption: [
    在电场中氢原子能级的分裂
  ],
)

钠双黄线问题中所用微扰近似：简并微扰

#figure(
  image("pic/2024-06-19-13-40-49.png", width: 80%),
  numbering: none,
)

#figure(
  image("pic/2024-06-19-13-44-51.png", width: 80%),
  numbering: none,
)

#figure(
  image("pic/2024-06-19-13-47-30.png", width: 80%),
  numbering: none,
)

#figure(
  image("pic/2024-06-19-13-47-44.png", width: 80%),
  numbering: none,
)

= 散射理论

== 散射波函数

如果入射波是波矢为$arrow(k)$的三维平面波，散射势能中心在原点，可以证明在$r→∞$时，散射波具有球面波的形式，综合起来：
$
lim_(r->oo) psi_(arrow(r)) = e^(i arrow(k) dot arrow(r)) + f(theta, phi)  e^(i arrow(k) dot arrow(r))/r
$
类似一维散射问题，这里只关心散射部分相对入射部分的概率，所以不必关心平面波的正交归一化系数（这里直接设为1）。

入射波概率流密度：
$
arrow(J)_i = rho arrow(v) = abs(psi_i)^2 (hbar arrow(k)) /m = (hbar arrow(k)) /m
$
散射波概率流密度：
$
arrow(J)_s = (i hbar)/(2 m) (psi_s grad psi_s^* - psi_s^* grad psi_s) = (hbar k)/( r^2 m) |f(theta, phi)|^2 arrow(e)_r
$

== 散射截面

散射相对入射的大小：
$
J_s/J_i = abs(f(theta, phi))^2/r^2
$
依赖于$r$。应该考虑在$(theta, phi)$附近立体角内的概率密度流，应计算的比例关系是
$
(J_s r^2 dd(Omega))/J_i = abs(f(theta, phi))^2 dd(Omega)
$
上式计算的是散射后*单位时间*内通过立体角$dd(Omega)$对应的球面微元$r^2 dd(Omega)$的*概率*，相对入射波*单位时间*内通过*单位横截面积*的概率的大小。可惜这两个量的量纲又不同。

解决方案：设入射波通过横截面微元$dd(sigma)$的概率，经散射后全部通过立体角$dd(Omega)$对应的球面微元$r^2 dd(Omega)$流出，则有
$
J_s r^2 dd(Omega) = J_i dd(sigma)
$
于是：
$
dd(sigma) = (J_s r^2)/J_i dd(Omega) = abs(f(theta, phi))^2dd(Omega) = sigma(theta, phi) dd(Omega)
$
其中$sigma(theta, phi) = abs(f(theta, phi))^2$称为称为*微分散射横截面积*(其中$f(theta,phi)$为*微分散射振幅*)，它的积分给出*总散射截面*：
$
sigma_"tot" = integral sigma(theta, phi) dd(Omega)
$
跟一维问题不同，在三维散射问题中，我们用*散射截面来作为散射发生强烈程度的度量*。散射势场对入射波的散射越强烈，散射截面就越大，即在入射波概率流密度保持不变的情况下，有更高的概率被势场散射。单位时间内入射波被散射的总概率为
$
J_i sigma_"tot" ("总散射速率")
$
单位时间内入射波被散射到$(theta, phi)$附近单位立体角内的概率为
$
J_i sigma(theta, phi) ("微分散射速率")
$

#newpara()

设无散射微扰势能$V$时系统定态方程(自由粒子)为：
$
hat(H)_0 ket(psi_0) = E ket(psi_0)
$
加上微扰$V$后：
$
(hat(H)_0 + hat(V)) ket(psi) = E ket(psi)\
cases(
  (E - hat(H)_0) ket(psi_0) = 0,
  (E - hat(H)_0) ket(psi) = hat(V) ket(psi)
)
$
两式相减得：
$
(E - hat(H)_0) (ket(psi) - ket(psi_0)) = hat(V) ket(psi)
$
可以形式解出$ket(psi)$*（Lippman-Schwinger方程）*：
$
ket(psi) = ket(psi_0) + (E - hat(H)_0)^(-1) V ket(psi) = ket(psi_0) +  hat(G) V ket(psi)
$
其中$hat(G) = (E - hat(H)_0)^(-1)$是*Green函数*。

第一项$ket(psi_0)$代表无微扰时的0级波函数，第二项$hat(G) V ket(psi)$代表微扰修正，$hat(G) = (E - hat(H)_0)^(-1)$为与传播子相关的*格林算符*。

可以用迭代法求级数解，即方程右边的$ket(psi)$用0级近似$ket(psi_0)$代替，求得$ket(psi)$后再代入方程的右边，如此循环往复得：
$
ket(psi) &= ket(psi_0) +  hat(G) V ket(psi_0) +  hat(G) V hat(G) V ket(psi_0) + ...\
&= (1 + hat(G)hat(T)_s) ket(psi_0), hat(T)_s = V + V hat(G) V + ...
$
如果$hat(T)_s$取1级近似（*波恩近似*），则
$
ket(psi) = (1 + hat(G) V) ket(psi_0)
$
$
braket(arrow(r),psi) = braket(arrow(r),psi_0) + braket(arrow(r),hat(G) V, psi_0) = braket(arrow(r),psi_0) + integral  dd(""^3 r') braket(arrow(r),hat(G), arrow(r')) braket(arrow(r'),V, psi_0)
$
其中$braket(arrow(r'),V, psi_0)$表示粒子在$arrow(r)'$处被散射，Green函数$braket(arrow(r),hat(G), arrow(r'))$表示粒子从$arrow(r')$传播到$arrow(r)$。下面求坐标表象中的格林函数：
$
braket(arrow(r),hat(G), arrow(r')) = integral  dd(""^3 k') braket(arrow(r),1/(E-hat(H)_0), arrow(k')) braket(arrow(k'), arrow(r')) = 1/(2 pi)^3 integral  dd(""^3 k') e^(i arrow(k') dot (arrow(r) - arrow(r'))) / (E - (hbar^2 k'^2) /( 2 m))
$
$
braket(arrow(r),hat(G), arrow(r')) = - (2m)/((2 pi)^3 hbar^2) integral  dd(""^3 k')   1 / ((k' + k)(k' - k)) e^(i arrow(k') dot (arrow(r) - arrow(r')))
$
其中定义
$
E = (hbar^2 k^2) /(2 m)
$
$k$为入射粒子的波矢。因为被积函数在$k = k'$处有奇点，利用留数定理，可以得到
$
braket(arrow(r),hat(G), arrow(r')) = - (m)/(2 pi^2 hbar^2) 1/abs(arrow(r) - arrow(r') )e^(i k abs(arrow(r) - arrow(r')))
$
在这里还要做一个近似：因为观察点离开散射中心的距离$r$通常远大于散射势场本身的尺度$r'$（即$r >> r'$），所以
$
1/abs(arrow(r) - arrow(r')) approx 1/r \
abs(arrow(r) - arrow(r')) = sqrt(r^2 - 2 r arrow(e)_r dot arrow(r') + r'^2) approx r sqrt(1 - 2 arrow(e)_r dot arrow(r')/r ) approx r (1 - arrow(e)_r dot arrow(r')/r) = r - arrow(e)_r dot arrow(r')
$
于是
$
braket(arrow(r),hat(G), arrow(r')) = - (m)/(2 pi^2 hbar^2) e^(i k e)/r e^(-i k arrow(e)_r dot arrow(r')) = - (m)/(2 pi^2 hbar^2)e^(i k e)/r e^(- i arrow(k') dot arrow(r')), arrow(k') = k arrow(e)_r
$
另外：
$
braket(arrow(r'),V(hat(arrow(r))), psi_0) = V(r') braket(arrow(r'), psi_0) = V(r') e^(i arrow(k) dot arrow(r'))
$
于是$braket(arrow(r),psi)  = braket(arrow(r),psi_0) + integral  dd(""^3 r') braket(arrow(r),hat(G), arrow(r')) braket(arrow(r'),V, psi_0)$就是：
$
psi(arrow(r)) = e^(i arrow(k) dot arrow(r)) - (m)/(2 pi^2 hbar^2)  e^(i k e)/r integral  dd(""^3 r') e^(- i arrow(k') dot arrow(r')) V(arrow(r')) e^(i arrow(k) dot arrow(r'))
$
对比可知
$
f(theta, phi) = - m/(2 pi^2 hbar^2) integral  dd(""^3 r) e^(- i arrow(k') dot arrow(r)) V(arrow(r)) e^(i arrow(k) dot arrow(r))
$
物理含义：微分散射振幅$f$正比于粒子在$V$的作用下从波矢为$arrow(k)$的平面波散射为波矢为$arrow(k')$的平面波的*跃迁矩阵元*（$abs(arrow(k)) = abs(arrow(k'))$）。

#figure(
  image("pic/2024-06-19-16-17-41.png", width: 20%),
  numbering: none,
)

如果$V$为中心势场：$V(arrow(r)) = V(r)$，则积分可以简化$arrow(q) = arrow(k')- arrow(k)$
$
integral e^( - i arrow(q) dot arrow(r)) V(r) dd(""^3x) &= integral e^(-i q r cos theta) V(r) r^2 dd(r) sin theta dd(theta) dd(phi)\
&= integral 2 pi (2 sin (q r) )/(q r) V(r) r^2 dd(r)\
&= (4 pi)/q integral_0^oo r V(r) sin(q r) dd(r)
$
于是
$
sigma(theta ) = (4 m^2)/(hbar^4 q^2) abs(integral_0^oo r V(r) sin(q r) dd(r))^2
$

== 散射问题定态解渐进形式的讨论

$
psi(arrow(r)) = e^(i arrow(k) dot arrow(r)) + f(theta, phi) e^(i arrow(k) dot arrow(r))/r
$

渐近形式解的几点讨论：

1. 只在$r >> r'$时有意义。
2. 此形式不是自由粒子波函数的精确解，只是近似解。
3. 如果选定z轴正向沿$arrow(k)$的方向，则一般说来$f$和$phi$无关(见后面讨论)。
4. 粒子概率流密度守恒对$f(theta, phi)$的形式作出了限制：光学定理。作为光学定理的推论，在$θ=0$的方向入射和散射这两部分的波有相消干涉。

== 卢瑟福散射的波恩近似

卢瑟福散射即α离子轰击中性原子，原子核受内层电子遮挡产生屏蔽势
$
V(r) = (Z Z' e_s^2)/r e^(- r/a), e_s = e/sqrt(4 pi epsilon_0)
$
其中$a$为原子半径，$Z'e$为入射粒子电量，代入中心力场情况下微分截面的波恩近似公式
$
sigma(theta) &= (4 m^2)/(hbar^4 q^2) abs(integral_0^oo r V(r) sin(q r) dd(r))^2\
&= (4 m^2 Z^2 Z'^2 e_s^4)/(hbar^4 q^2) abs(integral_0^oo e^(- r/a) sin(q r) dd(r))^2\
&= (4 m^2 Z^2 Z'^2 e_s^4)/(hbar^4 (q^2 + 1/a^2)^2)\
$

如果入射粒子能量很高，其德布罗意波长远小于散射势场半径(这里即原子半径$a$)，同时散射角$θ$不是特别小的情况下，我们有
$
q a =  2 k a sin theta/2 >> 1
$
这时散射截面就是
$
sigma(theta) approx (4 m^2 Z^2 Z'^2 e_s^4)/(hbar^4 q^2) = (Z^2 Z'^2 e_s^4)/(4 m^2 v^4 sin^4 theta/2)
$
这就是卢瑟福微分散射截面公式，由卢瑟福不考虑屏蔽作用(公式中不出现$a$)的情况下用经典力学方法计算库仑势散射得出。在高能情况下粒子的粒子性较强，波动性较弱，把其当作经典粒子进行处理是适合的(公式中不出现$hbar$)

== 散射问题中的角动量守恒

中心力场问题（束缚或非束缚）的定态薛定鄂方程为
$
(- hbar^2/(2 mu r^2) partial/(partial r) (r^2 partial/(partial r)) +V(r) + hat(L)^2/(2 mu r^2)) psi = E psi
$
在势场$V$与时间和角度(中心力场)无关的情况下，可分离变量求解，其中径向波函数满足
$
(1/r^2 dd("")/dd(r) (r^2 dd("")/dd(r)) + (2 mu)/(hbar^2) (E- V(r)  - (l(l+1)hbar^2)/(2 mu r^2)) ) R(r) = 0
$
同时，其角度部分波函数解即是球谐函数$Y_(l m) (theta, phi)$ ，而且角动量算符和哈密顿对易：
$
[hat(L)^2, hat(H)] = 0, [hat(L)_z, hat(H)] = 0
$
这也就意味着*角动量在散射过程中是守恒量*。如果把散射看成含时微扰，比如初态是$l=1$态，那么散射过后末态也必为$l=1$态。

设入射粒子沿$z$轴正向，于是其角动量$arrow(L) = arrow(r) crossproduct arrow(p)$没有$z$方向分量，只有$x、y$方向的分量，也即意味着球谐函数的磁量子数$m=0$，于是被散射粒子的角度分布不显含$φ$（但可以是不同$l$的球谐函数的叠加）。所以对于中心力场，如果设入射粒子沿z轴正向，则散射问题解的渐近形式为
$
lim_(r -> oo) psi_k (r, theta) = e^(i k z) + f(theta) e^(i k z)/r
$
其中微分散射振幅不再与$φ$有关系，而只与$θ$有关。

== 全同粒子散射

不考虑全同性时双粒子平面波：
$
psi(arrow(r)_1, arrow(r)_2) = 1/(2 pi)^3 e^(i (arrow(k)_1 dot arrow(r)_1 + arrow(k)_2 dot arrow(r)_2))
$
引入质心系坐标$arrow(R) = (arrow(r)_1 + arrow(r)_2)/2$和相对坐标$arrow(r) = arrow(r)_1 - arrow(r)_2$，总波矢$arrow(K) = arrow(k)_1 + arrow(k)_2$，相对波矢$arrow(k) = arrow(k)_1 - arrow(k)_2$，则
$
psi(arrow(R), arrow(r)) = 1/(2 pi)^3 e^(i (arrow(K) dot arrow(R) + arrow(k) dot arrow(r)))
$
对称化或反对称化之后的波函数：
$
psi_(plus.minus) (arrow(R), arrow(r)) = 1/(2 pi)^3 e^(i arrow(K) dot arrow(R)) 1/sqrt(2) (e^(i arrow(k) dot arrow(r)) plus.minus e^(- i arrow(k) dot arrow(r))) 
$
#newpara()
对于单粒子在固定势场中的散射，散射问题定态方程的解为
$
psi(arrow(r)) = e^(i arrow(k) dot arrow(r)) + f(theta, phi) e^(i arrow(k) dot arrow(r))/r
$
如果是双粒子散射，则考虑二者的质心坐标系，约化质量和相对坐标。这时相对坐标为
$
arrow(r) = arrow(r)_1 - arrow(r)_2
$
但是对于全同双粒子散射，$psi(arrow(r))$需要进行对称或反对称化：
$
psi(arrow(r)) = e^(i arrow(k) dot arrow(r)) plus.minus e^(- i arrow(k) dot arrow(r)) + (f(theta, phi) plus.minus f(pi - theta, phi + pi)) e^(i arrow(k) dot arrow(r))/r
$
如果两个粒子分别沿正负z轴方向入射，则
$
psi(arrow(r)) = e^(i arrow(k) dot arrow(r)) plus.minus e^(- i arrow(k) dot arrow(r)) + (f(theta) plus.minus f(pi - theta)) e^(i arrow(k) dot arrow(r))/r
$
注意这里没有因子$1/sqrt(2)$，因为我们仍旧只关心散射部分相对一个入射波的大小，同时两个粒子的散射都观测，有截面增大的效果。

#figure(
  image("pic/2024-06-19-17-15-55.png", width: 80%),
  numbering: none,
)

这时系统的微分散射截面为
$
sigma(theta) = abs(f(theta) plus.minus f(pi - theta))^2 = abs(f(theta))^2 + abs(f(pi - theta))^2 plus.minus 2 Re(f^*(theta) f(pi - theta))
$
其中最后一项为干涉项，是基于全同性原理的量子力学效应。

如果是非全同粒子，则在$θ$方向观测到两种粒子*任意一个*的总散射截面，应该是这样的非相干叠加：
$
sigma(theta) =abs(f(theta))^2 + abs(f(pi - theta))^2
$

#figure(
  image("pic/2024-06-19-17-19-11.png", width: 80%),
  numbering: none,
)

设电子空间散射不翻转自旋，两个自旋态为$ket(↑↑)$的电子的微分散射截面是
$
abs(f(theta) - f(pi - theta))^2
$
这是因为Feimi子是交换反对称的，但是自旋部分交换对称，所以空间部分是反对称的，因此散射截面也满足反对称性。两个自旋态为$ket(00)$的电子的微分散射截面是
$
abs(f(theta) + f(pi - theta))^2
$
这是因为$ket(00) = 1/sqrt(2) (ket(↑↓) - ket(↓↑))$，自旋部分交换反对称，所以空间部分是对称的，因此散射截面也满足对称性。两个自旋态为$ket(↑↓)$的电子的微分散射截面是
$
abs(f(theta))^2 + abs(f(pi - theta))^2
$
这是因为$ket(↑↓)$，没有交换对称性，不是全同粒子，所以散射截面是两个单独的散射截面之和。

- 如果粒子波函数里还包括自旋分量、偏振这些分立的指标，在做全同粒子波函数对称化或反对称化操作时，这些指标要随粒子坐标一起进行交换。
- 如果在一个物理过程中两个全同粒子有一个不相同的量子数（比如$S_z$），同时此量子数在过程中守恒，则原则上可以用这个量子数来区分这两个粒子，它们不再是全同的。
- 在上述情况下，是否对此量子数做对称化或反对称化的操作，物理结果是相同的（即不存在全同粒子交换效应）

_例：两个处于下列叠加态的电子散射，求微分散射截面。_
$
psi = 1/sqrt(2) (ket(1","0) + ket(0","0))
$
_方法一：根据角动量守恒_，自旋$ket(1","0)$态散射后仍为$ket(1","0)$态，自旋$ket(0","0)$态散射后仍为$ket(0","0)$态，所以散射截面可以分别计算后相加

- $ket(1","0)$态自旋对称，所以空间部分反对称，这部分微分散射截面为
$
1/2 abs(f(theta) - f(pi - theta))^2
$
- $ket(0","0)$态自旋对称，所以空间部分对称，这部分微分散射截面为
$
1/2 abs(f(theta) + f(pi - theta))^2
$
于是总微分散射截面为
$
sigma(theta) = 1/2 abs(f(theta) - f(pi - theta))^2 + 1/2 abs(f(theta) + f(pi - theta))^2 = abs(f(theta))^2 + abs(f(pi - theta))^2
$
_方法二：从耦合表象换为非耦合表象_
$
psi = 1/sqrt(2) (ket(1","0) + ket(0","0)) =  ket(arrow.t arrow.b)
$
也就是说两个电子处于不同的自旋本征态上，是非全同粒子，在$θ$方向上观察到任意电子的总截面应是两个分截面的非相干叠加：
$
sigma(theta) = abs(f(theta))^2 + abs(f(pi - theta))^2
$
#newpara()

_例：求两个总自旋为1的全同粒子散射的非极化微分散射截面_

极化的意思：一个粒子的总自旋量子数$>0$，则在其任意自旋态下，总能在空间中找到一个方向，该粒子在此方向上自旋投影的平均值为最大值且非$0$。此量子态是一个*纯态*。

非极化的意思：非极化是这样一种量子态，粒子自旋在空间任意方向投影的平均值为$0$。如果粒子总自旋非$0$，那么在纯态量子空间是找不到这种态的，而只能在*混合态*中找。比如给定一组全同粒子，一半自旋向上，一半自旋向下。这种几率混合不同于量子力学中的几率振幅的叠加，而就是一种纯粹*统计上*的纯量子态的混合，混合结果是该组粒子平均自旋为$0$。

两个粒子的总自旋$S$为$0、1、2$（单位为$hbar$），分别是一、三、五重态。根据CG系数的公式：
$
ket(j m j_1 j_2) = sum_(m_1+m_2=m) C(j m ";"j_1 j_2 m_1 m_2) ket(j_1 m_1) ket(j_2 m_2), (j_1 = j_2 =1)\
C(j m";"j_1 j_2 m_1 m_2) = braket(j_1 m_1 j_2 m_2, j m) = (-1)^(j_1 - j_2 + j) braket(j_2 m_2 j_1 m_1, j m)
$
即交换两个粒子后，自旋部分波函数的符号变为
$
(-1)^(j_1 - j_2 + j) = (-1)^(2 - j)
$
对$S=1$三重态来说此符号为负，其它态为正。

所以三重态自旋波函数反对称，所以空间部分也应该反对称（玻色子总体对称），相应的微分散射截面为
$
abs(f(theta) - f(pi - theta))^2
$
$S=0$和$S=2$的自旋态自旋波函数对称，所以空间部分也应该对称，相应的微分散射截面为
$
abs(f(theta) + f(pi - theta))^2
$
因为总自旋非极化，所以是混合态，粒子在这$9$个纯态上的统计概率都是$1/9$，所以最后的总微分截面为
$
sigma(theta) &= 3/9 abs(f(theta) - f(pi - theta))^2 + 6/9 abs(f(theta) + f(pi - theta))^2 \
&= abs(f(theta))^2 + abs(f(pi - theta))^2 + 2/3 Re(f^*(theta) f(pi - theta))
$

#newpara()

_例：两个电子散射，求非极化的微分散射截面_

#figure(
  image("pic/2024-06-19-23-57-24.png", width: 80%),
  numbering: none,
)

= 含时微扰

== 量子态跃迁

无微扰时系统从初态$ket(phi_i)$经时间$t$后跃迁到末态$ket(phi_f)$的*概率幅*为
$
A_(f i) = braket(phi_f, e^(- i/hbar hat(H) t), phi_i) = e^(- i/hbar E_i t)braket(phi_f, phi_i) = e^(- i/hbar E_i t) delta_(f i)
$
如果定义$ket(phi_i)$,$ket(phi_f)$所处表象的算符(如自旋)与$hat(H)$不对易，那么就会出现量子态跃迁的情况(参见自旋例题)。

如果系统有微扰(通常与原$hat(H)$不对易)，那么微扰项将会产生跃
$
hat(H) = hat(H)_0 + hat(H)'
$
其中$hat(H)_0$为原哈密顿算符(如自由粒子)，$hat(H)'$为微扰算符。Shrödinger方程为
$
i hbar partial/(partial t) ket(psi(t)) = (hat(H)_0 + hat(H)') ket(psi)
$
$ket(psi(t))$通常很难解析求出，所以用微扰近似，用无微扰时的波函数来展开波函数的修正项。

非微扰哈密顿算符定态本征值及本征函数为
$
hat(H)_0 ket(phi_n) = E_n ket(phi_n)
$
加上微扰后的薛定鄂方程
$
i hbar partial/(partial t) ket(psi(t)) = (hat(H)_0 + V(arrow(x), t)) ket(psi(t)) , V =H'
$
根据$hat(H)_0$的本征函数的完备性，方程任一解可以展开为：
$
psi = sum_n a_n (t) phi_n e^(- i/hbar E_n t) 
$
代入薛定鄂方程得：
$
i hbar sum dd(a_n)/dd(t) phi_n e^(- i/hbar E_n t) = sum_n a_n (E_n + V) phi_n e^(- i/hbar E_n t)
$
两边乘以末态波函数$phi_f^*$并对空间积分得：
$
i hbar dd(a_f)/dd(t) = sum_n a_n (t) integral dd(""^3 x) phi_f^* V phi_n e^(- i/hbar (E_n - E_f) t)
$
这并非求得了$a_f$的解，因为方程右方求和中仍有$a_f$项。

假设$t=-T/2$的初始时刻系统处于初始态$phi_i$：
$
cases(
  a_i (-T/2) = 1,
  a_n (-T/2) = 0 "for" n != i
)
$
$V$很小时，可用上面这套$a$的初始值代入右边作为*一级近似*得：
$
dd(a_f)/dd(t) = - i/hbar integral dd(""^3 x) phi_f^* V phi_i e^(- i/hbar (E_i - E_f) t)
$
于是在时间间隔$T$内$i→f$的*跃迁振幅*（记为$i T_(f i)$）为：
$
i T_(f i) = a_f (T/2) = 1/(i hbar) integral_(-T/2)^(T/2) dd(t) e^(- i/hbar (E_i - E_f) t)integral dd(""^3 x) phi_f^* V phi_i 
$

== 有限时常微扰

设$V$只在[−T/2, T/2]时间内起作用，且在产生作用的此时间窗口内不随时间变化（有限时常微扰），则有
$
V_(f i) = integral dd(""^3 x) phi_f^* V phi_i, i T_(f i) = V_(f i) /(i hbar) integral_(-T/2)^(T/2) dd(t) e^(- i/hbar (E_i - E_f) t)\
i T_(f i) = - 2 pi i (sin ((Delta E_(f i) T)/(2 hbar)))/(pi Delta E_(f i)) V_(f i)\
$
$T$足够大时
$
i T_(f i) = - 2 pi i delta(E_f - E_i) V_(f i)
$
跃迁速率（单位时间内从$i$跃迁到$f$的几率）：
$
W_(f i) = abs(T_(f i))^2/T = 4 abs(V_(f i))^2/T abs((sin ((Delta E_(f i) T)/(2 hbar)))/(Delta E_(f i)))^2
$
如果$T$足够大，可以利用渐进公式
$
lim_(t -> oo) (sin^2 (x t))/x^2 = pi t delta(x)
$
我们有
$
W_(f i) = abs(T_(f i))^2/T = 4 abs(V_(f i))^2/T pi (T/(2)) delta(E_f - E_i) = (2 pi)/hbar  abs(V_(f i))^2 delta(E_f - E_i)
$
实际上公式中的$E_f$还应包括所有可能的简并末态，所以
$
W_(f i) = sum_k (2 pi)/hbar  abs(V_(f i))^2 delta(E_(f k) - E_k)
$
其中$k$表征能量为$E_f$的所有简并末态。如果$E_f$为连续谱，则应对$E_f$求积分
$
W_(f i) = integral (2 pi)/hbar  abs(V_(f i))^2 delta(E_f - E_i) rho(E_f) dd(E_f)
$
其中$ρ(E_f)$为在$E_f→E_f+dd(E_f)$ 能量区间内的末态态密度（简并态密度，即单位能量间隔内的简并度），最后跃迁速率可表示为
$
W_(f i) = (2 pi)/hbar  abs(V_(f i))^2 rho(E_i)
$
这就是*费米黄金定则*：跃迁速率——单位时间内从初态$ket(phi_i)$跃迁到末态$ket(phi_f)$的几率。

费米黄金定则的物理意义：
- $V_n$：跃迁矩阵元(matrix element)，由微扰势能函数决定
- $rho(E)$：末态态密度，或末态相空间(phase space)
系统反应(或衰变)速率由*矩阵元和相空间*共同决定

_例：求箱归一化条件下的自由系统态密度$ρ$_

系统波函数：
$
psi = L^(-3/2) e^(i/hbar arrow(p) dot arrow(r))
$
动量本征值：
$
p_(x y z) = (2 pi n_(x y z) hbar)/L
$
$
rho(E) dd(E) = (4 pi p^2 dd(p))/((2 pi hbar)/L)^3 = (L/(2 pi hbar))^3 4 pi m sqrt(2 m E) dd(E)
$
如果只限于立体角$Omega(theta , phi)$附近，则
$
rho(E, Omega) dd(E) dd(Omega) = (L/(2 pi hbar))^3 m sqrt(2 m E) dd(E) dd(Omega)
$
有限时常微扰的实例——散射问题的处理。

== 散射问题的含时微扰处理

想象粒子沿$z$轴入射(波矢$arrow(k)$)，在原点处被固定势能散射，粒子态跃迁到$arrow(k)'$态，粒子能量不变(弹性散射)。

可用波包代表粒子，固定势也有作用范围(如半径$a$)，在波包和固定势场没有重合时，粒子处于自由运动状态，在两者重合时散射势$hat(H)'$产生作用，散射后势能又不起作用，这正是在$[-T/2, T/2]$时间内的常微扰问题，下面用*有限时常微扰方法*分析。注意：在波包尺度$→∞$时，相当于微扰时间$T→∞$，结果将和定态微扰一致。

初态入射粒子平面波（箱归一化）：
$
phi_i =L^(-3/2) e^(i/hbar arrow(p) dot arrow(r)) = L^(-3/2) e^(i arrow(k) dot arrow(r)) 
$
末态散射粒子平面波（箱归一化）：
$
phi_f = L^(-3/2) e^(i/hbar arrow(p)_f dot arrow(r)) = L^(-3/2) e^(i arrow(k)_f dot arrow(r))
$
于是：
$
V_(f i) &= integral phi^*_f V(arrow(r)) phi_i dd(""^3 x)\
&= L^(-3) integral e^(- i arrow(k)_f dot arrow(r)) V(arrow(r)) e^(i arrow(k) dot arrow(r)) dd(""^3 x)\
&= L^(-3) integral e^(- i arrow(q) dot arrow(r)) V(arrow(r)) dd(""^3 x)
$
其中
$
arrow(q) = arrow(k)_f - arrow(k)
$
跃迁速率：
$
W_(f i)=integral (2 pi)/hbar abs(V_(f i))^2 rho(E_i, Omega) dd(Omega)
$
末态态密度：
$
rho(E, Omega) = (L/(2 pi hbar))^3 m sqrt(2 m E)
$
代入跃迁速率公式中得：
$
W = (L^3 m)/(4 pi^2 hbar^4) integral abs(V_(f i))^2 sqrt(2 m E_i) dd(Omega)
$
于是在立体角$Omega(theta, phi)$附近的微分越迁速率$dd(W) = W(θ, φ) dd(Omega)$为
$
W(θ, φ) dd(Omega)  = (L^3 m)/(4 pi^2 hbar^4) abs(V_(f i))^2 p dd(Omega) = j_s r^2 dd(Omega)
$
于是
$
W(θ, φ) = (L^3 m)/(4 pi^2 hbar^4) abs(V_(f i))^2 p
$
入射粒子概率流密度为：
$
j_("in") = rho v = (1/L^3) (p/m)
$
最后，微分散射截面为：
$
sigma(θ, φ) &= W(θ, φ)/j_("in") = ((L^6 M^2)/(4 pi^2 hbar^4)) abs(V_(f i))^2\
&= (L^6 M^2)/(4 pi^2 hbar^4) abs(L^(-3) integral e^(- i arrow(q) dot arrow(r)) V(arrow(r)) dd(""^3 x))^2\
&= m^2/(4 pi^2 hbar^4) abs(integral e^(- i arrow(q) dot arrow(r)) V(arrow(r)) dd(""^3 x))^2
$
这和按定态微扰处理的波恩近似结果一致，且$L$被消掉了。

== 有限时周期微扰

有限时周期微扰，即微扰的加入是在一个有限的时间段内，但是在这段时间内又是呈周期性变化(如交变电磁场)，这时微扰引起的跃迁通常会伴随着系统与外界的能量、质量等交换，因为这时系统是非孤立系统。

显然，这时问题的处理只能用含时微扰的方法，而不能用定态微扰法，一般来说，我们引入的微扰具有下面的形式
$
hat(H) = hat(H)_0 + hat(H)'(t), hat(H)'(t) = hat(F) sin (omega t)
$
即$H'$有简谐微扰的形式，简谐震动角频率为$ω$，一般简谐振幅$F$与时间无关。

跃迁振幅公式为（注意时间积分改为从0开始了）
$
i T_(m k) = a_m (t) = 1/(i hbar) integral_0^t dd(t) e^(i/hbar (E_m - E_k) t) integral dd(""^3 x) phi_m^* hat(H)'(t) phi_k
$
矩阵元$H'_(m k)$也有类似时间依赖关系：
$
H'_(m k) (t) = F_(m k) sin (omega t), F_(m k) = integral  phi_m^* hat(F) phi_k dd(tau)
$
$
a_(k->m) &= F_(m k) 1/(i hbar) integral_0^t  sin(omega t') e^(i omega_(m k) t') dd(t')\
&= F_(m k) -1/(2 hbar) integral_0^t  (e^(i omega t') - e^(- i omega t')) e^(i omega_(m k) t') dd(t')\
&= - F_(m k) 1/(2 i hbar)( (e^(i(omega+omega_(m k))t)-1)/(omega + omega_(m k)) + (e^(-i(omega-omega_(m k))t)-1)/(omega - omega_(m k)))
$
其中
$
omega_(m k) = (E_m - E_k)/hbar
$
跃迁几率成为
$
P_(k->m) (t) &= abs(F_(m k))^2/(4 hbar^2) abs( (e^(i(omega+omega_(m k))t)-1)/(omega + omega_(m k)) + (e^(-i(omega-omega_(m k))t)-1)/(omega - omega_(m k)))^2\
&= abs(F_(m k))^2/(2 hbar^2)((1 - cos((omega + omega_(m k))t))/(omega + omega_(m k))^2 + (1 - cos((omega - omega_(m k))t))/(omega - omega_(m k))^2 + (2 cos (omega t) (cos(omega t) - cos(omega_(m k) t)))/(omega^2 - omega_(m k)^2))
$
只考虑时间间隔$t$足够长的情形，利用渐近公式
$
lim_(t -> oo) (1 - cos(x t))/x^2 = pi t delta(x)
$
$
P_(k->m) (t) ->  abs(F_(m k))^2/(2 hbar^2)pi t (delta(omega + omega_(m k)) + delta(omega - omega_(m k)))
$
简谐扰动引起的跃迁的若干重要特征：
- 跃迁几率包含两个$δ$函数，只在
  $
  omega_(m k) = ± omega
  $
  时跃迁几率才显著地不为0，而其它的跃迁几率都可以忽略不计。这种情况称为共振跃迁:
  $
  E_m - E_k = ± hbar omega
  $
  在$E_m=E_k+hbar ω$时称为共振吸收，在$E_m=E_k-hbar ω$时称为共振发射。

  原子或分子对光的共振吸收形成它的特征暗线光谱，而共振发射形成特征明线光谱。“核磁共振”也是共振跃迁的重要例子。
- 在发生共振跃迁的时候，跃迁几率与时间成正比，所以它的跃迁速率是常数：
  $
  W_(k->m) = (P_(k->m) (t))/t = abs(F_(m k))^2/(2 hbar^2)pi (delta(omega + omega_(m k)) + delta(omega - omega_(m k)))
  $
  严格来讲，等式右边的$δ$函数只在$t→∞$才成立，但只要$t$足够大，$δ$函数就已经是很好的近似了。$t$足够大的判据是$t ≫ 1/omega_min$，其中$omega_min$是系统最小的$omega_(m k)$。$1/omega_min$被称为系统的*特征时间*。
- 利用F的厄密性可以证明
  $
  abs(F_(m k)) = abs(F_(k m))
  $
  于是
  $
  W_(k->m) = W_(m->k)
  $
  这称为*细致平衡原理*，在统计力学里有重要的应用。

== 跃迁选择定则

在跃迁几率的表达式中包含有矩阵元
$
H'_(m k) (t) = integral phi_m^* hat(F) phi_k dd(tau)
$
对于某些$hat(H)'(t), psi_m,psi_n$，矩阵元$H'_(m k)$可能为0，这时跃迁几率也为0，这就是*选择定则*。

在$hat(H)'_(m k) != 0$的时跃迁是允许的，在$hat(H)'_(m k) = 0$的时跃迁是禁止的。

选择定则的存在通常是由于某些守恒定律，如动量守恒、能量守恒、角动量守恒、电荷守恒、宇称守恒等等。

== 原子发射和吸收光子

真空中传播的电磁波的能量密度为
$
E = 1/2 (epsilon_0 E^2 + B^2/mu_0)
$
其中电场能与磁场能各占一半，即平均来说
$
epsilon_0 macron(E^2) = 1/mu_0 macron(B^2) => macron(B)/macron(E) = sqrt(mu_0/epsilon_0) = c
$
结合洛仑兹力和库仑力公式
$
arrow(F) = q (arrow(E) + arrow(v) crossproduct arrow(B))
$
粒子所受这两种力的比值为
$
(v macron(B))/macron(E) = v/c << 1
$
所以说对低速带电粒子来说，其所感受到的电磁波中的洛仑兹力远小于库伦力，洛仑兹力可以忽略不计。

在目前的理论框架下，我们仍然把光看作是经典的电磁波，也就是波动的电磁场，由于电场远比磁场的作用显著，所以微扰近似为：
$
arrow(E)(arrow(r), t) = arrow(E)_0 sin (omega t - arrow(k) dot arrow(r))
$
如果光的波长≫一个原子的尺度，可以认为在一个原子的尺度内电场是*空间均匀*地随时间而振荡的，这样的话
$
arrow(E) = arrow(E)_0 sin (omega t)
$
这种近似称为*长波近似*，它引起的跃迁称为*电偶极跃迁*。
$
hat(H)' (arrow(r), t) = e arrow(r) dot arrow(E)(arrow(r), t) = hat(F) sin (omega t), hat(F) = e arrow(r) dot arrow(E)_0
$
相应的跃迁矩阵元是
$
F_(m k) = e integral phi_m^* arrow(r) dot arrow(E)_0 phi_k dd(""^3 x)
$
其中的波函数是
$
psi_(n l m) = R_(n l) (r) Y_(l m) (theta, phi)
$
随着$arrow(E)_0$方向的不同，$arrow(r) dot arrow(E)_0$有不同的表达式，但可能的分量有三个：$E_0 x, E_0 y, E_0 z$，所现在先假设电场沿$x$方向，看系统吸收光子，那么
$
W_(k->m) = (pi e^2 E_0^2)/(2 hbar^2) abs(x_(m k))^2 delta(omega - omega_(m k))
$
电磁波的平均能量密度为
$
I = 1/2 expval(epsilon_0 E^2 + B^2/mu_0)
$
其中左右尖括号表示在一个时间周期内求平均，再利用
$
expval(B^2/mu_0) = expval(epsilon_0 E^2) = epsilon_0 E_0^2 omega/(2 pi) integral_0^(2 pi/omega) sin^2 (omega t) dd(t) = 1/2 epsilon_0 E_0^2
$
所以
$
I = 1/2 epsilon_0 E_0^2
$
$
W_(k->m) = (pi e^2)/(epsilon_0 hbar^2) I abs(x_(m k))^2 delta(omega - omega_(m k)) = (4 pi^2 e_s^2)/hbar^2 I abs(x_(m k))^2 delta(omega - omega_(m k))
$
上面讨论的是入射光为单色偏振光的情形，一般来说，光能量密度随频率有一个分布，设角频率在$ω→ω+dd(ω)$之间的入射光（在单位角频率内的）能量密度为$I(ω)$，则此区间的能量为
$
I(omega) dd(omega)
$
于是
$
W_(k->m) = (4 pi^2 e_s^2)/hbar^2 abs(x_(m k))^2 integral I(omega) dd(omega) delta(omega - omega_(m k)) = (4 pi^2 e_s^2)/hbar^2 I(omega_(m k)) abs(x_(m k))^2
$
目前为止，我们一直假设入射光是沿$x$方向偏振的。实际上入射光波矢沿各个方向都有，偏振也是各个方向都有，均匀分布。因此，我们要把$y, z$方向的偏振贡献也加进来，然后除以$3$求平均
$
W_(k->m) = (4 pi^2 e_s^2)/(3 hbar^2) I(omega_(m k)) (abs(x_(m k))^2 + abs(y_(m k))^2 + abs(z_(m k))^2) = (4 pi^2 e_s^2)/(3 hbar^2) I(omega_(m k)) abs(arrow(r)_(m k))^2
$
公式中出现了电子电偶极矩$-e arrow(r)$的矩阵元，所以这种跃迁又称为*电偶极跃迁*。这个跃迁速率的表达式与费米黄金定则很相似，但两者存在根本不同。这里原子初末态能量不守恒，相差$hbar omega_(m k)$，而费米黄金定则适用于初末态能量守恒的过程。

系统能否发生跃迁由矩阵元$arrow(r)_(m k)$决定，其角度部分的积分牵涉到下列三角函数在球坐标系中的矩阵元

#figure(
  image("pic/2024-06-20-11-28-18.png", width: 80%),
  numbering: none,
)

再利用球谐函数的正交性，看出只有在
$
l' - l = plus.minus 1 , m' - m = 0, plus.minus 1
$
的时候，矩阵元才不为0，所以得到结论：电偶极跃迁的选择定则是
$
Delta l = plus.minus 1, Delta m = 0, plus.minus 1
$
从物理的角度来看，这是由于角动量守恒，因为光子的总自旋量子数是1。当然，在其它的过程中还会有类似的选择定则。

== 正常塞曼效应再探讨 —— 自旋的影响

#figure(
  image("pic/2024-06-20-11-30-41.png", width: 80%),
  numbering: none,
)

#figure(
  image("pic/2024-06-20-11-30-58.png", width: 80%),
  numbering: none,
)

#figure(
  image("pic/2024-06-20-11-41-59.png", width: 80%),
  numbering: none,
)

#figure(
  image("pic/2024-06-20-11-46-43.png", width: 80%),
  numbering: none,
)

== 自发辐射的爱因斯坦理论

在非相对论量子力学框架内无法解释原子的自发辐射问题，而只能近似处理受激辐射和吸收的问题。但爱因斯坦提出了一个半唯象理论把它们之间的关系找了出来。

设能级$E_k$小于能级$E_m$，定义三个系数：
- $B_(m k)$：受激发射系数
- $B_(k m)$：吸收系数
- $A_(m k)$：自发辐射系数

设在热平衡条件下，处于这两个能级的原子数分别为$N_k, N_m$，那么
$
N_k B_(k m) I(omega_(m k)) = N_m (A_(m k) + B_(m k) I(omega_(m k)))
$
即单位时间内从$k$跃迁至$m$（吸收）的原子数，和从$m$跃迁至$k$（自发或受激）的原子数相等。据此可解出
$
I(omega_(m k)) = A_(m k)/(N_k/N_m B_(k m) - B_(m k))
$
又根据麦克斯韦-波尔兹曼分布知
$
N_k = C e^(-E_k/(k T)), N_m = C e^(-E_m/(k T))
$
所以
$
N_k / N_m = e^((E_m - E_k)/(k T)) = e^((hbar omega_(m k))/(k T))
$
代入$I$的表达式得
$
I(omega_(m k)) = A_(m k)/(e^((hbar omega_(m k))/(k T)) B_(k m) - B_(m k))
$
或者说
$
I(nu_(m k)) = (2pi A_(m k))/(e^((hbar nu_(m k))/(k T)) B_(k m) - B_(m k))
$
又根据黑体辐射的公式
$
I(nu) = (8 pi h nu^3)/(c^3 (e^((h nu)/(k T)) - 1))
$
跟前面$I$的表达式对照得出：
$
A_(m k) = (4 h nu_(m k)^3)/(c^3) B_(k m) = (hbar omega_(m k)^3)/(pi^2 c^3) B_(k m)\
B_(m k) = B_(k m)
$
即受激辐射和吸收系数相等，而且得出了自发辐射系数与受激辐射系数之间的关系。由这些关系进一步得到
$
A_(m k) /(B_(m k) I(omega_(m k))) = e^((hbar omega_(m k))/(k T)) - 1
$
此即热平衡时自发辐射和受激辐射速率之比。若要二者相等，则在$T=300k$时需要
$
omega_(m k) = (k T)/hbar ln(2) = 2.74 times 10^13 s^(-1)
$
这个角频率对应的波长大约为$69$µm，远大于可见光波长，因此在可见光波段，原子的自发辐射占绝对主导地位。

前面已求出吸收跃迁速率
$
W_(k->m) = (4 pi^2 e_s^2)/(3 hbar^2) I(omega_(m k)) abs(arrow(r)_(m k))^2
$
这应等于单个原子吸收跃迁速率$B_(k m) I(omega_(m k))$，所以
$
B_(k m) = B_(m k) = (4 pi^2 e_s^2)/(3 hbar^2) abs(arrow(r)_(m k))^2
$
$A_(m k)$就代表单个原子*自发辐射跃迁速率*，处于$m$能级原子的平均寿命就可表示为
$
tau_(m k)= 1/A_(m k)
$
如果原子能从m态跃迁到一系列能量更低的k态，则总平均寿命为
$
tau_m = sum_k 1/A_(m k)
$
受激辐射产生的光单色性、相干性好（激光），相比之下自发辐射则是随机的、相干性差。