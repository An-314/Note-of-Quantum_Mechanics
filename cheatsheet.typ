#import "@preview/physica:0.9.2": *
#import "@preview/tablex:0.0.6": tablex, hlinex
#import "@preview/tablem:0.1.0": tablem

#let three-line-table = tablem.with(
  render: (columns: auto, ..args) => {
    tablex(
      columns: columns,
      auto-lines: false,
      align: center + horizon,
      hlinex(y: 0),
      hlinex(y: 1),
      ..args,
      hlinex(),
    )
  }
)


#set page(
  paper: "a4",
  margin: (x: 0.5cm, y: 0.5cm),
)


#show: rest => columns(3,gutter: 6pt, rest)

#let song = ("Linux Libertine", "SimSun")
#let hei = ("Linux Libertine", "SIMHEI")
#let kai = ("Linux Libertine", "KaiTi",)
#let xbsong = "FZXiaoBiaoSong-B05"
#let code = "Consolas"
#let title-font = hei
#let author-font = kai
#let body-font = song
#let heading-font = hei
#let caption-font = kai
#let header-font = kai
#let strong-font = hei
#let emph-font = kai
#let raw-font = code

#set block(spacing: 4pt)
#show strong: set text(font: strong-font)
#show strong: set block(spacing: 6pt)
#show emph: set text(font: emph-font)

#set text(size: 8pt, font: body-font)

#set par(leading: 4pt)

*量子力学基本公设*

#set text(size: 6pt)
+ 微观体系的状态由波函数描述，波函数满足单值、有限、连续条件
+ 波函数的动力学演化满足薛定鄂方程
+ 力学量用*厄密算符*表示，且有组成*完备集的本征函数系*
+ 公任一波函数可以展开为力学量算符本征函数的线性叠加，测得力学量为本征值$lambda_n$的几率为展开式中对应本征函数系数的模方$abs(c_n)^2$
+ 全同粒子相互调换不改变体系状态
#set text(size: 8pt)
  

*表象变换* 

$S$是*幺正变换*，$S^dagger = S^(-1)$

_基底_：$bra(u)^* = bra(u')^* S$，$S_(m n) = braket(u^'_m,u_n)$

即：$ket(u') = S^* ket(u)$, $ket(u) = S^TT ket(u')$

_态矢_：$ket(psi') = S ket(psi)$，$S_(m n) = braket(psi^'_m,psi_n)$

即：$ket(psi') = S ket(psi)$，$ket(psi) = S^dagger ket(psi')$

_算符_：$A' = S A S^dagger$ 

*Dirac符号*

基底$\{ket(n)\}$_正交归一_：$braket(m,n) = delta_(m n)$

_完备性_：$sum_n ketbra(n) = I$

上面中的一项被称作*投影算符*：$P_n = ket(n) bra(n)$，有性质$P_n^2 = P_n$

- 态矢量在具体表象中的表示

基矢量$ket(k)$，态矢量$ket(psi)$用基矢量展开：
$
ket(psi) = sum_k ket(k) braket(k, psi) = sum_k a_k ket(k)
$
其中$a_k = braket(k, psi)$是态矢量在$ket(k)$表象中的分量。

${a_k} = {braket(k, psi)}$是态矢量$ket(psi)$在$ket(k)$表象中的表示
$
mat(
    a_1, a_2, dots;
)^TT=
mat(
    braket(1, psi), braket(2, psi), dots;
)^TT
$

- 算符在具体表象中的表示

$
bra(j) hat(L) ket(psi) = sum_k bra(j) hat(L) ket(k) braket(k, psi) = braket(j, phi)
$
即：$sum_k L_(j k) a_k = b_j$
其中
$
L_(j k) = bra(j) hat(L) ket(k)
$
就是算符$hat(L)$在$F$表象中的*矩阵元*。

算符$hat(L)$的狄拉克符号表示为：
$
hat(L) = sum_(j k) L_(j k) ket(j) bra(k) = sum_(j k) ket(j) bra(j) hat(L) ket(k) bra(k)
$
算符$hat(F)$在其自身$F$表象中的表示为：
$
F_(m n) = bra(m) hat(F) ket(n) =  bra(m) f_n ket(n) = f_n delta_(m n)\
hat(F) = sum_n f_n ket(n) bra(n)
$
其中$f_n$是$hat(F)$在$ket(n)$表象中的本征值。*任何算符在其自身表象中是对角化的*。

*中心力场中的运动*

中心力场中两体问题的定态Schrödinger方程：
#set text(size: 6pt)
$
(-hbar^2/(2 m_1) nabla^2_1 -hbar^2/(2 m_2) nabla^2_2+ U(abs(arrow(r)_1 - arrow(r)_2))) Psi(arrow(r)_1, arrow(r)_2) = E_"tot" Psi(arrow(r)_1, arrow(r)_2)
$
#set text(size: 8pt)
总能量$E_"tot"$可分为整体*平动*动能和*相对运动*能量（相对动能+势能）两部分。坐标转化：
$
cases(
  arrow(R) = (m_1 arrow(r)_1 + m_2 arrow(r)_2)/(m_1 + m_2) "（质心系坐标）"\ 
arrow(r) = arrow(r)_1 - arrow(r)_2 "（相对坐标）"
), arrow(r)_1 , arrow(r)_2 => arrow(R), arrow(r)
$
#set text(size: 5pt)
$
1/m_1 nabla^2_1 + 1/m_2 nabla^2_2 = 1/M nabla^2_R + 1/mu nabla^2_r , mu = (m_1 m_2) /M
$
#set text(size: 8pt)
分离变量带入：
#set text(size: 5pt)
$
- hbar^2/(2 M) nabla^2_R phi(arrow(R)) = E_"cm" phi(arrow(R)) , (- hbar^2/(2 mu) nabla^2_r + U(r)) psi(arrow(r)) = (E_"tot" - E_"cm") psi(arrow(r))
$
#set text(size: 8pt)

第一个方程的解即自由粒子平面波：
$
phi(arrow(R)) = c_1 e^(i/hbar arrow(P) dot arrow(R)) + c_2 e^(-i/hbar arrow(P) dot arrow(R))
$
其中$P = sqrt(2 M E_"cm")$是质心动量。

第二个方程的解是中心势场问题的解：
#set text(size: 7pt)
$
(- hbar^2/(2 mu) nabla^2_r + U(r)) psi(arrow(r)) = E psi(arrow(r))
$
#set text(size: 8pt)
若$U(r -> oo) > E$，这个方程的能量本征值分立。

#set text(size: 6pt)
$
nabla^2_r = 1/r^2 partial/(partial r) (r^2 partial/(partial r)) + 1/(r^2 sin(theta)) partial/(partial theta) (sin(theta) partial/(partial theta)) + 1/(r^2 sin^2(theta)) partial^2/(partial phi^2)\
- hbar^2/(2 mu) nabla^2 = - hbar^2/(2 mu r^2) partial/(partial r) (r^2 partial/(partial r)) + 1/(2 mu  r^2) hat(L)^2
$
#set text(size: 8pt)
其中$hat(L)^2$是角动量平方算符。方程化为
#set text(size: 6pt)
$
(- hbar^2/(2 mu r^2) partial/(partial r) (r^2 partial/(partial r)) + 1/(2 mu  r^2) hat(L)^2 + U(r)) psi(arrow(r)) = E psi(arrow(r))
$
#set text(size: 8pt)


在球坐标系中分离变量：
$
psi(arrow(r)) = R(r) Y(theta, phi)\
hat(L)^2 Y_(l m) (theta, phi) = l (l + 1) hbar^2 Y_(l m) (theta, phi)
$
径向方程：
#set text(size: 6pt)
$
(- hbar^2/(2 mu r^2) dd("")/(dd(r)) (r^2 dd("")/(dd(r))) + l (l + 1) hbar^2/(2 mu r^2) + U(r)) R_l (r) = E R_l (r)
$
#set text(size: 8pt)
设$R(r) = u(r)/r$，则有
#set text(size: 5pt)
$
- hbar^2 / (2 mu) dd(""^2u) /dd(r^2)  + ((l (l + 1) hbar^2) / (2 mu r^2) + U(r)) u = E u
$
有$2l+1$重简并。边界条件：
$
lim_(r -> 0) r R_l (r) = lim_(r ->0) u(r) = 0
$
当$r→0$时，上式的渐进式是：
$
dd(""^2R)/dd(r^2) + 2/r dd(R)/dd(r) - (l(l+1))/r^2 R= 0
$
在正则奇点$r=0$的邻域内，设$R∝ r^s$，指标方程
$
s(s-1) r^(s-2) + 2 s r^(s-1) - l(l+1) r^(s-2) = 0
$
得到
#set text(size: 8pt)
$
u(r) = r R(r) prop r^(l + 1)
$
#parbreak()

*氢原子*
#set text(size: 5pt)
$
U(r) = - 1/(4 pi epsilon_0) (Z e^2) / r
,
dd(""^2 u)/dd(r^2) + ((2 mu) / hbar^2( E + k_1 (Z e^2) / r )- (l (l + 1)) /  (r^2)) u = 0
$
#set text(size: 8pt)
对于束缚态，$E<0$。定义无量纲量
#set text(size: 6pt)
$
rho = alpha r , alpha = sqrt(8 mu abs(E)) / (hbar), beta = (2 mu k_1 Z e^2)/(alpha hbar^2) = (k_1 Z e^2) / hbar sqrt(mu/(2 abs(E)))
$
#set text(size: 8pt)
于是方程化为：
#set text(size: 6pt)
$
dd(""^2 u)/dd(rho^2) + (- 1/4 + beta/rho - (l (l + 1)) / rho^2 ) u = 0
$
渐进分析假设$u(rho) = rho^(l+1) e^(-rho) nu(rho)$，得到超流几何方程：
$
 dd(""^2 nu)/dd(rho^2) + (2(l + 1)/rho - 1) dd(nu)/dd(rho) + (beta - (l + 1))/rho nu = 0
$
有解条件$n_r = beta - l - 1 = 0, 1, 2, ...$
#set text(size: 8pt)
得到$n^2$度简并的能级：
$
E_n = (mu k_1^2 Z^2 e^4)/(2 hbar^2 n^2) = - (k_1 Z e^2)/(2 a/Z n^2) , a = hbar^2/(mu k_1 e^2)
$
对应的波函数是：
#set text(size: 6pt)
$
psi_(n l m) (r, theta, phi) = R_(n l) (r) Y_(l m) (theta, phi)\,R_(n l) (r) = (u_(n l) (r))/r
$
其中$u_(n l) (r)$是缔合Laguerre多项式。
#set text(size: 8pt)

电子在$(r, r+dd(r))$*球壳中的概率*为：
$
|R_(n l) (r)|^2 r^2 dd(r) = |u_(n l) (r)|^2 dd(r) 
$
$|u_(n l) (r)|^2$的极大值位置称为*最可几半径*。

在$(theta , phi)$的*立体角$dd(Omega)$中的概率密度*和$phi$无无关：
$
|Y_(l m) (theta, phi)|^2 dd(Omega) prop | P_l^m (cos theta) |^2 sin theta dd(theta) dd(phi)
$

#set text(size: 6pt)
设氢原子中电子$psi_(n l m)$态电流密度：
$
arrow(j) = - e (i hbar)/(2 mu) (psi grad psi^* - psi^* grad psi)
$
而
$
psi_(n l m) tilde R_(n l) (r) P^m_l (cos theta) e^(i m phi)
$
其中$R,P$是实函数，从而$j_r = j_theta =0$。
#set text(size: 5pt)
$
j_phi &= (i e hbar)/(2 mu) 1/(r sin theta) (psi^* partial/(partial phi) psi - psi partial/(partial phi) psi^*)= - (e hbar m)/(mu r sin theta) |R_(n l) (r) P^m_l (cos theta)|^2
$
#set text(size: 6pt)
电流密度$j_phi$对应的磁矩为：
#set text(size: 4pt)
$
arrow(M) = 1/2 integral arrow(r) crossproduct arrow(j) dd(V),
M_z = 1/2 integral r sin theta j_phi dd(V) = - (e hbar m)/(2 mu) |R_(n l) (r) P^m_l (cos theta)|^2 dd(r) dd(theta) dd(phi) = - (e hbar m)/(2 mu)
$
#set text(size: 6pt)
$
M_z = - mu_B m,mu_B = (e hbar) / (2 mu), g = M_z / L_z = (- m mu_B)/(m hbar) = - e/(2 mu)
$
#set text(size: 8pt)
#parbreak()

*三维各向同性谐振子*

_直角坐标系下_
#set text(size: 5pt)
$
V(r) = 1/2 mu omega^2 r^2 = 1/2 mu omega^2 (x^2 + y^2 + z^2),H = sum_i H_i , H_i = - hbar^2/(2 mu) nabla_i^2 + 1/2 mu omega^2 r_i^2
$
#set text(size: 8pt)
系统波函数分离变量：
$
psi(x, y, z) = psi_(n_x) (x) psi_(n_y) (y) psi_(n_z) (z)
$
其中$psi_n$为*一维谐振子*与量子数$n$的本征函数。
$
E = (n_x + n_y + n_z + 3/2) hbar omega = (N + 3/2) hbar omega
$
#set text(size: 6pt)
Virial定理：$macron(T) = macron(V) = 1/2 macron(E)$；对于给定的$N$，其简并度为
$((N + 1)(N + 2))/2$
#set text(size: 8pt)

_极坐标系下_

#set text(size: 6pt)
$
(dd("")/dd(r) + (2 mu)/hbar^2 (E - 1/2 mu omega^2 r^2) -( l (l + 1) hbar )/(2 mu r^2)) u = 0 
$
$
rho = alpha r, alpha = sqrt(mu omega / hbar)\,lambda = 2 E / (hbar omega)
,
R(rho) = e^(-rho^2/2) rho^l v(rho)
$
代入原方程后，再做变量代换$ξ = ρ^2$得到*合流超几何方程*：
$
dd(""^2v)/dd(ξ^2)  + ((2l + 3)/(2 ξ) - 1) dd(v)/dd(ξ) + (lambda - 2l -3)/(4 eta) v = 0
$
有解条件为$n_r = (lambda - 2l - 3)/4 = 0, 1, 2, ...$
解得
$
lambda = 2N + 3 ,N = 2n_r + l
,
E_N = (N + 3/2) hbar omega
$
在$N$给定以后，$l$可以取值$l = N , N - 2, ..., 0"或"1$
最后系统径向波函数为
$
R(r) = C L_(n_r)^(l+1/2) (rho^2) rho^l e^(-rho^2/2) (rho = sqrt((mu omega) / hbar) r)
$
其中$L_(n_r)^(l+1/2)$是缔合Laguerre多项式。
#set text(size: 8pt)

_表象变换_ $(hat(H), hat(L)^2, hat(L)_z)$与$(hat(H)_x, hat(H)_y, hat(H)_z)$

#set text(size: 7pt)
$
Psi_000 = (2 alpha^(3/2))pi^(1/4) e^(-alpha^2 r^2/2) Y_00 = psi_0 (x) psi_0 (y) psi_0 (z) = Phi_000
$
$
mat(
  Psi_011; Psi_010 ; Psi_(01-1)
) = (sqrt(2) alpha^(5/2)) / pi^(3/4) e^(-alpha^2 r^2/2) mat(
  -1/sqrt(2) r sin theta e^(i phi); r cos theta; -1/sqrt(2) - r sin theta e^(-i phi)
) prop mat(
   -1/sqrt(2)(x + i y); x; 1/sqrt(2)(x - i y)
)
$
$
mat(
  Phi_100; Phi_010; Phi_001
) = (sqrt(2) alpha^(5/2)) / pi^(3/4) e^(-alpha^2 r^2/2) mat(
  x; y; z
)
Phi = S^* Psi, S = mat(
  - 1/sqrt(2), 0 ,1/sqrt(2); - i/sqrt(2), 0, -i/sqrt(2); 0, 1, 0
)
$

#set text(size: 8pt)

*带电粒子在电磁场中的运动*

带电粒子在外电磁场作用下的Hamilton算符：
$
hat(H) &= 1/(2 mu) (hat(p) - q hat(A))^2 + q Phi(arrow(r)) =  1/(2 mu) (- i hbar grad - q arrow(A))^2 + q Phi(arrow(r))
$
规范（gauge）变换不改变系统物理学性质
$
hat(H)' psi' = i hbar partial/(partial t) psi'
cases(
  psi -> psi' = e^(i theta) psi\
hat(A) -> hat(A)' = hat(A) + hbar/q grad theta\
Phi -> Phi' = Phi - hbar/q (partial theta) / (partial t)
)
$

*Zeeman效应(奇)：碱金属原子的能级在强磁场中分裂的现象称为正常Zeeman效应。【磁-轨耦合】*

#set text(size: 6pt)
$
1/(2m) (-hbar^2 nabla^2 psi + i hbar q ((grad dot arrow(A)) psi +2 arrow(A) dot grad psi) + q^2 arrow(A)^2 psi) = (E - q Phi) psi\
$
#set text(size: 8pt)
取*Coulomb规范*：$div arrow(A) = 1/2 div (arrow(B) crossproduct arrow(r)) = 0 $则
#set text(size: 6pt)
$
1/(2m) (-hbar^2 nabla^2 psi + 2i hbar q (arrow(A) dot grad psi) + q^2 arrow(A)^2 psi) = (E - q Phi) psi\
1/(2m) (-hbar^2 nabla^2 psi + i hbar q  (arrow(B) crossproduct arrow(r)) dot grad psi + 1/4q^2  (arrow(B) crossproduct arrow(r))^2 psi) = (E - q Phi) psi\
$
第二项（的变化量）来说，第三项可以忽略不计：
$
(- hbar^2/(2m) nabla^2 - q arrow(B) dot hat(arrow(L)) + q Phi) psi = E psi
$
如果选择$z$轴的方向为$arrow(B)$的方向，电子电荷$q = - e$，所以
$
(- hbar^2/(2m) nabla^2 - (q B)/(2 m) hat(L)_z + q Phi) psi = E psi
$
$
(- hbar^2/(2 mu) nabla^2 + (e B)/(2 mu) hat(L)_z - e Phi) psi = E psi
$
未加磁场$(B=0)$时碱金属原子的能级与波函数
$
E_(n l) , psi_(n l) (r, theta, phi) = R_(n l) (r) Y_(l m) (theta, phi)
$
#set text(size: 8pt)
每一个能级是$(2l+1)$度简并的。那么加上外磁场后(相当于$hat(H)' = hat(H) + (e B)/(2 mu)hat(L)_z$)，本征波函数不变，本征值发生改变，简并将被打破：
$
E_(n l m) = E_(n l) + (e B)/(2 mu)hbar m 
$

*Landau能级*

*对称规范*约定$arrow(A) = 1/2 arrow(B) crossproduct arrow(r)$，取$B$沿z轴方向。

进行规范变换，新的$arrow(A)$仍满足库伦规范，$arrow(B)$不变：
$
arrow(A) -> arrow(A)' = arrow(A) + grad f , f = 1/2 B x y
$
这时磁矢势变为*Landau规范*
$
arrow(A) = mat(-1/2 B y, 1/2 B x, 0)^TT -> arrow(A)' = mat(0, B x, 0)^TT
$
电子运动限制在$x-y$平面内（二维电子气模型）：
$
1/(2 m) (hat(p)_x^2 + (hat(p)_y + e B x)^2) psi = E psi
$
$[hat(p)_y, hat(H)] = 0$，分离：$psi(x, y) = e^(i k_y y) phi(x)$
$
(- hbar^2/(2 m) dd(""^2)/dd(x^2) + 1/2 m omega_c^2 (x + x_0)^2 ) phi(x) = E phi(x)
$

$
omega_c = (e B) / m, x_0 = k_y l_c^2 , l_c = sqrt(hbar /(m omega_c)) = 1/alpha , alpha^2 = (e B) / hbar
$
$ω_c$是回旋角频率，$l_c$是最小回旋半径
$
(m v)/R = e B ,
T = (2 pi R)/v = (2 pi m) / (e B), 
omega_c = (2 pi) / T =  (e B)/m
$
$
2 pi l_c = lambda = h/p,
e B = p/l_c
=> l_c = sqrt(hbar / (e B)) = sqrt(hbar / (m omega_c))
$
解是坐标平移了$x_0$是一维谐振子方程的解:
$
phi(x) = phi_n (x+x_0), 
E_n = (n + 1/2) hbar omega_c "Landau能级"
$

*粒子能量是转动产生的磁矩*与磁场的相互作用能
#set text(size: 6pt)
$
E = (n + 1/2) hbar omega_c = (n + 1/2) (e hbar) / m B = -mu_z B
,
mu_z = - (e hbar) / (2 m) = - mu_B
$
#set text(size: 8pt)
即磁矩方向与磁场方向相反——*Landau抗磁性*。

#set text(size: 6pt)
波函数$e^(i k_y y) phi_n (x+x_0)$，简并度是无穷大的：对于每个能级$E_n$，对应波函数中$k_y$可以任意取值。考虑电子气限于$L_x$宽的长条，有
$
0 < x_0 < L_x => 0 < k_y < L_x alpha^2
$
$y$向周期性边界条件：长条内每$L_y$长度内有一个电子（箱归一化）
$
k_y = (2 pi N)/L_y , N = 0, 1, 2, ...
,
0 < N < (L_x L_y alpha^2)/(2 pi) = (e B L_x L_y) / (h)
$
于是*单位面积内的能级简并度为*
$
g = (e B) / (h)
$

#set text(size: 8pt)

*电子自旋及其描述*

角动量算符的一般定义：$[hat(J)_i, hat(J)_j] = i hbar epsilon_(i j k) hat(J)_k$

$
[hat(J)^2 , hat(J)_i] = 0, hat(J)^2 = hat(J)_x^2 + hat(J)_y^2 + hat(J)_z^2
$
*角动量本征态*是$hat(J)^2$和$hat(J)_z$的共同本征态
$
hat(J)^2 ket(eta"," m) = eta hbar^2 ket(eta","  m),hat(J)_z ket(eta","  m) = m hbar ket(eta","  m)
$
*阶梯算符*：$hat(J)_(plus.minus) = hat(J)_x ± i hat(J)_y$

$
[hat(J)_z , hat(J)_(plus.minus)] = [hat(J)_z, hat(J)_x] ± i [hat(J)_z, hat(J)_y] = ± hbar hat(J)_(plus.minus)
$
#set text(size: 6pt)
$
hat(J)_z hat(J)_(plus.minus) ket(eta"," m) = (m ± 1) hbar hat(J)_(plus.minus) ket(eta"," m)
,
hat(J)_(plus.minus) ket(eta"," m) = c ket(eta'"," m ± 1)
$
#set text(size: 8pt)

$
[hat(J)^2 , hat(J)_(plus.minus)] = [hat(J)^2 , hat(J)_x] ± i [hat(J)^2 , hat(J)_y] = 0
$
#set text(size: 5pt)
$
hat(J)^2 hat(J)_(plus.minus) ket(eta"," m) = hat(J)_(plus.minus) hat(J)^2 ket(eta"," m) = eta hbar^2 hat(J)_(plus.minus) ket(eta"," m)
,
hat(J)^2 c ket(eta'"," m ± 1) = c eta' hbar^2 ket(eta'"," m ± 1) => eta = eta'
$
#set text(size: 8pt)

$m$极大值为$j$，$eta=j(j+1)$，用$ket(j "," m)$表示本征态：
$
hat(J)_(plus.minus) ket(j "," m) = hbar sqrt(j(j+1) - m(m ± 1)) ket(j "," m ± 1)
$

在$(hat(L)^2,hat(L)_z)$表象，基底$ket(11),ket(10),ket(1-1)$，
$
hat(J)_z = hbar mat(
  1, 0, 0; 0, 0, 0; 0, 0, -1
), hat(J)_x = (hat(J)_+ + hat(J)_-)/2 = hbar/sqrt(2) mat(
  0, 1, 0; 1, 0, 1; 0, 1, 0
)
$
久期方程得到$hat(J)_x$本征态与概率幅：
#set text(size: 6pt)
$
1/2 mat(1;sqrt(2);1) , 1/sqrt(2) mat(1;0;-1) , 1/2 mat(1;-sqrt(2);1); 1/4, 1/2, 1/4
$
#set text(size: 8pt)

*电子自旋角动量*
- 电子自旋任何方向投影只取$plus.minus hbar/2$，$hat(S)_i$本征值$± hbar/2$，
  $
  hat(S)_x^2 = hat(S)_y^2 = hat(S)_z^2 = (hbar/2)^2,hat(S)^2 = hat(S)_x^2 + hat(S)_y^2 + hat(S)_z^2 = 3/4 hbar^2
  $
- 自旋角动量导致电子有*自旋磁矩*，其$z$轴投影为
$
mu_z / S_z = - e/m_e , mu_z = minus.plus mu_B , hat(arrow(mu)) = - e /m_e hat(arrow(S))
$
在$(hat(S)^2, hat(S)_z)$表象下，利用升降算符，得Pauli矩阵
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
$
hat(S)_x = hbar/2 sigma_x, hat(S)_y = hbar/2 sigma_y, hat(S)_z = hbar/2 sigma_z\,hat(arrow(S)) = hbar/2 hat(arrow(sigma))
$
#set text(size: 6pt)
*Clifford代数*
$sigma_i sigma_j = delta_(i j) I + i epsilon_(i j k) sigma_k$,
$(arrow(a) dot arrow(sigma)) (arrow(b) dot arrow(sigma)) = arrow(a) dot arrow(b) I + i arrow(a) crossproduct arrow(b) dot arrow(sigma)
$
#set text(size: 8pt)

*二分量波函数*，*二分量旋量*：
$
Psi(arrow(r), t) = Psi_1 (arrow(r), t) v_+ + Psi_2 (arrow(r), t) v_- = mat(Psi_1 (arrow(r), t), Psi_2 (arrow(r), t))^TT
$
$
Psi_1 (arrow(r), t) = Psi(arrow(r) , t, S_z = hbar/2), Psi_2 (arrow(r), t) = Psi(arrow(r) , t, S_z = -hbar/2)
$

粒子轨道和自旋*无耦合*时，直积态：$Psi(arrow(r), t) = psi(arrow(r), t) ket(chi)$，$ket(chi) = c_1 ket(arrow.t) + c_2 ket(arrow.b) = mat(c_1; c_2)$，否则*耦合*。

算符的平均值是：
$
macron(A) = integral Psi^dagger hat(A) Psi dd(arrow(r)) = integral (mat(Psi_1^*, Psi_2^*) hat(A) mat(Psi_1, Psi_2)^TT) dd(arrow(r))
$
算符$hat(A)$和自旋无关则在自旋表象中对角化：$hat(A) = hat(A) I$

*静磁场*只考虑内禀自旋，$hat(H) = - g hat(arrow(S)) dot arrow(B) = - (g hbar B)/2 arrow(sigma) dot arrow(e)_B$
#set text(size: 6pt)
$
arrow(sigma) dot arrow(e)_B = mat(
  cos theta, sin theta e^(-i phi);
  sin theta e^(i phi), - cos theta
), nu_1 =  mat(
  cos theta/2 e^(-i phi/2);
  sin theta/2 e^(i phi/2)
), nu_(-1) = mat(
  sin theta/2 e^(-i phi/2);
  - cos theta/2 e^(i phi/2)
)
$
#set text(size: 8pt)

在静磁场中两个相邻本征能级之差为
$
Delta E = abs(g hbar B) = abs(hbar omega_L),omega_L = - g B "Lamor频率"
$
*时间演化与量子跃迁*：
$
ket(chi(t)) &= e^(- i/hbar t hat(H)) ket(chi(0))
= e^(- i (w_L t)/2 arrow(sigma) dot arrow(e)_B) mat(a_0; b_0) sigma_z"表象"\
&= (cos(w_L/2 t) - i sin(w_L/2 t) arrow(sigma) dot arrow(e)_B ) mat(a_0; b_0)\
$
$e^(- i/hbar t hat(H))$有非0非对角矩阵元，则可能自旋“*跃迁*”。

周期性跃迁*振荡*：$sigma_z$和$hat(H)$不对易，$ket(chi)$非定态。

*角动量的合成* $arrow(J) = arrow(J)_1 + arrow(J)_2$也是角动量算符

$hat(arrow(J))_1$$hat(arrow(J))_2$独立：$[hat(J)_i, hat(J)_j] = i hbar epsilon_(i j k) hat(J)_k$，$[hat(J)_(1 i), hat(J)_(2 j)] = 0$；

*未耦合*$(hat(J)_1^2, hat(J)_2^2, hat(J)_1^2, hat(J)_2^2)$
共同本征态是直积，并矢
$
ket(j_1","m_1","j_2","m_2) = ket(j_1","m_1) ket(j_2","m_2) [(2j_1 + 1)(2j_2 + 1)"维"]
$
*耦合*$(hat(J)^2, hat(J)_z, hat(J)_1^2, hat(J)_2^2)$共同本征态记做$ket(j","m";"j_1","j_2)$

#set text(size: 5.5pt)
用$hat(J)_z = hat(J)_(1z) + hat(J)_(2z)$的本征态展开：$(m - m_1 - m_2) C(j ,m; j_1, m_1, j_2, m_2) = 0$
#set text(size: 6pt)
*Clebsch-Gordan系数*
$
ket(j","m";"j_1","j_2) = sum_(m = m_1 + m_2) C(j ,m; j_1, m_1, j_2, m_2) ket(j_1","m_1) ket(j_2","m_2)
$
#set text(size: 8pt)

*最大投影态*：直积与耦合的本征态相同$ket(j_1 j_1) ket(j_2 j_2) = ket(j = j_1 + j_2","m = j_1 + j_2","j_1","j_2)$
$
hat(J)^2  = hat(arrow(J))_1^2 + hat(arrow(J))_2^2 + hat(J)_(1 +) hat(J)_(2 -) + hat(J)_(1 -) hat(J)_(2 +) + 2 hat(J)_(1 z) hat(J)_(2 z)
$
*三角形法则*：$abs(j_1 - j_2) <= j <= j_1 + j_2,m = m_1 + m_2$
- *电子的旋-轨耦合总角动量*
#set text(size: 6pt)
$
hat(arrow(J)) = hat(arrow(L)) + hat(arrow(S))
,
j_1 = l = 0, 1, 2, ...,
j_2 = s = 1/2
,
j = l+1/2 或 l-1/2
$
#set text(size: 8pt)

*CG系数的特点及符号约定*：CG系数都为*实数*，同时在$m=j,m_1=j_1$时，系数为*非负实数*。
#set text(size: 6pt)
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
#set text(size: 8pt)
L-S耦合表象本征态不是$hat(S)_z$的，$hat(S)_z$与$hat(J)_z, hat(L)^2, hat(S)^2$对易
#set text(size: 5pt)
$
[hat(J)^2, hat(S)_z] &= [hat(L)^2 + hat(S)^2 + hat(L)_+ hat(S)_- + hat(L)_- hat(S)_+ + 2 hat(L)_z hat(S)_z, hat(S)_z] 
 = hat(L)_+ [ hat(S)_-, hat(S)_z] + hat(L)_- [ hat(S)_+, hat(S)_z]\
& =^([hat(J)_plus.minus , hat(J)_z] = minus.plus hbar hat(J)_plus.minus) hbar( hat(L)_+ hat(S)_- - hat(L)_- hat(S)_+)!= 0
$
$
macron(S)_z &= hbar/(2(2l+1)) integral mat(
  sqrt(j + m) Y_(l,m-1/2);
  sqrt(j - m) Y_(l,m+1/2)
)^dagger sigma_z mat(
  sqrt(j + m) Y_(l,m-1/2);
  sqrt(j - m) Y_(l,m+1/2)
) dd(tau)= (m hbar)/(2l+1)
$
#set text(size: 6pt)
非耦合表象下$hat(S)_z$的矩阵形式是$
hbar/2 mat(
  1, 0;
  0, -1
)$，耦合表象下是$
U hbar/2 mat(
  1, 0;
  0, -1
) U^dagger = hbar/(2(2l+1)) mat(
  2m,-sqrt((2l+1)^2 - 4m^2);
  -sqrt((2l+1)^2 - 4m^2), -2m
)
$
#set text(size: 7pt)

$hat(arrow(L)) dot hat(arrow(S))= 1/2(hat(J)^2 - hat(L)^2 - hat(S)^2)$的本征态*是L-S耦合表象的基底*。

#set text(size: 8pt)

- *碱金属原子光谱双线结构*【*旋-轨耦合*】

$
hat(H) = hat(p)^2/(2 mu) + V(r) + xi(r) arrow(L) dot arrow(S)
,
xi(r) = 1/(2 m^2 c^2) dd(V)/dd(r)
$
#set text(size: 6pt)
$
Delta E &= braket(n","j","m_j","l","1/2 , xi(r) arrow(L) dot arrow(S) , n","j","m_j","l","1/2) \
&= braket(n l, xi(r), n l) braket(j","m_j","l","1/2 , arrow(L) dot arrow(S) , j","m_j","l","1/2)\
&= xi_(n l) hbar^2/2 (j(j+1) - l(l+1) - 3/4) = cases(
  1/2 hbar^2 xi_(n l) 当(j= l+1/2) , - (l + 1)/2 hbar^2 xi_(n l)当 (j= l-1/2)
)
$
#set text(size: 8pt)
- *反常Zeeman效应(偶)*【*旋-轨耦合*和*磁-轨耦合*】
#set text(size: 6pt)
$
hat(H) &= hat(p)^2/(2 mu) + V(r) + xi(r) arrow(L) dot arrow(S) + (e B)/(2 mu) (hat(L)_z + 2 hat(S)_z)\
&= hat(p)^2/(2 mu) + V(r) + xi(r)/2 (hat(J)^2 - hat(L)^2 - hat(S)^2) + (e B)/(2 mu) hat(J)_z + (e B)/(2 mu) hat(S)_z
$
#set text(size: 8pt)
无外磁场，旋-轨耦合表象$E_(n l j) , ket(n "," j "," m_j "," l "," 1/2)$

加入$(e B)/(2 mu) hat(J)_z$项，与原哈密顿算符对易，量子态不变：
$
E_(n l j m_j) = E_(n l j) + (e B)/(2 mu) m_j hbar [$(2j+1)$"重简并被完全消除"]
$
$(e B)/(2 mu) hat(S)_z$*微扰*，修正为
#set text(size: 6pt)
$
Delta E = braket(n "," j "," m_j "," l "," 1/2 , (e B)/(2 mu) hat(S)_z , n "," j "," m_j "," l "," 1/2)
= plus.minus (e B hbar)/(2 mu (2l + 1) ) m_j
$
#set text(size: 8pt)
于是最后修正后的能级为
$
E = cases(
  E_(n l j) + (e B)/(2 mu) (1  + 1/(2l + 1)) hbar m_j 当(j = l + 1/2) , E_(n l j) + (e B)/(2 mu) (1  - 1/(2l + 1)) hbar m_j 当(j = l - 1/2)
)
$

- *两个电子自旋的合成*
#set text(size: 6pt)

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


$hat(S)_+ = hat(S)_(1+) + hat(S)_(2+)$在旧表象$mat(
  ket(arrow.t "," arrow.t),
  ket(arrow.b "," arrow.b),
  ket(arrow.t "," arrow.b),
  ket(arrow.b "," arrow.t)
)^TT$下的矩阵和新表象$mat(
  ket(1","1),
  ket(1","-1),
  ket(1","0),
  ket(0","0)
)^TT
$下矩阵为

$
U mat(
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
#set text(size: 8pt)

*定态微扰论*

#set text(size: 6pt)
$
(hat(H)^((0)) + hat(H)^') (psi_n^((0)) + psi_n^((1)) + psi_n^((2)) + ...) = (E_n^((0)) + E_n^((1)) + E_n^((2)) + ...) (psi_n^((0)) + psi_n^((1)) + psi_n^((2)) + ...)\
(hat(H)^((0)) - E_n^((0)) + H' - E_n^((1)) - E_n^((2)) + ...) (psi_n^((0)) + psi_n^((1)) + psi_n^((2)) + ...) = 0
$
#set text(size: 8pt)

*零级方程*就是无微扰时$hat(H)^((0))$的本征方程
$
(hat(H)^((0)) - E_n^((0))) psi_n^((0)) = 0
$
*一级方程*$(hat(H)^((0)) - E_n^((0))) psi_n^((1)) = (E_n^((1)) - hat(H)^') psi_n^((0))$正交性
$
a_(n m )^((1)) (E_m^((0)) - E_n^((0)) )= - braket(psi_m^((0)) , hat(H)^' , psi_n^((0))) + E_n^((1)) delta_(m n)\
E_n^((1)) = braket(psi_n^((0)) , hat(H)^' , psi_n^((0))) = H_(n n)^',psi_n^((1)) = sum_m a_(n m)^((1)) psi_m^((0))
$
$
a_(n m)^((1)) = - braket(psi_m^((0)) , hat(H)^' , psi_n^((0))) / (E_m^((0)) - E_n^((0))),
psi_n^((1)) = sum_(m != n) H_(m n)^'  / (E_n^((0)) - E_m^((0))) psi_m^((0))
$

$
braket(psi_n^((0)), psi_n^((1))) + braket(psi_n^((1)), psi_n^((0))) = 0 => a^((1))_(n n) i  = a_n^((1)) =^"相位因子"_"不妨" 0
$

*二级方程*
#set text(size: 7pt)
$(hat(H)^((0)) - E_n^((0))) psi_n^((2)) = - (hat(H)^' - E_n^((1))) psi_n^((1)) + E_n^((2)) psi_n^((0))$
#set text(size: 6pt)
$
E_n^((2)) = sum_(m != n) H_(m n)^' / (E_n^((0)) - E_m^((0))) braket(psi_n^((0)) , hat(H)' ,psi_m^((0))) = sum_(m != n) (H_(m n)^'  H_(m n)^')/ (E_n^((0)) - E_m^((0))) = sum_(m != n) abs(H_(m n)^')^2 / (E_n^((0)) - E_m^((0)))
$
#set text(size: 8pt)
$
"非简并情形"
cases(
psi_n = psi_n^((0)) + sum_(m != n) H_(m n)^' / (E_n^((0)) - E_m^((0))) psi_m^((0)) + ...\
E_n = E_n^((0)) + H_(n n)^' + sum_(m != n) abs(H_(m n)^')^2 / (E_n^((0)) - E_m^((0))) + ...
)
$

*在静电场中的一维谐振子*
#set text(size: 6.5pt)
$hat(H)^((0)) = hat(p)^2/(2 mu) + 1/2 mu omega^2 hat(x)^2,hat(H)^' = - q E hat(x)$
#set text(size: 7pt)
利用递推关系$hat(x) psi_n^((0)) = sqrt(hbar/(2 mu omega)) (sqrt(n+1) psi_(n+1)^((0)) + sqrt(n) psi_(n-1)^((0)))$

二级修正$E_n = (n + 1/2) hbar omega - (q^2 E^2) / (2 mu omega^2)$
#set text(size: 8pt)
$
"简并情形" (hat(H)^((0)) - E_n^((0))) psi_(n)^((1)) = (E_n^((1)) - hat(H)^') sum_i c_(n i)^((0)) psi_(n i)^((0))\
sum_i c_(n i)^((0)) ( H_(j i)^' - E_n^((1)) delta_(j i) ) = 0 => det(H' - E_n^((1))I) = 0
$
*Stark效应*【静电场】$hat(H)^' = e E z = e E r cos theta$
$H' = braket(psi_(2 l' m')^((0)), hat(H)^', psi_(2 l m)^((0))) $利用球谐递推
#set text(size: 5pt)
$
cos theta Y_(l,m) = sqrt(((l+1)^2 - m^2)/((2l + 1)(2l+3))) Y_(l+1,m) + sqrt((l^2 - m^2)/((2l + 1)(2l-1))) Y_(l-1,m)
$

$
det
  mat(
     - E_2^((1)) , -3e E a_0, 0 , 0 ;
    -3e E a_0 ,  - E_2^((1)) , 0 ,0;
    0 , 0 ,  - E_2^((1)) , 0 ;
    0 , 0, 0 ,  - E_2^((1))
  )
= 0
 => 
E_2^((1)) = 0, plus.minus 3 e E a_0
$

#set text(size: 8pt)

*散射理论*

*散射波函数*$lim_(r->oo) psi_(arrow(r)) = e^(i arrow(k) dot arrow(r)) + f(theta, phi)  e^(i arrow(k) dot arrow(r))/r$
#set text(size: 6pt)
$
arrow(J)_i = rho arrow(v) = abs(psi_i)^2 (hbar arrow(k)) /m = (hbar arrow(k)) /m
,
arrow(J)_s = (i hbar)/(2 m) (psi_s grad psi_s^* - psi_s^* grad psi_s) = (hbar k)/( r^2 m) |f(theta, phi)|^2 arrow(e)_r
$
#set text(size: 8pt)
*散射截面*$dd(sigma) = (J_s r^2)/J_i dd(Omega) = abs(f(theta, phi))^2dd(Omega) = sigma(theta, phi) dd(Omega)$
单位时间内入射波被散射概率$J_i sigma ("散射速率")$
*Lippman-Schwinger方程* $(hat(H)_0 + hat(V)) ket(psi) = E ket(psi)$
#set text(size: 6pt)
$
cases(
  (E - hat(H)_0) ket(psi_0) = 0,
  (E - hat(H)_0) ket(psi) = hat(V) ket(psi)
) =>
ket(psi) = ket(psi_0) + (E - hat(H)_0)^(-1) V ket(psi) = ket(psi_0) +  hat(G) V ket(psi)
$
#set text(size: 8pt)
*Born近似*
#set text(size: 5pt)
$ket(psi) = (1 + hat(G) V) ket(psi_0)$$braket(arrow(r),psi)= braket(arrow(r),psi_0) + integral  dd(""^3 r') braket(arrow(r),hat(G), arrow(r')) braket(arrow(r'),V, psi_0)$
$
braket(arrow(r),hat(G), arrow(r')) = 1/(2 pi)^3 integral  dd(""^3 k') e^(i arrow(k') dot (arrow(r) - arrow(r'))) / (E - (hbar^2 k'^2) /( 2 m)) = - (m)/(2 pi^2 hbar^2) 1/abs(arrow(r) - arrow(r') )e^(i k abs(arrow(r) - arrow(r'))) = - (m)/(2 pi^2 hbar^2)e^(i k e)/r e^(- i arrow(k') dot arrow(r')), arrow(k') = k arrow(e)_r\
braket(arrow(r'),V(hat(arrow(r))), psi_0) = V(r') braket(arrow(r'), psi_0) = V(r') e^(i arrow(k) dot arrow(r')), E= (hbar^2 k^2)/(2 m)\
f(theta, phi) = - m/(2 pi^2 hbar^2) integral  dd(""^3 r) e^(- i arrow(k') dot arrow(r)) V(arrow(r)) e^(i arrow(k) dot arrow(r)),abs(arrow(k)) = abs(arrow(k'))
$
#set text(size: 8pt)
$"中心势场" 
sigma(theta ) = (4 m^2)/(hbar^4 q^2) abs(integral_0^oo r V(r) sin(q r) dd(r))^2
$
*Rutherford散射* $V(r) = (Z Z' e_s^2)/r e^(- r/a), e_s = e/sqrt(4 pi epsilon_0)$
#set text(size: 6pt)
$
sigma(theta) = (4 m^2 Z^2 Z'^2 e_s^4)/(hbar^4 (q^2 + 1/a^2)^2) approx^("忽略"a) (Z^2 Z'^2 e_s^4)/(4 m^2 v^4 sin^4 theta/2)
$
#set text(size: 8pt)
*角动量在散射过程中是守恒量*

*全同粒子散射*#set text(size: 7pt)
$psi(arrow(r)) = e^(i arrow(k) dot arrow(r)) plus.minus e^(- i arrow(k) dot arrow(r)) + (f(theta) plus.minus f(pi - theta)) e^(i arrow(k) dot arrow(r))/r$
#set text(size: 8pt)

#set text(size: 5pt)
#three-line-table[
 | 空间部分交换对称 $ket(00)$ | 空间部分交换反对称 $ket(↑↑)$ | 自旋部分不全同 $ket(↑↓)$ |
  | ---| ---| --- |
  | $abs(f(theta)+f(pi - theta))^2$ | $abs(f(theta)-f(pi - theta))^2$ | $abs(f(theta))^2 + abs(f(pi - theta))^2$ |
]
#set text(size: 7pt)
交换两个粒子后，自旋部分波函数的符号变为$(-1)^(j_1 - j_2 + j)$
#set text(size: 8pt)

*含时微扰*

*概率幅*
#set text(size: 7pt)
$A_(f i) = braket(phi_f, e^(- i/hbar hat(H) t), phi_i) = e^(- i/hbar E_i t)braket(phi_f, phi_i) = e^(- i/hbar E_i t) delta_(f i)$
$
i hbar partial/(partial t) psi = (hat(H)_0 + V(arrow(x), t)) psi , V =H',psi = sum_n a_n (t) phi_n e^(- i/hbar E_n t) 
$
$
i hbar dd(a_f)/dd(t) = sum_n a_n (t) integral dd(""^3 x) phi_f^* V phi_n e^(- i/hbar (E_n - E_f) t)
$

$
cases(
  a_i (-T/2) = 1,
  a_n (-T/2) = 0 "for" n != i
)
=>^"一级近似"

dd(a_f)/dd(t) = - i/hbar integral dd(""^3 x) phi_f^* V phi_i e^(- i/hbar (E_i - E_f) t)
$
$
i T_(f i) = a_f (T/2) = 1/(i hbar) integral_(-T/2)^(T/2) dd(t) e^(- i/hbar (E_i - E_f) t)integral dd(""^3 x) phi_f^* V phi_i "跃迁振幅"
$
#set text(size: 8pt)
*有限时常微扰*
#set text(size: 6pt)
$i T_(f i) = - 2 pi i (sin ((Delta E_(f i) T)/(2 hbar)))/(pi Delta E_(f i)) V_(f i) -> - 2 pi i delta(E_f - E_i) V_(f i)$
#set text(size: 8pt)
$
W_(f i) = abs(T_(f i))^2/T  = (2 pi)/hbar  abs(V_(f i))^2 delta(E_f - E_i) "跃迁速率"
$
密度积分得*Fermi黄金定则*$W_(f i) = (2 pi)/hbar  abs(V_(f i))^2 rho(E_i)$

*箱归一化*自由系统$psi = L^(-3/2) e^(i/hbar arrow(p) dot arrow(r))$,
$p_(x y z) = (2 pi n_(x y z) hbar)/L$
#set text(size: 5pt)
$
rho(E) dd(E) = (4 pi p^2 dd(p))/((2 pi hbar)/L)^3 = (L/(2 pi hbar))^3 4 pi m sqrt(2 m E) dd(E),
rho(E, Omega) dd(E) dd(Omega) = (L/(2 pi hbar))^3 m sqrt(2 m E) dd(E) dd(Omega)
$
#set text(size: 8pt)
*散射含时微扰*$phi_i =L^(-3/2) e^(i/hbar arrow(p) dot arrow(r)) = L^(-3/2) e^(i arrow(k) dot arrow(r)) $
$
V_(f i) &= integral phi^*_f V(arrow(r)) phi_i dd(""^3 x)= L^(-3) integral e^(- i arrow(q) dot arrow(r)) V(arrow(r)) dd(""^3 x)
$
#set text(size: 6pt)
$
W = (L^3 m)/(4 pi^2 hbar^4) integral abs(V_(f i))^2 sqrt(2 m E_i) dd(Omega), W(θ, φ) dd(Omega)  = (L^3 m)/(4 pi^2 hbar^4) abs(V_(f i))^2 p dd(Omega) = j_s r^2 dd(Omega)
$
$
j_("in") = rho v = (1/L^3) (p/m)
,
sigma(θ, φ) &= W(θ, φ)/j_("in")= m^2/(4 pi^2 hbar^4) abs(integral e^(- i arrow(q) dot arrow(r)) V(arrow(r)) dd(""^3 x))^2
$
#set text(size: 8pt)
*有限时周期微扰*$hat(H) = hat(H)_0 + hat(H)'(t), hat(H)'(t) = hat(F) sin (omega t)$
$
H'_(m k) (t) = F_(m k) sin (omega t), F_(m k) = integral  phi_m^* hat(F) phi_k dd(tau)
$
$
a_(k->m) &= F_(m k) 1/(i hbar) integral_0^t  sin(omega t') e^(i omega_(m k) t') dd(t'),
omega_(m k) = (E_m - E_k)/hbar\
&= - F_(m k) 1/(2 i hbar)( (e^(i(omega+omega_(m k))t)-1)/(omega + omega_(m k)) + (e^(-i(omega-omega_(m k))t)-1)/(omega - omega_(m k)))
$
#set text(size: 7pt)
$
P_(k->m) (t) &= abs(F_(m k))^2/(4 hbar^2) abs( (e^(i(omega+omega_(m k))t)-1)/(omega + omega_(m k)) + (e^(-i(omega-omega_(m k))t)-1)/(omega - omega_(m k)))^2\
&->  abs(F_(m k))^2/(2 hbar^2)pi t (delta(omega + omega_(m k)) + delta(omega - omega_(m k)))
$
#set text(size: 5pt)
$E_m=E_k plus.minus hbar ω$共振吸收/发射，$1/omega_min$特征时间， $W_(k->m) = W_(m->k)$*细致平衡原理*
#set text(size: 8pt)

*电偶极跃迁*$hat(F) = e arrow(r) dot arrow(E)_0$，$F_(m k) = e integral phi_m^* arrow(r) dot arrow(E)_0 phi_k dd(""^3 x)$

#set text(size: 7pt)
$
W_(k->m) = (pi e^2 E_0^2)/(2 hbar^2) abs(x_(m k))^2 delta(omega - omega_(m k)) 
= (4 pi^2 e_s^2)/hbar^2 I abs(x_(m k))^2 delta(omega - omega_(m k))
$
#set text(size: 8pt)
$
W_(k->m) = (4 pi^2 e_s^2)/(3 hbar^2) I(omega_(m k)) abs(arrow(r)_(m k))^2,
Delta l = plus.minus 1, Delta m = 0, plus.minus 1
$

*Einstein自发辐射*
#set text(size: 4pt)
$E_k<E_m$，$B_(m k)$受激发射系数，$B_(k m)$吸收系数，$A_(m k)$自发辐射系数
#set text(size: 6pt)
$
N_k B_(k m) I(omega_(m k)) = N_m (A_(m k) + B_(m k) I(omega_(m k))) 
=>
I(omega_(m k)) = A_(m k)/(e^((hbar omega_(m k))/(k T)) B_(k m) - B_(m k))
$
$
I(nu_(m k)) = (2pi A_(m k))/(e^((hbar nu_(m k))/(k T)) B_(k m) - B_(m k))
,
I(nu) = (8 pi h nu^3)/(c^3 (e^((h nu)/(k T)) - 1))
=>
A_(m k)  = (hbar omega_(m k)^3)/(pi^2 c^3) B_(k m)\

W_(k->m) = (4 pi^2 e_s^2)/(3 hbar^2) I(omega_(m k)) abs(arrow(r)_(m k))^2 = B_(k m) I(omega_(m k))
,
B_(k m) = B_(m k) = (4 pi^2 e_s^2)/(3 hbar^2) abs(arrow(r)_(m k))^2
$
#set text(size: 5pt)
$A_(m k)$单个原子*自发辐射跃迁速率*，平均寿命$tau_(m k)= 1/A_(m k)$总平均寿命为$tau_m = sum_k 1/A_(m k)$
#set text(size: 6pt)

- 无限深方势阱$E_n = p^2/(2m) = (n^2 pi^2 hbar^2) / (8 m a^2), psi_n (x) = sqrt(1/a) sin(n pi (x+a)/(2a))$
- 谐振子$E_n = (n + 1/2) hbar omega, psi_n (x) = N_n e^(-xi^2/2) H_n (xi) = N_n H_n (alpha x) e^(-alpha^2 x^2/2)$
- *Virial定理*：$macron(T) = (1/2 macron(sum_i hat(x)_i partial/(partial x_i) hat(V)))$ 
- *Heisenberg*：$Delta hat(F) = hat(F) - macron(F) , macron(Delta hat(F)^2) =macron(hat(F)^2) - macron(F)^2, macron(Delta hat(F)^2) macron(Delta hat(G)^2) >= 1/4 |macron(hat(C))|^2$
-  $hat(x)$：$x$,$i hbar dd("")/dd(p)$,$delta(x - x')$,$1/sqrt(2 pi hbar) e^(i/hbar p x)$；$hat(p)$：$- i hbar dd("")/dd(x)$,$p$,$1/sqrt(2 pi hbar) e^(-i/hbar p x)$ ,$delta(p - p')$
-  Bose:自旋整数，对称，光介子；Fermi:自旋半整数，反对称，电质中子

$
psi_+ (arrow(R), arrow(r)) = 1/(2 pi)^3 e^(i arrow(K) dot arrow(R))1/sqrt(2) (e^(i arrow(k) dot arrow(r)) + e^(i arrow(k) dot arrow(r))) = 1/(2 pi)^3 e^(i arrow(K) dot arrow(R)) sqrt(2) cos(arrow(k) dot arrow(r))
$
