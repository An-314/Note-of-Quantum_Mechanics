#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第8次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [9.1])[
  两个粒子可占据三个互不相同的单粒子态。分别在可分辨粒子、全同费米子和全同玻色子的情况下，求两粒子态的数目。
]
#solution[
  可分辨粒子有 $3^2=9$ 个直积态。无自旋全同费米子不能占据同一单粒子态，故有 $binom(3,2)=3$ 个反对称态。全同玻色子允许重复占据，故有
  $
    binom(3+2-1,2)=6
  $
  个对称态。
]

#exercise(subname: [9.2])[
  三个粒子可占据三个单粒子态。分别求可分辨粒子、全同费米子和全同玻色子的三粒子态数目。
]
#solution[
  可分辨粒子有 $3^3=27$ 个态。费米子必须每个单粒子态各占一个，只有 $binom(3,3)=1$ 个态。玻色子态数为
  $
    binom(3+3-1,3)=10.
  $
]

#exercise(subname: [9.3])[
  两个无相互作用的全同粒子处于 $-a<x<a$ 的无限深势阱中。分别求玻色子和无自旋费米子的基态能量与波函数。
]
#solution[
  单粒子能量为 $epsilon_n=n^2 pi^2 hbar^2/(8 m a^2)$。玻色子可同时占据 $n=1$，故
  $
    E_B=2 epsilon_1=pi^2 hbar^2/(4m a^2), quad
    Psi_B(x_1,x_2)=phi_1(x_1)phi_1(x_2).
  $
  无自旋费米子必须占据 $n=1,2$，故
  $
    E_F=epsilon_1+epsilon_2=5 pi^2 hbar^2/(8 m a^2),
  $
  $
    Psi_F=1/sqrt(2)[phi_1(x_1)phi_2(x_2)-phi_2(x_1)phi_1(x_2)].
  $
]

#exercise(subname: [9.4])[
  写出坐标表象中 $x,p,H$ 的矩阵元（积分核）。
]
#solution[
  $
    〈 x'|hat(x)|x'' 〉=x' delta(x'-x''),
  $
  $
    〈 x'|hat(p)|x'' 〉=-i hbar pdv(,x')delta(x'-x''),
  $
  $
    〈 x'|H|x'' 〉=
    [-hbar^2/(2m)frac(partial^2,partial x'^2)+V(x')]delta(x'-x'').
  $
]

#exercise(subname: [9.5])[
  写出动量表象中 $x,p,H$ 的矩阵元。
]
#solution[
  $
    〈 p'|hat(x)|p'' 〉=i hbar pdv(,p')delta(p'-p''),
  $
  $
    〈 p'|hat(p)|p'' 〉=p' delta(p'-p''),
  $
  对解析势可形式地写为
  $
    〈 p'|H|p'' 〉=
    [p'^2/(2m)+V(i hbar pdv(,p'))]delta(p'-p'').
  $
]

#exercise(subname: [9.6])[
  利用坐标与动量表象间的幺正变换推导上一题的结果。
]
#proof[
  插入坐标完备关系并使用
  $
    〈 x|p 〉=(2 pi hbar)^(-1/2)exp(i p x/hbar),
  $
  例如
  $
    〈 p'|x|p'' 〉
    =integral 〈 p'|x 〉 x 〈 x|p'' 〉 dif x
    =i hbar pdv(,p')delta(p'-p'').
  $
  对 $p$ 及 $H$ 作同样计算即可。
]

#exercise(subname: [9.7])[
  写出角动量在坐标与动量表象中的形式，并说明其对易关系是否改变。
]
#solution[
  坐标表象中
  $
    bold(L)=-i hbar bold(r) times nabla_r.
  $
  动量表象中 $bold(r)=i hbar nabla_p$，故
  $
    bold(L)=-i hbar bold(p) times nabla_p.
  $
  表象变换是幺正变换，$L_i'=U L_i U^dagger$，因此
  $
    [L_i',L_j']=i hbar epsilon_(i j k)L_k'
  $
  与原表象完全相同。
]
