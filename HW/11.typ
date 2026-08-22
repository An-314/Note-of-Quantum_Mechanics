#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第11次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [13.1])[
  对 $l=1$，在 $L_z$ 表象中求 $L_x$ 的本征值、本征矢，并求系统处于 $|1,0 〉$ 时测得各 $L_x$ 本征值的概率。
]
#solution[
  $
    L_x=hbar/sqrt(2)mat(0,1,0;1,0,1;0,1,0).
  $
  本征值为 $hbar,0,-hbar$，可取归一化本征矢
  $
    v_+=1/2 mat(1;sqrt(2);1), quad
    v_0=1/sqrt(2)mat(1;0;-1), quad
    v_-=1/2 mat(1;-sqrt(2);1).
  $
  $|1,0 〉=mat(0;1;0)$，故三个概率依次为 $1/2,0,1/2$。
]

#exercise(subname: [13.2])[
  求三个泡利矩阵的本征值和归一化本征矢，并说明自旋算符的本征值。
]
#solution[
  $
    sigma_x=mat(0,1;1,0), quad sigma_y=mat(0,-i;i,0), quad sigma_z=mat(1,0;0,-1).
  $
  三者本征值均为 $plus.minus 1$。可分别取
  $
    sigma_x: 1/sqrt(2)mat(1;plus.minus 1),
  $
  $
    sigma_y: 1/sqrt(2)mat(1;plus.minus i), quad
    sigma_z: mat(1;0),mat(0;1).
  $
  因 $S_i=hbar sigma_i/2$，自旋分量的本征值为 $plus.minus hbar/2$。
]

#exercise(subname: [13.3])[
  对单位矢量 $bold(n)=(sin theta cos phi,sin theta sin phi,cos theta)$，求 $bold(n) dot bold(sigma)$ 的本征值与本征矢。
]
#solution[
  $
    bold(n) dot bold(sigma)=mat(cos theta,exp(-i phi)sin theta;exp(i phi)sin theta,-cos theta).
  $
  本征值为 $plus.minus 1$，一组归一化本征矢为
  $
    chi_+=mat(cos(theta/2);exp(i phi)sin(theta/2)),
  $
  $
    chi_-=mat(-exp(-i phi)sin(theta/2);cos(theta/2)).
  $
]

#exercise(subname: [13.4])[
  系统处于 $S_z=hbar/2$ 的本征态，求 $S_x,S_y$ 的平均值与不确定度，并检验不确定关系。
]
#solution[
  在 $chi_z^+=mat(1;0)$ 中，$〈 S_x 〉=〈 S_y 〉=0$，且 $S_x^2=S_y^2=hbar^2/4$，故
  $
    Delta S_x=Delta S_y=hbar/2.
  $
  从而 $Delta S_x Delta S_y=hbar^2/4$，恰等于 $hbar abs(〈 S_z 〉)/2$。
]

#exercise(subname: [13.5])[
  自旋 $1/2$ 粒子处于沿方向 $(theta,phi)$ 的“向上”态。求测量 $S_x,S_y,S_z$ 得到正、负本征值的概率和平均值。
]
#solution[
  对任意单位方向 $bold(a)$，有
  $
    P_a(plus.minus)=1/2(1 plus.minus bold(a) dot bold(n)), quad
    〈 S_a 〉=hbar/2 bold(a) dot bold(n).
  $
  因而分别代入 $n_x=sin theta cos phi$、$n_y=sin theta sin phi$、$n_z=cos theta$ 即得各结果。
]

#exercise(subname: [13.6])[
  自旋 $1/2$ 粒子在沿 $z$ 轴的恒定磁场中运动。若初态为 $S_x=hbar/2$ 的本征态，求随时间演化的态以及三个自旋分量的平均值。
]
#solution[
  令 $H=hbar omega sigma_z/2$，则
  $
    |psi(t)〉=1/sqrt(2)mat(exp(-i omega t/2);exp(i omega t/2)).
  $
  直接计算得
  $
    〈 S_x 〉=hbar/2 cos(omega t), quad
    〈 S_y 〉=-hbar/2 sin(omega t), quad
    〈 S_z 〉=0.
  $
  符号随 $H$ 与旋磁比的约定而改变，但物理上是绕磁场方向的拉莫尔进动。
]

#exercise(subname: [13.7])[
  利用角动量代数证明 $J_+|j,m 〉 prop |j,m+1 〉$、$J_-|j,m 〉 prop |j,m-1 〉$，并求比例系数。
]
#proof[
  由 $[J_z,J_ plus.minus]=plus.minus hbar J_ plus.minus$ 可知升降后的磁量子数改变 $plus.minus 1$；又因 $[J^2,J_ plus.minus]=0$，$j$ 不变。范数为
  $
    norm(J_ plus.minus|j m 〉)^2
    =hbar^2[j(j+1)-m(m plus.minus 1)].
  $
  因此
  $
    J_+|j m 〉=hbar sqrt((j-m)(j+m+1))|j,m+1 〉,
  $
  $
    J_-|j m 〉=hbar sqrt((j+m)(j-m+1))|j,m-1 〉.
  $
]
