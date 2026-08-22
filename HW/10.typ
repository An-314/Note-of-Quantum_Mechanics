#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第10次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [12.1])[
  在三维谐振子球坐标基底 $|011 〉,|010 〉,|01{-1} 〉$ 中，写出笛卡尔激发态 $|100 〉,|010 〉,|001 〉$ 的列矩阵。
]
#solution[
  相位取通常的球基矢约定，则
  $
    |100 〉=1/sqrt(2)(-|011 〉+|01{-1} 〉),
  $
  $
    |010 〉=i/sqrt(2)(|011 〉+|01{-1} 〉), quad
    |001 〉=|010 〉.
  $
  因而依题给基底顺序，三者的列矩阵分别为
  $
    1/sqrt(2) mat(-1;0;1), quad
    i/sqrt(2) mat(1;0;1), quad mat(0;1;0).
  $
]

#exercise(subname: [12.2])[
  已知 $l=1$ 时 $L_x$ 的三个归一化本征矢（在 $|1,1 〉,|1,0 〉,|1,-1 〉$ 基底中）分别对应本征值 $hbar,0,-hbar$，求 $L_x$ 的矩阵。
]
#solution[
  把三个本征矢作为列组成幺正矩阵 $S$，则
  $
    L_x=S "diag"(hbar,0,-hbar)S^dagger
    =hbar/sqrt(2) mat(0,1,0;1,0,1;0,1,0).
  $
]

#exercise(subname: [12.3])[
  在 $|100 〉,|010 〉,|001 〉$ 基底中，系统态为 $|psi 〉=mat(cos theta;exp(i phi)sin theta;0)$。求测得 $L_z=hbar,0,-hbar$ 的概率。
]
#solution[
  将笛卡尔基底变换到球基底，得到
  $
    c_+=(-cos theta-i exp(i phi)sin theta)/sqrt(2), quad
    c_0=0,
  $
  $
    c_-=(cos theta-i exp(i phi)sin theta)/sqrt(2).
  $
  因而
  $
    P(hbar)=(1+sin(2theta)sin phi)/2, quad P(0)=0,
  $
  $
    P(-hbar)=(1-sin(2theta)sin phi)/2.
  $
]

#exercise(subname: [12.4])[
  自由电子在沿 $z$ 轴的均匀磁场中运动，并被限制在 $x-y$ 平面。朗道规范下定态波函数的 $y$ 因子为平面波，求定态中 $y$ 方向的平均机械速度。
]
#solution[
  取电子电荷为 $-e$、朗道规范 $bold(A)=(0,B x,0)$，并令 $psi=exp(i k_y y)phi_n(x-x_0)$。机械动量为 $Pi_y=p_y+e B x$，谐振子中心满足 $x_0=-hbar k_y/(e B)$，故
  $
    〈 v_y 〉=1/m(hbar k_y+e B 〈 x 〉)
    =1/m(hbar k_y+e B x_0)=0.
  $
  同理 $〈 v_x 〉=0$。
]

#exercise(subname: [12.5])[
  分别在对称规范和朗道规范下求电子机械动量 $Pi_x,Pi_y$ 的对易关系，并说明其是否依赖规范。
]
#solution[
  对电荷 $q=-e$，$bold(Pi)=bold(p)-q bold(A)=bold(p)+e bold(A)$。利用 $partial_x A_y-partial_y A_x=B$，
  $
    [Pi_x,Pi_y]=-i hbar e B.
  $
  对称规范 $bold(A)=(-B y/2,B x/2,0)$ 与朗道规范 $bold(A)=(0,B x,0)$ 均给出同一结果。它只依赖磁场 $bold(B)=nabla times bold(A)$，因而与规范选择无关。
]

#exercise(subname: [12.6])[
  在上题磁场中再加入沿 $x$ 轴负向的均匀电场 $cal(E)$。在朗道规范下求定态能级、波函数以及平均机械速度。
]
#solution[
  势能为 $e cal(E)x$。令 $psi=exp(i k_y y)phi(x)$，哈密顿量化为平移谐振子。以 $omega_c=e B/m$，中心
  $
    x_c=-hbar k_y/(e B)-e cal(E)/(m omega_c^2),
  $
  则本征函数为 $exp(i k_y y)phi_n(x-x_c)$，能量为
  $
    E_(n,k_y)=(n+1/2)hbar omega_c
    -e cal(E)hbar k_y/(e B)-e^2 cal(E)^2/(2m omega_c^2).
  $
  平均速度满足
  $
    〈 v_x 〉=0, quad 〈 v_y 〉=-cal(E)/B,
  $
  即产生与电荷符号无关的 $bold(E) times bold(B)$ 漂移（方向由题设坐标确定）。
]
