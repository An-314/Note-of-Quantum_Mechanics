#import "@preview/scripst:1.1.2": *

#show: scripst.with(
  title: [量子力学第1次作业],
  author: "Anzreww",
  time: "2024年",
  matheq-depth: 1,
  cb-counter-depth: 1,
)

#exercise(subname: [1.1])[
  按照玻尔兹曼分布分别计算频率为 $nu$ 的经典谐振子与能量量子化谐振子的平均能量和热容量，并讨论高、低温极限。
]

#solution[
  经典谐振子的配分函数为 $Z = integral_0^infinity exp(-beta E) dif E = 1/beta$，故
  $
    bar(E)=-pdv(ln Z,beta)=k_B T, quad C_V=k_B.
  $
  量子谐振子的能级为 $E_n=(n+1/2)hbar omega$，其中 $omega=2 pi nu$。于是
  $
    Z=exp(-beta hbar omega/2)/(1-exp(-beta hbar omega)),
  $
  $
    bar(E)=hbar omega/2+(hbar omega)/(exp(beta hbar omega)-1),
  $
  $
    C_V=k_B (beta hbar omega)^2
    exp(beta hbar omega)/(exp(beta hbar omega)-1)^2.
  $
  当 $k_B T << hbar omega$ 时 $C_V -> 0$；当 $k_B T >> hbar omega$ 时 $C_V -> k_B$，回到经典结果。
]

#exercise(subname: [1.2])[
  由普朗克黑体辐射公式求低频和高频极限，并求温度为 $2.7 "K"$ 的黑体辐射中能谱最强处所对应的波长和频率。
]

#solution[
  按频率表示的能量密度为
  $
    u(nu,T)=8 pi h nu^3/c^3 dot 1/(exp(h nu/(k_B T))-1).
  $
  低频时 $h nu << k_B T$，展开指数得到瑞利—金斯公式
  $
    u(nu,T) approx 8 pi k_B T nu^2/c^3.
  $
  高频时 $h nu >> k_B T$，得到维恩公式
  $
    u(nu,T) approx 8 pi h nu^3/c^3 exp(-h nu/(k_B T)).
  $
  按波长的峰值满足 $lambda_m T=2.898 times 10^(-3) "m K"$，故
  $
    lambda_m approx 1.07 "mm".
  $
  按频率的峰值满足 $h nu_m/(k_B T)=2.82144$，故
  $
    nu_m approx 1.59 times 10^11 "Hz".
  $
  两个峰值不能用 $nu_m=c/lambda_m$ 相互换算，因为 $u_nu dif nu=u_lambda abs(dif lambda)$ 含有雅可比因子。
]

#exercise(subname: [1.3])[
  从普朗克公式导出维恩位移定律，并求出其中的常数。
]

#solution[
  波长谱为
  $
    u_lambda(lambda,T)=8 pi h c/lambda^5 dot 1/(exp(h c/(lambda k_B T))-1).
  $
  令 $x=h c/(lambda k_B T)$，对 $lambda$ 求极值得
  $
    5(1-exp(-x))=x.
  $
  非零根为 $x=4.96511$，所以
  $
    lambda_m T=(h c)/(4.96511 k_B)=2.898 times 10^(-3) "m K".
  $
]

#exercise(subname: [1.4])[
  推导康普顿散射的波长改变量，并说明为什么可见光的康普顿效应难以观察。
]

#solution[
  初始电子静止。由能量和动量守恒
  $
    h nu+m c^2=h nu'+E_e, quad bold(p)_gamma=bold(p)'_gamma+bold(p)_e.
  $
  消去电子的能量、动量并用 $p_gamma=h/lambda$，得
  $
    lambda'-lambda=h/(m c)(1-cos theta)=lambda_C(1-cos theta).
  $
  电子的康普顿波长 $lambda_C=2.426 times 10^(-12) "m"$。可见光波长约为 $10^(-7) "m"$，相对改变量最多仅约 $10^(-5)$，因而难以分辨；对 X 射线则容易观察。
]

#exercise(subname: [1.5])[
  入射 X 射线与运动电子发生背向散射。分别讨论电子初速度与入射光同向和反向的情形，求散射光波长，并说明何时有一条谱线的波长保持不变。
]

#solution[
  设电子初始四动量为 $(E/c,p)$，入射、背散射光子的动量分别为 $h/lambda$ 与 $-h/lambda'$。由四动量守恒并利用 $E^2-p^2c^2=m^2c^4$，可得
  $
    lambda'=lambda (E-p c)/(E+p c)+2 h c/(E+p c).
  $
  这里 $p$ 取沿入射方向为正；电子反向运动时令 $p<0$ 即可。写成 $E=gamma m c^2$、$p=gamma m v$，则
  $
    lambda'=lambda (1-beta)/(1+beta)+2 lambda_C/(gamma(1+beta)).
  $
  若要求 $lambda'=lambda$，则
  $
    lambda=lambda_C/(gamma beta)=h/(gamma m v).
  $
  因而在电子与光同向运动且入射波长等于电子的德布罗意波长时，会出现不发生波长移动的情形。
]
