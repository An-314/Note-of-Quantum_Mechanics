#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第9次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [11.1])[
  矩阵 $A,B$ 满足 $A^2=0$、$A A^dagger+A^dagger A=1$、$B=A^dagger A$。（a）证明 $B^2=B$；（b）设 $B$ 的本征值无简并，在 $B$ 表象中求 $A$ 的矩阵表示。
]
#solution[
  利用 $A A^dagger=1-A^dagger A$，
  $
    B^2=A^dagger A A^dagger A=A^dagger(1-A^dagger A)A=A^dagger A=B.
  $
  因而 $B$ 的本征值只能是 $0,1$。按 $B|0 〉=0$、$B|1 〉=|1 〉$ 排列基底，$A$ 把 $|1 〉$ 映到 $|0 〉$，相位可吸收到基矢中，故可取
  $
    B=mat(0,0;0,1), quad A=mat(0,1;0,0).
  $
  更一般地右上元可为模长等于一的相位因子。
]

#exercise(subname: [11.2])[
  厄米算符 $A,B$ 的本征值均不简并，且 $A^2=B^2=1$、$A B+B A=0$。在 $A$ 自身表象中给出 $A,B$ 的矩阵。
]
#solution[
  $A$ 的本征值为 $plus.minus 1$，故 $A=mat(1,0;0,-1)$。设 $B=mat(a,b;c,d)$，反对易关系给出 $a=d=0$；厄米性给出 $c=b^*$；$B^2=1$ 给出 $abs(b)=1$。因此
  $
    B=mat(0,exp(i theta);exp(-i theta),0).
  $
  适当改变基矢相位可令 $theta=0$。
]

#exercise(subname: [11.3])[
  证明中心力场中 $r -> 0$ 时，正则径向波函数满足 $R_l(r) prop r^l$。
]
#proof[
  径向方程在原点附近的主导项为
  $
    R''+2/r R'-l(l+1)/r^2 R=0.
  $
  设 $R prop r^s$，得到 $s(s+1)=l(l+1)$，即 $s=l$ 或 $s=-l-1$。后一解在原点发散且不可接受，故 $R_l(r) prop r^l$。
]

#exercise(subname: [11.4])[
  利用类氢能级公式讨论电子偶素、缪原子和缪子偶素的能谱。
]
#solution[
  两体库仑问题只需把电子质量换成约化质量 $mu=m_1m_2/(m_1+m_2)$：
  $
    E_n=-mu Z^2 e^4/(2(4 pi epsilon_0)^2 hbar^2 n^2).
  $
  电子偶素中 $mu=m_e/2$，故能级间隔是氢的约一半；缪原子中 $mu=m_mu M/(m_mu+M)$，轨道尺度约缩小为电子原子的 $m_e/mu$；缪子偶素 $(mu^+e^-)$ 中 $mu=m_e m_mu/(m_e+m_mu)$，能谱接近氢原子但有可测的约化质量修正。
]

#exercise(subname: [11.5])[
  求氢原子基态中径向坐标和动量的不确定度。
]
#solution[
  对 $psi_(100)=1/sqrt(pi a_0^3) exp(-r/a_0)$，
  $
    〈 r 〉=3a_0/2, quad 〈 r^2 〉=3a_0^2,
  $
  故 $Delta r=sqrt(3)a_0/2$。由球对称性 $〈 bold(p) 〉=0$；用维里定理或直接作用拉普拉斯算符得
  $
    〈 bold(p)^2 〉=hbar^2/a_0^2, quad Delta p=hbar/a_0.
  $
  因而 $Delta r Delta p=sqrt(3)hbar/2$。
]

#exercise(subname: [11.6])[
  对类氢离子的圆轨道态 $l=n-1$，其径向函数满足 $R_(n,n-1)(r) prop r^(n-1)exp(-Z r/(n a_0))$。求最概然半径、平均半径和径向涨落。
]
#solution[
  径向概率密度 $P(r)=r^2 abs(R)^2 prop r^(2 n)exp(-2 Z r/(n a_0))$，故
  $
    r_("mp")=n^2 a_0/Z.
  $
  利用伽马积分得到
  $
    〈 r 〉=(n^2+n/2)a_0/Z,
  $
  $
    Delta r=(n a_0)/(2Z)sqrt(2n+1).
  $
]
