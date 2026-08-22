#import "@preview/scripst:1.1.2": *
#show: scripst.with(title: [量子力学第5次作业], author: "Anzreww", time: "2024年", matheq-depth: 1, cb-counter-depth: 1)

#exercise(subname: [5.1])[
  在平面波基底中验证动量算符 $hat(p)=-i hbar dif_x$ 的厄米性。
]
#proof[
  取 $〈 x|p 〉=(2 pi hbar)^(-1/2)exp(i p x/hbar)$，则
  $
    〈 p'|hat(p)|p 〉
    =p delta(p-p').
  $
  另一方面
  $
    〈 p'|hat(p)|p 〉^*=p' delta(p-p').
  $
  因为 $(p-p')delta(p-p')=0$，两者相等，故动量算符在广义归一化意义下为厄米算符。
]

#exercise(subname: [5.2])[
  对平方可积的束缚态波函数证明 $hat(p)=-i hbar dif_x$ 是厄米算符。
]
#proof[
  分部积分得
  $
    integral phi^*(-i hbar psi') dif x
    =[-i hbar phi^* psi]_{-infinity}^{infinity}
    +integral (-i hbar phi')^* psi dif x.
  $
  束缚态在无穷远处趋于零，边界项消失，因此
  $
    〈 phi|hat(p) psi 〉=〈 hat(p) phi|psi 〉.
  $
]

#exercise(subname: [5.3])[
  证明伴随运算的恒等式
  $
    (A^dagger)^dagger=A, quad (A B)^dagger=B^dagger A^dagger,
    quad (A^(-1))^dagger=(A^dagger)^(-1).
  $
]
#proof[
  由伴随的定义 $〈 phi|A psi 〉=〈 A^dagger phi|psi 〉$，连续使用两次即得第一式。对乘积有
  $
    〈 phi|A B psi 〉
    =〈 A^dagger phi|B psi 〉
    =〈 B^dagger A^dagger phi|psi 〉,
  $
  故顺序反转。最后由 $A A^(-1)=1$ 取伴随，得 $(A^(-1))^dagger A^dagger=1$，从而得到第三式。
]

#exercise(subname: [5.4])[
  由 $[L_i,L_j]=i hbar epsilon_(i j k)L_k$ 证明 $[L^2,L_i]=0$。
]
#proof[
  以 $i=x$ 为例，
  $
    [L^2,L_x]=[L_y^2,L_x]+[L_z^2,L_x]
    =-i hbar(L_y L_z+L_z L_y)
    +i hbar(L_z L_y+L_y L_z)=0.
  $
  循环置换即可得到其余两式。
]

#exercise(subname: [5.5])[
  粒子处于 $0<x<a$ 的无限深势阱中，初态为 $psi(x)=A x(a-x)$。求归一化常数以及测得各能级的概率。
]
#solution[
  由 $integral_0^a abs(psi)^2 dif x=1$ 得
  $
    A=sqrt(30/a^5).
  $
  用本征函数 $phi_n=sqrt(2/a)sin(n pi x/a)$ 展开，
  $
    c_n=integral_0^a phi_n(x) psi(x) dif x
    =2 sqrt(60)/(n^3 pi^3)(1-(-1)^n).
  $
  因而
  $
    P_n=abs(c_n)^2=240/(n^6 pi^6)(1-cos(n pi))^2.
  $
  偶数 $n$ 的概率为零；基态概率 $P_1=960/pi^6$。
]

#exercise(subname: [5.6])[
  对可作幂级数展开的函数 $F(x,p)$，证明
  $
    [x,F]=i hbar pdv(F,p), quad [p,F]=-i hbar pdv(F,x).
  $
]
#proof[
  由 $[x,p^n]=sum_(r=0)^(n-1)p^r[x,p]p^(n-1-r)=i hbar n p^(n-1)$，对每个单项式逐项应用乘积对易关系即可得到第一式。类似地 $[p,x^n]=-i hbar n x^(n-1)$，从而得到第二式。
]

#exercise(subname: [5.7])[
  定义径向动量
  $
    p_r=1/2(bold(e)_r dot bold(p)+bold(p) dot bold(e)_r).
  $
  求其坐标表示、$[r,p_r]$、$p_r^2$，并给出它与总动量平方的关系。
]
#solution[
  注意到 $nabla dot bold(e)_r=2/r$，故
  $
    p_r=-i hbar (pdv(,r)+1/r), quad [r,p_r]=i hbar.
  $
  再作用一次得到
  $
    p_r^2=-hbar^2 (frac(partial^2,partial r^2)+2/r pdv(,r)).
  $
  与球坐标中的拉普拉斯算符比较可知
  $
    bold(p)^2=p_r^2+L^2/r^2.
  $
]
