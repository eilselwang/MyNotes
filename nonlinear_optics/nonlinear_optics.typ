#import "@preview/physica:0.9.8": *

#set heading(numbering: "1.1")
#set par(first-line-indent: (amount: 2em, all: false), justify: true, leading: 1em)

// reset counter at each chapter
// if you want to change the number of displayed
// section numbers, change the level there
#show heading.where(level: 1): it => {
  counter(math.equation).update(0)
  it
}

#set math.equation(
  numbering: n => {
    numbering("(1.1)", counter(heading).get().first(), n)
    // if you want change the number of number of displayed
    // section numbers, modify it this way:
    /*
    let count = counter(heading).get()
    let h1 = count.first()
    let h2 = count.at(1, default: 0)
    numbering("(1.1.1)", h1, h2, n)
    */
  },
  supplement: [方程],
)

#let indent = h(2em)

= 理论

== Maxwell方程到波动方程(Wave equation)
$
   div vb(D) & = rho \
   div vb(B) & = 0 \
  curl vb(E) & = - (partial vb(B)) / (partial t) \
  curl vb(H) & = vb(J) + (partial vb(D)) / (partial t)
$<Maxwell>
#indent @Maxwell 所有的场振幅和相位依赖于空间$vb(r)=(x, y)$，时间$t$，和传播方向$z$。电位移矢量的表达式为：
$
  hat(vb(D))(vb(r),omega,z)=epsilon_0 epsilon(omega)hat(vb(E))(vb(r),omega,z)+hat(vb(P))(vb(r),omega,z)
$<Displacement>
对@Maxwell 中的第三式取旋度，并带入@Displacement 化简得：
$
  laplacian vb(E) - grad(div vb(E)) - 1/c^2 pdv(, t, 2) integral_(-oo)^t epsilon(t-t^prime)vb(E)(vb(r), t^prime,z))dd(t^prime)=mu_0 (pdv(vb(J), t)+pdv(vb(P), t, 2))
$<time_domain>
@time_domain 中$vb(E)$，$vb(J)$，$vb(P)$依赖于$(vb(r),t,z)$，变换到空间-频率域：
$
  laplacian vb(hat(E)) - grad(div vb(hat(E))) + (omega^2 n^2(omega))/c^2 vb(hat(E))=mu_0 (-i omega vb(hat(J))-omega^2 vb(hat(P)))
$<fre_domain>
@fre_domain 中$vb(hat(E))$，$vb(hat(J))$，$vb(hat(P))$依赖于$(vb(r),omega,z)$，介电常数$epsilon(omega)=n^2(omega)$，表示材料的线性折射率。

== 标量波动方程
#indent 假设电场强度垂直于传播方向线性偏振，偏振方向为$vb(e)_s$。那么，$vb(E)=E vb(e)_s$，$vb(P)=P vb(e)_s$，$vb(J)=J vb(e)_s$。@fre_domain 中的$grad(div vb(hat(E)))$因此可以被忽略。@fre_domain 可以写成标量形式：
$
  (partial_z^2+nabla_perp^2) hat(E)(vb(r),omega,z) + k^2(omega) hat(E)(vb(r),omega,z)=-mu_0 omega^2 hat(P)(vb(r),omega,z)
$<scalar_form>
其中，$k(omega)=n(omega)omega\/c$。在需要加入电流面密度$hat(J)$时，可以将极化强度$hat(P)$替换为$hat(P)+i hat(J)\/omega$。

== 一些变换
=== 从实验室参考系变换到脉冲参考系
#indent 实验室参考系$(z,t)$变换到脉冲参考系$(zeta,tau)$：
$
  zeta & =z \
   tau & =t-z\/v_g
$
对于某一依赖关系$pdv(A(z,t), z)=pdv(A(zeta,tau), zeta) pdv(zeta, z)+pdv(A(zeta,tau), tau) pdv(tau, z)$，因此：
$
  partial_z & = partial_zeta-(1\/v_g) partial_tau \
  partial_t & = partial_tau
$<frame>
通常选脉冲的群速度作为$v_g$。在频域中，$partial_z & = partial_zeta+i(omega\/v_g)$

=== 傅里叶变换
==== 时间傅里叶变换
#indent 时间$t$和频率$omega$的傅里叶变换为：
$
  hat(A)(omega)=integral_(-oo)^(+oo)A(t)exp(i omega t)dd(t)
  approx
  T/N sum_(n=0)^(N-1)A(t_n)exp(i omega t_n)
$
逆变换为：
$
  A(t)=1/(2pi)integral_(-oo)^(+oo)hat(A)(omega)exp(-i omega t)dd(omega)
  approx
  1/T sum_(n=0)^(N-1)hat(A)(omega_n)exp(-i omega_n t)
$
其中_T_为时间窗口，_N_为采样点数。

对于连续时间傅里叶变换的卷积定理：
$
  scr(F)[f(t)*g(t)] = scr(F)[f(t)]scr(F)[g(t)]
$
其中，$scr(F)$表示连续时间傅里叶变换。对于归一化的DFT（numpy或者matlab）：
$
  scr(F)_D [f(t)*g(t)] = T scr(F)_D [f(t)] scr(F)_D [g(t)]
$
对于非归一化的DFT（FFTW或者cuFFT）：
$
  scr(F)_D [f(t)*g(t)] = T/N scr(F)_D [f(t)] scr(F)_D [g(t)]
$
由于DFT默认时间零点从索引0开始，一个函数的时间零点要从索引0开始，另一个函数的时间零点在序列中心。
==== 空间二维傅里叶变换
#indent 空间$(x,y)$和空频域$(k_x,k_y)$的傅里叶变换为：
$
  hat(A)(k_x, k_y) = integral.double_(-oo)^(+oo) A(x,y) exp(-i (k_x x + k_y y))dd(x)dd(y) approx \
  X/M Y/N sum_(m=0)^(M-1)sum_(n=0)^(N-1)A(x_m,y_n)exp(-i (k_x x_m + k_y y_n))
$
逆变换为：
$
  A(x, y) = 1/(2pi)^2 integral.double_(-oo)^(+oo) hat(A)(k_x,k_y) exp(i(k_x x + k_y y))dd(k_x)dd(k_y) approx \
  1/(X Y)sum_(m=0)^(M-1)sum_(n=0)^(N-1)hat(A)(k_x_m,k_y_n)exp(i (k_x_m x + k_y_n y))
$
其中，_X_和_Y_为两个方向的窗口长度，_M_和_N_为两个方向的采样点数。数值实现里注意fft和ifft的使用。

== 载波分辨的单向传输方程
=== 标量波动方程分解<factor_method>
#indent 采用分解@scalar_form 的方式推导Forward Maxwell Equation(FME)@Feit:88，将传播算子分解为前向和后向：
$
  (partial_z + i k(omega))(partial_z - i k(omega))hat(E) = -nabla_perp^2 hat(E)-mu_0 omega^2 hat(P)(vb(r),omega,z)
$<factor>
若@factor 不含右手边表示衍射和非线性极化的项，方程的解为前向传输的波与反向传输的波的叠加：
$
  hat(E)(omega,z) = hat(A)_+(omega)exp(i k(omega) z)+hat(A)_-(omega)exp(-i k(omega) z)
$
假设前向传输的分量远远强于后向传输的分量，$|hat(A)_-|<<|hat(A)_+|$，那么近似$partial_z + i k(omega) approx 2i k(omega)$成立，@factor 近似为FME：
$
  pdv(hat(E), z) = i k(omega) hat(E) + i/(2k(omega))nabla_perp^2 hat(E) + i/n(omega) omega/c hat(P)/(2epsilon_0)
$<FME>

#indent 利用@frame 将@FME 变换到脉冲参考系，利用傅里叶变换，得到FME在空频-时频域$(vb(k)_perp,omega)$的表示：
$
  pdv(tilde(E), zeta) = i (k(omega) - k_perp^2/(2k(omega)) - omega/v_g) tilde(E) + i omega/(c n(omega)) tilde(P)/(2epsilon_0)
$<FME_local>

=== 标准化形式
#indent 在空频-时频域，不仅仅是FME，任何单向传输方程(Unidirectional Propagation Equations)都可以化为标准形式：
$
  pdv(tilde(E), z) = i K_z (vb(k)_perp,omega) tilde(E) + i Q(vb(k)_perp,omega) tilde(P)/(2epsilon_0)
$<canonical>
在脉冲坐标下，$Q$的形式不变，$K_z$变换为$K_z - omega/v_g$：
$
  pdv(tilde(E), zeta) = i (K_z (vb(k)_perp,omega) - omega/v_g) tilde(E) + i Q(vb(k)_perp,omega) tilde(P)/(2epsilon_0)
$<canonical_local>
例如，FME的$K_z$和$Q$分别为：
$
  K_z^(("FME"))(vb(k)_perp,omega) & = k(omega) - k_perp^2/(2k(omega)) \
    Q^(("FME"))(vb(k)_perp,omega) & = omega/(c n(omega))
$

=== Slowly Evolving Wave Approximation
#indent @frame 直接带入到@scalar_form 并变换到时频域：
$
  pdv(hat(E), zeta, 2) + 2i omega/v_g pdv(hat(E), zeta) = -nabla_perp^2 hat(E) - (k^2(omega)-omega^2/v_g^2)hat(E) - omega^2/c^2 hat(P)/epsilon_0
$<scalar_form_local>
引入Slowly Evolving Wave Approximation(SEWA)，忽略@scalar_form_local 左手边$hat(E)$关于$zeta$的二阶导，即认为
$abs(pdv(hat(E), zeta, 2))<<abs(2 omega/v_g pdv(hat(E), zeta))$，或者等效的：
$
  abs(pdv(hat(E), zeta))<<abs(omega/v_g hat(E))
$<SEWA>
@SEWA 的物理意义是$hat(E)(vb(r),omega)$发生明显变化的长度尺寸远远大于$v_g/omega$。

将SEWA应用于@scalar_form_local 得到Forward Wave Equation(FWE)，实验室坐标下FWE的标准形式$K_z^(("FWE"))(vb(k)_perp,omega)$和$Q^(("FWE"))(vb(k)_perp,omega)$为：
$
  K_z^(("FWE"))(vb(k)_perp,omega) & = k(omega) + v_g/(2omega) ((k(omega)-omega/v_g)^2-k_perp^2) \
    Q^(("FWE"))(vb(k)_perp,omega) & = (omega v_g)/c^2
$<FWE>

#indent 比较FME和FWE的线性和非线性算符，可以得到以下的关系式：
$
  K_z^(("FME")) - K_z^(("FWE")) & = -omega/(2v_g) (1-(n(omega)v_g)/c)^2-k_perp^2/(2k(omega))(1-(n(omega)v_g)/c) \
      Q^(("FME")) - Q^(("FWE")) & = omega/(n(omega)c)(1-(n(omega)v_g)/c)
$<compare>

=== Unidirectional Pulse Propagation Equation
#indent @scalar_form 变换到空频-时频域：
$
  (partial_z^2 + (k^2(omega) - k_perp^2))tilde(E) = -mu_0 omega^2 tilde(P)
$
引入$K_z^(("UPPE"))=sqrt(k^2(omega) - k_perp^2)$，上面方程变成：
$
  (partial_z + i K_z)(partial_z - i K_z)tilde(E) = -mu_0 omega^2 tilde(P)
$
通过@factor_method 中类似的分解方法，可以得到Unidirectional Pulse Propagation Equation(UPPE)：
$
  pdv(tilde(E), z) = i K_z^(("UPPE"))(vb(k_perp), omega)tilde(E) + i omega^2/(c^2 K_z^(("UPPE"))(vb(k_perp), omega)) tilde(P)/(2epsilon_0)
$<UPPE>
@UPPE 取标准形式时，线性和非线性算符为：
$
  K_z^(("UPPE"))(vb(k_perp), omega) & = sqrt(k^2(omega) - k_perp^2) \
  Q_z^(("UPPE"))(vb(k_perp), omega) & = omega^2/(c^2 sqrt(k^2(omega) - k_perp^2))
$<UPPE_c>
当考虑傍轴近似时，$abs(k_perp)<<abs(k(omega))$，UPPE可以近似为FME。以FME为主。

== 基于包络的传输方程
#indent 电场强度被描述为脉冲包络$cal(E)$和一个载波频率$omega_0$的叠加：$E(vb(r_perp),t,z)=cal(E)(vb(r_perp),t,z)exp(i(k_0 z - omega_0 t))$。在脉冲坐标下：$E(vb(r_perp),tau,zeta)=cal(E)(vb(r_perp),tau,zeta)exp(i(k_0 - omega_0\/v_g)zeta - i omega_0 tau))$。

在实验室坐标下：
$
  partial_z E & = exp(i k_0 z-i omega_0 t)(partial_z + i k_0)cal(E) \
  partial_t E & = exp(i k_0 z-i omega_0 t)(partial_t - i omega_0)cal(E)
$<lab_frame>
在脉冲坐标下：
$
  partial_zeta E & = exp(i(k_0 - omega_0\/v_g)zeta - i omega_0 tau)(partial_zeta + i (k_0 - omega_0\/v_g))cal(E) \
   partial_tau E & = exp(i(k_0 - omega_0\/v_g)zeta - i omega_0 tau)(partial_tau - i omega_0)cal(E)
$<pulse_frame>
包络传输方程也有如下的标准形式：
$
  pdv(tilde(cal(E)), zeta) = i cal(K)(vb(k)_perp,Omega)tilde(cal(E)) + i cal(Q)(vb(k)_perp,Omega) tilde(cal(P))/(2epsilon_0)
$<env_ca>
其中，$Omega=omega-omega_0$。

=== Forward Envelope Equation
#indent 将@frame 和@pulse_frame 带入@canonical 可得：
$
  cal(K)(vb(k)_perp,Omega) & = K_z (vb(k)_perp, omega=Omega+omega_0)-kappa(omega=Omega+omega_0) \
  cal(Q)(vb(k)_perp,Omega) & = Q(vb(k)_perp,omega=Omega+omega_0)
$
其中，$kappa(omega=Omega+omega_0) = k_0 + (omega-omega_0)\/v_g$。

== 至此，主要的文本已经完成，还需要总结和非线性极化强度的具体形式。主要的参考文献@couairon_practitioners_2011

#bibliography("参考文献.bib", title: [参考文献])
