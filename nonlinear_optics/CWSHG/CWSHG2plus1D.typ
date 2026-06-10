#import "@preview/physica:0.9.8": *
#import "@preview/cetz:0.4.2": *

#set document(title: [Continuous Wave Second Harmonics Generation], author: "Wei")
#set page(header: align(right, context document.title), paper: "a4", numbering: "1")
#set par(justify: true, first-line-indent: 2em)
#set text(size: 10pt)

#place(top + center, float: true, scope: "parent", clearance: 2em)[
  #show title: set text(size: 17pt)
  #show title: set block(below: 1.2em)
  #title()
  #align(center)[Wei]
]

#show heading.where(level: 1): set align(center)
#show heading.where(level: 1): set text(
  size: 12pt,
  weight: "regular",
)
#show heading.where(level: 1): smallcaps
#show heading.where(level: 2): set text(size: 10pt, style: "italic", weight: "regular")
#show heading.where(level: 2): it => { it.body + [.] }
#set math.equation(numbering: "(1)")

= Theory
== Coupled wave equations
Considering diffraction, walk off and nonlinearity, the coupled wave equations are:
$
  pdv(A_("SH"), z) &= i/(2k_("SH"))nabla_perp^2A_("SH") + tan rho pdv(A_("SH"), y) + (i omega_("SH"))/(2n_("SH")c)d_("eff")A_"F"^2exp(-i Delta k z) \
  pdv(A_("F"), z) &= i/(2k_("F"))nabla_perp^2A_("F") + (i omega_("F"))/(n_("F")c)d_("eff") A_"SH" A_"F"^* exp(i Delta k z) \
  Delta k &= k_("SH")-2k_"F"
$<CWE>
For LBO in XY plane, $d_"eff"=d_(32) cos phi$, where $d_(32)=-0.98$ pm/V.

The input gaussian beam is:
$
  A(x,y,z) = A_0 omega_0/omega(z) exp(-(x^2+y^2)/(omega^2(z)))exp(i k(z+(x^2+y^2)/(2R(z)))-i arctan(z/z_R))
$<gaussian_beam>
where $omega_0$ and $omega(z)$ is beam radius at $1\/e^2$ intensity and z is distance from the beam waist. If the average
power is $P$ and the refractive index of the medium is $n$, other parameters are:
$
         k & = (2pi)/lambda \
  omega(z) & = omega_0 sqrt(1+(z/z_R)^2) \
      R(z) & = z+z_R^2/z \
       z_R & = (pi omega_0^2)/lambda \
       A_0 & = sqrt((4P)/(pi omega_0^2 n c epsilon_0))
$

#figure(image("../fig/beam.svg"), caption: [Gaussian beam is focused to the nonlinear crystal])

The definition of D4σ radius is:
$
  omega_x = 2sqrt((integral x^2 I(x,y)dd(x)dd(y))/(integral I(x,y)dd(x)dd(y)))
$

== Numerical algorithm
@CWE is solved by using RK4IP method@balac_embedded_2013, but without adaptive step size control. @CWE can be organized as follows:
$
  pdv(A_S, z) & = cal(D)_S A_S + cal(N)_S (A_S, A_F) \
  pdv(A_F, z) & = cal(D)_F A_S + cal(N)_F (A_S, A_F)
$
where $cal(D)$ represents linear effects and $cal(N)$ represents nonlinear coupling. Then transform $A_S$ and $A_F$ to the
interaction picture:
$
  A_S^"ip" = exp(-(z-z^prime)cal(D)_S)A_S \
  A_F^"ip" = exp(-(z-z^prime)cal(D)_F)A_F
$
The coupled equations in the interaction picture are:
$
  pdv(A_S^"ip", z) & = exp(-(z-z^prime)cal(D)_S) cal(N)_S (exp((z-z^prime)cal(D)_S)A_S^"ip", exp((z-z^prime)cal(D)_F)A_F^"ip") \
  pdv(A_F^"ip", z) & = exp(-(z-z^prime)cal(D)_F) cal(N)_F (exp((z-z^prime)cal(D)_S)A_S^"ip", exp((z-z^prime)cal(D)_F)A_F^"ip")
$<IP_CWE>
By choosing $z^prime = z+h/2$, the @IP_CWE can be solved by RK4 method:
$
  alpha_"1S" = exp(h/2 cal(D)_S) cal(N)_S (A_S, A_F) \
  alpha_"1F" = exp(h/2 cal(D)_F) cal(N)_F (A_S, A_F) \
  alpha_"2S" = cal(N)_S (A_S^"ip"+h/2 alpha_"1S", A_F^"ip"+h/2 alpha_"1F") \
  alpha_"2F" = cal(N)_F (A_S^"ip"+h/2 alpha_"1S", A_F^"ip"+h/2 alpha_"1F") \
  alpha_"3S" = cal(N)_S (A_S^"ip"+h/2 alpha_"2S", A_F^"ip"+h/2 alpha_"2F") \
  alpha_"3F" = cal(N)_F (A_S^"ip"+h/2 alpha_"2S", A_F^"ip"+h/2 alpha_"2F") \
  alpha_"4S" = cal(N)_S (exp(h/2 cal(D)_S)(A_S^"ip"+h alpha_"3S"), exp(h/2 cal(D)_F)(A_F^"ip"+h alpha_"3F")) \
  alpha_"4F" = cal(N)_F (exp(h/2 cal(D)_S)(A_S^"ip"+h alpha_"3S"), exp(h/2 cal(D)_F)(A_F^"ip"+h alpha_"3F")) \
  A_S (z+h) = exp(h/2 cal(D)_S) (A_S^"ip"+h/6 (alpha_"1S"+2 alpha_"2S"+2 alpha_"3S"))+h/6 alpha_"4S" \
  A_F (z+h) = exp(h/2 cal(D)_F) (A_F^"ip"+h/6 (alpha_"1F"+2 alpha_"2F"+2 alpha_"3F"))+h/6 alpha_"4F"
$
One should use the form of @IP_CWE and finally apply $z^prime = z+h/2$ when deriving the above relations and remember to include
the phase matching term when doing actual numerical calculation.


#bibliography("ref.bib", style: "american-physics-society")
