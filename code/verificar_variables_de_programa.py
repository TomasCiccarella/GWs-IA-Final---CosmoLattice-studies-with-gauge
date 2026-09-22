"""Verifica con sympy las cuentas de la base teorica sobre variables de
programa, condiciones iniciales y restriccion de Gauss usadas por
CosmoLattice (ver ../bases-teoricas-modelos-gauge.html).

Cada bloque reproduce, de forma simbolica, una cuenta que en el codigo de
CosmoLattice (Final/CosmoLattice/include/..., Final/CosmoLattice/models/...)
aparece hecha "a mano" o implicita en los valores por defecto de los .in, y
chequea que sea consistente con las formulas publicadas en Figueroa et al.
2006.15122 ("Art I") y 2102.01031 ("User Manual"). No corre CosmoLattice ni
lee sus fuentes: es una verificacion algebraica independiente de esas
formulas, hecha a partir de lo que ambas fuentes (paper y codigo) dicen.
"""
import sympy as sp


def check(nombre, expr_deberia_ser_cero, mostrar=None):
    resultado = sp.simplify(expr_deberia_ser_cero)
    if resultado != 0:
        # sp.simplify no siempre cancela exponenciales complejas: se reintenta
        # pasando a senos/cosenos (expand_complex) antes de declarar FALLA.
        resultado = sp.simplify(sp.expand_complex(resultado))
    ok = resultado == 0
    estado = "OK" if ok else "FALLA"
    print(f"[{estado}] {nombre}")
    if mostrar is not None:
        print(f"       {mostrar}")
    if not ok:
        raise AssertionError(f"{nombre}: residuo no nulo -> {resultado}")
    return ok


print("=" * 78)
print("1) Reescalado del potencial lphi4: V = (lambda/4) phi^4")
print("=" * 78)

lam, fstar, phi_pr = sp.symbols('lambda f_star varphi_pr', positive=True)
omega_star = sp.sqrt(lam) * fstar          # models/lphi4.h: omegaStar = sqrt(lambda)*fStar
phi_fisico = fstar * phi_pr                # phi = fStar * phi_pr  (definicion de variable de programa)

V_fisico = sp.Rational(1, 4) * lam * phi_fisico**4
V_programa = sp.simplify(V_fisico / (fstar**2 * omega_star**2))

check(
    "V_programa == varphi_pr^4 / 4  (potentialTerms(0) del codigo: 0.25*fldS(0)^4)",
    V_programa - phi_pr**4 / 4,
    mostrar=f"V_programa calculado = {V_programa}",
)

print()
print("=" * 78)
print("2) Reescalado del acoplamiento hijo: q = g^2/lambda")
print("=" * 78)

q, chi_pr = sp.symbols('q varphi_pr_hijo', positive=True)
g = sp.sqrt(q * lam)                        # models/lphi4.h linea 97: g = sqrt(q*lambda)
chi_fisico = fstar * chi_pr

Vint_fisico = sp.Rational(1, 2) * g**2 * phi_fisico**2 * chi_fisico**2
Vint_programa = sp.simplify(Vint_fisico / (fstar**2 * omega_star**2))

check(
    "Vint_programa == (q/2) * varphi_pr^2 * chi_pr^2  (potentialTerms: 0.5*q*(fldS0*fldS1)^2)",
    Vint_programa - sp.Rational(1, 2) * q * phi_pr**2 * chi_pr**2,
    mostrar=f"Vint_programa calculado = {Vint_programa}",
)

print()
print("=" * 78)
print("3) Derivadas del potencial de programa (lo que codifican potDeriv/potDeriv2)")
print("=" * 78)

V_total_pr = phi_pr**4 / 4 + sp.Rational(1, 2) * q * phi_pr**2 * chi_pr**2

dV_dphi = sp.diff(V_total_pr, phi_pr)
dV_dchi = sp.diff(V_total_pr, chi_pr)
d2V_dphi2 = sp.diff(V_total_pr, phi_pr, 2)
d2V_dchi2 = sp.diff(V_total_pr, chi_pr, 2)

print(f"       dV_pr/d(varphi_pr)       = {dV_dphi}")
print(f"       dV_pr/d(chi_pr)          = {dV_dchi}")
print(f"       d^2V_pr/d(varphi_pr)^2   = {d2V_dphi2}")
print(f"       d^2V_pr/d(chi_pr)^2      = {d2V_dchi2}")
check(
    "d^2V_pr/d(varphi_pr)^2 en el punto homogeneo (chi_pr=0) define la masa m_pr^2 usada en omega_k",
    d2V_dphi2.subs(chi_pr, 0) - 3 * phi_pr**2,
)

print()
print("=" * 78)
print("4) Ecuacion de Friedmann en variables de programa (regla de la cadena)")
print("=" * 78)

alpha, omega_star_s, fstar_s, Mpl, rho_pr, a, H = sp.symbols(
    'alpha omega_star f_star M_pl rho_pr a H', positive=True
)

# Definicion del tiempo de programa: d(eta)/dt = a^{-alpha} * omega_star
# => a'(eta) = (da/dt) * (dt/d eta) = adot * a^{alpha}/omega_star
a_prima_sobre_a = H * a**alpha / omega_star_s      # esto ES a'/a, por la regla de la cadena de arriba

friedmann_fisica = H**2 - rho_pr * omega_star_s**2 * fstar_s**2 / (3 * Mpl**2)
# rho_fisica = omega_star^2 * fstar^2 * rho_pr  (reescalado de energias, Energies::kineticS etc.)

H2_en_terminos_de_rho_pr = sp.solve(friedmann_fisica, H**2)[0]
a_prima_sobre_a_cuadrado = sp.simplify(
    (a_prima_sobre_a**2).subs(H**2, H2_en_terminos_de_rho_pr)
)

formula_codigo = a**(2 * alpha) * fstar_s**2 * rho_pr / (3 * Mpl**2)

check(
    "(a'/a)^2 == a^(2 alpha) * (f_star/M_pl)^2 * rho_pr/3  (scalefactorinitializer.h / hubbleconstraint.h)",
    a_prima_sobre_a_cuadrado - formula_codigo,
    mostrar=f"(a'/a)^2 calculado = {a_prima_sobre_a_cuadrado}",
)

print()
print("=" * 78)
print("5) Normalizacion de las fluctuaciones de vacio iniciales")
print("=" * 78)

omega_star_n, fstar_n, Lbox, dx, Nlat, omega_k = sp.symbols(
    'omega_star f_star L dx N omega_k', positive=True
)

# fluctuationsgenerator.h: Norm(k) es el rms de Re(phi_k) (y, por separado, de Im(phi_k))
Norm_k = (omega_star_n / fstar_n) * (Lbox / dx**2) ** sp.Rational(3, 2) \
    * (1 / sp.sqrt(2 * omega_k)) * (1 / sp.sqrt(2))

# <|phi_k|^2> = <Re^2> + <Im^2> = 2*Norm(k)^2  (Re e Im independientes, misma varianza)
varianza_codigo = sp.simplify(2 * Norm_k**2)

# Art I eq. 435-436 / User Manual eq. 62, evaluada en a=1 (a_* = 1 por convencion) y con L = N*dx
varianza_paper = (omega_star_n / fstar_n) ** 2 * (Lbox / dx**2) ** 3 / (2 * omega_k)

check(
    "2*Norm(k)^2 == (omega_star/f_star)^2 (L/dx^2)^3 / (2 omega_k)  [formula publicada, a=1]",
    varianza_codigo - varianza_paper,
    mostrar=f"2*Norm(k)^2 = {varianza_codigo}",
)

# Chequeo adicional: con L = N*dx, (L/dx^2)^3 == (N/dx)^3, tal como se escribe en el paper
check(
    "(L/dx^2)^3 == (N/dx)^3 cuando L = N*dx",
    sp.simplify((Lbox / dx**2) ** 3 - (Nlat / dx) ** 3).subs(Lbox, Nlat * dx),
)

print()
print("=" * 78)
print("6) Condicion inicial del campo de gauge: solucion de la ley de Gauss discreta")
print("=" * 78)

kx, ky, kz = sp.symbols('k_x k_y k_z', real=True)
j0 = sp.symbols('j_0')  # corriente J_0 en espacio de Fourier (tratada como escalar simbolico)

# Derivada discreta "hacia atras" en espacio de Fourier: Delta_i^- f(n) <-> (1 - e^{-i k_i}) * f(k)
keff = [1 - sp.exp(-sp.I * kx), 1 - sp.exp(-sp.I * ky), 1 - sp.exp(-sp.I * kz)]

keff2 = sp.nsimplify(sp.simplify(sum(sp.conjugate(k) * k for k in keff)))
keff2_forma_cerrada = sp.simplify(
    keff2.rewrite(sp.cos).expand()
)
print(f"       |k_eff|^2 en forma cerrada = {keff2_forma_cerrada}")

# La condicion inicial (u1initializer.h) asigna pi_i(k) = conj(keff_i)/|keff|^2 * j0(k).
pi_asignado = [sp.conjugate(k) / keff2 * j0 for k in keff]

# La ley de Gauss discreta en espacio de Fourier es sum_i keff_i * pi_i(k) == j0(k).
lado_izquierdo = sp.simplify(sum(k * p for k, p in zip(keff, pi_asignado)))

check(
    "sum_i k_eff_i * pi_i(k) == j_0(k)  (la asignacion de u1initializer.h resuelve Gauss por construccion)",
    sp.simplify(lado_izquierdo - j0),
)

print()
print("Todas las verificaciones terminaron OK.")
