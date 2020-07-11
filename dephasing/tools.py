import numpy as np
import qutip as qt
from typing import NamedTuple


class Params(NamedTuple):
    wc: float
    wa: float

    g: float
    chi_e: float
    chi_f: float
    kappa: float
    gamma: float
    N: float
    n_th_a: float
    n_th_c: float

    tlist: np.ndarray

    cav0: qt.qobj.Qobj
    qub0: qt.qobj.Qobj
    psi0: qt.qobj.Qobj
    psi_e: qt.qobj.Qobj

    a: qt.qobj.Qobj
    sm: qt.qobj.Qobj

    H: qt.qobj.Qobj

    c_ops: list


wc = .1 * 2*np.pi  # cavity frequency
wa = 0 * 2*np.pi   # atom frequency

g = 0.1 * 2*np.pi  # coupling strength
chi_e = 2*np.pi*0.05
chi_f = 2*np.pi*0.005
kappa = 1e-6       # cavity dissipation rate
gamma = 0.004      # atom dissipation rate
N = 5              # number of cavity fock states
qub_lvls = 3       # Total number of qubit levels
n_th_a = 0.02      # avg number of atom bath excitation
n_th_c = 0.01      # avg number of cavity bath excitation

tlist = np.linspace(0, 70/chi_e, 1000)  # MUST: max(tlist) >> 1/chi

# intial state
cav0 = (qt.basis(N, 1) + qt.basis(N, 0))/np.sqrt(2)
qub0 = qt.basis(qub_lvls, 0)
psi0 = qt.ket2dm(qt.tensor(cav0, qub0))
psi_e = qt.ket2dm(qt.tensor(cav0, qt.basis(qub_lvls, 1)))

# operators
a = qt.tensor(qt.destroy(N), qt.qeye(qub_lvls))
sm = qt.tensor(qt.qeye(N), qt.destroy(qub_lvls))

# Hamiltonian (no RWA) jayness cummings
# H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag()*sm + a*sm.dag())
# g = qt.tensor(qt.qeye(N), qt.basis(qub_lvls, 0))
e = qt.tensor(qt.qeye(N), qt.basis(qub_lvls, 1))
f = qt.tensor(qt.qeye(N), qt.basis(qub_lvls, 2))
H = wc * a.dag() * a + wa * sm.dag() * sm + a.dag() * \
    a*(chi_e*e*e.dag() + chi_f*f*f.dag())

# collapse operators
c_ops = [
    np.sqrt(0.5*kappa*(1+n_th_c)) * a,     # cavity relaxation (Thermal)
    np.sqrt(0.5*kappa*n_th_c) * a.dag(),   # cavity excitation (Thermal)
    np.sqrt(0.5*gamma*(1+n_th_a)) * sm,    # qubit relaxation  (Thermal)
    np.sqrt(0.5*gamma*n_th_a) * sm.dag(),  # qubit excitation  (Thermal)
    #     np.sqrt(0.5*gamma)            * sm,        # qubit relaxation
    #     np.sqrt(0.5*gamma*rel)        * sm
]

params = Params(
    wc,
    wa,

    g,
    chi_e,
    chi_f,
    kappa,
    gamma,
    N,
    n_th_a,
    n_th_c,

    tlist,

    cav0,
    qub0,
    psi0,
    psi_e,

    a,
    sm,

    H,

    c_ops
)


def proj(state, cavity_levels, mu=0, rnd=False):
    """
    Projecting the qubit into the ground state while keeping the cavity state the same

    State(qobj) - The cavity-transmon state to prject
    mu(float)   - The error rate (0-no errors  1-all errors)
    rnd(bool)   - Weather to use density matrix probabilites for the errors or do them randomly (True=do randomly)
    """
    g = qt.basis(qub_lvls, 0)
    e = qt.basis(qub_lvls, 1)
    N = qt.qeye(cavity_levels)
    N_gg = qt.tensor(N, g*g.dag())  # |N> X |g><g|
    N_ge = qt.tensor(N, g*e.dag())  # |N> X |g><e|
    N_eg = qt.tensor(N, e*g.dag())  # |N> X |e><g|
    N_ee = qt.tensor(N, e*e.dag())  # |N> X |e><e|

    if rnd:
        if np.random.random() > mu:
            return N_gg*state*N_gg + N_ge*state*N_eg
        else:
            return N_eg*state*N_ge + N_ee*state*N_ee

    # (1-u)*(|g><g| S |g><g|)
    return (1 - mu)*(N_gg*state*N_gg + N_ge*state*N_eg) + mu*(N_eg*state*N_ge + N_ee*state*N_ee)
