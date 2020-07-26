import numpy as np
import qutip as qt


def proj(state, cavity_levels, qub_lvls=3, mu=[1, 1, 1], rnd=False):
    '''
    Projecting the qubit into the ground state while keeping the cavity state the same

    State(qobj) - The cavity-transmon state to prject
    mu(float)   - The error rate (0-no errors  1-all errors)
    rnd(bool)   - Weather to use density matrix probabilites for the errors or do them randomly (True=do randomly)
    '''
    g = qt.basis(qub_lvls, 0)   # |g>
    e = qt.basis(qub_lvls, 1)   # |e>
    f = qt.basis(qub_lvls, 2)   # |f>
    N = qt.qeye(cavity_levels)  # I

    N_gg = qt.tensor(N, g*g.dag())  # |N> X |g><g|
    N_ge = qt.tensor(N, g*e.dag())  # |N> X |g><e|
    N_eg = qt.tensor(N, e*g.dag())  # |N> X |e><g|
    N_ee = qt.tensor(N, e*e.dag())  # |N> X |e><e|

    N_ff = qt.tensor(N, f*f.dag())  # |N> X |e><e|
    N_gf = qt.tensor(N, g*f.dag())  # |N> X |e><e|
    N_fg = qt.tensor(N, f*g.dag())  # |N> X |e><e|

    g_proj = mu[0]*N_gg*state*N_gg + (1 - mu[0])*N_eg*state*N_ge
    e_proj = mu[1]*N_ge*state*N_eg + (1 - mu[1])*N_ee*state*N_ee
    f_proj = mu[2]*N_gf*state*N_fg + (1 - mu[2])*N_ff*state*N_ff

    return g_proj + e_proj + f_proj
