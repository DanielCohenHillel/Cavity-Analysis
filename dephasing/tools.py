import numpy as np
import qutip as qt


def proj(state, cavity_levels, mu=0, rnd=False):
    """
    Projecting the qubit into the ground state while keeping the cavity state the same

    State(qobj) - The cavity-transmon state to prject
    mu(float)   - The error rate (0-no errors  1-all errors)
    rnd(bool)   - Weather to use density matrix probabilites for the errors or do them randomly (True=do randomly)
    """
    g = qt.basis(2, 0)
    e = qt.basis(2, 1)
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

    return (1 - mu)*(N_gg*state*N_gg + N_ge*state*N_eg) + mu*(N_eg*state*N_ge + N_ee*state*N_ee)
