import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.style.use('ggplot')


def I(x): return qt.qeye(x)


def sinexp(x, w, r, A, C, D):
    return A*(1+np.cos(w*x)+D)*np.exp(-x*r)+C


wc = .1 * 2*np.pi      # cavity frequency
wa = 0 * 2*np.pi       # atom frequency

g = 0.1 * 2*np.pi  # coupling strength
chi = 1*2*np.pi
kappa = 0.001         # cavity dissipation rate
gamma = 0.08          # atom dissipation rate
N = 5             # number of cavity fock states
n_th_a = 0.08e1        # avg number of atom bath excitation
n_th_c = 0.01          # avg number of cavity bath excitation

tlist = np.linspace(0, 200, 2000)

# intial state
psi0 = qt.ket2dm(
    qt.tensor((qt.basis(N, 1) + qt.basis(N, 0))/np.sqrt(2), qt.basis(2, 0)))
psi1 = qt.ket2dm(
    qt.tensor((-qt.basis(N, 1) + qt.basis(N, 0))/np.sqrt(2), qt.basis(2, 0)))


# operators
a = qt.tensor(qt.destroy(N), qt.qeye(2))
sm = qt.tensor(qt.qeye(N), qt.destroy(2))

# Hamiltonian (no RWA) jayness cummings
# H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag()*sm + a*sm.dag())
H = wc * a.dag() * a + wa * sm.dag() * sm + chi * a.dag()*a * sm.dag()*sm

rel = 10
# collapse operators
c_ops = [
    np.sqrt(0.5*kappa*(1+n_th_c)) * a,         # cavity relaxation (Thermal)
    np.sqrt(0.5*kappa*n_th_c) * a.dag(),   # cavity excitation (Thermal)
    np.sqrt(0.5*gamma*(1+n_th_a)) * sm,        # qubit relaxation  (Thermal)
    np.sqrt(0.5*gamma*n_th_a) * sm.dag(),  # qubit excitation  (Thermal)
    #     np.sqrt(0.5*gamma)            * sm,        # qubit relaxation
    #     np.sqrt(0.5*gamma*rel)        * sm
]


def proj(state):
    return qt.tensor(qt.qeye(N), qt.basis(2, 0)*qt.basis(2, 0).dag()) * state * qt.tensor(qt.qeye(N), qt.basis(2, 0)*qt.basis(2, 0).dag()) +\
        qt.tensor(qt.qeye(N), qt.basis(2, 0)*qt.basis(2, 1).dag()) * \
        state * qt.tensor(qt.qeye(N), qt.basis(2, 1)*qt.basis(2, 0).dag())


# Results. Heisenberg picture, project on photon number and "atom number"
e_ops = [psi0, psi1, sm.dag()*sm]
# output = qt.mesolve(H, psi0, tlist, c_ops, e_ops)
# outputi = qt.mesolve(H, psi0, tlist[:1000], c_ops)


life_time = []
meas_range = list(range(5, len(tlist)//3, 1))
for meas_num in meas_range:
    # meas_num = 50
    states = [psi0]
    t = 0
    while t < len(tlist):
        dt = meas_num
        # t = dt*i

        # print(f'{t} -> {t+dt}')

        output = qt.mesolve(
            H, proj(states[-1]), tlist[t:min(t+dt, len(tlist))], c_ops)
        states.extend(output.states)
        t += dt
    states = states[1:]
    fix = len(tlist)-len(states)
    states.extend([psi0]*fix)

    qstates = np.array([np.real(state.ptrace(1)) for state in states])
    cstates = np.array([np.real(state.ptrace(0)) for state in states])

    n = [[abs(np.trace(state*e_ops[j])) for state in states] for j in range(2)]

    popt, pcov = curve_fit(sinexp, tlist, n[0], maxfev=100000000,
                           bounds=([0, 0, 0, -0.1, -1], [10, 5, 1.2, 1.1, 1]), p0=[g, gamma, 1, 0, 0.5])

    print(f'{popt[1]:.2e}  ||  {fix}')
    if meas_num % 30 == 0:
        plt.figure()
        print(len(popt))
        plt.plot(tlist, n[0])
        fit = sinexp(tlist, *popt)
        plt.plot(tlist, fit)
        plt.plot(tlist, popt[2]*(2+popt[4])*np.exp(-tlist*popt[1])+popt[3])
    life_time.append(popt[1])
plt.figure()
plt.plot(meas_range, life_time)
plt.show()
