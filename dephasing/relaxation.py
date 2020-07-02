import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.style.use('ggplot')


def I(x): return qt.qeye(x)


def sinexp(x, w, r, A, C, D):
    return A*(1+np.cos(w*x)+D)*np.exp(-x*r)+C


wc = .1 * 2*np.pi  # cavity frequency
wa = 0 * 2*np.pi   # atom frequency

g = 0.1 * 2*np.pi  # coupling strength
chi = 1*2*np.pi
kappa = 0.001      # cavity dissipation rate
gamma = 0.08       # atom dissipation rate
N = 5              # number of cavity fock states
n_th_a = 0.08e1    # avg number of atom bath excitation
n_th_c = 0.01      # avg number of cavity bath excitation

tlist = np.linspace(0, 2, 10000)

# intial state
cav0 = (qt.basis(N, 1) + qt.basis(N, 0))/np.sqrt(2)
qub0 = qt.basis(2, 0)
psi0 = qt.ket2dm(qt.tensor(cav0, qub0))
psi_e = qt.ket2dm(qt.tensor(cav0, qt.basis(2,1)))

# operators
a = qt.tensor(qt.destroy(N), qt.qeye(2))
sm = qt.tensor(qt.qeye(N), qt.destroy(2))

# Hamiltonian (no RWA) jayness cummings
# H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag()*sm + a*sm.dag())
H = wc * a.dag() * a + wa * sm.dag() * sm + chi * a.dag()*a * sm.dag()*sm

rel = 10
# collapse operators
c_ops = [
    np.sqrt(0.5*kappa*(1+n_th_c)) * a,     # cavity relaxation (Thermal)
    np.sqrt(0.5*kappa*n_th_c) * a.dag(),   # cavity excitation (Thermal)
    np.sqrt(0.5*gamma*(1+n_th_a)) * sm,    # qubit relaxation  (Thermal)
    np.sqrt(0.5*gamma*n_th_a) * sm.dag(),  # qubit excitation  (Thermal)
    #     np.sqrt(0.5*gamma)            * sm,        # qubit relaxation
    #     np.sqrt(0.5*gamma*rel)        * sm
]

def proj(state):
    return qt.tensor(qt.qeye(N), qt.basis(2, 0)*qt.basis(2, 0).dag()) * state * qt.tensor(qt.qeye(N), qt.basis(2, 0)*qt.basis(2, 0).dag()) +\
        qt.tensor(qt.qeye(N), qt.basis(2, 0)*qt.basis(2, 1).dag()) * \
        state * qt.tensor(qt.qeye(N), qt.basis(2, 1)*qt.basis(2, 0).dag())


life_time = []
meas_range = list(range(1, len(tlist), 100))
for meas_num in meas_range:
    states = [psi0]
    t = 0
    while t < len(tlist):
        dt = meas_num

        output = qt.mesolve(H, proj(states[-1]), tlist[t:min(t+dt, len(tlist))], c_ops)
        states.extend(output.states)
        t += dt
    states = states[1:]
    # Some bandate bug fix - could be removed seftely
    fix = len(tlist)-len(states)
    states.extend([psi0]*fix)

    n = [abs(np.trace(state*(psi0+psi_e))) for state in states]

    popt, pcov = curve_fit(sinexp, tlist, n, maxfev=100000000,
                           bounds=([0, 0, 0, -0.1, -1], [10, 5, 1.2, 1.1, 1]), p0=[g, gamma, 1, 0, 0.5])

    print(f'{meas_num}/{meas_range[-1]} --> {popt[1]:.4e}', end='\r')
    # if meas_num % 30 != 0:
    #     plt.figure()
    #     plt.plot(tlist, n)
    #     fit = sinexp(tlist, *popt)
    #     plt.plot(tlist, fit)
    #     plt.plot(tlist, popt[2]*(2+popt[4])*np.exp(-tlist*popt[1])+popt[3])
    life_time.append(popt[1])
plt.figure(figsize=(20,10))
plt.plot(np.array(meas_range)*tlist[-1]/meas_range[-1], life_time)
# plt.axvline(chi, color='black')
plt.axvline(1/chi, color='black')

plt.figure(figsize=(20,10))
plt.plot(np.array(meas_range)*tlist[-1]/meas_range[-1], 1/np.array(life_time))
plt.show()
