import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tools import proj
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

tlist = np.linspace(0, 150/chi, 1000)  # MUST: max(tlist) >> 1/chi

# intial state
cav0 = (qt.basis(N, 1) + qt.basis(N, 0))/np.sqrt(2)
qub0 = qt.basis(2, 0)
psi0 = qt.ket2dm(qt.tensor(cav0, qub0))
psi_e = qt.ket2dm(qt.tensor(cav0, qt.basis(2, 1)))

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

# def proj(state):
#     return qt.tensor(qt.qeye(N), qt.basis(2, 0)*qt.basis(2, 0).dag()) * state * qt.tensor(qt.qeye(N), qt.basis(2, 0)*qt.basis(2, 0).dag()) +\
#         qt.tensor(qt.qeye(N), qt.basis(2, 0)*qt.basis(2, 1).dag()) * \
#         state * qt.tensor(qt.qeye(N), qt.basis(2, 1)*qt.basis(2, 0).dag())

for mu in [1]:
    t_max = 30*np.argmax(-abs(tlist-1/chi))
    life_time = []
    meas_range = list(range(2, t_max, 10))
    for meas_num in meas_range:
        states = [psi0]
        t = 0
        while t < len(tlist):
            dt = meas_num

            output = qt.mesolve(
                H, proj(states[-1], N, mu), tlist[t:min(t+dt, len(tlist))], c_ops)
            states.extend(output.states)
            t += dt
        states = states[1:]
        # Some bandate bug fix - could be removed seftely
        fix = len(tlist)-len(states)
        states.extend([psi0]*fix)

        n = [abs(np.trace(state*(psi0+psi_e))) for state in states]

        popt, pcov = curve_fit(sinexp, tlist, n, maxfev=100000000,
                               bounds=([0, 0, 0, -0.1, -1], [10, 5, 1.2, 1.1, 1]), p0=[g, gamma, 1, 0, 0.5])

        print(
            f'{mu} | {100*meas_num/meas_range[-1]:.2f}% --> {popt[1]:.4e}', end='\r')
        life_time.append(popt[1])
        # if meas_num % 6 != 0:
        # plt.figure()
        # plt.plot(tlist, n)
        # fit = sinexp(tlist, *popt)
        # plt.plot(tlist, fit)
        # plt.plot(tlist, popt[2]*(2+popt[4])*np.exp(-tlist*popt[1])+popt[3])
    fig, (ax1, ax2) = plt.subplots(2, figsize=(30, 20))
    ax1.set_title('Dephasing rate vs $dt$ (measure spacing)')
    ax1.set_xlabel(r'dt [us]')
    ax1.set_ylabel('Dephasing rate')
    ax1.plot(np.linspace(tlist[0], tlist[t_max], len(
        life_time)), life_time, label='Decay rate')
# plt.axvline(chi, color='black')
    ax1.axvline(1/chi, color='black', label=r'$\frac{1}{\chi}$')
    ax1.legend()

    ax2.set_title(r'Dephasing time ($\frac{1}{r}$) vs $dt$ (measure spacing)')
    ax2.set_xlabel(r'dt [us]')
    ax2.set_ylabel('Dephasing rate')
    ax2.plot(np.linspace(tlist[0], tlist[t_max], len(
        life_time)), 1/np.array(life_time), label='Life time')
    ax2.axvline(1/chi, color='black', label=r'$\frac{1}{\chi}$')
    ax2.legend()

    plt.subplots_adjust(hspace=0.2)
plt.show()
