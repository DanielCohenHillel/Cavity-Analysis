import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
import tools
plt.style.use('ggplot')


def sinexp(x, w, r, A, C, D):
    return A*(1+np.cos(w*x)+D)*np.exp(-x*r)+C


p = tools.params

mu = 0.001

life_time = []
t_max = 10*np.argmax(-abs(p.tlist-1/p.chi_e))
meas_range = list(range(2, t_max, 5))
for meas_num in meas_range:
    states = [p.psi0]
    t = 0
    while t < len(p.tlist):
        dt = meas_num

        output = qt.mesolve(
            p.H, tools.proj(states[-1], p.N, mu), p.tlist[t:min(t+dt, len(p.tlist))], p.c_ops)
        states.extend(output.states)
        t += dt
    states = states[1:]
    # Some bandate bug fix - could be removed seftely
    fix = len(p.tlist)-len(states)
    states.extend([p.psi0]*fix)

    n = [abs(np.trace(state*(p.psi0+p.psi_e))) for state in states]

    popt, pcov = curve_fit(sinexp, p.tlist, n, maxfev=100000000,
                           bounds=([0, 0, 0, -0.1, -1], [10, 5, 1.2, 1.1, 1]), p0=[p.g, p.gamma, 1, 0, 0.5])

    print(
        f'{mu} | {100*meas_num/meas_range[-1]:.2f}% --> {popt[1]:.4e}', end='\r')
    life_time.append(popt[1])
    # if meas_num % 6 != 0:
    plt.figure()
    plt.plot(p.tlist, n)
    fit = sinexp(p.tlist, *popt)
    plt.plot(p.tlist, fit)
    plt.plot(p.tlist, popt[2]*(2+popt[4])*np.exp(-p.tlist*popt[1])+popt[3])
fig, (ax1, ax2) = plt.subplots(2, figsize=(30, 20))
ax1.set_title('Dephasing rate vs $dt$ (measure spacing)')
ax1.set_xlabel(r'dt [us]')
ax1.set_ylabel('Dephasing rate')
ax1.plot(np.linspace(p.tlist[0], p.tlist[t_max], len(
    life_time)), life_time, label='Decay rate')
# ax1.xlim([0, 15])
# plt.axvline(chi, color='black')
ax1.axvline(1/p.chi_e, color='black', label=r'$\frac{1}{\chi_e}$')
ax1.legend()

ax2.set_title(r'Dephasing time ($\frac{1}{r}$) vs $dt$ (measure spacing)')
ax2.set_xlabel(r'dt [us]')
ax2.set_ylabel('Dephasing rate')
ax2.plot(np.linspace(p.tlist[0], p.tlist[t_max], len(
    life_time)), 1/np.array(life_time), label='Life time')
ax2.axvline(1/p.chi_e, color='black', label=r'$\frac{1}{\chi_e}$')
ax2.legend()

plt.subplots_adjust(hspace=0.2)
plt.show()
