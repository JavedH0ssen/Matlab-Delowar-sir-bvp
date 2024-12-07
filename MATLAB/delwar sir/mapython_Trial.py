import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt


def main(aa, bb, cc, name):
    sol = jobs(aa)
    sol1 = jobs(bb)
    sol2 = jobs(cc)
    fig, axs = plt.subplots(3, 2)
    plotting(2, "Pimary Velocity(F')", sol, sol1, sol2, aa, bb, cc, name, axs[0, 0])
    plotting(4, "Secondary Velocity(G)", sol, sol1, sol2, aa, bb, cc, name, axs[0, 1])
    plotting(6, "Concentation 1 (\u039E)", sol, sol1, sol2, aa, bb, cc, name, axs[1, 0])
    plotting(8, "Temparature(\u03B8)", sol, sol1, sol2, aa, bb, cc, name, axs[1, 1])
    axs[1, 1].set_xlim([0, 3.2])
    plotting(
        10,
        "Non Partical Volume Fraction(\u03C6)",
        sol,
        sol1,
        sol2,
        aa,
        bb,
        cc,
        name,
        axs[2, 0],
    )
    axs[2, 0].set_xlim([0, 3.5])
    plt.show()
    plt.close()


def plotting(line, titles, sol, sol1, sol2, aa, bb, cc, name, ax):
    ax.plot(sol.x, sol.y[line, :], "r", label=f"{name}={aa}")
    ax.plot(sol1.x, sol1.y[line, :], "g", label=f"{name}={bb}")
    ax.plot(sol2.x, sol2.y[line, :], "b", label=f"{name}={cc}")
    ax.axhline(0, color="black")
    ax.axvline(0, color="black")
    ax.legend()
    ax.set_title(titles, weight='bold')  
    ax.set_xlabel("eta(\u03B7)")
    ax.set_ylabel(titles)
    # Adjust the spacing between subplots
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    ax.tick_params(axis='both', which='major', labelsize='small', pad=1)

    


def jobs(VARABLE):
    RC = 0.1
    FW = 0.5
    M = 0.5
    BI = 1.5
    BE = 1.5
    AE = 1 + BE * BI
    R = 0.6
    GR = 8
    GM = 6
    PR = VARABLE
    EC = 0.3
    S = 0.5
    S0 = 1
    K = 0.5
    SC = 0.6
    D = 2
    GAMMA = 1
    EPS = 0.8
    ST = 0.5
    STT = 0.5
    NEBLA = 2
    LAMB = 0.8
    short = AE**2 + BE**2
    XX = 1 + D

    def bvp2d(x, y):
        yy3 = -(
            D * y[6]
            + GR * y[7]
            + GM * y[9]
            + y[0] * y[2] / (EPS**2)
            - (K + M * AE / short) * y[1]
            + (R - (M * BE / short)) * y[3]
            - GAMMA * y[1] ** 2
        ) / ((1 + D) / EPS)
        yy5 = -(
            y[0] * y[4] / EPS**2
            - (K + M * AE / short) * y[3]
            - (R - M * BE / short) * y[1]
            - GAMMA * y[3] ** 2
        ) / ((1 + D) / EPS)
        yy7 = -(y[5] * y[1] + y[0] * y[6] - 2 * LAMB * y[5] - LAMB * y[2]) / NEBLA
        Holder = -(
            (1 + D) * PR * EC * (y[2] ** 2 + y[4] ** 2)
            + EC * M * PR * (y[1] ** 2 + y[3] ** 2) / short
            + PR * y[0] * y[8]
            - ST * PR * y[1]
        )
        yy9 = Holder
        yy11 = -(
            y[0] * SC * y[10] - STT * SC * y[1] + S0 * SC * Holder - RC * SC * y[9]
        )
        return np.vstack(
            [y[1], y[2], yy3, y[4], yy5, y[6], yy7, y[8], yy9, y[10], yy11]
        )

    def bc2d(y0, yinf):
        return np.array(
            [
                y0[1] - 1,
                y0[0] - FW,
                y0[3],
                y0[5] + y0[2] / 2,
                y0[7] - 1 + ST / 2,
                y0[9] - 1 + STT / 2,
                yinf[1],
                yinf[3],
                yinf[5],
                yinf[7],
                yinf[9],
            ]
        )

    x = np.linspace(0, 6, 100)
    y = np.zeros((11, x.size))
    sol = solve_bvp(bvp2d, bc2d, x, y)
    return sol


main(0.71, 1, 1.5, "B")
