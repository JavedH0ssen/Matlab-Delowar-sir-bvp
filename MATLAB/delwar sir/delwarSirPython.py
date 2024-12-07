import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

def bvp2d(x, y, PR):
    AE = 1.4; BE = 0.5; M = 1.2; G = 0.2; EC = 0.1; R = 0.6; FW = 1.5; GR = 10.0
    X = (G + M * AE / (AE**2 + BE**2))
    Z = (R - M * BE / (AE**2 + BE**2))
    yy1 = -(y[0]*y[2] - X*y[1] + Z*y[3] + GR*y[5])
    yy2 = -(y[0]*y[4] - X*y[3] - Z*y[1])
    yy3 = -(PR*y[0]*y[6] - 2*PR*y[1]*y[5] + PR*EC*(y[2]**2 + y[4]**2) + PR*EC*M*(y[1]**2 + y[3]**2) / (AE**2 + BE**2))
    return np.vstack((y[1], y[2], yy1, y[4], yy2, y[6], yy3))

def bc2d(ya, yb):
    FW = 1.5
    return np.array([ya[1] - 1, ya[0] - FW, ya[3], ya[5] - 1, yb[1], yb[3], yb[5]])

def jobs(PR, ROW, LENGTH):
    AE = 1.4; BE = 0.5; M = 1.2; G = 0.2; EC = 0.1; R = 0.6; FW = 1.5; GR = 10.0
    x = np.linspace(0, LENGTH, 100)
    y = np.zeros((7, x.size))
    sol = solve_bvp(lambda x, y: bvp2d(x, y, PR), bc2d, x, y)
    plt.plot(sol.x, sol.y[ROW, :])

def sbplot(ROW, LENGTH):
    jobs(1.38, ROW, LENGTH)
    jobs(1, ROW, LENGTH)
    jobs(0.7, ROW, LENGTH)
    plt.legend(["1.38", "1", "0.7"])

plt.figure()
plt.subplot(2, 2, 1)
sbplot(2, 3)
plt.title("Primary Velocity(F)")
plt.subplot(2, 2, 2)
sbplot(4, 3.5)
plt.title("Secondary Velocity(G)")
plt.subplot(2, 2, 3)
sbplot(6, 2.5)
plt.title("Temperature (θ)")
plt.show()