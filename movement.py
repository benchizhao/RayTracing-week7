import numpy as np
import input
import Radiation_Force
from decimal import Decimal
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def mass():
    r = input.radius * 1E-6
    V = 4 / 3 * np.pi * r ** 3
    m = input.density * V
    return m

def gravity():
    N = mass()*input.g
    return np.array([-N, float(0.)])

def total_force():
    g_f = gravity()
    rad_f = Radiation_Force.radiation_force()
    return g_f + rad_f

def acc():
    return total_force()/mass()

def state():
    pos = input.droplet_pos[0:2]
    v = input.droplet_pos[2:4]
    iteration = 10
    delt = 0.01
    # plt.figure()
    x = []
    y = []
    plt.figure()
    for i in range(iteration):
        print(input.droplet_pos)
        print(acc(), input.y_displacement)
        dvx = acc()[0] * delt
        dvy = acc()[1] * delt
        dx = input.droplet_pos[2]*delt
        dy = input.droplet_pos[3]*delt

        input.droplet_pos = [input.droplet_pos[0] + dx, input.droplet_pos[1] + dy, input.droplet_pos[2] + dvx, input.droplet_pos[3] + dvy, i]
        if input.droplet_pos[2]>0:
            input.power -= abs(input.droplet_pos[2])*0.03
        elif input.droplet_pos[2] < 0:
            input.power += abs(input.droplet_pos[2])*0.03

        # if input.droplet_pos[3]>0:
        #     input.y_displacement += abs(input.droplet_pos[3])*0.2
        # elif input.droplet_pos[3]<0:
        #     input.y_displacement -= abs(input.droplet_pos[3])*0.2
        x.append(input.droplet_pos[0])
        y.append(input.droplet_pos[1])

        # plt.plot(input.droplet_pos[0],input.droplet_pos[1],'k.')
        # plt.plot(i,input.droplet_pos[0],'k.')
        # plt.plot(i, input.droplet_pos[1], 'b.')
    plt.plot(x,y)
    plt.show()
    return x,y


#
#     for i in range(iteration):
#         a = total_force()/mass()
#         y_new = x[-1]*delt
#         y.append(y_new)
#         v_new = a*delt
#         v.append((v_new))
#         if y[-1] < 0


if __name__ == '__main__':
    # print(gravity())
    # state()



    # x, y = state()
    #
    # def update_points(num):
    #     point_ani.set_data(x[num], y[num])
    #     return point_ani,
    #
    # fig = plt.figure(tight_layout=True)
    # plt.plot(x, y)
    # point_ani, = plt.plot(x[0], y[0], "ro")
    # plt.grid(ls="--")
    # # 开始制作动画
    # ani = animation.FuncAnimation(fig, update_points, np.arange(0, len(x)), interval=100, blit=False)
    #
    # ani.save('position.gif', writer='imagemagick', fps=10)
    # plt.show()




    print(acc())
    print(gravity())
    print(Radiation_Force.radiation_force())
    print(total_force())


    # g = 9.8
    # l = 1
    #
    #
    # def diff2(d_list, t):
    #     omega, theta = d_list
    #     return np.array([-g / l * theta, omega])
    #
    #
    # t = np.linspace(0, 20, 200)
    # result = odeint(diff2, [0, 35 / 180 * np.pi], t)
    # # 结果是一个两列的矩阵， odeint中第二个是初始单摆角度35度
    # plt.plot(t, result[:, 0],'.')  # 输出omega随时变化曲线
    # plt.plot(t, result[:, 1],'.')  # 输出theta随时变化曲线，即方程解
    # plt.show()
