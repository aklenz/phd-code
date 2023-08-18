#todo: linting and commenting

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# Length bending point tube to spur
len_neck = 7.815 # mm

# Start coordinates bending point tube
coord_bend_ini = np.array((0, -len_neck))

# Start angle tube bending from horizontal
alpha_ini = 90 # degree

# Start angle lid bending between lever arms
beta_ini = 90 # degree

# Angular steps for simulation
step = 1 # degree


def calc_angle(angle_ini, angular_change):
    angle = angle_ini - angular_change
    return angle

def calc_coord_spur(alpha, len_neck, coord_bend_ini):
    # Calculate r
    radius = len_neck
    # Calculate theta in rad
    theta = alpha * np.pi / 180

    # Coordinates spur
    coord_spur = np.array((coord_bend_ini[0] + np.cos(theta)*radius,  # xs
                          coord_bend_ini[1] + np.sin(theta)*radius)).T # ys
    return coord_spur

def calc_coord_lid(alpha, beta, len_neck, coord_bend_ini):
    # Calculate r
    radius = 2 * np.sin(beta/2  * np.pi / 180.) * len_neck
    # Calculate theta in rad
    theta = (alpha * np.pi / 180) - np.arccos(radius/2 / len_neck)

    # Coordinates midlid
    coord_lid = np.array((coord_bend_ini[0] + np.cos(theta)*radius,  # xs
                          coord_bend_ini[1] + np.sin(theta)*radius)).T # ys
    return coord_lid

def cal_vectors(coord_lid):

    x_vec = np.zeros((len(coord_lid)))
    y_vec = np.zeros((len(coord_lid)))
    for i in range(1, len(coord_lid)-1):
        if i == 0:
            x_vec[i] = 1
            y_vec[i] = 0
        else:
            x_vec[i] = coord_lid[i+1, 1] - coord_lid[i-1, 1]
            y_vec[i] = coord_lid[i+1, 0] - coord_lid[i-1, 0]

    percent_v = 100* y_vec**2 / (y_vec**2 + x_vec**2)
    x_vec = abs(coord_lid[0, 1] - coord_lid[:, 1])

    return x_vec, percent_v

def plot_movement(coord_bend, coord_spur, coord_lid, ax, step, col):
    plt.sca(ax)
    plt.xlim(-0.1, 1.45*len_neck)
    plt.axis('equal')
    for i in range(0, len(coord_lid), 15):
        xs = [coord_bend[0], coord_spur[i, 0], coord_lid[i, 0]]
        ys = [coord_bend[1], coord_spur[i, 1], coord_lid[i, 1]]
        alpha = 1-0.8*step*i/90
        plt.plot(xs, ys, color='#777777', alpha = alpha)

    plt.plot(coord_lid[:, 0], coord_lid[:, 1], color=col)
    xs = [coord_bend[0], coord_spur[0, 0]]
    ys = [coord_bend[1], coord_spur[0, 1]]
    plt.scatter(xs, ys, marker='x', color='black')
    plt.scatter(coord_lid[0, 0], coord_lid[0, 1], marker='o', color='black')
    xs = [coord_bend[0], coord_spur[0, 0], coord_lid[0, 0]]
    ys = [coord_bend[1], coord_spur[0, 1], coord_lid[0, 1]]
    plt.plot(xs, ys, color='black')
    plt.xlabel("x in mm")
    return

def plot_vector(x_vec, y_vec, ax, col):
    plt.sca(ax)
    plt.plot(y_vec, -x_vec, color=col)
    for i in range(0, len(x_vec), 15):
        plt.scatter(y_vec[i-1], -x_vec[i-1], color=col, marker='+', s=100)
    return

def case(ax, ax_vec, part, col):

    angular_change_alpha = np.linspace(0, part*alpha_ini, int(alpha_ini/step) +1)
    angular_change_beta = np.linspace(0, (1-part)*alpha_ini, int(alpha_ini/step) +1)
    alpha = calc_angle(alpha_ini, angular_change_alpha)
    beta = calc_angle(beta_ini, angular_change_beta)
    coord_bend = coord_bend_ini
    coord_spur = calc_coord_spur(alpha, len_neck, coord_bend_ini)
    coord_lid = calc_coord_lid(alpha, beta, len_neck, coord_bend_ini)

    plot_movement(coord_bend, coord_spur, coord_lid, ax, step, col)

    x_vec, y_vec = cal_vectors(coord_lid)
    plot_vector(x_vec, y_vec, ax_vec, col)
    return


def main():
    simulation_no =3
    parts = np.linspace(0, 1, simulation_no)

    col = ['#666666', '#dba400', '#0062a8', '#71094f', '#578e0b', '#000000', '#999999', '#ffcf33', '#29a6ff', '#f471c8', '#96ed1c']

    # Plot all cases
    gs = gridspec.GridSpec(2, int((simulation_no + 5)/2))

    fig = plt.figure()
    fig.set_size_inches(40, 7)
    ax0 = fig.add_subplot(gs[:, :2])
    plt.sca(ax0)
    plt.ylabel("Downward movement in mm")
    plt.xlabel("Ratio of\nforward/backward movement to\ntotal movement in %")

    for i in range(len(parts)):
        if i < simulation_no/2:
            ax1 = fig.add_subplot(gs[0, i + 2], sharey=ax0)
        else:
            ax1 = fig.add_subplot(gs[1, (i + 1)-int(simulation_no/2)], sharey=ax0)
        case(ax1, ax0, parts[i], col[i])

    plt.show()

    fig.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    figname = 'Kinematic_Model.pdf'
    fig.savefig(figname)
    plt.close()


if __name__ == '__main__':
    main()
