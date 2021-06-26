import matplotlib.pyplot as plt
from random import randint, uniform
from math import sqrt, log10, log2, pi, cos, sin
from numpy import random


BW = 10E6
TEMP = 27 + 273.15
BASE_P = 23 - 30
UE_P = 0 - 30
TX_G = 0
RX_G = 0
B_H = 51.5
UE_H = 1.5
BOLTZ_CONST = 1.38E-23
EX = 4
SQRT_3 = sqrt(3)
SQRT_3_div_2 = (sqrt(3) / 2)
NEG_SQRT_3 = (-1) * sqrt(3)
NEG_SQRT_3_div_2 = (-1) * (sqrt(3) / 2)
UE_NUM = 75
SCALE = 250 / SQRT_3_div_2
PAIR_DIS = 200


class Map():
    def __init__(self):
        self._cluster = []

    def add_cluster(self, cluster):
        self._cluster.append(cluster)

    @property
    def cluster(self):
        return self._cluster


class Cluster():
    def __init__(self, loc_x, loc_y, ma):
        self._bs = []
        self._map = ma
        self._bs.append(Bs(0 + loc_x, 0 + loc_y, self))
        self._bs.append(Bs(0 + loc_x, 500 + loc_y, self))
        self._bs.append(Bs(0 + loc_x, -500 + loc_y, self))
        self._bs.append(Bs(250 * NEG_SQRT_3 + loc_x, 250 + loc_y, self))
        self._bs.append(Bs(250 * NEG_SQRT_3 + loc_x, -250 + loc_y, self))
        self._bs.append(Bs(250 * SQRT_3 + loc_x, 250 + loc_y, self))
        self._bs.append(Bs(250 * SQRT_3 + loc_x, -250 + loc_y, self))

    @property
    def bs(self):
        return self._bs

    @property
    def map(self):
        return self._map

    def plot_map(self, filename):
        x = []
        y = []
        ue_x = []
        ue_y = []
        ue_x_recv = []
        ue_y_recv = []
        for i in range(len(self.bs)):
            x.append(self.bs[i].x)
            y.append(self.bs[i].y)
            for j in range(len(self.bs[i].ue)):
                if self.bs[i].ue[j]._recv == 0:
                    ue_x.append(self.bs[i].ue[j].x)
                    ue_y.append(self.bs[i].ue[j].y)
                else:
                    ue_x_recv.append(self.bs[i].ue[j].x)
                    ue_y_recv.append(self.bs[i].ue[j].y)
        plt.scatter(x, y, color='r', marker='^')
        plt.scatter(ue_x, ue_y, marker='.', color='b')
        plt.scatter(ue_x_recv, ue_y_recv, marker='.', color='g')
        plt.axis('square')
        plt.title('UE Map')
        plt.xlabel("X(m)")
        plt.ylabel("Y(m)")
        plt.savefig(filename)
        plt.close()


class Bs():
    def __init__(self, x, y, cluster):
        self._loc_x = x
        self._loc_y = y
        self._ue = []
        self._cluster = cluster
        self._map = self._cluster.map

    @property
    def ue(self):
        return self._ue

    @property
    def x(self):
        return self._loc_x

    @property
    def y(self):
        return self._loc_y

    def add_ue(self, ue):
        self._ue.append(ue)

    def gen_ue(self):
        for i in range(UE_NUM):
            ue = Ue(self)
            self._ue.append(ue)
            ue2 = Ue(self, 1, ue, ue.x, ue.y)
            self._ue.append(ue2)


class Ue():
    def __init__(self, bs, recv=0, tx=None, gen_x=None, gen_y=None):
        if gen_x is None or gen_y is None:
            self._x, self._y = gen_loc()
        else:
            self._x, self._y = gen_loc_with_initial(gen_x, gen_y)
        self._dis = sqrt(self._x ** 2 + self._y ** 2)
        self._x += bs.x
        self._y += bs.y
        self._bs = bs
        self._recv = recv
        self._tx = tx

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def bs(self):
        return self._bs


def gen_loc():
    while True:
        x = uniform(-1, 1)
        y = uniform(-1, 1)
        if (y <= SQRT_3_div_2) and \
           (y >= NEG_SQRT_3_div_2) and \
           (SQRT_3 * x + y <= SQRT_3) and \
           (SQRT_3 * x + y >= NEG_SQRT_3) and \
           (NEG_SQRT_3 * x + y <= SQRT_3) and \
           (NEG_SQRT_3 * x + y >= NEG_SQRT_3):
            return x * SCALE, y * SCALE


def gen_loc_with_initial(gen_x, gen_y):
    while True:
        theta = uniform(0, 2 * pi)
        dis = uniform(0, PAIR_DIS)
        x = gen_x / SCALE + dis / SCALE * cos(theta)
        y = gen_y / SCALE + dis / SCALE * sin(theta)
        if (y <= SQRT_3_div_2) and \
           (y >= NEG_SQRT_3_div_2) and \
           (SQRT_3 * x + y <= SQRT_3) and \
           (SQRT_3 * x + y >= NEG_SQRT_3) and \
           (NEG_SQRT_3 * x + y <= SQRT_3) and \
           (NEG_SQRT_3 * x + y >= NEG_SQRT_3):
            return x * SCALE, y * SCALE


def db_to_int(n):
    return 10 ** (n / 10)


def down_rxp(dis):
    g = (B_H * UE_H) ** 2 / (dis ** EX)
    g_db = 10 * log10(g)
    rx_p = g_db + BASE_P + TX_G + RX_G
    return rx_p


def up_rxp(dis):
    g = (B_H * UE_H) ** 2 / (dis ** EX)
    g_db = 10 * log10(g)
    rx_p = g_db + UE_P + TX_G + RX_G
    return rx_p


def up_rxp_d2d(dis):
    g = (UE_H * UE_H) ** 2 / (dis ** EX)
    g_db = 10 * log10(g)
    rx_p = g_db + UE_P + TX_G + RX_G
    return rx_p


def Sinr(power_db, inf):
    noise = BOLTZ_CONST * TEMP * BW / UE_NUM
    p = db_to_int(power_db)
    s = p / (noise + inf)
    return 10 * log10(s)


def shannon(sinr):
    return BW * log2(1 + db_to_int(sinr))


if __name__ == "__main__":
    ################1-1################
    ma = Map()
    clus = Cluster(0, 0, ma)
    cent_bs = clus.bs[0]
    cent_bs.gen_ue()
    clus.plot_map('./fig1_1.jpg')
    print('problem 1-1')
    ################1-2-a################

    sinr_lst = []
    count_lst = []
    count = 0
    for ue in cent_bs.ue:
        up_p = up_rxp(ue._dis)
        sinr_lst.append(Sinr(up_p, 0))
        count_lst.append(count)
        count += 1
    sinr_lst.sort()
    plt.scatter(sinr_lst, count_lst, marker='.')
    plt.savefig('./fig1_2_a.jpg')
    plt.close()
    print('problem 1-2a')

    ################1-2-b################

    sinr_lst = []
    count_lst = []
    count = 0
    for ue in cent_bs.ue:
        down_p = down_rxp(ue._dis)
        inf_p = 0
        for j in range(1, len(clus.bs)):
            dis = sqrt((clus.bs[j].x - ue.x) ** 2 + (clus.bs[j].y - ue.y) ** 2)
            inf_p += db_to_int(down_rxp(dis))
        sinr_lst.append(Sinr(down_p, inf_p))
        count_lst.append(count)
        count += 1
    sinr_lst.sort()
    plt.scatter(sinr_lst, count_lst, marker='.')
    plt.savefig('./fig1_2_b.jpg')
    plt.close()
    print('problem 1-2b')

    ################1-3################

    rate = 0
    for sinr in sinr_lst:
        rate += shannon(sinr)
    print(f'Throughput of Downlink : {rate}')
    print('problem 1-3')

    ################1-4################
    sinr_lst = []
    count_lst = []
    count = 0
    for ue in cent_bs.ue:
        if ue._recv == 1:
            d2d_p = 0
            for tx_ue in cent_bs.ue:
                if tx_ue._recv == 0:
                    dis = sqrt((ue.x - tx_ue.x) ** 2 + (ue.y - tx_ue.y) ** 2)
                    d2d_p += db_to_int(up_rxp_d2d(dis))
            dis = sqrt((ue.x - ue._tx.x) ** 2 + (ue.y - ue._tx.y) ** 2)
            up_p = up_rxp_d2d(dis)
            sinr_lst.append(Sinr(up_p, d2d_p - db_to_int(up_p)))
            count_lst.append(count)
            count += 1
    sinr_lst.sort()

    plt.scatter(sinr_lst, count_lst, marker='.')
    plt.savefig('./fig1_4.jpg')
    plt.close()
    print('problem 1-4')

    ################1-5################

    rate = 0
    for sinr in sinr_lst:
        rate += shannon(sinr)
    print(f'Throughput of D2D systems : {rate}')
    print('problem 1-5')

    ################1-6################
    rate_lst = []
    ue_num_lst = []

    for i in range(10):
        UE_NUM += 25
        ma = Map()
        clus = Cluster(0, 0, ma)
        cent_bs = clus.bs[0]
        cent_bs.gen_ue()
        sinr_lst = []
        count_lst = []
        count = 0
        for ue in cent_bs.ue:
            if ue._recv == 1:
                d2d_p = 0
                for tx_ue in cent_bs.ue:
                    if tx_ue._recv == 0:
                        dis = sqrt((ue.x - tx_ue.x) ** 2 + (ue.y - tx_ue.y) ** 2)
                        d2d_p += db_to_int(up_rxp_d2d(dis))
                dis = sqrt((ue.x - ue._tx.x) ** 2 + (ue.y - ue._tx.y) ** 2)
                up_p = up_rxp_d2d(dis)
                sinr_lst.append(Sinr(up_p, (d2d_p - db_to_int(up_p))))
                count_lst.append(count)
                count += 1
        sinr_lst.sort()
        rate = 0
        for sinr in sinr_lst:
            rate += shannon(sinr)
        rate_lst.append(rate)
        ue_num_lst.append(UE_NUM)
    plt.scatter(ue_num_lst, rate_lst, marker='.')
    plt.savefig('./fig1_6.jpg')
    plt.close()
    print('problem 1-6')
