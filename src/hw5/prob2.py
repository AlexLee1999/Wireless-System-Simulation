import matplotlib.pyplot as plt
from random import randint, uniform
from math import sqrt, log10, log2
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

BS_BUFFER_SIZE = 15E6
UE_BUFFER_SIZE = 0.5E6

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
        for i in range(len(self.bs)):
            x.append(self.bs[i].x)
            y.append(self.bs[i].y)
            for j in range(len(self.bs[i].ue)):
                ue_x.append(self.bs[i].ue[j].x)
                ue_y.append(self.bs[i].ue[j].y)
        plt.scatter(x, y, color='r', marker='^')
        plt.scatter(ue_x, ue_y, marker='.')
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

    @property
    def adj(self):
        return self._adj

    def add_ue(self, ue):
        self._ue.append(ue)

    def gen_ue(self):
        for i in range(UE_NUM):
            ue = Ue(self)
            self._ue.append(ue)
            ue2 = Ue(self, 1, ue)
            while sqrt((ue.x - ue2.x) ** 2 + (ue.y - ue2.y) ** 2) > 200:
                ue2 = Ue(self, 1, ue) #add Tx pair
            self._ue.append(ue2)
            ue.set_pair(ue2)

    def cal_buffer(self):
        tot = 0
        for ue in self.ue:
            if ue._is_recv == 1:
                tot += ue._tx_buffer
        return tot


class Ue():
    def __init__(self, bs, is_recv=0, pair=None):
        self._x, self._y = gen_loc()
        self._x += bs.x
        self._y += bs.y
        self._bs = bs
        self._rate = None
        self._is_recv = is_recv
        self._pair = pair
        self._tx_buffer = 0
        self._buffer = 0
        self._d2d_rate = 0

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def bs(self):
        return self._bs

    @property
    def rate(self):
        return self._rate

    def set_tx_buffer(self, da):
        self._tx_buffer = da

    def set_buffer(self, da):
        self._buffer = da

    @rate.setter
    def rate(self, r):
        self._rate = r

    def set_pair(self, ue):
        self._pair = ue


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
    ma = Map()
    clus = Cluster(0, 0, ma)
    cent_bs = clus.bs[0]
    cent_bs.gen_ue()
    d2d_p = 0
    for ue in cent_bs.ue:
        if ue._is_recv == 1:
            dis = sqrt((ue.x - ue._pair.x)**2+(ue.y - ue._pair.y)**2)
            d2d_p += db_to_int(up_rxp_d2d(dis))
    for ue in cent_bs.ue:
        dis = sqrt((ue.x-ue._pair.x)**2+(ue.y-ue._pair.y)**2)
        up_p = up_rxp_d2d(dis)
        ue._d2d_rate = shannon((Sinr(up_p, d2d_p - db_to_int(up_p))))
        print(ue._d2d_rate)
    # Calculate rate for downlink and uplink

    for ue in cent_bs.ue:
        if ue._is_recv == 1:
            inf_p = 0
            for j in range(1, len(clus.bs)):
                dis = sqrt((clus.bs[j].x - ue.x) ** 2 + (clus.bs[j].y - ue.y) ** 2)
                inf_p += db_to_int(down_rxp(dis))
            sinr = Sinr(down_rxp(sqrt((cent_bs.x - ue.x) ** 2 + (cent_bs.y - ue.y) ** 2)), inf_p)
            shan = shannon(sinr)
            ue.rate = shan
        else:
            sinr = Sinr(up_rxp(sqrt((cent_bs.x - ue.x) ** 2 + (cent_bs.y - ue.y) ** 2)), 0)
            shan = shannon(sinr)
            ue.rate = shan



    for rate in range(100000, 2000000, 200000):
        loss_data = 0
        total_data = 0
        for _ in range(1000):
            for ue in cent_bs.ue:
                if ue._is_recv == 0: # Sending UE
                    data = random.poisson(lam=rate)
                    total_data += data
                    ready_data = data + ue._buffer
                    ue.set_buffer(max((ready_data - ue.rate), 0))
                    if ue._buffer > UE_BUFFER_SIZE:
                        loss_data += ue._buffer - UE_BUFFER_SIZE
                        ue.set_buffer(UE_BUFFER_SIZE)
                    if ready_data - ue.rate < 0:
                        ue._pair._tx_buffer = ready_data
                    else:
                        ue._pair._tx_buffer = ue.rate

            # BS bit loss
            for ue in cent_bs.ue:
                if ue._is_recv == 1:
                    ue.set_tx_buffer(max(ue._tx_buffer - ue.rate, 0))
            current = cent_bs.cal_tx_buffer()
            # print(current)
            current_loss_data = current - BS_BUFFER_SIZE
            if current_loss_data > 0:
                loss_data += current_loss_data
                for ue in cent_bs.ue:
                    if ue._is_recv == 1:
                        ue.set_tx_buffer(ue._tx_buffer * BS_BUFFER_SIZE / current)
        # print(loss_data, total_data)


    for ue in cent_bs.ue:
        ue.set_buffer(0)
    for rate in range(100000, 2000000, 200000):
        loss_data = 0
        total_data = 0
        for _ in range(1000):
            for ue in cent_bs.ue:
                if ue._is_recv == 0: # Sending UE
                    data = random.poisson(lam=rate)
                    total_data += data
                    ready_data = data + ue._buffer
                    ue.set_buffer(max((ready_data - ue._d2d_rate), 0))
                    if ue._buffer > UE_BUFFER_SIZE:
                        loss_data += ue._buffer - UE_BUFFER_SIZE
                        ue.set_buffer(UE_BUFFER_SIZE)
        print(loss_data/total_data)
