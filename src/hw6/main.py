import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from shapely.geometry import Point, Polygon


SQRT3 = np.sqrt(3)


class Env():
    temperature = T = 300           # K
    k           = 1.38064852E-23    # Bolzmann constant
    msNum       = 50

    def __init__(self, cells):
        self.cells = cells
        self.MSs = None
        for cell in self.cells:
            cell.env = self

    def clear(self):
        for cell in self.cells:
            cell.clear()
        self.MSs = None

    def __str__(self):
        coords = str([[int(c) for c in cell.coord] for cell in self.cells])[2:-2]
        coords = coords.replace("], [", ")\n(")
        return "({})".format(coords)

    def twoRayGround(ht_m, hr_m, d_m):
        return 10 * np.log10(
            (ht_m * hr_m)**2 / d_m**4
        )

    def dBm2dB(p_dBm):
        return p_dBm - 30

    def dB2ratio(p_dB):
        return 10 ** (p_dB * 0.1)

    def ratio2dB(p_ratio):
        return 10 * np.log10(p_ratio)

    def plot(self):
        fig, ax = plt.subplots(1)
        ax.set_aspect("equal")
        colors = [["yellow"]] + 6 * [["orange"]] + 12 * [["purple"]]
        coords = np.array([cell.getCoord() for cell in self.cells])
        idx = 0
        for x, y, c in zip(coords[:, 0], coords[:, 1], colors):
            color = c[0].lower()
            cell = RegularPolygon(
                (x, y),
                numVertices = 6,
                radius = Cell.R,
                orientation = np.radians(30),
                facecolor = color,
                alpha = 0.2,
                edgecolor = "k"
            )
            ax.add_patch(cell)
            ax.scatter(x, y, color="k", marker="^")
            idx += 1
            plt.text(x, y, str(idx), fontsize=10)
        plt.xlim(-4 * Cell.R, 4 * Cell.R)
        plt.ylim(-2.5 * Cell.ISD, 2.5 * Cell.ISD)
        plt.savefig("./img/Cell_ID.jpg")
        plt.close()



class Cell():
    cellNum = 0
    interSiteDistance = ISD = 500   # m
    radius = R = ISD / SQRT3        # m

    def __init__(self, coord):
        self.env = None
        Cell.cellNum += 1
        self.cellID = Cell.cellNum
        self.coord = np.array(coord)
        self.bounds = + np.array([
            [Cell.R / 2, Cell.ISD / 2],
            [Cell.R, 0],
            [Cell.R / 2, -Cell.ISD / 2],
            [-Cell.R / 2, -Cell.ISD / 2],
            [-Cell.R, 0],
            [-Cell.R / 2, Cell.ISD / 2]
        ]) + self.coord
        self.hexagon = Polygon(self.bounds)
        self.MSs = []
        self.msNum = 0
        self.adjUp = None
        self.adjUR = None
        self.adjUL = None
        self.adjLo = None
        self.adjLR = None
        self.adjLL = None

    def clear(self):
        self.MSs = []
        self.msNum = 0

    def __str__(self):
        return "({}, {})".format(
            round(self.coord[0]),
            round(self.coord[1])
        )

    def getCoord(self):
        coord = self.coord
        return coord

    def getCellID(self):
        return int(self.cellID)

    def addMS(self, ms):
        ms.link(self)
        self.MSs.append(ms)
        self.msNum += 1

    def removeMS(self, ms):
        ms.unlink(self)
        self.MS.remove(ms)
        self.msNum -= 1

    def getMS(self):
        MSs = ""
        for ms in self.MSs:
            MSs += "{}\n".format(str(ms))
        return MSs[:-1]

    def addRandomMSs(self, num):
        count = 0
        while count < num:
            x = np.random.uniform(
                -2 * Cell.R,
                2 * Cell.R
            ) + self.coord[0]
            y = np.random.uniform(
                -2 * Cell.ISD,
                2 * Cell.ISD
            ) + self.coord[1]
            if Point(x, y).within(self.hexagon):
                ms = MS((x, y))
                ms.env = self.env
                self.addMS(ms)
                count += 1

    def plot(self):
        fig, ax = plt.subplots(1)
        ax.set_aspect("equal")
        cell = RegularPolygon(
            self.coord,
            numVertices = 6,
            radius = Cell.R,
            orientation = np.radians(30),
            facecolor = "blue",
            alpha = 0.2,
            edgecolor = "k"
        )
        ax.add_patch(cell)
        ax.scatter(self.coord[0], self.coord[1], color="k", marker="^")
        for ms in self.MSs:
            coord = ms.getCoord()
            ax.scatter(
                coord[0], coord[1], alpha = 0.5, marker = "."
            )
        plt.xlim(self.coord[0] - Cell.R, self.coord[0] + Cell.R)
        plt.ylim(self.coord[1] - .5 * Cell.ISD, self.coord[1] + .5 * Cell.ISD)





class BS(Cell):
    power = 33              # dBm
    bandwidth = BW = 10E6   # Hz
    gainTx = 14             # dB
    gainRx = 14             # dB
    height = H = 1.5 + 50   # m
    powerNoise = N = 10 * np.log10(Env.k * Env.T * BW)  #dB

    def __init__(self, coord):
        super().__init__(coord)
        self.ID = Cell.cellNum
        self.bufferSize = 6E6
        self.cbr_l, self.cbr_m, self.cbr_h = 1E6, 2E6, 3E6
        self._lambda_l, self._lambda_m, self._lambda_h = 1E6, 2E6, 3E6

    def __str__(self):
        return "({}, {})".format(
            round(self.coord[0]),
            round(self.coord[1])
        )

    def getCoord(self):
        coord = self.coord
        return coord


class MS():
    power = 0               # dBm
    gainTx = 14             # dB
    gainRx = 14             # dB
    height = H = 1.5        # m

    def __init__(self, coord):
        self.env = None
        self.coord = np.array(coord)
        self.bs = None
        self.BW = 0
        self.N = 0
        self.capacity = 0

    def __str__(self):
        return "({}, {})".format(
            round(self.coord[0]),
            round(self.coord[1])
        )

    def getCoord(self):
        return self.coord

    def link(self, bs):
        self.bs = bs

    def unlink(self):
        self.bs = None

    def getBandwidth(self):
        self.BW = self.bs.BW / self.bs.msNum
        bw = self.BW
        return bw

    def getPowerNoise(self):
        self.N = 10 * np.log10(Env.k * Env.T * self.getBandwidth())
        noise = self.N
        return noise

    def powerRx(self, bs, dist=None):
        if dist is None:
            dist = bs.coord - self.coord
        return Env.twoRayGround(
            MS.H, BS.H, np.linalg.norm(dist)
        ) + Env.dBm2dB(BS.power) + BS.gainTx + MS.gainRx

    def constsinr(self, bs):
        powerSig = sum([Env.dB2ratio(self.powerRx(cell)) for cell in self.env.cells[:7]])
        powerInf = sum([Env.dB2ratio(self.powerRx(cell)) for cell in self.env.cells[7:]])

        return Env.ratio2dB(powerSig) - Env.ratio2dB(powerInf + Env.dB2ratio(self.getPowerNoise()))

    def sinr(self, bs):
        powerSig = self.powerRx(bs)
        powerInf = sum(
            [Env.dB2ratio(self.powerRx(cell)) for cell in self.env.cells]
        ) - Env.dB2ratio(powerSig)
        return powerSig - Env.ratio2dB(
            powerInf + Env.dB2ratio(self.getPowerNoise())
        )

    def shannonCapacity(self, bs, sinr=0):
        if sinr == 0:
            self.capacity = self.getBandwidth() * np.log2(1 + Env.dB2ratio(self.sinr(bs)))
        else:
            self.capacity = self.getBandwidth() * np.log2(1 + Env.dB2ratio(sinr))
        return self.capacity


# The coordinates of the base stations
coordBSs = np.array(
     [
         [0, 0],
         [0, Cell.ISD],
         [1.5 * Cell.R, 0.5 * Cell.ISD],
         [1.5 * Cell.R, -0.5 * Cell.ISD],
         [0, -Cell.ISD],
         [-1.5 * Cell.R, -0.5 * Cell.ISD],
         [-1.5 * Cell.R, 0.5 * Cell.ISD],
         [0, 2 * Cell.ISD],
         [1.5 * Cell.R, 1.5 * Cell.ISD],
         [3 * Cell.R, Cell.ISD],
         [3 * Cell.R, 0],
         [3 * Cell.R, -Cell.ISD],
         [1.5 * Cell.R, -1.5 * Cell.ISD],
         [0, -2 * Cell.ISD],
         [-1.5 * Cell.R, -1.5 * Cell.ISD],
         [-3 * Cell.R, -Cell.ISD],
         [-3 * Cell.R, 0],
         [-3 * Cell.R, Cell.ISD],
         [-1.5 * Cell.R, 1.5 * Cell.ISD],
     ]
)
BSs = np.array(
    [BS(coord) for coord in coordBSs]
)

# Construct the 19-cell environment
env0 = Env(BSs)
env0.cells[0].adjUp = env0.cells[1]
env0.cells[0].adjUL = env0.cells[6]
env0.cells[0].adjUR = env0.cells[2]
env0.cells[0].adjLo = env0.cells[4]
env0.cells[0].adjLL = env0.cells[5]
env0.cells[0].adjLR = env0.cells[3]

env0.cells[1].adjUp = env0.cells[7]
env0.cells[1].adjUL = env0.cells[18]
env0.cells[1].adjUR = env0.cells[8]
env0.cells[1].adjLo = env0.cells[0]
env0.cells[1].adjLL = env0.cells[6]
env0.cells[1].adjLR = env0.cells[2]

env0.cells[2].adjUp = env0.cells[8]
env0.cells[2].adjUL = env0.cells[1]
env0.cells[2].adjUR = env0.cells[9]
env0.cells[2].adjLo = env0.cells[3]
env0.cells[2].adjLL = env0.cells[0]
env0.cells[2].adjLR = env0.cells[10]

env0.cells[3].adjUp = env0.cells[2]
env0.cells[3].adjUL = env0.cells[0]
env0.cells[3].adjUR = env0.cells[10]
env0.cells[3].adjLo = env0.cells[12]
env0.cells[3].adjLL = env0.cells[4]
env0.cells[3].adjLR = env0.cells[11]

env0.cells[4].adjUp = env0.cells[0]
env0.cells[4].adjUL = env0.cells[5]
env0.cells[4].adjUR = env0.cells[3]
env0.cells[4].adjLo = env0.cells[13]
env0.cells[4].adjLL = env0.cells[14]
env0.cells[4].adjLR = env0.cells[12]

env0.cells[5].adjUp = env0.cells[6]
env0.cells[5].adjUL = env0.cells[16]
env0.cells[5].adjUR = env0.cells[0]
env0.cells[5].adjLo = env0.cells[14]
env0.cells[5].adjLL = env0.cells[15]
env0.cells[5].adjLR = env0.cells[4]

env0.cells[6].adjUp = env0.cells[18]
env0.cells[6].adjUL = env0.cells[17]
env0.cells[6].adjUR = env0.cells[1]
env0.cells[6].adjLo = env0.cells[5]
env0.cells[6].adjLL = env0.cells[16]
env0.cells[6].adjLR = env0.cells[0]

env0.cells[7].adjUp = env0.cells[11]
env0.cells[7].adjUL = env0.cells[12]
env0.cells[7].adjUR = env0.cells[15]
env0.cells[7].adjLo = env0.cells[1]
env0.cells[7].adjLL = env0.cells[18]
env0.cells[7].adjLR = env0.cells[8]

env0.cells[8].adjUp = env0.cells[15]
env0.cells[8].adjUL = env0.cells[7]
env0.cells[8].adjUR = env0.cells[14]
env0.cells[8].adjLo = env0.cells[2]
env0.cells[8].adjLL = env0.cells[1]
env0.cells[8].adjLR = env0.cells[9]

env0.cells[9].adjUp = env0.cells[14]
env0.cells[9].adjUL = env0.cells[8]
env0.cells[9].adjUR = env0.cells[13]
env0.cells[9].adjLo = env0.cells[10]
env0.cells[9].adjLL = env0.cells[2]
env0.cells[9].adjLR = env0.cells[17]

env0.cells[10].adjUp = env0.cells[9]
env0.cells[10].adjUL = env0.cells[2]
env0.cells[10].adjUR = env0.cells[17]
env0.cells[10].adjLo = env0.cells[11]
env0.cells[10].adjLL = env0.cells[3]
env0.cells[10].adjLR = env0.cells[16]

env0.cells[11].adjUp = env0.cells[10]
env0.cells[11].adjUL = env0.cells[3]
env0.cells[11].adjUR = env0.cells[16]
env0.cells[11].adjLo = env0.cells[7]
env0.cells[11].adjLL = env0.cells[12]
env0.cells[11].adjLR = env0.cells[15]

env0.cells[12].adjUp = env0.cells[3]
env0.cells[12].adjUL = env0.cells[4]
env0.cells[12].adjUR = env0.cells[11]
env0.cells[12].adjLo = env0.cells[18]
env0.cells[12].adjLL = env0.cells[13]
env0.cells[12].adjLR = env0.cells[7]

env0.cells[13].adjUp = env0.cells[4]
env0.cells[13].adjUL = env0.cells[14]
env0.cells[13].adjUR = env0.cells[12]
env0.cells[13].adjLo = env0.cells[17]
env0.cells[13].adjLL = env0.cells[9]
env0.cells[13].adjLR = env0.cells[18]

env0.cells[14].adjUp = env0.cells[5]
env0.cells[14].adjUL = env0.cells[15]
env0.cells[14].adjUR = env0.cells[4]
env0.cells[14].adjLo = env0.cells[9]
env0.cells[14].adjLL = env0.cells[8]
env0.cells[14].adjLR = env0.cells[13]

env0.cells[15].adjUp = env0.cells[16]
env0.cells[15].adjUL = env0.cells[11]
env0.cells[15].adjUR = env0.cells[5]
env0.cells[15].adjLo = env0.cells[8]
env0.cells[15].adjLL = env0.cells[7]
env0.cells[15].adjLR = env0.cells[14]

env0.cells[16].adjUp = env0.cells[17]
env0.cells[16].adjUL = env0.cells[10]
env0.cells[16].adjUR = env0.cells[6]
env0.cells[16].adjLo = env0.cells[15]
env0.cells[16].adjLL = env0.cells[11]
env0.cells[16].adjLR = env0.cells[5]

env0.cells[17].adjUp = env0.cells[13]
env0.cells[17].adjUL = env0.cells[9]
env0.cells[17].adjUR = env0.cells[18]
env0.cells[17].adjLo = env0.cells[16]
env0.cells[17].adjLL = env0.cells[10]
env0.cells[17].adjLR = env0.cells[6]

env0.cells[18].adjUp = env0.cells[12]
env0.cells[18].adjUL = env0.cells[13]
env0.cells[18].adjUR = env0.cells[7]
env0.cells[18].adjLo = env0.cells[6]
env0.cells[18].adjLL = env0.cells[17]
env0.cells[18].adjLR = env0.cells[1]







if __name__ == "__main__":
    simIt = 100

    ################################ 1-1 ################################
    env0.plot()

    ################################ 1-2 ################################
    innerCells = [cell for cell in env0.cells[:7]]
    dataPoints_1, dataPoints_2 = [], []
    resourceEffs_1, resourceEffs_2 = [], []
    for idx in range(simIt):
        for cell in innerCells:
            cell.clear()
            cell.addRandomMSs(np.random.randint(5, 16))
        poorestSINRs_1, poorestSINRs_2 = [], []
        for cell in innerCells:
            poorestSINRs_1.append(
                min([ms.sinr(cell) for ms in cell.MSs])
            )
            poorestSINRs_2.append(
                min([ms.constsinr(cell) for ms in cell.MSs])
            )
        dataPoints_1.append(poorestSINRs_1)
        dataPoints_2.append(min(poorestSINRs_2))
        resourceEff_1 = []; i = 0
        msCount = 0
        for cell in innerCells:
            shnCp = cell.MSs[0].shannonCapacity(cell, sinr=dataPoints_1[-1][i])
            resourceEff_1.append(
                shnCp / cell.MSs[0].getBandwidth() * cell.msNum
            )
            msCount += cell.msNum
            i += 1
        resourceEffs_1.append(resourceEff_1)
        shnCp = innerCells[0].MSs[0].shannonCapacity(innerCells[0], sinr=dataPoints_2[-1])
        resourceEffs_2.append(
            shnCp / innerCells[0].MSs[0].getBandwidth() * msCount
        )

    dataPoints_1 = np.array(dataPoints_1)
    for idx in range(7):
        dataPoint_1 = dataPoints_1[:, idx]
        dataPoint_1.sort()
        count = 0; Y = []
        for point in dataPoint_1:
            count += 1
            Y.append(count/len(dataPoint_1))
        plt.plot(dataPoint_1, Y)
        plt.savefig(f"./img/Cell_{idx+1}.jpg")
        plt.close()

    ################################ 1-3 ################################
    resourceEffs_1 = np.array(resourceEffs_1)
    resourceEffs_1 = resourceEffs_1.sum(axis=0) / simIt
    aveDataRate_1 = []
    for idx in range(simIt):
        poorestSINRs_1 = dataPoints_1[idx]
        dataRates = []
        for poorestSINR_1 in poorestSINRs_1:
            shnCp = innerCells[0].MSs[0].shannonCapacity(innerCells[0], sinr=poorestSINR_1)
            dataRates.append(shnCp)
        aveDataRate_1.append(dataRates)
    aveDataRate_1 = np.array(aveDataRate_1)
    aveDataRate_1 = aveDataRate_1.sum(axis=0) / simIt
    print(aveDataRate_1)

    ################################ 2-1 ################################
    dataPoints_2 = np.array(dataPoints_2)
    dataPoints_2.sort()
    plt.plot(dataPoints_2, Y)
    plt.show()
    plt.savefig("./img/2-1.jpg")
    plt.close()

    ################################ 2-2 ################################
    resourceEffs_2 = np.array(resourceEffs_2) / simIt
    aveDataRate_2 = sum(
        [innerCells[0].MSs[0].shannonCapacity(innerCells[0], sinr=dataPoints_2[idx]) for idx in range(simIt)]
    ) / simIt
    print(aveDataRate_2)
