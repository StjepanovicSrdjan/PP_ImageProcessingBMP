from matplotlib import pyplot as plt
import numpy as np
import os
from interpolation import lagrange_interpolation


def getDataByFilterSize(directory):
    A = []
    Bserial = []
    Bparallel = []

    for filename in os.scandir(directory):
        tokens = filename.name.split("-")
        if tokens[0].split("_")[1] == "350":
            f = open(filename, "r")
            lines = f.readlines()
            A.append(float(tokens[1].split("_")[1][:-3]))
            Bserial.append(float(lines[0].split(":")[1].strip()))
            Bparallel.append(float(lines[1].split(":")[1].strip()))

            f.close()

    return A, Bserial, Bparallel


def getDataByCutOff(directory, con):
    A = []
    B = []
    for filename in os.scandir(directory):
        tokens = filename.name.split("-")
        if tokens[1].split("_")[1].strip()[:-4] == con:
            f = open(filename, "r")
            lines = f.readlines()
            A.append(float(tokens[0].split("_")[1].strip()))
            B.append(float(lines[1].split(":")[1].strip()))

            f.close()

    return A, B


def filterSize_plot():
    lists = getDataByFilterSize("../test_data/Image2Prewitt")
    A = np.array(lists[0])
    Bserial = np.array(lists[1])
    Bparallel = np.array(lists[2])

    p_serial = lagrange_interpolation(A, Bserial)
    p_parallel = lagrange_interpolation(A, Bparallel)

    x = np.linspace(np.min(A), np.max(A), 100)
    yS = np.polyval(p_serial, x)
    yP = np.polyval(p_parallel, x)

    fig, ax = plt.subplots()
    ax.plot(x, yS, 'red', label='PrewittSerial')
    ax.plot(x, yP, 'b', label='PrewittParallel')
    leg = ax.legend()
    plt.show()


def cutOffPlot(directory, con):
    lists = getDataByCutOff(directory, con)
    A = np.array(lists[0])
    A = np.sort(A)
    B = np.array(lists[1])
    B = np.sort(B)

    print(A)
    print(B)


    fig, ax = plt.subplots()
    ax.plot(A, B, 'red')
    plt.show()


if __name__ == '__main__':
    filterSize_plot()
    cutOffPlot("../test_data/Prewitt", "5")
    cutOffPlot("../test_data/EdgeDetection", "1")
