#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from scipy import stats

def mrplot(df):
    mlim = [0.5, 3.0]
    rlim = [8, 15.0]

    x1 = df.r1
    y1 = df.m1

    xmin1 = x1.min()
    xmax1 = x1.max()
    ymin1 = y1.min()
    ymax1 = y1.max()

    X1, Y1 = np.mgrid[xmin1:xmax1:200j, ymin1:ymax1:200j]
    positions = np.vstack([X1.ravel(), Y1.ravel()])
    values = np.vstack([x1, y1])
    kernel = stats.gaussian_kde(values)
    Z1 = np.reshape(kernel(positions).T, X1.shape)

    x2 = df.r2
    y2 = df.m2

    xmin2 = x2.min()
    xmax2 = x2.max()
    ymin2 = y2.min()
    ymax2 = y2.max()

    X2, Y2 = np.mgrid[xmin2:xmax2:200j, ymin2:ymax2:200j]
    positions = np.vstack([X2.ravel(), Y2.ravel()])
    values = np.vstack([x2, y2])
    kernel = stats.gaussian_kde(values)
    Z2 = np.reshape(kernel(positions).T, X2.shape)

    fig, ax = plt.subplots()
    # ax.scatter(x1, y1, alpha=0.1, color="red")
    ax.contourf(X1, Y1, Z1, [4.0e-1, Z1.max()], alpha=0.60, colors=["yellow"])
    # ax.scatter(x2, y2, alpha=0.1, color="blue")
    ax.contourf(X2, Y2, Z2, [6.0e-1, Z2.max()], alpha=0.60, colors=["tab:olive"])

    ax.text(9.0, 1.7, "GW170817", color="indigo")

    ax.set_xlabel("R (km)")
    ax.set_ylabel("M (M⊙)")
    ax.set_xlim(rlim)
    ax.set_ylim(mlim)

    plt.savefig("universalrel_mr.png")

def mlambdaplot(df):
    x1 = df.m1
    y1 = df.lambda1
    y1log = np.log10(df.lambda1)

    xmin1 = x1.min()
    xmax1 = x1.max()
    ymin1 = y1log.min()
    ymax1 = y1log.max()

    # TODO: pra plotar tem que ser utilizando a escala log, vou ter que arrumar um jeito
    # de pegar os valores de x e y que correspondem a uma densidade no log e então fazer
    # talvez um poligono com esses ou algo do tipo

    X1, Y1 = np.mgrid[xmin1:xmax1:200j, ymin1:ymax1:200j]
    positions = np.vstack([X1.ravel(), Y1.ravel()])
    values = np.vstack([x1, y1log])
    kernel = stats.gaussian_kde(values)
    Z1 = np.reshape(kernel(positions).T, X1.shape)

    x2 = df.m2
    y2 = df.lambda2
    y2log = np.log10(df.lambda1)

    xmin2 = x2.min()
    xmax2 = x2.max()
    ymin2 = y2log.min()
    ymax2 = y2log.max()

    X2, Y2 = np.mgrid[xmin2:xmax2:200j, ymin2:ymax2:200j]
    positions = np.vstack([X2.ravel(), Y2.ravel()])
    values = np.vstack([x2, y2log])
    kernel = stats.gaussian_kde(values)
    Z2 = np.reshape(kernel(positions).T, X2.shape)

    lambdalim = np.log10([2.0, 2000])
    mlim = [1.0, 3.5]

    fig, ax = plt.subplots()

    print(Z1.min(), Z1.max(), Z1.mean())
    print(Z2.min(), Z2.max(), Z1.mean())

    # ax.scatter(x1, y1log, color="red", alpha=0.1, s=2)
    ax.contourf(X1, Y1, Z1, [0.3, Z1.max()], colors=["yellow"], alpha=0.60)
    # ax.scatter(x2, y2log, color="blue", alpha=0.1, s=2)
    ax.contourf(X2, Y2, Z2, [0.5, Z2.max()], colors=["tab:olive"], alpha=0.60)

    ax.set_yscale("log")
    ax.set_yticks(np.log10([5.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 2000.0]))
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("M (M⊙)")
    ax.set_ylabel(r"$\log \Lambda$")
    ax.set_xlim(mlim)
    # ax.set_ylim(lambdalim)

    plt.savefig("universalrel_lambdam.png")

def main():
    df = pd.read_csv("EoS-insensitive_posterior_samples.dat", sep="\s+", skiprows=1,
                     usecols=[0, 1, 2, 3, 4, 5], names=["m1", "m2", "lambda1", "lambda2", "r1", "r2"])

    mrplot(df)
    mlambdaplot(df)

if __name__ == "__main__":
    main()
