# -*- coding: utf-8 -*- 

import numpy as np
import click
from matplotlib import pyplot as plt

def plot_result(org_csv, inter_csv, bsp_csv):
    """
    スプライン補間比較用
    """
    org = np.loadtxt(org_csv, delimiter=",", skiprows=1)
    inter = np.loadtxt(inter_csv, delimiter=",", skiprows=1)
    bsp = np.loadtxt(bsp_csv, delimiter=",", skiprows=1)

    plt.subplot(2,1,1)
    plt.title("Cubic Spline")
    plt.plot(org[:,0], org[:,1], linewidth=0, marker="o", label="src")
    plt.plot(inter[:,0], inter[:,1], label="Vec")
    plt.plot(inter[:,0], inter[:,2], linestyle="dashed", label="Nalgebra")
    plt.plot(inter[:,0], inter[:,3], linestyle="dotted", label="Linear")
    plt.legend()

    plt.subplot(2,1,2)
    plt.title("B Spline")
    plt.plot(org[:,0], org[:,1], linewidth=0, marker="o", label="src")
    plt.plot(bsp[:,0], bsp[:,1], label="Bspline(1d)")
    plt.plot(bsp[:,2], bsp[:,3], label="Bspline(2d)")
    plt.plot(bsp[:,4], bsp[:,5], label="Bspline(3d)")
    plt.legend()
    plt.show()


@click.command()
@click.option("--org", default="org.csv", help="original csv file")
@click.option("--cmp", default="out.csv", help="interpolated csv file")
@click.option("--bsp", default="bsp.csv", help="b-spline csv file")
def main(org, cmp, bsp):
    plot_result(org, cmp, bsp)

if __name__ == "__main__":
    main()