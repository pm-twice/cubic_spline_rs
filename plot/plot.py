# -*- coding: utf-8 -*- 

import numpy as np
import click
from matplotlib import pyplot as plt

def plot_result(org_csv, inter_csv):
    """
    スプライン補間比較用
    """
    org = np.loadtxt(org_csv, delimiter=",", skiprows=1)
    inter = np.loadtxt(inter_csv, delimiter=",", skiprows=1)

    plt.plot(org[:,0], org[:,1], linewidth=0, marker="o", label="src")
    plt.plot(inter[:,0], inter[:,1], label="Vec")
    plt.plot(inter[:,0], inter[:,2], linestyle="dashed", label="Nalgebra")
    plt.plot(inter[:,0], inter[:,3], linestyle="dotted", label="Linear")
    plt.legend()
    plt.show()


@click.command()
@click.option("--org", default="org.csv", help="original csv file")
@click.option("--cmp", default="out.csv", help="interpolated csv file")
def main(org, cmp):
    plot_result(org, cmp)

if __name__ == "__main__":
    main()