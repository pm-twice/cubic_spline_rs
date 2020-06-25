# -*- coding: utf-8 -*- 

import numpy as np
import click
from matplotlib import pyplot as plt

def plot_result(bas1_csv, bas2_csv, bas3_csv):
    """
    スプライン補間比較用
    """
    bas1 = np.loadtxt(bas1_csv, delimiter=",", skiprows=0)
    n = bas1.shape[1] -1

    plt.subplot(3,1,1)
    plt.xlim([0,100])
    plt.ylim([0,1])
    for i in range(n):
        plt.plot(bas1[:,0]*100/8, bas1[:,i+1])

    bas2 = np.loadtxt(bas2_csv, delimiter=",", skiprows=0)
    n = bas1.shape[1] -1

    plt.subplot(3,1,2)
    plt.xlim([0,100])
    plt.ylim([0,1])
    for i in range(n):
        plt.plot(bas2[:,0]*100/8, bas2[:,i+1])

    bas3 = np.loadtxt(bas3_csv, delimiter=",", skiprows=0)
    n = bas1.shape[1] -1

    plt.subplot(3,1,3)
    plt.xlim([0,100])
    plt.ylim([0,1])
    for i in range(n):
        plt.plot(bas3[:,0]*100/8, bas3[:,i+1])

    plt.show()


@click.command()
@click.option("--bas1", default="bas1.csv", help="original csv file")
@click.option("--bas2", default="bas2.csv", help="original csv file")
@click.option("--bas3", default="bas3.csv", help="original csv file")
def main(bas1, bas2, bas3):
    plot_result(bas1, bas2, bas3)

if __name__ == "__main__":
    main()