#this program compares the values for the binding energies in two files and prints any disagreements to the screen

import numpy as np
import sys

def main():
    file1 = np.loadtxt('6grainbe_Ilse.inp')
    file2 = np.loadtxt('6grainbe_CNOiso.inp')

    print file1.shape
    sys.exit()


main()
