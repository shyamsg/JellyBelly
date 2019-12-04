import sys
import numpy as np


def read_data(filehandle):
    data = np.loadtxt(filehandle,delimiter="\t")
    return data

def main():
    infile = sys.argv[1]
    data = read_data(open(infile,"r"))
    #sys.stderr.write(str(data.shape[0]) + '\t' + str(data.shape[1]) + '\n')
    if len(data.shape) != 2:
        print("ERROR: reading file")
        return -1
    print(np.linalg.norm(data[0]-data[1]))

if __name__ == "__main__":
    main()
