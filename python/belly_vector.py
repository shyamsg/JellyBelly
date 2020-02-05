import sys
import numpy as np

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("ERROR: Not enough parameters.\n")
        usage()

    jellyvec = JELLYVECS(sys.argv[1])

    for vec in jellyvec.belly_loopvec():
        print(vec.shape)

    for vec in jellyvec.belly_loadvec():
        print(vec.shape, vec[0])



def belly_readhead(fp):
        fp.seek(0,0)
        header = fp.read(26)
        if len(header) != 26:
            belly_headerr(fp)
            return -1,-1,-1,-1
        #read 7 0 bytes
        for i in header[:7]:
            if i != 0:
                belly_headerr(fp)
                return -1,-1,-1,-1
        #read int smer length
        smerlength = int.from_bytes(header[7:11],
                                    byteorder=sys.byteorder,
                                    signed=False)
        #read 3 0 bytes
        for i in header[11:14]:
            if i != 0:
                belly_headerr(fp)
                return -1,-1,-1,-1
        #read int kmer length
        kmerlength = int.from_bytes(header[14:18],
                                    byteorder=sys.byteorder,
                                    signed=False)
        #read 2 0 bytes
        for i in header[18:20]:
            if i != 0:
                belly_headerr(fp)
                return -1,-1,-1,-1

        if header[20] == 1:
            h = 1
            fmt = np.int32
        else:
            h = 0
            fmt = np.float32
        #read 5 bytes check for  values 5,4,3,2,1
        signature = [74,69,76,76,89]
        for i in range(1,6):
            if header[-i] != signature[-i]:
                belly_headerr(fp)

        sys.stderr.write("INFO:\tkmer length: " + str(kmerlength) + "\n")
        sys.stderr.write("\tspaced kmer length: " + str(smerlength) + "\n")
        sys.stderr.write("\tvector length: " + str(4**smerlength) + "\n")
        sys.stderr.write("\tfile format: " + "raw" if h else "scaled" + "\n")
        return smerlength, kmerlength, 4**smerlength, fmt


def belly_readtail(fp):
    try:
        fp.seek(-17, 2)
    except:
        belly_tailerr(fp)

    tail = fp.read(17)
    if len(tail) != 17:
        belly_tailerr(fp)
    for i in tail[0:4]:
        if i != 0:
            belly_tailerr(fp)

    numsamples = int.from_bytes(tail[4:12], byteorder=sys.byteorder, signed=False)

    signature = [66,69,76,76,89]
    for i in range(0,5):
        if tail[12:][i] != signature[i]:
            belly_tailerr(fp)

    sys.stderr.write("INFO:\t" + str(numsamples) + " sample(s) in jellyfile.\n")
    return numsamples


def belly_headerr(fp):
    sys.stderr.write("ERROR: jellyfile header misformated.\n")
    raise


def belly_tailerr(fp):
    sys.stderr.write("ERROR: jellyfile tail could not be read.\n")
    raise


def belly_usage():
    sys.stderr.write("Usage:\n\tpython bellyparse.py <binfile>\n")
    sys.exit(-1)



class JELLYVECS:
    def __init__(self, filename):
        sys.stderr.write("<<<<>>>>\n")
        sys.stderr.write("binary or text?\n")
        try:
            self.fp = open(filename, "rb")
            self.fp.seek(0,0)
            header = fp.read(7)
            self.fp.close()
            if len(header) != 7:
                belly_headerr(fp)
        #read 7 0 bytes
        self.fformat = "binary"
        for i in header:
            if i != 0:
                sys.stderr.write("\t7 0 bytes not found, assuming text format.\n")
                self.fformat = "text"
        if self.fformat == "binary":
            sys.stderr.write("\t7 0 bytes found, assuming binary format.\n")
            try:
                self.fp = open(filename, "rb")
                self.setvars()
                self.errflag = False
            except:
                sys.stderr.write("ERROR: Could not open input file: " + filename +  " \n")
                self.errflag = True
        sys.stderr.write("<<<<>>>>\n")
        if self.errflag:
            raise Exception


    def setvars(self):
        self.smerlength, self.kmerlength, self.veclength, self.fmt = belly_readhead(self.fp)
        self.numsamples = belly_readtail(self.fp)


    def belly_loopvec(self):
        self.fp.seek(26,0)
        for i in range(self.numsamples):
            buff = self.fp.read(4*self.veclength)
            vecarray = np.frombuffer(buff, dtype=self.fmt)
            yield vecarray


    def belly_loadvec(self, numsamples=100):
        if self.fformat == "text":
            #figure out how to yield numsamples vectors at a time
        elif self.fformat == "binary":
            self.fp.seek(0,2)
            vecdatacap = self.fp.tell() - 26 - 17
            totalvecs = int(vecdatacap / (self.veclength*4))
            if numsamples == -1:
                numsamples = totalvecs
            numcycles = float(vecdatacap/(numsamples*self.veclength*4))
            if numcycles < 1.0:
                numcycles = 1
                datalength = vecdatacap
                shape = (totalvecs,self.veclength)
            else:
                numcycles = int(numcycles)
                datalength = numsamples*self.veclength*4
                shape = (numsamples, self.veclength)

            self.fp.seek(26,0)
            for i in range(numcycles):
                totalvecs -= numsamples
                buff = self.fp.read(datalength)
                vecarray = np.frombuffer(buff, dtype=self.fmt)
                vecarray = np.reshape(vecarray, shape)
                yield vecarray
            if totalvecs > 0:
                datalength = totalvecs*self.veclength*4
                shape = (totalvecs, self.veclength)
                buff = self.fp.read(datalength)
                vecarray = np.frombuffer(buff, dtype=self.fmt)
                vecarray = np.reshape(vecarray, shape)
                yield vecarray


if __name__ == "__main__":
    main()
