
def readheader(header):
    for i in header[:7]:
        if i != 0:
            return -1,-1
    smerlength = int.from_bytes(header[7:11], byteorder=sys.byteorder, signed=False)
    for i in header[11:14]:
        if i != 0:
            return -1,-1
    kmerlength = int.from_bytes(header[14:18], byteorder=sys.byteorder, signed=False)
    for i in header[18:21]:
        if i != 0:
            return -1,-1
    for i in range(5):
        if header[-i] != i:
            return -1,-1
    return smerlength, kmerlength

def readtail(tail):
    for i in tail[0:4]:
        if i != 0:
            return -1
    numsamples = int.from_bytes(tail[4:8], byteorder=sys.byteorder, signed=False)
    for i in tail[8:12]:
        if i != 0:
            return -1
    if tail[-1] != 1:
        return -1
    return numsamples

def readvecs(fp, numsampls, veclength):
    buff = fp.read(4*veclength*numsampls)
    vecarray = np.frombuffer(buff, dtype=np.float32)
    vecarray = np.reshape(vecarray, (numsampls,veclength))
    yield vecarray
