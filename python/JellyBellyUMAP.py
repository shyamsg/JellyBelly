###############################################################################
# Julian Regalado julian.perez@sund.ku.dk
# Run UMAP on a JellyBelly output file
#
# MIT License
#
#Copyright (c) 2019 Julian Regalado Perez
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
###############################################################################
import sys
import getopt
components=(0,1,2)

def computeUMAP(matrix, dims):
    import umap
    data_fit = umap.UMAP(n_components=dims, n_neighbors=5)
    try:
        model = data_fit.fit(matrix)
    except Exception as err:
        sys.stdout.write("\nERROR:\n\t" + str(err) + "\n")
        #TODO excit gracefully
        sys.exit(-1)
    projection = data_fit.transform(matrix)
    return model, projection


def prepareScatter(model,
                   projection,
                   sample_info,
                   labelby,
                   outprefix,
                   components=components):
    """
    Prepares a plotly instance for figure output. A scatterplot will be created
    based on the PCA coordinates and coloring will be based on labels provided
    by the sample file. ONLY THE FIRST TWO COMPONENTS ARE PLOTED
    Taken from:
        https://plot.ly/ipython-notebooks/principal-component-analysis/
        with minor modifications to adapt for a specific format
    """
    
    import plotly.graph_objs as go
    import numpy as np
    import colorlover as cl
    traces = []
    labels = np.array([i[labelby] for i in sample_info])
    sampleNames = np.array([i[0] for i in sample_info])

    xp = projection[:, components[0]]
    yp = projection[:, components[1]]
    # For Site
    colors = cl.scales['11']['div']['RdBu']

    # loop over labels
    unique_labels = list(sorted(set(labels)))
    print(len(unique_labels),"unique labels.")

    for i in range(len(unique_labels)):
        label = unique_labels[i]
        print(label, colors[i])
        if label == 'NA':
            showlegend = False
        else:
            showlegend = False
        # Create traces
        x1 = projection[labels == label, components[0]]
        y1 = projection[labels == label, components[1]]
        z1 = projection[labels == label, 2]
        trace = go.Scattergl(x=x1,#*scalex,
                             y=y1,#*scaley,
                             #z=z1,
                             mode='markers',
                             text=sampleNames[labels == label],
                             name=label[0:6],
                             marker=dict(size=35,
                                         symbol='circle',
                                         line=dict(color=colors[i],
                                                     width=0.5),
                                         color=colors[i],
                                         opacity=1))
        traces.append(trace)
    data = traces

    layout = go.Layout(title="",
                       font=dict(size=15,),
                       paper_bgcolor='rgba(0,0,0,0)',
                       plot_bgcolor='rgba(0,0,0,0)',
                       titlefont=dict(size=50,
                                      family='arial'),
                       margin=dict(b=70,
                                   t=0,
                                   l=90,
                                   r=10),
                       xaxis=dict(title="",
                                  showline=True,
                                  visible=True,
                                  zeroline=True,
                                  tickfont=dict(size=50)),

                       yaxis=dict(title="",
                                  showline=True,
                                  zeroline=True,
                                  tickfont=dict(size=50)),

                       #shapes=loading_vecs,
                       legend=dict(x=0,
                                   y=-.2,
                                   orientation="h",
                                   bgcolor='#E2E2E2'))
    fig = go.Figure(data=data, layout=layout)
    return fig


def parseInfo(sampleFile):
    info = []
    for line in sampleFile:
        if line[0] == '#':
            continue
        data = line.rstrip().split('\t')
        info.append(data)
    return info


def usage():
    sys.stderr.write("Usage:\n\tpython3 JellyBellyUMAP.py [optons] -f <bellyfile> " + \
                                                    "-s <samplefile> -c <column #>\n")
    sys.stderr.write("\nRequired options:\n\n")
    sys.stderr.write("\t-f <file>\tInput JellyFile. Either a binary JellyFile, "\
                     "or a tab delimited\n\t\t\ttextfile with one vector per line. Format "\
                     "is detected automatically.\n")
    sys.stderr.write("\n\t-s <file>\tInput sample information. Tab delimited text file with"\
                     " sample metadata.\n\t\t\tSamples are stored row-wise in the "\
                     "following format:\n\n\t\t\tSample1ID\tvar1\tvar2\t...\tvarN\n\t\t"\
                     "\tSample2ID \tvar1\tvar2\t...\tvarN\n\t\t\t...\n\t\t\tSampleMID\t"\
                     "var1\tvar2\t...\tvarN\n")
    sys.stderr.write("\n\t-c <int>\tVariable used for plotting samples. Specifies the "\
                     "column in metadata\n\t\t\tfile used to color individual datapoints.\n")
    sys.stderr.write("\nOther options:\n\n")
    sys.stderr.write("\t-n <int>\tNumber of samples in input file to load. "\
                     "(Default: -1 - all)\n")
    sys.stderr.write("\n\t-d <int>\tNumber of dimention UMAP should output. (Default: 10)\n")
    sys.stderr.write("\n")


def readopts(optstr):
    try:
        opts, args = getopt.getopt(optstr, "f:s:c:n:d:h", ["help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        sys.stderr.write("ERROR:\n\t" + str(err) + "\n")
        usage()
        sys.exit(-1)
    inputfile = None
    samplefile = None
    columnnumber = None
    numsamples = ("numsamples", -1)
    dims = ("dims", 10)
    for o, a in opts:
        if o == "-f":
            inputfile = ("inputfile",a)
        elif o == "-s":
            samplefile = ("samplefile",a)
        elif o == "-c":
            try:
                columnnumber = ("column",int(a))
            except ValueError as err:
                sys.stderr.write("ERROR:\n\t" + str(err) + "\n")
                sys.exit(-1)
        elif o == "-n":
            try:
                numsamples = ("numsamples",int(a))
            except ValueError as err:
                sys.stderr.write(str(err) + "\n")
                sys.exit(-1)
        elif o == "-d":
            try:
                dims = ("dims",int(a))
            except ValueError as err:
                sys.stderr.write(str(err) + "\n")
                sys.exit(-1)
        elif o in ("-h", "--help"):
            usage()
            sys.exit(-1)
        else:
            assert False, "unhandled option"

    if inputfile is None or samplefile is None or columnnumber is None:
        print("ERROR:\n\tPlease specify valid parameters for -f, -s, and -c\n")
        usage()
        sys.exit(-1)
    return inputfile, samplefile, columnnumber, numsamples, dims


def jellyload(filename, mode = "bin"):
    from belly_vector import JELLYVECS
    if mode == "bin":
        try:
            jellyvecs = JELLYVECS(filename)
        except:
            return None
        return jellyvecs


def loadmatrix(opts):
    jellyvecs = jellyload(opts["inputfile"])
    if jellyvecs is None:
        print("ERROR")
        sys.exit(-1)
    matrix = next(jellyvecs.belly_loadvec(opts["numsamples"]))
    return matrix


def displayopts(opts):
    sys.stdout.write("Runing JellyBellyUMAP.py\n")
    for key in ["inputfile", "samplefile", "column", "numsamples", "dims"]:
        if key == "inputfile":
            sys.stderr.write("\t" + "-f " + opts[key] + "\n")
        if key == "samplefile":
            sys.stderr.write("\t" + "-s " + opts[key] + "\n")
        if key == "column":
            sys.stderr.write("\t" + "-c " + str(opts[key]) + "\n")
        if key == "numsamples":
            sys.stderr.write("\t" + "-n " + str(opts[key]) + "\n")
        if key == "dims":
            sys.stderr.write("\t" + "-d " + str(opts[key]) + "\n")


def main():
    try:
        assert sys.version_info >= (3, 6)
    except AssertionError:
        sys.stderr.write("PLEASE USE python >= 3.6\n")
        sys.exit(-1)
    opts = readopts(sys.argv[1:])
    opts = dict(opts)
    displayopts(opts)

    #Load jellyfile TODO: Arange for text mode
    matrix = loadmatrix(opts)
    model, projection = computeUMAP(matrix, opts["dims"])

    samplefile = open(opts["samplefile"])
    labelby = int(opts["column"]) - 1
    #sys.stderr.write("Loading data:")
    #sys.stderr.flush()
    #matrix = np.loadtxt(matrixfile)

    print("Matrix has shape: ", matrix.shape)
    sample_info = parseInfo(samplefile)
    # run  automatic PCA
    #model, projection = computeUMAP(matrix)

    fig = prepareScatter(model,
                         projection,
                         sample_info,
                         labelby,
                         "UMAPout")
    from plotly.offline import plot
    plot(fig, filename="UMAPout.html", auto_open=False, image='svg',
         image_width=1600, image_height=1600, output_type='file')


if __name__ == "__main__":
    main()
