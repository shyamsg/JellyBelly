###############################################################################
# Julian Regalado - julian.regalado@tuebingen.mpg.de
# PCA on a matrix
#
# Create PCA plot of an input matrix of dimensions (nxm) == (samplesxfeatures)
###############################################################################
import sys
#from sklearn.decomposition import PCA
import umap
from plotly.offline import plot
import plotly.graph_objs as go
import numpy as np
import colorlover as cl
components=(0,1,2)

def computeUMAP(matrix):
    data_fit = umap.UMAP(n_components=10, n_neighbors=5)
    model = data_fit.fit(matrix)
    projection = data_fit.transform(matrix)
    print("^^", projection.shape)
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


def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Not enough parameters\n")
        sys.exit(1)
    matrixfile = sys.stdin
    samplefile = open(sys.argv[1],'r')
    lableby = int(sys.argv[2])

    # prepare data
    sys.stderr.write("Loading data:")
    sys.stderr.flush()
    matrix = np.loadtxt(matrixfile)

    print("Matrix has shape: ", matrix.shape)
    sample_info = parseInfo(samplefile)
    # run  automatic PCA
    model, projection = computeUMAP(matrix)

    fig = prepareScatter(model,
                         projection,
                         sample_info,
                         lableby,
                         "UMAPout")
    plot(fig, filename="UMAPout.html", auto_open=False, image='svg',
         image_width=1600, image_height=1600, output_type='file')


if __name__ == "__main__":
    main()
