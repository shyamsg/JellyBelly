# JellyBelly's python routines

Here you will find python scripts that load and manipulate JellyBelly's output for downstream analysis. In [example](https://www.example.com) you will find a complete pipeline from sequencing reads to plots.
  
### Dependencies

These scripts assume python 3.4 or later. In addition to the other python dependencies shown in the main README.

# Loading JellyBelly's binary files into numpy arrays

belly_vector.py contains a simple interface to load JellyBelly's binary output files (binfiles) into numpy arrays. JellyBelly's binfiles contain all the information necessary for this.

Start by loading the sys module as well as the routines in belly_vector.py

```python
import sys
# This is safe as all routines start with "belly_*" 
from belly_vector import *
```
Make sure belly_vector.py is in either your python module path or in the same folder you are running python from.

Next, instatiate a JELLYVECS class with the path to the binfile you want to manipulate.

```python
jellyvec = JELLYVECS("binfilename")
```

If successfull, no errors should be raised.

Instances of JELLYVECS have a couple of variables and methods: 

Variables:
        
        JELLYVECS.smerlength: Length of spaced kmer
        
        JELLYVECS.kmerlength: Length of kmer
        
        JELLYVECS.fmt: Data format (scaled or raw)
       
        JELLYVECS.numsamples: Number of samples in binfile

Methods:
        
        JELLYVECS.belly_loopvec(): Generator of vectors one a time.
        
        JELLYVECS.belly_loadvec(numsamples = 100): Generator of vectors in batches (DEFAULT 100 vectors per iteration) 

There are two ways to manipulate JellyBelly's vectors. With JELLYVECS.belly_loopvec and with JELLYVECS.belly_loadvec

1) JELLYVECS.belly_loopvec
This is a generator that outputs a one dimentional vector per iteration.

```python
for vec in jellyvec.belly_loopvec():
        print(vec.shape)
        #Do some work on the vector
```
Vectors can be accessed in only one order unless you store the data somewhere else.

2) JELLYVECS.belly_loadvec
This is a generator that outputs two dimentional vectors in batches of a specific size. If you would like to load more than one vector at a time use this function. Any number of vectors can be loaded. DEAULT: 100

```python
for vec in jellyvec.belly_loadvec(1000):
        print(vec.shape)
        #Do some work with these vectors
```


Future routines such as random vector access will be implemented. 

## Run

  Run these utilities without parameters to print help messages
  
## Dependencies

  python scripts require plotly, numpy, colorlover, umap and sklearn. These can be installed with:
  
    pip install plotly numpy colorlover sklearn umap

  These scripts create html files that you can open and interact with in a web browser.
