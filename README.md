# Jplots

#============================================================
SUMMARY
#============================================================


We can separate structures in an image using dendrograms
and then use the principal moment of inertia of these structures to
classify their shapes.

This code is able to separate centrally concentrated structures
(cores), elongated structures (filaments) and hollow circular
structures (bubbles) from the main population of ‘slightly
irregular blobs’ that make up most astronomical images.

This can be applied to any greyscale image (single
wavelength/tracer or column density).


#============================================================
USAGE
#============================================================

The main file is jplots.py
This reads in the specified parameters file (see example), 
builds a dendrogram, and analyses all the structures. It outputs
an interactive J plot and image of the data, and if you are 
analysing multiple files it creates a plot of all the structures 
from all the files.

It should be pretty simple to run this code and get the plots 
out, but if you want to do anything more complicated please
contact the authors on sarah.jaffa@astro.cf.ac.uk as the
code is still a bit of a mess.

Thanks!


#============================================================
PAPER
#============================================================

A full description of this algorithm, the proof of concept 
tests, and some example stronomical applications are described 
in the paper. A PDF is included in this repository or can be 
found on arXiv at https://arxiv.org/abs/1803.01640
