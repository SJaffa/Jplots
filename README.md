## J plots and J3D

We can use the principal moments of inertia of a shape to
classify astronomically ineresting structures in 2D (J plots) and 3D (J3D).


# J plots (2D)

This code is able to separate centrally concentrated structures
(cores), elongated structures (filaments) and hollow circular
structures (bubbles) from the main population of ‘slightly
irregular blobs’ that make up most astronomical images.

This can be applied to any 2D greyscale image (single
wavelength/tracer or column density).

The main file is jplots.py which uses the example_Params.txt file
to set up analyse the HiGAL tile descibed in the paper.

This reads in the specified parameters file, builds a 
dendrogram, and analyses all the structures. It outputs
an interactive J plot and image of the data, and if you are 
analysing multiple files it creates a plot of all the structures 
from all the files.

A full description of this algorithm, the proof of concept 
tests, and some example astronomical applications are described 
in the paper. A PDF is included in this repository or can be 
found on arXiv at https://arxiv.org/abs/1803.01640

# J3D (3D)

This code is able to separate centrally concentrated structures, 
elongated structures (filaments), hollow structures, and prolate/oblate 
spheroids.

This can be applied to any 3D greyscale data cube (single
wavelength/tracer or column density in PPP or PPV).

It should be pretty simple to run this code and get the plots 
out, but if you want to do anything more complicated please
contact the authors on s.jaffa@herts.ac.uk.

A full description of this algorithm, the proof of concept 
tests, and some example astronomical applications are described 
in the paper which is currently in preparation.