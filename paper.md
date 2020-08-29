---
title: 'J3D: using moments to quantify the shape of molecular clouds in 3D'
tags:
  - Python
  - astronomy
  - molecular clouds
  - simulations
  - observations
authors:
  - name: Sarah E. Jaffa
    orcid: 0000-0002-6711-6345
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Anthony P. Whitworth
    affiliation: 2
  - name: Seamus D. Clarke
    orcid: 0000-0001-9751-4603
    affiliation: 3
affiliations:
 - name: Centre for Astrophysics Research, University of Hertfordshire
   index: 1
 - name: Cardiff University
   index: 2
 - name: I. Physikalisches Institut, Universit&auml;t zu K&ouml;ln
   index: 3
date: 11 February 2020
bibliography: paper.bib
---
# Todo

Open source license\\
Acknowledge financial support\\
Installation instructions and dependencies\\
3D example usage for testing\\

Checklist: https://joss.readthedocs.io/en/latest/review_criteria.html

# Summary

Molecular clouds are the birthplace of stars. In these chaotic environments, the interplay of gravity, pressure, turbulence, chemistry, radiation, heating and cooling, and feedback from young stars creates a wide variety of structures. We can observe these projected on the sky in two dimensions (position-position, PP), and, if we use line radiation to infer radial velocity from the Doppler Effect, in three dimensions (position-position-velocity, PPV). Line emission at different radial velocities is often presumed to come from physically separate structures along the line of sight, so used as a proxy for the third spatial dimension. We can also simulate many of these processes and create theoretical molecular clouds in two (PP) or three dimensions (PPP), and then use our understanding of the physics and chemistry of the interstellar medium to create synthetic observations of these clouds (PPV), which can then be compared with real observations. In all of these studies, the key is to be able to compare pixelated greyscale images in two or three spatial and velocity dimensions.

``J plots`` is a tool to quantify the shape of 2D (PP) greyscale images using their moment of inertia. Written in Python, this allows users to import their data set from any format Python can understand (including the commonly used FITS format using the [``pyfits`` package](https://pyfits.readthedocs.io/en/latest/)) to make it into a numpy array and segment it using a method of their choice into meaningful structures (dendrograms are used in the example scripts, relying on the ``astrodendro`` package, but other techniques such as ``clumpfind`` [@clumpfind] or simple thresholding could equally be used). These structures can be individually or collectively analysed to reveal trends in their shape with environment, observational or simulated constraints, or segmentation method. This allows astronomers to quantify these chaotic structures and compare them in a statistical sense.

# The J plots method

The basics of the J plots method, as described in @Jaffa:2018, takes as input set of real-valued pixels representing a part of a greyscale image (called a 'structure'). For each structure, we calculate the area, $A$ (number of pixels), mass $M$ (sum of pixel values), centre of mass and principal moments of inertia, $I_{1...D}, where $D$ is the number of dimensions. We also calculate what the moments would be for a reference shape of the same mass and area, $I_{0}$. The J moments are then defined as
$$J_{i} = \frac{I_{0} + I_{i}}{I_{0} - I_{i}}, \hspace{2cm} i=1, 2, 3.$$

## The two-dimensional case

In two dimensions the reference shape is a filled circle of constant surface density, so $I_{0} = AM/4\pi$. If the shape is centrally concentrated, such as a collapsing core, both principle moments will be smaller than $I_{0}$ so both J values will be positive. For a hollow ring shape such as bubble blown in the cloud by stellar feedback, both principal moments will be greater than $I_{0}$ so both J values will be negative. For elongated shapes such as filaments, which are prevalent in molecular clouds, one moment will be larger and one smaller, so $J_{1}$ will be positive and $J_{2}$ will be negative. This gives us a simple diagnostic of these common shapes and allows us to place quantitative restrictions on shapes that fall between these categories.

![Proof of concept of 2D J plots. The J values of several simple shapes are plotted, representing morphologies observed in molecular clouds. This demonstrates that distinct categories of shape are distributed in different regions on the J plot.](tests/proof.pdf)

``J plots`` has been used in 2D to analyse the shape of structures within simulated filaments. @Clarke:2018 examined the J values of sub-filaments (small elongated structures formed inside the main filament identified in 3D PPP and projected into 2D PP) and fibres (elongated structures inside the main filament identified observationally in PPV, then projected into the PP) and found that the PPV detected fibres did not represent the same gas as the PPP identified sub-filaments. Their distributions of J values showed that the shape of the structures recovered was changed by the nature of the observations, so observed velocity coherent structures should not be taken to represent physically separate objects.

## J3D

``J3D`` now also allows analysis of 3D objects without projection. This is mainly intended for the analysis of simulations in their native space, but can also be used for observations including velocity, as will be discussed later. In PPP space, the reference shape analogous to a uniform surface density circle becomes a uniform volume density sphere. For the moment of inertia, we care about the projection on a uniform density sphere into two dimensions. Therefore in PPP,
$$I_{0} = \frac{2}{5}*M*(\frac{3V}{4\pi})^{2/3}$$,
where M is the mass of the structure (sum of pixel values) and V is the volume (number of pixels).

Again, hollow shapes where the mass is further from the centre than this reference shape will have all negative J values, centrally concentrated shapes will have all positive J values, and shapes that are elongated on one or two dimensions will have one or two positive J values and the others negative. In this case, the delineation between astrophysically relevant shapes is not as clear because most of the shapes we often discuss are defined in observations and therefore 2D, and the prevalence of filaments, oblate or prolate spheroids, and thin sheet-like structures is only discussed theoretically or in simulations. Algorithms to identify these kinds of shapes and studies of their relevance to astrophysical phenomena (such as turbulence, feedback, gravitational collapse, etc) are not well studied. We believe J3D meets an important need for robust quantification of 3D structures in simulations to enable statistical comparison of different studies.


### Velocity dimension 
Real observations of molecular lines can provide information on the velocity along the line of sight of the gas, and this can help us understand the motion of the interstellar medium. Simulations also have velocity information, so can theoretically be analysed in 6-dimensional PPPVVV, or in projected PPV space to be comparable to observations. In such a case we caution against the blind interpretation of velocity as another dimension to be treated the same as spatial dimensions. One axis of the moment analysis can be confined to be along the velocity axis to keep the types of dimensions separate, but this option is left up to the user. Combining dimensions may provide some interesting insights, but will need much more thoughtful analysis to draw physically meaningful conclusions. 

# Concluding remarks

We present ``J3D``, a code to quantify the three-dimensional structure of astrophysical objects in simulations or observations, based on the two-dimensional J plots. Comparing the moment of inertia to a suitable reference shape separated hollow, centrally concentrated, elongated or elliptical shapes which may relate to important physical processes. Quantifying these shapes allows us to analyse and compare them in a statistical sense. The code is available on [GitHub](https://github.com/SJaffa/Jplots). This release is also updated to be Python 3 compatible.

# Acknowledgements
This research made use of astrodendro, a Python package to compute dendrograms of Astronomical data (http://www.dendrograms.org/) and Astropy [@astropy:2013; @astropy:2019]. SEJ gratefully acknowledges support from the STFC grant number ST/R000905/1. APW gratefully acknowledges the support of the consolidated grant ST/K00926/1 from the UK Science and Technology Facilities Council, and of the EU-funded VIALACTEA NetworkFP7-SPACE-607380. SDC acknowledges support from the ERC starting grant no. 679852 ‘RADFEEDBACK’.


# References
