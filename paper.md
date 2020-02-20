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

Check spelling
reference clumpfind
Reference oter shape-finding codes
Open source license
Restructure code
Acknowledge financial support
Ongoing projects with J3D
Installation instructions and dependencies
3D example usage for testing
Community guidelines: Are there clear guidelines for third parties wishing to 1) Contribute to the software 2) Report issues or problems with the software 3) Seek support

Checklist: https://joss.readthedocs.io/en/latest/review_criteria.html

# Summary

Molecular clouds are the birthplace of stars. In these chaotic environments, the interplay of gravity, pressure, turbulence, chemistry, radiation, heating and cooling, and feedback from young stars creates a wide variety of structures. We can observe these projected on the sky in two dimensions (position-position, PP), and also use chemical tracers to infer a velocity, which is often used as a third dimension to separate structures which are moving together and therefore assumed to be spatially related (position-position-velocity, PPV). We can also simulate many of these processes and create theoretical molecular clouds in two (PP) or three dimensions (PPP), and then use ur understanding of the physics and chemistry of the interstellar medium to create synthetic obvservations of these clouds, which can then be compared with real observations. In all of these studies the key is to be able to compare pixelated greyscale images in two or three spatial and velocity dimensions.

``J plots`` is a tool to quantify the shape of 2D (PP) greyscale images using the moment of inertia. Written in Python, this allows users to import their data set from any format Python can understand (including the commonly used FITS format using the ``pyfits`` package) to make it into a numpy array and segment it using a method of their choice into meaningful structures (dendrograms are used in the example scripts, relying on the ``astrodendro`` package, but other techniques such as ``clumpfind`` or simple thresholding could equally be used). These structures can be individually or collectively analysed to reveal trends in their shape with environment, obervational or simulated constraints, or segmentation method. This allows astronomers to quantify these chaotic structures and compare them in a statistical sense.

# The J plots method

The basics of the J plots method, as described in @Jaffa:2018, takes as input set of real-valued pixels representing a part of a greyscale image (called a 'structure'). For each structure, we calculate the area, $A$ (number of pixels), mass $M$ (sum of pixel values), centre of mass and principle moments of inertia, $I_{1..D}$, where $D$ is the number of dimensions. We also calculate what the moments would be for a reference shape of the same mass and area, $I_{0}$. The J momets are then defined as

$J_{1...D} = I_{0} + I_{1...D}/I_{0} - I_{1...D}$

## The two-dimensional case

In two dimensions the reference shape is a filled circle of constant surface density, so $I_{0} = AM/4\pi$. If the shape is centrally concentrated, such as a collapsing core, both principle moments will be smaller than $I_{0}$ so both J values will be positive. For a hollow ring shape such as bubble blown in the cloud by stellar feedback, both rinciple moments will be greater than $I_{0}$ so both J values will be negative. For elongated shapes such as filaments, which are prevalent in molecular clouds, one moment will be larger and one smaller, so $J_{1}$ will be positive and $J_{2}$ will be negative. This gives us a simple diagnostic of these common shapes and allows us to place quantitative restrictions on shapes that fall between these categories.

\begin{figure}
\centering
\includegraphics[width=0.5\linewidth]{proof}
\caption{Proof of concept of 2D J plots. The J values of several simple shapes are plotted, representing morphologies observed in molecular clouds. This demonstrates that distinct categories of shape are distributed in different regions on the J plot.}
\label{fig:foo}
\end{figure}

![Proof of concept of 2D J plots. The J values of several simple shapes are plotted, representing morphologies observed in molecular clouds. This demonstrates that distinct categories of shape are distributed in different regions on the J plot.](proof.pdf)

``J plots`` has been used in 2D to analyse the shape of structures within simulated filaments. @Clarke:2018 examined the J values of sub-filaments (small elongated structures formed inside the main filament identified in 3D PPP and projected into 2D PP) and fibres (elongated structures inside the main filament identified observationally in PPV, then projected into the PP) and found that the PPV detected fibres did not represent the same gas as the PPP identified sub-filaments. Their distributions of J values showed that the shape of the structures recovered was changed by the nature of the observations, so observed velocity coherent structures should not be taken to represent physically separate objects.

## J3D

The new release of ``J plots`` (called ``J3D`` now also allows analysis of 3D objects without projection. This is mainly inteded for the analysis of simulations in their native space, but can also be used for observations including velocity, as will be discussed later. In PPP space, the reference shape analogous to a uniform surface density circle becomes a uniform volume density sphere. From the point of view of moment of inertia we care about the projection on a uniform density sphere into two dimensions. Therefore in PPP,

$I_{0} = (2/5)*M*(3*V/4*np.pi)^{2/3}$,

where M is the mass of the structure (sum of pixel values) and V is the volume (number of pixels).

Again, hollow shapes where the mas sis further from the centre than this reference shape will have all negative J values, centrally concentrated shapes will have all positive J values, and shapes that are elongated on one or two dimensions will have one or two positive J values and the others negative. In this case the deliniation between astrophysically relevant shapes is not as clear because most of hte shapes we often discuss are defined in obervations and therefore 2D, and the prevalence of filaments, oblate or prolate spheroids, and thin sheet-like structures is only discussed theoretically or in simulations. Algorithms to identify these kind of shapes and studies of their relevance to astrophysical phenomena (such as turbulence, feedback, gravitational collapse, etc) are not well studied. We believe J3D meets an important need for robust quantification of 3D structures in simulations to enable statistical comparrison of different studies.


### Velocity dimension 
Real observations of molecular lines can provide information on the velocity along the line of sight of the emitting gas, and this can be analysed to understand the motion of the interstellar medium. Simulations also have velocity information,  In such a case we caution agains the blind interpretation of velocity as another dimension to be treated the same as spatial dimensions. One axis of the moment analysis can be confined to be along the velocity axis to keep the types f dimensions separate, but this option is left up to the user as combining dimensions may provide some interesting insights, although will need much more thoughtful analysis to draw physically meaningful conclusions. 

# Concluding remarks

This release is also updated to be Python 3 compatible.

# Figures



# Acknowledgements



# References
