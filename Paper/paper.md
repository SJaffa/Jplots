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

# Abstract

# Introduction

Molecular clouds are the birthplace of stars. In these chaotic environments, the interplay of gravity, pressure, turbulence, chemistry, radiation, heating and cooling, nad feedback from young stars creates a wide variety of structures. We can observe these projected on the sky in two dimensions, and also use chemical tracers to infer a velocity, which is often used as a third dimension to separate structures of different velocities. We can also simulate many of these processes and create theoretical molecular clouds in three dimensions, and then create synthetic ovservations of these clouds to compare with observations. In all of these studies the key is to be able to robustly compare pixelated greyscale images in two or three spatial and velocity dimensions.

``J plots`` is a tool to quantify the shape of 2D or 3D greyscale images using the moment of inertia. Written in Python, this allows users to import their data set from any format Python can understand (including the commonly used FITS format using the ``pyfits`` package) to make it into a numpy array and segment it using a method of their choice into interesting structures (dendrograms are used in the example scripts, relying on the ``astrodendro`` package). These structures can be individually or collectively analysed to reveal trends in the shape of structures with environment, obervational or simulated constraints, or segmentation method. This allows astronomers to quantify these chaotic structures and compare them in a statistical sense.

# The J plots method

The basics of the J plots method, as described in `@Jaffa, S. E.:2018`, segments the grayscale image into 'structures' using dendrograms. For each structure, we calculate the area, $A$ (number of pixels), mass $M$ (sum of pixel values), centre of mass and principle moments of inertia, $I_{1,2}$. We also calculate what the moments would be for a flat cirle of the same mass and area, $I_{0} = AM/4\pi$. The J momets are then defined as

$J_{1,2} = I_{0} + I_{1,2}/I_{0} - I_{1,2}

If the shape is centrally concentrated, such as a collapsing core, both principle moments will be smaller than $I_{0}$ so both J values will be positive. For a hollow ring shape such as bubble blown in the cloud by stellar feedback, both rinciple moments will be greater than $I_{0}$ so both J values will be negative. For elongated shapes such as filaments, which are prevalent in molecular clouds, one moment will be larger and one smaller, so $J_{1}$ will be positive and $J_{2}$ will be negative. This gives us a simple diagnostic of these common shapes and allows us to place quantitative restrictions on shapes that fall between these categories.

![Proof of concept of 2D J plots. Sereral simple shapes commonly identified by eye in molecular clouds are distributed in different regions on the J plot.](proof.pdf)

``J plots`` has been used in 2D to analyse the shape of structures within simulated filaments. `@Clarke, S. D.:2018` examined the J values of sub-filaments (small elongated structures formed inside the main filament identified in 3D and projected into 2D) and fibres (elongated structures inside the main filament identified observationally in projected physical space with velocity as a third dimension, then projected into the spatial 2D) and found that the observtionally detected fibres did not represent the same gas as the 3D spatially identified sub-filaments. Their distributions of J values showed that hte shape of hte structures recovered was changed by the nature of the observations, so observed velocity coherent structures should not be taken to represent physically separate objects.

# J plots 3

The latest release of ``J plots`` now also allows analysis of 3D objects without projection. This is mainly inteded for the analysis of simulations in hteir native space, but can also be used for observations including velocity. In such a case we caution agains the blind interpretation of velocity as another dimension to be treated the same as spatial dimensions. One axis of the moment analysis can be confined to be along hte velocity axis to keep the types f dimensions separate, but this option is left up to the user as combining dimensions may provide some interesting insights, although will need much more thoughtful analysis to draw physically meaningful conclusions. This release is also updated to be Python 3 compatible.

# Figures



# Acknowledgements



# References
