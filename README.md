# Coherent Phi Electroproduction Event Generator

Monte Carlo event generator for coherent exclusive phi-meson electroproduction from helium-4:

**e + ⁴He → e′ + ⁴He′ + φ**, with **φ → K⁺ + K⁻**.

This code supports dissertation studies of the momentum-transfer distribution and diffractive structure of this reaction using CLAS12 and ALERT at Jefferson Lab.

## Author

Md Masud  
New Mexico State University

## Research context

Dissertation: *Exploring the diffractive minima in coherent φ electroproduction off helium-4 using CLAS12 and ALERT at Jefferson Lab.*

The study develops a simulation and reconstruction framework for a future measurement of this coherent channel. Generated events provide input to detector simulation and reconstruction studies. Model-weighted distributions are simulation predictions and should not be interpreted as measured cross sections or an extraction of the nuclear gluon distribution.

## Repository contents

- `generator.cc`: C++ source code for the event generator.
- `README.md`: project overview and usage documentation.

## Software requirements

The generator uses C++ and ROOT. Exact dependencies and the ROOT version used for the shared source should be recorded here after checking `generator.cc`.

Detector simulation and reconstruction are separate steps; this repository description does not imply that GEMC, CLAS12, or ALERT reconstruction software is included.

## Building and running

Build and run instructions are pending verification against the uploaded source. Check the source for its entry point, required input files, configuration parameters, and output filenames before execution.

## Reproducibility

When reporting results produced with this generator, record:

- The Git commit identifying the source version.
- Software versions, beam energy, and generation settings.
- Phase-space limits and selection cuts.
- The event-weight model and normalization convention.
- The random seed and number of generated events.

## Scope

This repository is shared for dissertation committee review and scientific documentation. Detector acceptance, reconstruction efficiency, background subtraction, and cross-section extraction require additional analysis beyond event generation.
