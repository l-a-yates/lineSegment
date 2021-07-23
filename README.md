## Spatial pattern analysis of line-segment data in ecology

## About

This repository provides code to reproduce the analyses in the manuscript "Spatial pattern analysis of line-segment data in ecology" by Yates, Luke; Brook, Barry;   Buettel, Jessie C (2021).

## Details

The series of four R scripts 00_lineSegment...R, etc. reproduce the full set of results for the manuscript:

  * 00 defines two functions: one to estimate the summary functions S(r), K(r), and L(r), and a second to generate restricted randomisations. Both functions act on a supplied line-segment pattern.
  * 01 performs the Monte Carlo tests and generates the summary plots (Fig 3).
  * 02 defines and fits all parametric models, and computes the AIC estimates.
  * 03 plots the treefall data (Fig 1) and the example realisations of line-segment models (Fig 2).

The time-consuming computation of the summary function estimates for the Monte Carlo analyses have been saved as an R data object in the folder 'output'. Alternatively, the analyses can be re-run in their entirety.
