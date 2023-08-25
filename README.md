# Reactor_Sizing_Optimization
# Chemical Reaction Engineering 

This repository contains MATLAB code for reactor sizing in chemical reaction engineering. The code calculates the minimum required volume for different reactor configurations based on given reaction data.

## Introduction

The purpose of this project is to determine the minimum volume of various reactor configurations required to achieve specific conversion and production rates. The provided MATLAB code offers a versatile tool for chemical engineers to explore different reactor types and find optimal operating conditions.

## Code Structure

The code is organized into sections based on different reactor configurations :

- **PFR (Plug Flow Reactor):** Determines the minimum volume for a single PFR by fitting a curve to given input points and calculating the integral to find the area under the curve.

- **CSTR (Continuous Stirred-Tank Reactor):** Calculates the minimum volume for a single CSTR by integrating the inverse rate function between the inlet and outlet concentrations.

- **Two CSTRs in Series:** Calculates the volumes of two CSTRs in series, where the first CSTR is fed with the inlet concentration and the second CSTR is fed with the outlet of the first.

- **PFR and MFR in Series:** Estimates the minimum volume of a PFR and a Mixed Flow Reactor (MFR) in series by finding the point of minima of the fitted inverse rate function.

- **PFR with Recycle:** Determines the minimum volume of a PFR with recycle by iteratively finding two concentrations that result in nearly equal areas under the curve.


## Dependencies

- MATLAB

