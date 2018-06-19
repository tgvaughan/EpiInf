---
title: Tutorial
layout: doc
---

Tutorial: Basic Analysis
========================

This tutorial is designed to briefly walk new EpiInf users through the
process of inferring epidemiological parameters and prevalence
trajectories from a set of serially-sampled genetic sequence data.

Loading an alignment
--------------------

For this tutorial, we will use a set of simulated sequence data which is
included with the EpiInf package.

Setting up the analysis
-----------------------

Setting up the EpiInf tree prior
--------------------------------

Running the analysis
--------------------

Run the analysis just as you would any other BEAST 2 analysis. That is,

1.  Start BEAST 2.
2.  Select the XML you produced in the previous section from the file
    selection dialog box.

Once BEAST is running, you should see output periodically printed to
standard out (if you're running BEAST from a terminal emulator) or the
output window. The analysis we've set up should take around half an hour
to complete on a modern computer.

Analyzing the results
---------------------

During the analysis results are written to several files which can
usually located in the same directory as the directory containing the
input XML. These are:

1.  The **log** file, which ends in the extension .log and contains
    sampled parameter values,
2.  The **tree** file, which ends in the extension .trees and contains
    sampled trees.
2.  The **traj** file, which ends in the extension .traj and contains
    sampled stochastic population trajectories.

### Parameter posteriors

### Tree posteriors

### Trajectory posteriors

Wrapping up
-----------

This completes the first tutorial. In a future tutorial we will
demonstrate how to use EpiInf to analyze a combination of genetic
and non-genetic data.
