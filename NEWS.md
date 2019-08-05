Version 0.2.0:  August 5, 2019
-------------------------------------------------------------------------------

Deprecated:

+ function `analyzeClones` is deprecated. The clonal analysis has been added 
  to the main function `defineClonesScoper` as an argument `analyze_clones`. 
+ the out put class is deprecated. Results would be reported as a list if 
  argument `analyze_clones` set to be true, otherwise a single dataframe is
  returned.
+ `plot_neighborhoods` from clonal analysis has been deprecated.
+ `neighborhoods` from clonal analysis has been deprecated.

General:

+ New models, `hierarchical` for hierarchical-clustering based, and `identical` 
  for clustering among identical junction sequences are added. 
+ New method for spectral-clustering based model has been added through the 
  `vj` in argument `method`.

Clonal analysis:

+ Switched the meaning of the "inter" and "intra" labels in 
  `calculateInterVsIntra` function. Now, "inter" is the label used to form 
  distances that mean between clones, and "intra" is the label used to form 
  distances that mean on the inside, within each clone.
+ Changed the `plotInterVsIntra` output from a density plot to a histogram.
+ Changed the way to calculate the effective threshold. Now, desntiy approach 
  is used.
    

Version 0.1.0:  October 4, 2018
-------------------------------------------------------------------------------

Initial public release.
