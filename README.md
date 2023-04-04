# CosmogenicTadpole
A landscape evolution model with cosmogenic nuclide conservation based on Ferrier and Perron (2020).

The 2D, MATLAB-based model computes transient changes in topography, soil thickness, soil mineral concentrations, and soil cosmogenic nuclide concentrations stemming
changes in rock uplift rate ("tectonic") or hillslope process efficiencies ("climatic"). 

1) Generate initial topography grid in another landscape evolution model like Tadpole (Richardson et al. 2020) [link]. This version of the model does not include dynamic
channels (e.g., stream-power incision). Extract channel network grid beginning at selected drainage area (e.g., 5000 m2).

2) Generate grids of soil thickness (H), mineral concentrations (X), and soil and bedrock nuclide concentrations (Ns & Nzb) using startup.m. 

3) Alter admwx.m parameters.

4) Run admwx()
