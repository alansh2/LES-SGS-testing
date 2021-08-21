# LES-SGS-testing
## Turbmat Tools
The contents of **@TurbulenceService** and **Matlab-Fast-SOAP-0.9.1** and scripts `getVelocity.m` and `getVector.m` are part of the John Hopkins Turbmat Tools library. These items are required to run `channel.m`
## Data Handlers
### channel.m
This parallel script collects DNS data at one instance in time specified by variable *t*. The parallel pool size is determined by the number of CPUs specified in the SLURM batch job script. The output is a save file `data.mat` containing variables *U*, *V*, *W*, *X*, *Y*, *Z*.
### ChannelResolvent.m
This serial script applies a Gaussian filter using `GaussianFilter.m` on the data contained in `data.mat` that filters out features with length less than *Delta*. The standard deviation for the Gaussian filter, *Gstd*, is calculated from the user-defined value for *Delta*. The outputs are the save files `properties.mat` containing variables *Delta*, *Gstd* and `resolved.mat` containing variables *Gu1*, *Gu2*, *Gu3*.
### stress.m
This serial script computes the exact subgrid-scale stress tensor. Resolved scales are read from `resolved.mat`. The remaining filtering operations are carried out by `GaussianFilter.m` using the standard deviation *Gstd* from `properties.mat`.
### strain.m
### SGSmodels.m
## Functions
### GaussianFilter.m
This filtering function applies MATLAB's built-in 2D filter to the uniform *x* and *z* axes and uses an iterating convolution for the nonuniform *y* axis.
### Smagorinsky.m
This function calculates eddy viscosity for the Smagorinsky model with additional wall modeling. The Van Driest solution is used for the near wall region with A+ = 26.
### WallAdapting.m
This function calculates eddy viscosity for the WALE model. The current implementation is believed to have bugs.
### closure.m
This function calculates the modeled SGS stress by eddy viscosity closure. The input model functions must have all input arguments self-contained (see `SGSmodels.m` for an example of this implementation). The outputed stress tensors are tagged with an additional field, *Name*.
### TensorPlots.m
This function flexibly plots tensor components of the exact and modeled SGS stresses. The user can ask for any components to be plotted as long as they exist in the input arguments *T* and *mod*. The outputed plots are 2 x 1 (row x col) subplots with the exact above and the modeled below. The subplot titles are dynamic, with the model name coming from the *Name* field in the *mod* input.
