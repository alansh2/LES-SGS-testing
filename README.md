# LES-SGS-testing
## Turbmat Tools
The contents of **@TurbulenceService** and **Matlab-Fast-SOAP-0.9.1** and scripts `getVelocity.m` and `getVector.m` are part of the John Hopkins Turbmat Tools library. These items are required to run `channel.m`
## Data Handlers
### channel.m
This parallel script collects DNS data at one instance in time specified by variable *t*. The parallel pool size is determined by the number of CPUs specified in the SLURM batch job script. The output is a save file `data.mat` containing variables *U*, *V*, *W*, *X*, *Y*, *Z*.
### ChannelResolvent.m
This serial script applies a Gaussian filter on the data contained in `data.mat` that filters out features with length less than *Delta*. The standard deviation for the Gaussian filter, *Gstd*, is calculated from the user-defined value for *Delta*. The outputs are the save files `properties.mat` containing variables *Delta*, *Gstd* and `resolved.mat` containing variables *Gu1*, *Gu2*, *Gu3*.
### stress.m
### strain.m
## Functions
### TensorPlots.m
