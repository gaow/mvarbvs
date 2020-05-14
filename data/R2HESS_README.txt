R2HESS Documentation
====================

Dependencies
------------
To run R2HESS you will need:

- the Gnu Scientific Library, GSL (if you have run previous versions of R2HESS you should have this installed already)
- gcc version 4.8.x or later
- the Boost Program Options library

If possible, please install Boost to its default location, ideally using the package manager for your version of Linux. The R2HESS configuration script doesn't always work with non-standard Boost installations.

Installation
------------
Hopefully, all you should have to do is:

$ R CMD INSTALL R2HESS_1.99-1.tar.gz

If GSL and/or Boost are in a non-standard location, you may need to specify that in the command line:

$ R CMD INSTALL R2HESS_1.99-1.tar.gz --configure-args='--with-gsl-prefix=/path/to/gsl --with-boost=/path/to/boost'

If you run into any installation problems, please let me know.

Running HESS
------------

1. Data:

The X and Y data sets should be loaded as matrices or similar in R. If you have column headings (for example SNP labels) then the post-processing will use these to label the outputs.

If you have any confounders, i.e. variables that should always be a part of the model, these must be the first columns of the X data.

If you want to do multi-tissue analysis, the tissues should be interleaved. For example, if tissue 1 has columns A1 A2 ... and tissue 2 columns B1 B2 ... then the Y matrix should have columns A1 B1 A2 B2 ... . You will need to do this manually for the moment, but I should have a script somewhere I can send you if it will help.

2. Output directory:

You will need to specify an output directory. All of the input and output data from your run will be stored here. Note that if you provide a relative path, things might break if you change R's working directory.

IMPORTANT WARNING:
- If you re-run an analysis, it will overwrite any existing files. If you wish to keep the old data, you should instead create a new configuration with a different output directory.

3. Create a configuration object:

The configuration object is a list containing various settings. The 'r' parameter is the number of tissues; this can be ignored unless you are doing multi-tissue analysis.

> library(R2HESS)
> config <- r2hess.makeConfig(Xdat, Ydat, output.dir, r=1)

Most of the defaults should be suitable for you, but you can modify various elements of the list. In particular, you might want to set the expected value and standard deviation of the individual model sizes. The default is 2 and 2, but you can set config$Egam (expected value) and config$Sgam (standard deviation) to values that better match your data.

If you have confounding variables, you need to set config$nConfounders to the number of confounders.

You can also modify the total number of iterations, config$nSweep, and the number of burn-in iterations, config$burn.in - this may be useful for short test runs.

4. Run HESS

> r2hess.run(config)

If you have a large data set, the run can take quite some time (hours, even days). To check the progress, look at the file 'log.txt' in the output directory. 

5. Post processing

After you run HESS, the configuration object and the data is saved in the output directory so you can close R and then analyse the output at a later date. Running load("config.Rda") is sufficient to run the post-processing.

There are several post-processing methods included in the R2HESS package. They all take the configuration object as a parameter, and most write output in the form of pdf files in the config's output directory.

r2hess.checkConvergence(config)
r2hess.plotParams(config)
r2hess.plotAssoc(config, threshold=0.5)
out <- r2hess.toptableAssoc(config, minthresh=0.5)
out <- r2hess.toptablePredictors(config, threshold=0.5)

If you have used previous versions of the package, these methods are more or less the same as before. Proper documentation, and hopefully some better functionality, will follow soon in an updated package.
