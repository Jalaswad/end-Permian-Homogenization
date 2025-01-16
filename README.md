# READ ME

This GitHub repository contains the code used in analysis of all datasets used in Physiologically driven homogenization of marine ecosystems after the end-Permian mass extinction.

by Jood Al Aswad, Justin Penn, Pedro Monarrez, Mohamad Bazzi, Curtis Deutsch and Jonathan Payne

Code written and maintained by Jood Al Aswad, Mohamad Bazzi and Justin Penn

Contact: jalaswad@stanford.edu

**All codes and data used in this paper are provided here for reproducibility.**

Please start with download_pbdb_new.R to download the raw dataset used in this study. If you do not wish to download the data, or you wish to simply look at the data utilized here, then simply open pbdb.csv which is provided here and shows the data used after filtering and cleaning for spatial and taxonomic errors.

The fossil occurrence analyes can be found and reproduced by opening index.html, which is accompanied by the source code index.qmd. All codes and files pertaining to the oceanographic models and ecophysiotype similarity measures are in the "Model Analyses" folder.

Access the study’s .Rdata object, containing all data and results, using: piggyback!

# This R code shows how to access a .Rdata file from a GitHub release using the piggyback package.
```
# 1. First install and load the piggyback package.
# 2. Create a temporary directory and download the default.RData file from the specified GitHub repository release version using the pb_download() function.
# 3. Finally, load the downloaded .Rdata file into the R environment using the load() function.

 install.packages(piggyback)
 require(piggyback)

 Create temporary directory and load .Rdata into R environment.
 pb_download(file = "default.RData",dest = tempdir(),repo = "jalaswad/end-Permian-Homogenization",tag = "v1.0.0")
 load(file = file.path(tempdir(),"default.RData"))

```

