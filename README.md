# Plastisphere resistome case (BIO328)

In this case you analyse the real pipeline output of Witsø et al. (2025), *Plastispheres as
reservoirs of antimicrobial resistance: Insights from metagenomic analyses across aquatic
environments*, PLOS One 20(9): e0330754, <https://doi.org/10.1371/journal.pone.0330754>.

## Getting started

1. Install R (4.3 or newer) and RStudio, then install the packages in R:
   `install.packages(c("tidyverse", "vegan", "rmarkdown"))`.
2. Download this repository (green **Code** button, **Download ZIP**, then unzip) or clone it.
   On Windows, keep the folder path short (for example `C:\Users\<you>\Documents\bio328`):
   Windows cannot read files whose full path is longer than 260 characters.
3. Open `plastpath_case.Rproj`, then `plastpath_case.Rmd`. Set your group number in the
   first chunk and follow the instructions in the document.

## Contents

| | |
|---|---|
| `plastpath_case.Rmd` | the case and your report: tasks, questions and answer boxes; knit it to HTML |
| `R/helpers.R` | helper functions, loaded in the first chunk |
| `metadata/` | sample metadata, group design and drug-class groups |
| `pipeline_output/` | the files written by the pipeline; `pipeline_output/README.md` explains every file and column |
| `docs/` | the paper and the KMA output specification |
| `figures/` | satellite image of the sampling sites |
