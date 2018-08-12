Tools for metatranscriptome analysis
====================================

mmtravis is an R-package to conveniently visualise and analyse NGS
transcriptomics data in different ways. **NOTE: This is still in
development, so don't expect everything to be working perfectly - yet!**

Installing mmtravis
-------------------

First, install [R (3.4.1 or later)](https://mirrors.dotsrc.org/cran/)
and
[RStudio](https://www.rstudio.com/products/rstudio/download/#download).
Windows users should also install
[RTools](https://mirrors.dotsrc.org/cran/bin/windows/Rtools/). Then open
RStudio as administrator (!) and run the commands below to install
mmtravis from the console:

    install.packages("remotes")
    remotes::install_github("TYMichaelsen/mmtravis")

Get started
-----------

Keep an eye out for a get started page, Comming soon! For now, go to the
[References](https://tymichaelsen.github.io/mmtravis/reference/index.html)
page to see a list of available functions.

Blog posts about mmtravis
-------------------------

Keep an eye out for upcomming blogposts at <http://albertsenlab.org/>
about mmtravis

Releases
--------

To install a specific (older) release of mmtravis use `"@release"` as
suffix, fx `remotes::install_github("TYMichaelsen/mmtravis@1.0")` to
install the first release (1.0) of mmtravis. The latest stable release
can be installed with
`remotes::install_github("TYMichaelsen/mmtravis@*release")`.
