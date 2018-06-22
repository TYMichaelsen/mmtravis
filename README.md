Tools for metatranscriptome analysis
====================================

mmtravis is an R-package to conveniently visualise and analyse NGS
transcriptomics data in different ways.

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

For a quick guide on how to use mmtravis go to the [Get Started](XXX)
page. Detailed documentation of all mmtravis functions can be found at
the [Functions](XXX) page.

Blog posts about mmtravis
-------------------------

Check out the blog posts at <http://albertsenlab.org/> about mmtravis:

-   [Introducing mmtravis: Our metatranscriptome analysis package](XXX)

Releases
--------

To install a specific (older) release of ampvis2 use `"@release"` as
suffix, fx `remotes::install_github("TYMichaelsen/travis@1.0")` to
install the first release (1.0) of mmtravis. The latest stable release
can be installed with
`remotes::install_github("TYMichaelsen/mmtravis@*release")`.
