
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the simulations of Gerard et al (2025)

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15784734.svg)](https://doi.org/10.5281/zenodo.15784734)
<!-- badges: end -->

This repo contains the analysis scripts needed to reproduce the
simulations of Gerard et al. (2025). These analyses were evaluated using
R version 4.5.0. The exact versions of all packages used are listed in
the “renv.lock” file.

1.  Edit the `nc` argument in the `Makefile` to be the number of cores
    you wish to use for the analysis.

2.  Open up R and restore the `renv`:

    ``` r
    renv::restore()
    ```

3.  Run `make` in the terminal:

    ``` bash
    make
    ```

4.  Get coffee:

    - [Aslin Coffee DC](https://maps.app.goo.gl/n8vVbjkwwrC9fiyy5)
    - [Doubles](https://maps.app.goo.gl/CXNaN1HpgVxZDk9h6)
    - [Bar Americano](https://maps.app.goo.gl/U6XJmTazJssadUS4A)
    - [BREATHE CO](https://maps.app.goo.gl/CpVTvioWjSbm8zWx5)
    - [Sad Coffee Co.](https://maps.app.goo.gl/KYKTVSi57dWizNTQA)

# References

Gerard, D, Ambrosano, GB, Pereira, GdS, & Garcia, AAF (2025). “Tests for
segregation distortion in higher ploidy F1 populations.” *bioRxiv*,
p. 1–20.
[bioRxiv:2025.06.23.661114](https://doi.org/10.1101/2025.06.23.661114)
