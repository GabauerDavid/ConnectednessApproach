# ConnectednessApproach

## Step 1: Install the devtools package

To install a R package, start by installing the devtools package. The best way to do this is from CRAN, by typing:

```r
install.packages("devtools")
```

## Step 2: Install the package of interest from GitHub

Install the package of interest from GitHub using the following code, where you need to remember to list both the author and the name of the package (in GitHub jargon, the package is the repo, which is short for repository). In this example, we are installing the ConnectednessApproach package created by GabauerDavid.

```r
library(devtools)
install_github("GabauerDavid/ConnectednessApproach")
```

## Step 3: Go through tutorial

[ConnectednessApproach Tutorial](https://gabauerdavid.github.io/ConnectednessApproach/#Dynamic_Connectedness_Measures)

## BibTeX Citation

If you use this package in a scientific publication, I would appreciate if you use the following citation:

```
@article{gabauer2022,
  title={Package ‘ConnectednessApproach’},
  author={Gabauer, David},
  year={2022}
}
```
