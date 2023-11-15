# <img src="inst/logo/logo.png" alt="CRISRPball" height="200"> ASOspotter

So you wanna look at anti-sense oligonucleotide target mutations, eh? Well, you've come to the right place.

## Installation

**ASOspotter** is available on Github and can be installed as follows:

```r
devtools::install_github("j-andrews7/ASOspotter")
```

## Usage

The app just requires a named list of BAM files and corresponding VCF file.

It will create a "tracks" directory in the working directory for internal use for data subsetting. 
Currently, this directory cannot be changed and is not deleted after the app is closed due to the internals of `shinyIGV`.

```r
library("ASOspotter")

# Using built-in example data.
bam <- system.file("extdata", "NA12878.sub.sorted.bam", package = "ASOspotter")
vcf <- system.file("extdata", "NA12878_HG001.hg38.benchmark.subset.vcf.gz", package = "ASOspotter")

bams <- list("NA12878" = bam, "NA12878_again" = bam)

ASOspotter(vcf = vcf, bams = bams, genome = "hg38")
```

The app is composed of an IGV.js instance and a table of variants located below the tracks. 
Clicking on a variant in the table will jump to that variant in the IGV.js instance. 
Clicking on a variant in the IGV.js instance will display additional info on the variant. 

The wing size can be adjusted to increase or decrease the size of the locus displayed in the IGV.js instance when a variant is clicked. 
The default is 50 bp on either side of the variant. 
Users are free to zoom in/out and scroll as they please.

## Adding Other Tracks

Optionally, users can also provide a named list of BED files that will be displayed.

```r
library("ASOspotter")

# Using built-in example data.
bam <- system.file("extdata", "NA12878.sub.sorted.bam", package = "ASOspotter")
vcf <- system.file("extdata", "NA12878_HG001.hg38.benchmark.subset.vcf.gz", package = "ASOspotter")
bed <- system.file("extdata", "hg38_cgi.bed", package = "ASOspotter")

beds <- list("CpG Islands" = bed)
bams <- list("NA12878" = bam, "NA12878_again" = bam)

ASOspotter(vcf = vcf, bams = bams, beds = beds, genome = "hg38")
```