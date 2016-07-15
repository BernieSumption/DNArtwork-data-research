# DNArtwork-data-research
R & D utils for the project of generating the list of SNPs that DNArtwork will use.

## List of tasks / components:

### Allele list

Status: TODO

Make a list of relatively rare, homogeneously geographically distributed SNP alleles.

Includes defining a standard format for representation of the allele list.

### Obtain genotype data for relatives

Status: TODO

Get genotype info for one more family groups. Multipel generations ideally.

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

```
CHS population:

HG00657               HG00656
   |                     |
   _______________________
       |             |
    HG00658       HG00702       HG00701
                     |             |
                     _______________
                            |
                         HG00703
```

Data download links: http://www.1000genomes.org/data-portal/sample/HG00658 **Note:** only the "HD Genotype Chip" results contain samples for the children **Note 2:** It's a massive file, so build a stream processor to split it into separate results files for each individual without unzipping it.

### Algorithm to take an individual's genome and produce 23 numbers

Status: TODO

As part of this, validate that this algorithm, when run on the relative data sets, produces the results expected

### Parsers for the file format of major suppliers

Status: TODO

Write parsers for 23andme, national genographic and ancestry.com.

Data sets [available here](https://my.pgp-hms.org/public_genetic_data?data_type=23andMe).

### Validate that all 3 gene data suppliers produce identical results for me

Status: TODO