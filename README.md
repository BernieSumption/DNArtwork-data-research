# DNArtwork-data-research
R & D utils for the project of generating the list of SNPs that DNArtwork will use.

### Obtain genotype data for relatives

Status: TODO

Get genotype info for one more family groups. Multiple generations ideally.

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


### Validate that all 3 gene data suppliers produce identical results for me

Status: TODO

Program to create list of markers
Filter rsIds not in 23andme
Filter rsIds not in genographic
Filter rsIds not in ancestry
Obtain genotype data for relatives (from 1000genomes project)
Program to generate 23 numbers from file in 1000genomes format
Validate that program produces sane results

Add an item…
