#!/bin/sh

sort markers-list-23andme.txt markers-list-ancestry.txt | uniq -d > tmp.txt

sort tmp.txt markers-list-genographic.txt | uniq -d > markers-list-intersection.txt

rm tmp.txt
