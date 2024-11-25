#!/bin/bash

#To merge Salmon quantification files ; make sure all the quant.sf files are in all folders


salmon quantmerge --quants quant/* --column NumReads --output merged_output.tsv


salmon quantmerge --quants salmon_quant/* --column NumReads --output salmon_merged_numreads.tsv