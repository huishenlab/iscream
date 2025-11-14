# Manuscript data availabiliy

## Datasets

Two datasets were used to produce the figure in the paper:

1.  snmC-seq2 scWGBS data ([Luo et al. 2018](#ref-luo2018a)) for
    single-cell benchmarks.
2.  bulk WGBS from a project in our lab

Both datasets were aligned with BISCUIT ([Zhou et al.
2024](#ref-zhou2024)). To produce Bismark BEDgraph files ([Krueger and
Andrews 2011](#ref-krueger2011)), we used a python script to convert the
beta and coverage value columns to percentage methylation, unmethylated
reads and methylated read columns.

``` py
import sys
import gzip

bedfile = sys.argv[1]

def bed2cov(line):
    chr, start, end, beta, cov = line.split('\t')
    percent, meth, unmeth = convert(float(beta), int(cov))
    return '\t'.join([chr, start, start, str(percent), str(meth), str(unmeth)])

def convert(beta, cov):
    percent = round(beta * 100)
    meth = round(beta * cov)
    unmeth = cov - meth
    return [percent, meth, unmeth]

with gzip.open(bedfile, 'rt') as bed:
    for line in bed:
        print(bed2cov(line))
```

Using GNU parallel:

``` sh
parallel python bed2cov.py <biscuit_path>/{} '>' <bismark_path>/{/.} ::: biscuit_path/*.bed.gz
```

Then using the `iscream.paper` package at
<https://github.com/huishenlab/iscream.paper> we ran the benchmarks and
produced the figures.

## References

Krueger, Felix, and Simon R. Andrews. 2011. “Bismark: A Flexible Aligner
and Methylation Caller for Bisulfite-Seq Applications.” *Bioinformatics*
27 (11): 1571–72. <https://doi.org/10.1093/bioinformatics/btr167>.

Luo, Chongyuan, Angeline Rivkin, Jingtian Zhou, Justin P. Sandoval,
Laurie Kurihara, Jacinta Lucero, Rosa Castanon, et al. 2018. “Robust
Single-Cell DNA Methylome Profiling with snmC-seq2.” *Nat Commun* 9 (1):
3824. <https://doi.org/10.1038/s41467-018-06355-2>.

Zhou, Wanding, Benjamin K Johnson, Jacob Morrison, Ian Beddows, James
Eapen, Efrat Katsman, Ayush Semwal, et al. 2024. “BISCUIT: An Efficient,
Standards-Compliant Tool Suite for Simultaneous Genetic and Epigenetic
Inference in Bulk and Single-Cell Studies.” *Nucleic Acids Research* 52
(6): gkae097. <https://doi.org/10.1093/nar/gkae097>.
