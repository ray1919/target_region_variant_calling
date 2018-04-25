# Target region variant calling
This is pipeline script I use for variant calling in amplicon based sequencing data, mainly gene exon sequencing.

**Importtant:** *I am a novice programmer and I bear no responsibility whatsoever if any of these scripts that I have written wipes your computer, destroys your data, crashes your car, or otherwise causes mayhem and destruction.  USE AT YOUR OWN RISK.*

## Dependency
* FASTQ preprocessor: [fastp](https://github.com/OpenGene/fastp)
* Primers cutter: [cutPrimers](https://github.com/ray1919/cutPrimers)
* reads aligner: [minimap2](https://github.com/lh3/minimap2)
* sam sort: [samtools](http://www.htslib.org/)
* coverage stats: [bamstats04](https://lindenb.github.io/jvarkit/BamStats04.html)
* variant caller: [strelka2](https://github.com/Illumina/strelka)
* Analyze, share, reproduce: [R markdown](https://rmarkdown.rstudio.com/)

## Usage
**Importtant** currently the pipeline is design for pair end data only, single end parameters temporarily not working. pair reads should be longer than 70bp, so PE100, PE150, PE300 should be ok.
### Input files
* **sample_file** one sample id/name per aligner
* **primer_file** interleaved 5p/3p primers of amplicon PCR design
* **gene_file** one gene symbol per line, should be official gene symbol
* **gene_gtf** ensembl gene GTF file downloaded in [ensembl website](http://asia.ensembl.org/info/data/ftp/index.html)
* **fastq_files** stored in **data_dir**, named as [sample id]_R1.fastq.gz, [sample_id]_R2.fastq.gz
* **ref_fasta** reference genome fasta file, must use a relative path format, also should have a minimap2 index file in the same location, with extension mmi. To create the minimap2 index use command ```minimap2 -d ref.mmi ref.fa```

### Parameters
```

```
