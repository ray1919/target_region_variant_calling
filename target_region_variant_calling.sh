#!/usr/bin/env bash

# ##################################################
# My Generic BASH script template
#
version="0.0.5"               # Sets version variable
#
scriptTemplateVersion="1.5.0" # Version of scriptTemplate.sh that this script is based on
#  v1.1.0 -  Added 'debug' option
#  v1.1.1 -  Moved all shared variables to Utils
#           -  Added $PASS variable when -p is passed
#  v1.2.0 -  Added 'checkDependencies' function to ensure needed
#           Bash packages are installed prior to execution
#  v1.3.0 -  Can now pass CLI without an option to $args
#  v1.4.0 -  checkDependencies now checks gems and mac apps via
#           Homebrew cask
#  v1.5.0 - Now has preferred IFS setting
#  - Preset flags now respect true/false
#  - Moved 'safeExit' function into template where it should
#   have been all along.
#
# HISTORY:
#
# * DATE - v1.0.0  - First Creation
#
# ##################################################

# Provide a variable with the location of this script.
scriptPath=/home/zhaorui/ct208/tool/shell-scripts
snpsift=/home/zhaorui/ct208/tool/SnpEff/snpEff/SnpSift.jar
DBSNP=/opt/data/db/snp/gatk/hg38/dbsnp_146.hg38.vcf

# Source Scripting Utilities
# -----------------------------------
# These shared utilities provide many functions which are needed to provide
# the functionality in this boilerplate. This script will fail if they can
# not be found.
# -----------------------------------

utilsLocation="${scriptPath}/lib/utils.sh" # Update this path to find the utilities.

if [ -f "${utilsLocation}" ]; then
  source "${utilsLocation}"
else
  echo "Please find the file util.sh and add a reference to it in this script. Exiting."
  exit 1
fi

# trapCleanup Function
# -----------------------------------
# Any actions that should be taken if the script is prematurely
# exited.  Always call this function at the top of your script.
# -----------------------------------
function trapCleanup() {
  echo ""
  # Delete temp files, if any
  if is_dir "${tmpDir}"; then
    rm -r "${tmpDir}"
  fi
  die "Exit trapped. In function: '${FUNCNAME[*]}'"
}

# safeExit
# -----------------------------------
# Non destructive exit for when script exits naturally.
# Usage: Add this function at the end of every script.
# -----------------------------------
function safeExit() {
  # Delete temp files, if any
  if is_dir "${tmpDir}"; then
    rm -r "${tmpDir}"
  fi
  trap - INT TERM EXIT
  exit
}

# Set Flags
# -----------------------------------
# Flags which can be overridden by user input.
# Default values are below
# -----------------------------------
data_dir=data
quiet=false
printLog=false
verbose=false
force=false
strict=false
debug=false
threads=4
min_length=15
args=()

# Set Temp Directory
# -----------------------------------
# Create temp directory with three random numbers and the process ID
# in the name.  This directory is removed automatically at exit.
# -----------------------------------
tmpDir="/tmp/${scriptName}.$RANDOM.$RANDOM.$RANDOM.$$"
(umask 077 && mkdir "${tmpDir}") || {
  die "Could not create temporary directory! Exiting."
}
outDir="output`date +%s`"

# Logging
# -----------------------------------
# Log is only used when the '-l' flag is set.
#
# To never save a logfile change variable to '/dev/null'
# Save to Desktop use: $HOME/Desktop/${scriptBasename}.log
# Save to standard user log location use: $HOME/Library/Logs/${scriptBasename}.log
# -----------------------------------
# logFile="$HOME/Library/Logs/${scriptBasename}.log"
logFile='/dev/null'

# Check for Dependencies
# -----------------------------------
# Arrays containing package dependencies needed to execute this script.
# The script will fail if dependencies are not installed.  For Mac users,
# most dependencies can be installed automatically using the package
# manager 'Homebrew'.  Mac applications will be installed using
# Homebrew Casks. Ruby and gems via RVM.
# -----------------------------------
homebrewDependencies=()
caskDependencies=()
gemDependencies=()

function mainScript() {
############## Begin Script Here ###################
####################################################

echo -n

##Step0: Input Verification
DIR=`dirname "$0"`

if [[ ! -e "${DIR}/reduced2report.R" ]]; then
    echo "${DIR}/reduced2report.R not found"
    exit 1
fi

if [[ ! -e ${DIR}/forcedGT2report.R ]]; then
    echo "${DIR}/forcedGT2report.R not found"
    exit 1
fi

if [[ -z ${sample_file+x} ]]; then
    usage >&2
    echo ""
fi

if [[ ! -e ${sample_file} ]]; then
    echo 'Invalid SAMPLE NAME: '${sample_file}
    exit 1
fi

if [[ ! -z ${primer_file+x} ]]; then
    if [[ ! -e ${primer_file} ]]; then
        echo 'Invalid primer file: '${primer_file}
        exit 1
    fi
fi

for s in `cat $sample_file`; do
    if [[ ! -e ${data_dir}/${s}_R1.fastq.gz ]]; then
        echo 'Invalid SAMPLE Fastq: '${data_dir}/${s}_R1.fastq.gz
        exit 1
    fi
    if [[ ! -e ${data_dir}/${s}_R2.fastq.gz ]]; then
        echo 'Invalid SAMPLE Fastq: '${data_dir}/${s}_R2.fastq.gz
        exit 1
    fi
done

if [[ -d ${outDir} ]]; then
    if ! ${force} ; then
        echo "outDir ${outDir} already exists."
        echo "existed files would be skipped."
        while true ; do
                read -p "Do you want to preceed? (Y/N): " confirm
                case ${confirm} in
                        Y|y)
                                echo "Confirmed. Initiating..."
                                echo
                                break
                                ;;
                        N|n)
                                echo "Terminated."
                                echo
                                exit 1
                                ;;
                        * )
                                echo "Please enter Y or N."
                                echo
                                ;;
                esac
        done
    fi
fi

if [[ -z ${adapter1+x} ]]; then
    adapter1=auto
fi
if [[ -z ${adapter2+x} ]]; then
    adapter2=${adapter1}
fi

if ${printLog}; then
    logFile=${outDir}/log.txt
fi

echo ###############################################
echo adapter1: $adapter1
echo adapter2: $adapter2
echo ###############################################

##Step1: Check tool version
mkdir -p $outDir
fastp --version | tee ${logFile}
cutPrimers.py -v | tee -a ${logFile}
echo "minimap2 version:" | tee -a ${logFile}
minimap2 --version | tee -a ${logFile}
samtools --version | tee -a ${logFile}
echo "bamstats05 version:" `bamstats05 --version` | tee -a ${logFile}
echo "strelka2 version:" | tee -a ${logFile}
configureStrelkaGermlineWorkflow.py --version | tee -a ${logFile}
Rscript --version | tee -a ${logFile}

##Step2: Raw data qc
for sample in `cat ${sample_file}`; do
    fastp_dir=$outDir/fastp/$sample
    mkdir -p $fastp_dir
    if [[ ! -e ${fastp_dir}/log.txt ]]; then
        fastp -z 3 -w ${threads} \
            -l $min_length \
            -i ${data_dir}/${sample}_R1.fastq.gz \
            -I ${data_dir}/${sample}_R2.fastq.gz \
            -o $fastp_dir/${sample}_R1.fastq.gz \
            -O $fastp_dir/${sample}_R2.fastq.gz \
            -j $fastp_dir/${sample}.json \
            -h $fastp_dir/${sample}.html \
            -R "${sample} QC report" \
            --trim_front1 2 \
            --trim_front2 2 \
            --adapter_sequence ${adapter1} \
            --adapter_sequence_r2 ${adapter2} \
            2>&1 | tee ${fastp_dir}/log.txt
        # if there is GC in both read1 and read2, trim them using following para
    fi
done

##Step3: cut primers
if [[ ! -z ${primer_file+x} ]]; then
    for sample in `cat ${sample_file}`; do
        CUTDIR=$outDir/cutprimers/${sample}
        mkdir -p ${CUTDIR}
        fastp_dir=$outDir/fastp/$sample
        if [[ ! -e ${CUTDIR}/log.txt ]]; then
            cutPrimers.py -r1 $fastp_dir/${sample}_R1.fastq.gz \
                -r2 $fastp_dir/${sample}_R2.fastq.gz \
                -pr ${primer_file} \
                -tr1 $CUTDIR/${sample}_R1.fastq.gz \
                -tr2 $CUTDIR/${sample}_R2.fastq.gz \
                -utr1 $CUTDIR/r1_untrimmed.fastq.gz \
                -utr2 $CUTDIR/r2_untrimmed.fastq.gz \
                --identify-dimers $CUTDIR/dimer.txt \
                -insa $CUTDIR/nsa.txt --error-number 3 \
                --primersStatistics $CUTDIR/ps.txt \
                --min-primer3-length 6 \
                --primer-location-buffer 0 \
                --primer3-absent -rnsa -notrim \
                --threads 5 2>&1 | tee $CUTDIR/log.txt
            #  -notrim -rnsa
        fi
    done
fi

##Step4: alignment, sorting & indexing sam files
if [[ ! -z ${primer_file+x} ]]; then
    for sample in `cat ${sample_file}`; do
        CUTDIR=$outDir/cutprimers/${sample}
        ALNDIR=$outDir/alignment
        mkdir -p ${ALNDIR}
        if [[ ! -e ${ALNDIR}/${sample}.aligned.sorted.bam ]]; then
            STRING=$(head -n 1 < <(zcat $CUTDIR/${sample}_R1.fastq.gz))
            # https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
            # {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
            RGPU=`cut -f3-4,10 -d: <<< $STRING | cut -f1 -d' '`
            # flowcell + lane + library
            RGID=`cut -f3-4 -d: <<< $STRING`.$sample
            #if [[ -z ${vcf_file+x} ]]; then
            #  index_file="${ref_fasta%.*}.mmi"
            #  minimap2 -ax sr \
            #    -R "@RG\tID:$RGID\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}\tPU:$RGPU" \
            #    $index_file \
            #    $CUTDIR/${sample}_R1.fastq.gz \
            #    $CUTDIR/${sample}_R2.fastq.gz \
            #    > $ALNDIR/${sample}.aligned.unsorted.sam
            #else
              bwa mem -t ${threads} \
                -R "@RG\tID:$RGID\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}\tPU:$RGPU" \
                ${ref_fasta} \
                $CUTDIR/${sample}_R1.fastq.gz \
                $CUTDIR/${sample}_R2.fastq.gz \
                > $ALNDIR/${sample}.aligned.unsorted.sam
            #fi
            samtools sort --threads ${threads} -m 2G \
                $ALNDIR/${sample}.aligned.unsorted.sam \
                -o $ALNDIR/${sample}.aligned.sorted.bam
            samtools index $ALNDIR/${sample}.aligned.sorted.bam
            rm -f $ALNDIR/${sample}.aligned.unsorted.sam
            # samtools rmdup --reference ${ref_fasta} $ALNDIR/${sample}.aligned.sorted.bam $ALNDIR/${sample}.dedup.bam
            # samtools index $ALNDIR/${sample}.dedup.bam
        fi
    done
else
    for sample in `cat ${sample_file}`; do
        fastp_dir=$outDir/fastp/$sample
        ALNDIR=$outDir/alignment
        mkdir -p ${ALNDIR}
        if [[ ! -e ${ALNDIR}/${sample}.aligned.sorted.bam ]]; then
            STRING=$(head -n 1 < <(zcat $fastp_dir/${sample}_R1.fastq.gz ))
            RGPU=`cut -f3-4,10 -d: <<< $STRING | cut -f1 -d' '`
            # flowcell + lane + library
            RGID=`cut -f3-4 -d: <<< $STRING`.$sample
            #if [[ -z ${vcf_file+x} ]]; then
            #  index_file="${ref_fasta%.*}.mmi"
            #  minimap2 -ax sr \
            #    -R "@RG\tID:$RGID\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}\tPU:$RGPU" \
            #    $index_file \
            #    $fastp_dir/${sample}_R1.fastq.gz \
            #    $fastp_dir/${sample}_R2.fastq.gz \
            #    > $ALNDIR/${sample}.aligned.unsorted.sam
            #else
              bwa mem -t ${threads} \
                -R "@RG\tID:$RGID\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}\tPU:$RGPU" \
                ${ref_fasta} \
                $fastp_dir/${sample}_R1.fastq.gz \
                $fastp_dir/${sample}_R2.fastq.gz \
                > $ALNDIR/${sample}.aligned.unsorted.sam
            #fi
            samtools sort --threads ${threads} -m 2G \
                $ALNDIR/${sample}.aligned.unsorted.sam \
                -o $ALNDIR/${sample}.aligned.sorted.bam
            samtools index $ALNDIR/${sample}.aligned.sorted.bam
            rm -f $ALNDIR/${sample}.aligned.unsorted.sam
            # samtools rmdup --reference ${ref_fasta} $ALNDIR/${sample}.aligned.sorted.bam $ALNDIR/${sample}.dedup.bam
            # samtools index $ALNDIR/${sample}.dedup.bam
        fi
    done
fi

##Step5: check gene exon coverage (optional)
#if [[ -e ${gene_file} && -e ${gene_gtf} ]]; then
#    # create exon bed file
#    bed_file=${outDir}/exon.bed.gz
#    if [[ ! -e ${bed_file} ]]; then
#        for gene in `cat ${gene_file}`; do
#            echo $gene
#            grep -w "gene_name \"${gene}\"" ${gene_gtf} | grep -w exon |cut -f1,4,5|sort -u >> ${outDir}/exon.bed
#        done
#        sed -i -s "s/^/chr/" ${outDir}/exon.bed
#        bgzip ${outDir}/exon.bed
#        tabix ${bed_file}
#    fi
#fi
if [[ -e ${cov_file} ]]; then
    if [[ ! -e ${outDir}/coverage_MAPQ0.txt ]]; then
        bamstats05 --filter "" \
            -B ${cov_file} \
            -o ${outDir}/coverage_MAPQ0.txt \
            $outDir/alignment/*.sorted.bam \
            2>&1 | tee -a ${logFile}
    fi
    if [[ ! -e ${outDir}/coverage_MAPQ20.txt ]]; then
        bamstats05 --filter "mapqlt(20) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()" \
            -B ${cov_file} \
            -o ${outDir}/coverage_MAPQ20.txt \
            $outDir/alignment/*.sorted.bam \
            2>&1 | tee -a ${logFile}
    fi
fi

##Step6: variant calling
BAMS=($outDir/alignment/*.sorted.bam)
function join_by { local IFS="$1"; shift; echo "$*"; }
BAMSP=(${BAMS[@]/#/--bam=})
SAMPLES=$(join_by ' ' ${BAMSP[@]})

if [[ ! -z ${bed_file+x} ]]; then
    # callRegions="--callRegions ${bed_file} --noCompress ${bed_file}"
    callRegions="--callRegions ${bed_file}"
fi
if [[ ! -z ${int_file+x} ]]; then
    intervalpaths="--intervalpaths ${int_file}"
fi

if [[ ! -z ${vcf_file+x} ]]; then
    forcedGT="--forcedGT ${vcf_file}"
fi

mkdir -p $outDir/pisces
if [[ $(ls $outDir/pisces/*.vcf |wc -l) < $(wc -l $sample_file|cut -f1 -d' ') ]];then
  dotnet ~/ct208/tool/pisces/Pisces_5.2.9.122/Pisces.dll \
    -g /opt/data/db/genome/pisces/hg38_chrom \
    --bam $outDir/alignment \
    --outfolder $outDir/pisces \
    --ploidy somatic \
    --MinimumFrequency 0.001 \
    $intervalpaths \
    --crushvcf false \
    --callmnvs false \
    --NoiseLevelForQModel 30 \
    --MinVariantQScore 20 \
    --minmapquality 0 \
    --minbasecallquality 30
    # --diploidgenotypeparameters 0.01,0.99,0.90 \
fi
#else
#  VARDIR=$outDir/variant
#  if [[ ! -e ${VARDIR}/runWorkflow.py ]]; then
#      cmd="configureStrelkaGermlineWorkflow.py --targeted \
#          --referenceFasta ${ref_fasta} \
#          ${callRegions} ${forcedGT} \
#          ${SAMPLES} \
#          --runDir ${VARDIR} 2>&1 | tee -a ${logFile}"
#      eval ${cmd}
#  fi
#
#  if [[ ! -e ${VARDIR}/log.txt ]]; then
#      ${VARDIR}/runWorkflow.py -j ${threads} -m local 2>&1 | tee ${VARDIR}/log.txt
#  fi
#fi

  echo STOP HERE `date`
  exit
##Step7: vcf annotate & convert
RESDIR=${VARDIR}/results/variants
if [[ -z ${vcf_file+x} ]]; then
    if [[ ! -e ${RESDIR}/annotated.vcf && -e ${RESDIR}/variants.vcf.gz ]]; then
        java -jar $snpsift annotate ${DBSNP} ${RESDIR}/variants.vcf.gz > ${RESDIR}/annotated.vcf
    fi

    if [[ ! -e ${RESDIR}/reduced.vcf && -e ${RESDIR}/annotated.vcf ]]; then
        ( grep -v '^##' ${RESDIR}/annotated.vcf | sed 's/^#//' ) > ${RESDIR}/reduced.vcf
        Rscript ${DIR}/reduced2report.R ${RESDIR}/reduced.vcf
    fi
else
    if [[ ! -e ${RESDIR}/reduced.vcf ]]; then
        ( zgrep -v '^##' ${RESDIR}/variants.vcf.gz | sed 's/^#//' ) > ${RESDIR}/reduced.vcf
        Rscript ${DIR}/forcedGT2report.R ${RESDIR}/reduced.vcf ${vcf_file}
    fi
fi

##Step8: QC report
if [[ ! -z ${primer_file+x} && ! -z ${gene_file+x} ]]; then
    Rscript -e "require( 'rmarkdown' ); rmarkdown::render('${DIR}/qcreport.Rmd',output_dir='${PWD}/${outDir}', params=list(outdir='${PWD}/${outDir}',primer_file='${PWD}/${primer_file}', gene_file='${PWD}/${gene_file}'))"
fi

####################################################
############### End Script Here ####################
}

############## Begin Options and Usage ###################


# Print usage
usage() {
  echo -n "${scriptName} [OPTION]... [FILE]...

This is pipeline script I use for variant calling in amplicon based sequencing data, mainly gene exon sequencing.

 ${bold}Options:${reset}
  -s, --sample_file sample list file, one sample per line
  -p, --primer_file interleaved 5p/3p primers of amplicon PCR design
  -d, --data_dir    fastq data dir, should contain [sample]_R1/2.fastq.gz, default: data
  -a, --adapter1    the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
  -A, --adapter2    the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter1> (string [=])
  -r, --ref_fasta   alignment referece fasta file, must use a relative path format, like ../broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
  -o, --out_dir     output dir name, default is a random name outputXXXX
  -V, --vcf_file    force genotype VCFs
  -b, --bed_file    bed file for calling variant in vcf file, must be bgzip/tabix
  -i, --int_file    interval file for calling variant using Pisces
  -B, --cov_file    bed file for coverage statistics, using bamstat05
  -g, --gene_file   one gene symbol per line, should be official gene symbol
  -L, --min_length  reads shorter than length_required will be discarded, default is 15. (int [=15])
  -t, --threads     threads to use, default: 4
  --force           Skip all user interaction.  Implied 'Yes' to all actions.
  -q, --quiet       Quiet (no output)
  -l, --log         Print log to file
  -s, --strict      Exit script with null variables.  i.e 'set -o nounset'
  -v, --verbose     Output more information. (Items echoed to 'verbose')
  -x, --debug       Runs script in BASH debug mode (set -x)
  -h, --help        Display this help and exit
      --version     Output version information and exit
"
# -f, --gene_gtf    ensembl gene GTF file downloaded in [ensembl website](http://asia.ensembl.org/info/data/ftp/index.html), ungziped file
}

# Iterate over options breaking -ab into -a -b when needed and --foo=bar into
# --foo bar
optstring=h
unset options
while (($#)); do
  case $1 in
    # If option is of type -ab
    -[!-]?*)
      # Loop over each character starting with the second
      for ((i=1; i < ${#1}; i++)); do
        c=${1:i:1}

        # Add current char to options
        options+=("-$c")

        # If option takes a required argument, and it's not the last char make
        # the rest of the string its argument
        if [[ $optstring = *"$c:"* && ${1:i+1} ]]; then
          options+=("${1:i+1}")
          break
        fi
      done
      ;;

    # If option is of type --foo=bar
    --?*=*) options+=("${1%%=*}" "${1#*=}") ;;
    # add --endopts for --
    --) options+=(--endopts) ;;
    # Otherwise, nothing special
    *) options+=("$1") ;;
  esac
  shift
done
set -- "${options[@]}"
unset options

# Print help if no arguments were passed.
# Uncomment to force arguments when invoking the script
# [[ $# -eq 0 ]] && set -- "--help"

# Read the options and set stuff
while [[ $1 = -?* ]]; do
  case $1 in
    -h|--help) usage >&2; safeExit ;;
    --version) echo "$(basename $0) ${version}"; safeExit ;;
    -s|--sample_file) shift; sample_file=${1} ;;
    -p|--primer_file) shift; primer_file=${1} ;;
    -d|--data_dir) shift; data_dir=${1} ;;
    -g|--gene_file) shift; gene_file=${1} ;;
#   -f|--gene_gtf) shift; gene_gtf=${1} ;;
    -r|--ref_fasta) shift; ref_fasta=${1} ;;
    -a|--adapter1) shift; adapter1=${1} ;;
    -A|--adapter2) shift; adapter2=${1} ;;
    -o|--out_dir) shift; outDir=${1} ;;
    -V|--vcf_file) shift; vcf_file=${1} ;;
    -b|--bed_file) shift; bed_file=${1} ;;
    -i|--int_file) shift; int_file=${1} ;;
    -B|--cov_file) shift; cov_file=${1} ;;
    -t|--threads) shift; threads=${1} ;;
    -L|--min_length) shift; min_length=${1} ;;
    -v|--verbose) verbose=true ;;
    -l|--log) printLog=true ;;
    -q|--quiet) quiet=true ;;
    -s|--strict) strict=true;;
    -b|--debug) debug=true;;
    --force) force=true ;;
    --endopts) shift; break ;;
    *) die "invalid option: '$1'." ;;
  esac
  shift
done

# Store the remaining part as arguments.
args+=("$@")

############## End Options and Usage ###################




# ############# ############# #############
# ##       TIME TO RUN THE SCRIPT        ##
# ##                                     ##
# ## You shouldn't need to edit anything ##
# ## beneath this line                   ##
# ##                                     ##
# ############# ############# #############

# Trap bad exits with your cleanup function
trap trapCleanup EXIT INT TERM

# Set IFS to preferred implementation
IFS=$'\n\t'

# Exit on error. Append '||true' when you run the script if you expect an error.
set -o errexit

# Run in debug mode, if set
if ${debug}; then set -x ; fi

# Exit on empty variable
if ${strict}; then set -o nounset ; fi

# Bash will remember & return the highest exitcode in a chain of pipes.
# This way you can catch the error in case mysqldump fails in `mysqldump |gzip`, for example.
set -o pipefail

# Invoke the checkDependenices function to test for Bash packages.  Uncomment if needed.
checkDependencies

# Run your script
mainScript

# Exit cleanlyd
safeExit
