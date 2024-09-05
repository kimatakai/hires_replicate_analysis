# The code and shell script for NGS analysis


This is NGS analysis pipeline from *_r1.fastq and *_r2.fastq of SRA to scHi-C and RNA-seq.

GEO : GSE223917
GSM : GSM6998595
SRA : SRR22522629	


## Download raw NGS data from Sequence Reads Archive

### Firstly, install SRA-toolkit.

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzvf sratoolkit.current-ubuntu64.tar.gz
vim ~/.bashrc
# add export PATH=$PATH:/home/{user.name}/sratoolkit.3.1.1-ubuntu64/bin
source ~/.bashrc
```

### Download raw fastq data from SRA.

```bash
fastq-dump --split-files SRR22522629
```


## Split raw fastq data to DNA fastq (without adapter seq) and RNA fastq (including adapter seq)

### Install cutadapt library.

```bash
pip isntall cutadapt
cutadapt --version
```

### Split raw data.

-G : determine adapter sequence. 'o=20' means minimum overlap.
-j : number of cpu cores

```bash
cutadapt -G "XGGTTGAGGTAGTATTGCGCAATG;o=20" -j 8 \
--untrimmed-output raw-data/SRR22522629.dna.r1.fastq \
--untrimmed-paired-output raw-data/SRR22522629.dna.r2.fastq \
-o raw-data/SRR22522629.rna.r1.fastq \
-p raw-data/SRR22522629.rna.r2.fastq \
raw-data/SRR22522629_1.fastq raw-data/SRR22522629_2.fastq
```

### Minimize genome contamination in RNA data, probably from exonuclease activity of reverse transcriptase.

```bash
cutadapt --action=none --discard-untrimmed \
-G "XNNNNNNNNTTTTTTTTTTTTTTT;o=18" -j 8 \
-o raw-data/SRR22522629.rna.clean.r1.fastq \
-p raw-data/SRR22522629.rna.clean.r2.fastq \
raw-data/SRR22522629.rna.r1.fastq raw-data/SRR22522629.rna.r2.fastq
```


## Analysis for RNA-seq

### Install UMI-tools.

```bash
pip install umi_tools
```

### Extract UMI (Unique Molecular Identifier) from cleaned rna fastq file.

using NNNNNNNN pattern (eight base pair UMI)

```bash
gzip raw-data/SRR22522629.rna.clean.r1.fastq
gzip raw-data/SRR22522629.rna.clean.r2.fastq
umi_tools extract -p "NNNNNNNN" \
-I raw-data/SRR22522629.rna.clean.r2.fastq.gz \
-S raw-data/umi.SRR22522629.rna.r2.fastq.gz \
--read2-in=raw-data/SRR22522629.rna.clean.r1.fastq.gz \
--read2-out=raw-data/umi.SRR22522629.rna.r1.fastq.gz
gunzip raw-data/umi.SRR22522629.rna.r1.fastq.gz
cp raw-data/umi.SRR22522629.rna.r1.fastq raw-data/umibycell.SRR22522629.rna.r1.fastq
cp raw-data/umibycell.SRR22522629.rna.r1.fastq raw-data/rnaAll.fastq
```

### Clean RNA data preprocessing.

```bash
cutadapt -a CTGTCTCTTATA raw-data/rnaAll.fastq -j 8 | sed 'N;N;N;/\\n\\n/d' > raw-data/rnaAll.clean.fastq
```

### Download reference genome (GRCm38) and GTF (gene feature format) file.

```bash
wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
mv Mus_musculus.GRCm38.dna.primary_assembly.fa.gz raw-data/
gunzip raw-data/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.chr.gtf.gz
mv Mus_musculus.GRCm38.99.chr.gtf.gz raw-data/
gunzip raw-data/Mus_musculus.GRCm38.99.chr.gtf.gz
```

### Install STAR library.

```bash
conda install -c bioconda -y star
STAR --version
STAR --help
```

### Create STAR index.

Create mapping index from reference genome. This processing takes long time.

```bash
mkdir star_index
STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir star_index \
--genomeFastaFiles raw-data/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--limitGenomeGenerateRAM 340000000000 \
--sjdbGTFfile raw-data/Mus_musculus.GRCm38.99.chr.gtf
```

### STAR mapping.

```bash
mkdir starOut
STAR --runThreadN 10 \
--genomeDir star_index \
--readFilesIn raw-data/rnaAll.clean.fastq \
--outFileNamePrefix starOut/star. \
--outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes NH HI NM MD
```

### Install featureCounts

"subread" includes featureCounts. 

```bash
conda install -c bioconda subread -y
featureCounts # show help
```

### Install samtools.

```bash
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -jxvf samtools-1.9.tar.bz2
cd samtools-1.9/
./configure
make
make install prefix=$HOME/.local
vim ~/.bashrc
# add export PATH=$PATH:/home/tynkkai/samtools-1.9
source ~/.bashrc
```


### Count UMI which represents gene expression.

```bash
mkdir rna-seq

# type gene
featureCounts -a raw-data/Mus_musculus.GRCm38.99.chr.gtf \
-o raw-data/gene_assigned \
-R BAM starOut/star.Aligned.sortedByCoord.out.bam \
-T 10 \
-Q 30 \
-t gene \
-g gene_name \
-M -O --fracti

samtools sort raw-data/star.Aligned.sortedByCoord.out.bam.featureCounts.gene.bam -o raw-data/samsort.gene.bam
samtools index raw-data/samsort.gene.bam

mkdir rna-seq
umi_tools count --per-gene --gene-tag=XT --wide-format-cell-counts --assigned-status-tag=XS \
-I raw-data/samsort.gene.bam \
-S rna-seq/counts.gene.SRR22522629.tsv

# type exon
featureCounts -a raw-data/Mus_musculus.GRCm38.99.chr.gtf \
-o raw-data/gene_assigned \
-R BAM starOut/star.Aligned.sortedByCoord.out.bam \
-T 10 \
-Q 30 \
-t exon \
-g gene_name \
-M -O --fracti

samtools sort raw-data/star.Aligned.sortedByCoord.out.bam.featureCounts.bam -o raw-data/samsort.exon.bam
samtools index raw-data/samsort.exon.bam

mkdir rna-seq
umi_tools count --per-gene --gene-tag=XT --wide-format-cell-counts --assigned-status-tag=XS \
-I raw-data/samsort.exon.bam \
-S rna-seq/counts.exon.SRR22522629.tsv
```

### Convert count format

fraction.py
```python
import sys
import pandas as pd

# check arg
if len(sys.argv) != 3:
    sys.exit("Usage: python fraction.py input.tsv output.tsv")

# input and output file path
input_file = sys.argv[1]
output_file = sys.argv[2]

# read tsv matrix
count_matrix = pd.read_csv(input_file, sep='\t')

if count_matrix.shape[1] < 2:
    count_matrix.to_csv(output_file, sep='\t', index=False)
else:
    # count 
    count_matrix['overlapTimes'] = count_matrix['gene'].str.count(",") + 1

    # normalixation
    dat_normed = count_matrix.drop(columns=['gene']).apply(lambda x: x / count_matrix['overlapTimes'])

    # add 
    dat_normed.insert(0, 'gene', count_matrix['gene'])

    dat_normed = dat_normed.assign(gene=dat_normed['gene'].str.split(",")).explode('gene')

    dat_normed = dat_normed.groupby('gene').sum().reset_index()

    dat_normed.to_csv(output_file, sep='\t', index=False)
```

```bash
# type gene
python3 fraction.py rna-seq/counts.gene.SRR22522629.tsv rna-seq/counts.gene.SRR22522629.format.tsv

# type exon
python3 fraction.py rna-seq/counts.exon.SRR22522629.tsv rna-seq/counts.exon.SRR22522629.format.tsv
```

