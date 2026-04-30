# pen-stack/spec:0.1.0
# Extends pen-stack/biophysics:0.1.0 with:
#   - BWA (Burrows-Wheeler Aligner) for S_Spec off-target genome scan
#   - samtools for SAM/BAM manipulation
#   - pysam Python bindings
#
# Build:
#   docker build -t pen-stack/spec:0.1.0 -t pen-stack/spec:latest \
#       ~/pen-stack/code/repos/pen-score/containers/spec
#
# Requires GRCh38 index on the VM BEFORE running S_Spec:
#   # Download (one-time, ~3.2 GB):
#   cd ~/pen-stack/data/genomes
#   wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/\
# GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/\
# GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#   gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#   ln -s GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fa
#
#   # Index (one-time, ~90 min, produces .amb .ann .bwt .pac .sa):
#   docker run --rm -v ~/pen-stack/data:/data pen-stack/spec:0.1.0 \
#       bwa index /data/genomes/GRCh38.fa
#
# Run S_Spec:
#   docker run --rm \
#       -v ~/pen-stack/data:/data \
#       -v ~/pen-stack/code/repos/pen-score:/pkg \
#       -w /pkg pen-stack/spec:0.1.0 \
#       bash -c "pip install -e '.[spec]' --quiet && python -m pen_score.axes.spec" \
#       2>&1 | tee ~/pen-stack/logs/pen-score/S_Spec_$(date +%Y%m%d).log

FROM pen-stack/biophysics:0.1.0

# System packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    bwa \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# Python packages
RUN pip install --no-cache-dir "pysam==0.22.1" "biopython==1.84"

# Smoke test
RUN bwa 2>&1 | head -1 && samtools --version | head -1 && \
    python -c "import pysam; print('pysam:', pysam.__version__, 'OK')" && \
    echo "spec smoke OK"

CMD ["/bin/bash"]
