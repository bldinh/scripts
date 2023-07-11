### Scripts
These scripts are for the imputation of a target dataset against a reference dataset. Imputation (with Minimac4) is performed in chunks in parallel before merging together. Additional scripts include a phasing script (using Eagle) and a shell script to convert the reference data to expected formats for the imputation script (using bcftools and Minimac3).

Credit to [Minhui Chen](https://github.com/Minhui-Chen) whose previous work some of this code is based on/extended from.

### Individual scripts
- **imputation_pipeline.py**
  - Impute a chromosome against a given reference
  - Expect target VCF to be QC'd prior to running
  - Breaks up target vcf into overlapping 25 Mb chunks (default)
  - After imputation on each chunk, stitches results into one final file


  | Argument | Description |
  | --- | --- |
  | `--tar` | Prefix for data to be imputed (in PLINK format). |
  | `--ref` | Phased vcf of reference panel. |
  | `--refm3` | Phased m3vcf of reference panel. |
  | `--chrom` | Chromosome to use (e.g. "chr22" for hg38). |
  | `--outdir` | Directory to write create and output to. |
  | `--workers` | Number of simultaneous threads. |
  | `--mem` |  Number of gb of RAM to limit to (e.g. 30 for 30gb). |
  | `--windowsize` | Window size for Minimac4 in bp (default: 500000). |
  | `--chunksize` | Chunk size to break chromosome into in bp (default: 25000000). |
  | `--overlap` | Overlap size of chunks in bp (default: 5000000). |
  | `--multicore` | Number of cores to use for phasing/imputation (default: 5). Suggest a factor of the number of `--workers`. |

- **phasing_in_chunks.py**
  - Phase a target chromosome in overlapping 25 Mb chunks (default)
  - If a reference vcf is provided, will phase chunks against the reference
  - After phasing on each chunk, stitches results into one final file


  | Argument | Description |
  | --- | --- |
  | `--tar` | VCF file to be phased. Chrom format should match ref. chrom format if phasing vs. reference. |
  | `--chrom` | Should be formatted the same way "chrom" field is in VCF, e.g. "chr22" or "22". |
  | `--outdir` | Output directory to use. Do not include trailing "/"! Directory will be created if it does not exist. |
  | `--workers` | Number of chunks to work on at same time (i.e. threads). |
  | `--ref` | Phased reference VCF. Not required. If not provided, target VCF will be phased w/o reference. |
  | `--chunksize` | Chunk size in bp (default: 25e6). |
  | `--overlap` | Overlap size in bp (default: 5e6). |
  | `--multicore` | Number of workers to use for when multithreading (default: 5). |

- **create_m3vcf_for_ref.slurm**
  - Creates m3vcf and bim from a phased vcf
