# Scripts
These scripts are for the creation of a reference dataset and imputation of a target dataset against it. Code uses some previous work by @mhc.

# Individual scripts
- **imputation_pipeline.py**
  - Impute a chromosome against a given reference
  - Breaks up target vcf into 25 Mb chunks (default, can be changed)
  - After imputation, stitches results into one final vcf
- **phasing_in_chunks.py**
  - Phase a target chromosome in chunks
  - If a reference is provided, will phase chunks against the reference
  - Stitches chunks into one final vcf for chromosome
- **create_m3vcf_for_ref.slurm**
  - Creates m3vcf and bim from a phased vcf
