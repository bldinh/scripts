### Scripts
These scripts are for the creation of a reference dataset and imputation of a target dataset against it.

Code is extension of previous work by [Minhui Chen](https://github.com/Minhui-Chen).

### Individual scripts
- **imputation_pipeline.py**
  - Impute a chromosome against a given reference
  - Breaks up target vcf into overlapping 25 Mb chunks (default)
  - After imputation on each chunk, stitches results into one final file
- **phasing_in_chunks.py**
  - Phase a target chromosome in overlapping 25 Mb chunks (default)
  - If a reference vcf is provided, will phase chunks against the reference
  - After phasing on each chunk, stitches results into one final file
- **create_m3vcf_for_ref.slurm**
  - Creates m3vcf and bim from a phased vcf
