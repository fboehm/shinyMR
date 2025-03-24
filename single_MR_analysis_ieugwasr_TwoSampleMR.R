##############################################
## Example: Single MR Analysis using 
##          ieugwasr + TwoSampleMR
##############################################

## 1. Load required packages (OpenGWAS token must be set)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)

## 2. Define exposure and outcome study IDs
exposure_id <- "ieu-a-2"   # Example: BMI GWAS
outcome_id  <- "ieu-a-7"   # Example: CHD/CAD GWAS

## 3. Extract top SNPs from the exposure GWAS
#    - pval=5e-8, clump=0: disable strict LD clumping
#    - You can use clump=1 for LD clumping if needed
exposure_snps <- tophits(
  id = exposure_id, 
  pval = 5e-8, 
  clump = 0
)

# Sort by p-value and keep top 50 hits (adjustable)
exposure_snps <- exposure_snps[order(exposure_snps$p), ]
exposure_snps <- head(exposure_snps, 50)

## 4. Get beta and SE for the same SNPs from the outcome GWAS
candidate_snps <- unique(exposure_snps$rsid)  # Remove duplicates
outcome_snps <- associations(
  variants = candidate_snps,
  id = outcome_id,
  proxies = 0
)

## 5. Format both datasets for TwoSampleMR
#    (Column names must match the expected format, e.g., beta.exposure)
exposure_snps_for_mr <- exposure_snps %>%
  rename(
    SNP                    = rsid,
    beta.exposure          = beta,
    se.exposure            = se,
    effect_allele.exposure = ea,
    other_allele.exposure  = nea,
    eaf.exposure           = eaf,
    pval.exposure          = p,
    samplesize.exposure    = n
  ) %>%
  mutate(id.exposure = exposure_id)

outcome_snps_for_mr <- outcome_snps %>%
  rename(
    SNP                    = rsid,
    beta.outcome           = beta,
    se.outcome             = se,
    effect_allele.outcome  = ea,
    other_allele.outcome   = nea,
    eaf.outcome            = eaf,
    pval.outcome           = p,
    samplesize.outcome     = n
  ) %>%
  mutate(id.outcome = outcome_id)

## 6. Harmonisation (align alleles, remove ambiguous or mismatched SNPs)
harmonised_data <- harmonise_data(
  exposure_dat = exposure_snps_for_mr,
  outcome_dat  = outcome_snps_for_mr
)

## 7. Run MR analysis
mr_result <- mr(harmonised_data)

## 8. Display results
print(mr_result)

