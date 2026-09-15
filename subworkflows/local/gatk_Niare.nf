// GATK4 Steps from Niare et al. workflow: nf-core modules and placeholders
// (at this point just for record keeping)

// step:gatk SamFormatConverter
include { GATK4_CLEANSAM                  } from '../modules/nf-core/gatk4/cleansam/main'
include { PICARD_SORTSAM                  } from '../modules/nf-core/picard/sortsam/main'
include { GATK4_MARKDUPLICATES            } from '../modules/nf-core/gatk4/markduplicates/main'
// step: gatk DepthOfCoverage
include { PICARD_COLLECTINSERTSIZEMETRICS } from '../modules/nf-core/picard/collectinsertsizemetrics/main'
include { GATK4_HAPLOTYPECALLER           } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT          } from '../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS             } from '../modules/nf-core/gatk4/genotypegvcfs/main'
// step: gatk GatherVCFs
include { 
  GATK4_VARIANTRECALIBRATOR as GATK_VARCALLINDELS      
          } from '../modules/nf-core/gatk4/variantrecalibrator/main'
include { 
  GATK4_APPLYVQSR as GATK4_VQSRINDELS
  } from '../modules/nf-core/gatk4/applyvqsr/main'
include { 
  GATK4_VARIANTRECALIBRATOR as GATK_VARCALLSNPS      
          } from '../modules/nf-core/gatk4/variantrecalibrator/main'
include { 
  GATK4_APPLYVQSR as GATK4_VQSRSNPS
  } from '../modules/nf-core/gatk4/applyvqsr/main'

workflow {
    
}