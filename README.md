# purple

performs purity and ploidy estimation

## Overview

## Dependencies

* [PURPLE](https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md)


## Usage

### Cromwell
```
java -jar cromwell.jar run purple.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumour_bam`|File|Input tumor file (bam).
`tumour_bai`|File|Input tumor file index (bai).
`normal_bam`|File|Input normal file (bam).
`normal_bai`|File|Input normal file index (bai).
`SV_vcf`|File|somatic Structural Variant File (.vcf) from gridss.
`smalls_vcf`|File|somatic small (SNV+indel) Variant File (.vcf) [tested with mutect2].
`normal_name`|String|Name of normal sample
`tumour_name`|String|Name of tumour sample


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`amber.amberScript`|String|"java -Xmx32G -cp $HMFTOOLS_ROOT/amber.jar com.hartwig.hmftools.amber.AmberApplication"|location of AMBER script
`amber.PON`|String|"$HMFTOOLS_DATA_ROOT/GermlineHetPon.38.vcf.gz"|Panel of Normal (PON) file, generated for AMBER
`amber.genomeVersion`|String|"V38"|genome version for AMBER, default set to V38
`amber.modules`|String|"argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"|Required environment modules
`amber.threads`|Int|8|Requested CPU threads
`amber.memory`|Int|32|Memory allocated for this job (GB)
`amber.timeout`|Int|100|Hours before task timeout
`cobalt.colbaltScript`|String|"java -Xmx8G -cp $HMFTOOLS_ROOT/cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication"|location of COBALT script
`cobalt.gcProfile`|String|"$HMFTOOLS_DATA_ROOT/GC_profile.1000bp.38.cnp"|GC profile, generated for COBALT
`cobalt.gamma`|String|100|gamma (penalty) value for segmenting
`cobalt.modules`|String|"argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"|Required environment modules
`cobalt.threads`|Int|8|Requested CPU threads
`cobalt.memory`|Int|32|Memory allocated for this job (GB)
`cobalt.timeout`|Int|100|Hours before task timeout
`filterSV.bcftoolsScript`|String|"$BCFTOOLS_ROOT/bin/bcftools view"|location of BCFTOOLS script, including view command
`filterSV.modules`|String|"bcftools"|Required environment modules
`filterSV.threads`|Int|8|Requested CPU threads
`filterSV.memory`|Int|32|Memory allocated for this job (GB)
`filterSV.timeout`|Int|100|Hours before task timeout
`filterSMALL.bcftoolsScript`|String|"$BCFTOOLS_ROOT/bin/bcftools view"|location of BCFTOOLS script, including view command
`filterSMALL.modules`|String|"bcftools"|Required environment modules
`filterSMALL.threads`|Int|8|Requested CPU threads
`filterSMALL.memory`|Int|32|Memory allocated for this job (GB)
`filterSMALL.timeout`|Int|100|Hours before task timeout
`runPURPLE.ensemblDir`|String|"$HMFTOOLS_DATA_ROOT/ensembl"|Directory of Ensembl data for PURPLE
`runPURPLE.refFasta`|String|"$HMFTOOLS_DATA_ROOT/hg38_random.fa"|fasta of reference genome
`runPURPLE.genomeVersion`|String|"V38"|genome version for AMBER, default set to V38
`runPURPLE.gcProfile`|String|"$HMFTOOLS_DATA_ROOT/GC_profile.1000bp.38.cnp"|GC profile, generated for COBALT
`runPURPLE.purpleScript`|String|"java -Xmx8G -jar $HMFTOOLS_ROOT/purple.jar"|location of PURPLE script
`runPURPLE.modules`|String|"argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"|Required environment modules
`runPURPLE.threads`|Int|8|Requested CPU threads
`runPURPLE.memory`|Int|32|Memory allocated for this job (GB)
`runPURPLE.timeout`|Int|100|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`purple_directory`|File|zipped directory containing all PURPLE output


## Commands
 This section lists commands run by the PURPLE workflow
 
 ### Run AMBER 
 
     ~{amberScript} \
       -reference ~{normal_name} -reference_bam ~{normal_bam} \
       -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
       -output_dir ~{tumour_name}.amber/ \
       -loci ~{PON} \
       -ref_genome_version ~{genomeVersion}
 
 ### Run COBALT 
 
     ~{colbaltScript} \
       -reference ~{normal_name} -reference_bam ~{normal_bam} \
       -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
       -output_dir ~{tumour_name}.cobalt/ \
       -gc_profile ~{gcProfile} \
       -pcf_gamma ~{gamma}
 
 ### Filter VCFs (both for SVs and for SNV+InDel) 
 
     ~{bcftoolsScript} -f 'PASS' ~{vcf}  >~{tumour_name}.PASS.vcf
 
 ### Run PURPLE 
 
     ~{purpleScript} \
       -no_charts \
       -ref_genome_version ~{genomeVersion} \
       -ref_genome ~{refFasta}  \
       -gc_profile ~{gcProfile} \
       -ensembl_data_dir ~{ensemblDir}  \
       -reference ~{normal_name} -tumor ~{tumour_name}  \
       -amber ~{tumour_name}.amber -cobalt ~{tumour_name}.cobalt \
       -somatic_vcf ~{smalls_vcf} \
       -structural_vcf ~{SV_vcf} \
       -output_dir ~{tumour_name}.purple/
 
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
