version 1.0

workflow hmftools {
  input {
    File tumorBam
    File tumorBai
    File normBam
    File normBai
    File svVCFfile
    File smallsVCFfiles
    String normName = basename("~{normBam}", ".filter.deduped.realigned.recalibrated.bam")
    String tumorName = basename("~{tumorBam}", ".filter.deduped.realigned.recalibrated.bam")
    Boolean runLinx = false
  }

  parameter_meta {
    tumorBam: "Input tumor file (bam)."
    tumorBai: "Input tumor file index (bai)."
    normBam: "Input normal file (bam)."
    normBai: "Input normal file index (bai)."
    svVCFfile: "somatic Structural Variant File (.vcf) from gridss."
    smallsVCFfiles: "somatic small (SNV+indel) Variant File (.vcf) [tested with mutect2]."
  }

  call amber {
    input:
      tumorBam = tumorBam,
      tumorBai = tumorBai,
      normBam = normBam,
      normBai = normBai,
      normName = normName,
      tumorName = tumorName
  }

  call cobalt {
    input:
      tumorBam = tumorBam,
      tumorBai = tumorBai,
      normBam = normBam,
      normBai = normBai,
      normName = normName,
      tumorName = tumorName
  }

  call purple {
    input:
      normName = normName,
      tumorName = tumorName,
      amberDir = amber.amberDir,
      cobaltDir = cobalt.cobaltDir,
      svVCFfile = svVCFfile,
      smallsVCFfiles = smallsVCFfiles
  }

  if(runLinx == true){
    call linx {
      input:
        tumorName = tumorName,
        purpleDir = purple.purpleDir
    }
  }

  meta {
    author: "Felix Beaudry"
    email: "fbeaudry@oicr.on.ca"
    description: "performs purity and ploidy estimation"
    dependencies: [
    {
      name: "PURPLE",
      url: "https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md"
    },
    {
      name: "LINX",
      url: "https://github.com/hartwigmedical/hmftools/blob/master/linx/README.md"
    }
    ]
  }

  output {
    File purpleDir = "~{tumorName}.purple.zip"
    File? linxDir = "~{tumorName}.linx.zip"
  }
}

task amber {
  input {
    File tumorBam
    File tumorBai
    File normBam
    File normBai
    String normName = basename("~{normBam}", ".filter.deduped.realigned.recalibrated.bam")
    String tumorName = basename("~{tumorBam}", ".filter.deduped.realigned.recalibrated.bam")
    String amberScript = "java -Xmx32G -cp $HMFTOOLS_ROOT/amber.jar com.hartwig.hmftools.amber.AmberApplication"
    String PON = "$HMFTOOLS_DATA_ROOT/GermlineHetPon.38.vcf.gz"
    String genomeVersion = "V38"
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    mkdir ~{tumorName}.amber  

    ~{amberScript} \
      -reference ~{normName} -reference_bam ~{normBam} \
      -tumor ~{tumorName} -tumor_bam ~{tumorBam} \
      -output_dir ~{tumorName}.amber/ \
      -loci ~{PON} \
      -ref_genome_version ~{genomeVersion}

    zip ~{tumorName}.amber/

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File amberDir = "~{tumorName}.amber.zip"
  }
}

task cobalt {
  input {
    File tumorBam
    File tumorBai
    File normBam
    File normBai
    String normName = basename("~{normBam}", ".filter.deduped.realigned.recalibrated.bam")
    String tumorName = basename("~{tumorBam}", ".filter.deduped.realigned.recalibrated.bam")
    String colbaltScript = "java -Xmx8G -cp $HMFTOOLS_ROOT/cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication"
    String gcProfile = "$HMFTOOLS_DATA_ROOT/GC_profile.1000bp.38.cnp"
    String gamma = 100
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    mkdir ~{tumorName}.cobalt 

    ~{colbaltScript} \
      -reference ~{normName} -reference_bam ~{normBam} \
      -tumor ~{tumorName} -tumor_bam ~{tumorBam} \
      -output_dir ~{tumorName}.cobalt/ \
      -gc_profile ~{gcProfile} \
      -pcf_gamma ~{gamma}

    zip ~{tumorName}.cobalt

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File cobaltDir = "~{tumorName}.cobalt.zip"
  }
}

task purple {
  input {
    File svVCFfile
    File smallsVCFfiles
    String normName
    String tumorName
    File amberDir
    File cobaltDir
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String refFasta = "$HMFTOOLS_DATA_ROOT/hg38_random.fa"
    String bcftoolsScript = "$BCFTOOLS_ROOT/bin/bcftools view"
    String purpleScript = "java -Xmx8G -jar $HMFTOOLS_ROOT/purple.jar"
    String gcProfile = "$HMFTOOLS_DATA_ROOT/GC_profile.1000bp.38.cnp"
    String ensemblDir = "$HMFTOOLS_DATA_ROOT/ensembl"
    String genomeVersion = "V38"
    String gamma = 100
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    unzip ~{amberDir} ~{cobaltDir} 

    ~{bcftoolsScript} -f 'PASS' \
      ~{smallsVCFfiles} \
      -s ~{tumorName} \
      >~{tumorName}.SMALLS.PASS.vcf

    ~{bcftoolsScript} -f 'PASS' \
      ~{svVCFfile} \
      -s ~{tumorName} \
      >~{tumorName}.SV.PASS.vcf

    mkdir ~{tumorName}.purple 

    ~{purpleScript} \
      -no_charts \
      -ref_genome_version ~{genomeVersion}
      -ref_genome ~{refFasta}  \
      -gc_profile ~{gcProfile} \
      -ensembl_data_dir ~{ensemblDir}  \
      -reference ~{normName} -tumor ~{tumorName}  \
      -amber ~{tumorName}.amber/ -cobalt ~{tumorName}.cobalt/ \
      -somatic_vcf ~{tumorName}.SMALLS.PASS.vcf \
      -structural_vcf ~{tumorName}.SV.PASS.vcf
      -output_dir ~{tumorName}.purple/

      zip ~{tumorName}.purple

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File purpleDir = "~{tumorName}.purple.zip"
  }
}


task linx {
  input {
    File purpleDir
    String tumorName
    String fragileSitesFile = "$HMFTOOLS_DATA_ROOT/fragile_sites_hmf.38.csv"
    String lineElementFile = "$HMFTOOLS_DATA_ROOT/line_elements.38.csv"
    String knownFusionFile = "$HMFTOOLS_DATA_ROOT/known_fusion_data.38.csv"
    String ensemblDir = "$HMFTOOLS_DATA_ROOT/ensembl"
    String genomeVersion = "38"
    String linxScript = "java -Xmx8G -jar $HMFTOOLS_ROOT/linx.jar"
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    Int threads = 8
    Int memory = 8
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    unzip ~{purpleDir}

    mkdir ~{tumorName}.linx 

     ~{linxScript} \
      -check_fusions -log_debug \
      -ref_genome_version ~{genomeVersion} \
      -fragile_site_file ~{fragileSitesFile} \
      -line_element_file ~{lineElementFile} \
      -known_fusion_file ~{knownFusionFile} \
      -ensembl_data_dir ~{ensemblDir}
      -sample ~{tumorName} \
      -sv_vcf ~{tumorName}.purple/~{tumorName}.purple.sv.vcf.gz \
      -purple_dir ~{tumorName}.purple/ \
      -output_dir ~{tumorName}.linx/

    zip ~{tumorName}.linx

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File linxDir = "~{tumorName}.linx.zip"
  }
}

