version 1.0

workflow hmftools {
  input {
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
  }

  parameter_meta {
    tumorBam: "Input tumor file (bam or sam)."
    normalBam: "Input normal file (bam or sam)."
  }

  meta {
    author: "Felix Beaudry"
    email: "afortuna@oicr.on.ca"
    description: "performs somatic genomic rearrangement detection and classification"
    dependencies: [
    {
      name: "GRIDSS",
      url: "https://github.com/PapenfussLab/gridss"
    },
    {
      name: "PURPLE",
      url: "https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md"
    },
    {
      name: "LINX",
      url: "https://github.com/hartwigmedical/hmftools/blob/master/linx/README.md"
    },
    {
      name: "VIRUSBreakend",
      url: "https://github.com/PapenfussLab/gridss/blob/master/VIRUSBreakend_Readme.md"
    }
    ]
  }

  call gridss {
    input:
      tumorBam = tumorBam
      normalBam = normalBam
  }

  call cobalt {
  
  }

  call amber {
  
  }



  call purple {
  
  }

  call linx {
  
  }

  call virusbreakend {

  }

  output {
  }
}

task gridss {
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String refFasta = "$HMFTOOLS_DATA_ROOT/hg38_random.fa"
    String gridssScript = "${GRIDSS_ROOT}/gridss --jar ${GRIDSS_ROOT}/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    File normBam
    File normBai
    File tumorBam
    File tumorBai
    Int threads = 8
    Int memory = 50
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    ~{gridssScript} \
    --reference ~{refFasta} \
    --output ./ \
    ~{normBam} ~{tumorBam}

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcf = "~{outputVcf}.allocated.vcf"
  }
}

task amber {
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String refFasta = "${HMFTOOLS_DATA_ROOT}/hg38_random.fa"
    String amberScript = "java -Xmx32G -cp ${HMFTOOLS_ROOT}/amber.jar com.hartwig.hmftools.amber.AmberApplication"
    String colbaltScript = "java -Xmx8G -cp $HMFTOOLS_ROOT/cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication"
    String bcftoolsScript = "$BCFTOOLS_ROOT/bin/bcftools view"
    String purpleScript = "java -Xmx8G -jar $HMFTOOLS_ROOT/purple.jar"
    File normBam
    File normBai
    String normName = basename("~{normBam}", ".filter.deduped.realigned.recalibrated.bam")
    File tumorBam
    File tumorBai
    String tumorName = basename("~{tumorBam}", ".filter.deduped.realigned.recalibrated.bam")
    String PON = "${HMFTOOLS_DATA_ROOT}/GermlineHetPon.38.vcf.gz"
    String gcProfile = "${HMFTOOLS_DATA_ROOT}/GC_profile.1000bp.38.cnp"
    String ensemblDir = "${HMFTOOLS_DATA_ROOT}/ensembl"
    String genomeVersion = "V38"
    String gamma = 100
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    mkdir amber  

    ~{amberScript} \
      -reference ~{normName} -reference_bam ~{normBam} \
      -tumor ~{tumorName} -tumor_bam ~{tumorBam} \
      -output_dir amber/ \
      -loci ~{PON} \
      -ref_genome_version ~{genomeVersion}

    zip amber/

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File amberDir = "amber.zip"
  }
}

task cobalt {
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String refFasta = "${HMFTOOLS_DATA_ROOT}/hg38_random.fa"
    String amberScript = "java -Xmx32G -cp ${HMFTOOLS_ROOT}/amber.jar com.hartwig.hmftools.amber.AmberApplication"
    String colbaltScript = "java -Xmx8G -cp $HMFTOOLS_ROOT/cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication"
    String bcftoolsScript = "$BCFTOOLS_ROOT/bin/bcftools view"
    String purpleScript = "java -Xmx8G -jar $HMFTOOLS_ROOT/purple.jar"
    File normBam
    File normBai
    String normName = basename("~{normBam}", ".filter.deduped.realigned.recalibrated.bam")
    File tumorBam
    File tumorBai
    String tumorName = basename("~{tumorBam}", ".filter.deduped.realigned.recalibrated.bam")
    String PON = "${HMFTOOLS_DATA_ROOT}/GermlineHetPon.38.vcf.gz"
    String gcProfile = "${HMFTOOLS_DATA_ROOT}/GC_profile.1000bp.38.cnp"
    String ensemblDir = "${HMFTOOLS_DATA_ROOT}/ensembl"
    String genomeVersion = "V38"
    String gamma = 100
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    mkdir cobalt 

    ~{colbaltScript} \
      -reference ~{normSample} -reference_bam ~{normBam} \
      -tumor ~{tumorSample} -tumor_bam ~{tumorBam} \
      -output_dir cobalt/ \
      -gc_profile ~{gcProfile} \
      -pcf_gamma ~{gamma}

    zip cobalt

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcf = "~{outputVcf}.allocated.vcf"
  }
}

task purple {
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String refFasta = "${HMFTOOLS_DATA_ROOT}/hg38_random.fa"
    String amberScript = "java -Xmx32G -cp ${HMFTOOLS_ROOT}/amber.jar com.hartwig.hmftools.amber.AmberApplication"
    String colbaltScript = "java -Xmx8G -cp $HMFTOOLS_ROOT/cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication"
    String bcftoolsScript = "$BCFTOOLS_ROOT/bin/bcftools view"
    String purpleScript = "java -Xmx8G -jar $HMFTOOLS_ROOT/purple.jar"
    File normBam
    File normBai
    String normName = basename("~{normBam}", ".filter.deduped.realigned.recalibrated.bam")
    File tumorBam
    File tumorBai
    String tumorName = basename("~{tumorBam}", ".filter.deduped.realigned.recalibrated.bam")
    String PON = "${HMFTOOLS_DATA_ROOT}/GermlineHetPon.38.vcf.gz"
    String gcProfile = "${HMFTOOLS_DATA_ROOT}/GC_profile.1000bp.38.cnp"
    String ensemblDir = "${HMFTOOLS_DATA_ROOT}/ensembl"
    String genomeVersion = "V38"
    String gamma = 100
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    unzip amber cobalt 

    ~{bcftoolsScript} -f 'PASS' \
      ~{smalls_vcf} \
      -s ~{tumorSample} \
      >~{tumorSample}.SMALLS.PASS.vcf

    ~{bcftoolsScript} -f 'PASS' \
      ~{SV_vcf} \
      -s ~{tumorSample} \
      >~{tumorSample}.SV.PASS.vcf

    mkdir purple 

    ~{purpleScript} \
      -no_charts \
      -ref_genome_version ~{genomeVersion}
      -ref_genome ~{refFasta}  \
      -gc_profile ~{gcProfile} \
      -ensembl_data_dir ~{ensemblDir}  \
      -reference ~{normSample} -tumor ~{tumorSample}  \
      -amber amber/ -cobalt cobalt/ \
      -somatic_vcf ~{tumorSample}.SMALLS.PASS.vcf \
      -structural_vcf ~{tumorSample}.SV.PASS.vcf
      -output_dir purple/

      zip purple

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcf = "~{outputVcf}.allocated.vcf"
  }
}


task linx {
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String refFasta = "$HMFTOOLS_DATA_ROOT/hg38_random.fa"
    File normBam
    File normBai
    File tumorBam
    File tumorBai
    string fragileSitesFile = ${HMFTOOLS_DATA_ROOT}/fragile_sites_hmf.38.csv
    string lineElementFile = ${HMFTOOLS_DATA_ROOT}/line_elements.38.csv
    string knownFusionFile = ${HMFTOOLS_DATA_ROOT}/known_fusion_data.38.csv
    String genomeVersion = "38"
    String gamma = 100
    Int threads = 8
    Int memory = 8
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    unzip purple

    mkdir linx 

    java -Xmx8G -jar $HMFTOOLS_ROOT/linx.jar \
      -check_fusions -log_debug \
      -ref_genome_version ~{genomeVersion} \
      -fragile_site_file ~{fragileSitesFile} \
      -line_element_file ~{lineElementFile} \
      -known_fusion_file ~{knownFusionFile} \
      -ensembl_data_dir ~{ensemblDir}
      -sample $tumorSample \
      -sv_vcf purple/~{tumorSample}.purple.sv.vcf.gz \
      -purple_dir purple/ \
      -output_dir linx/

    zip linx

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcf = "~{outputVcf}.allocated.vcf"
  }
}

task virusbreakend {
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String refFasta = "$HMFTOOLS_DATA_ROOT/hg38_random.fa"
    File normBam
    File normBai
    File tumorBam
    File tumorBai
    File gcProfile = ${HMFTOOLS_DATA_ROOT}/GC_profile.1000bp.38.cnp
    String gamma = 100
    Int threads = 8
    Int memory = 8
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    $GRIDSS_ROOT/virusbreakend \
      --kraken2db $VIRUSBREAKEND_DB_ROOT/ \
      --reference ${testDir}/hg38_random.fa \
      --output ${sample}.virusbreakend.vcf \
      ${sample}.bam

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcf = "~{outputVcf}.allocated.vcf"
  }
}