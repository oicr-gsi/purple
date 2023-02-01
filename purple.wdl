version 1.0

workflow purple {
  input {
    File tumour_bam
    File tumour_bai
    File normal_bam
    File normal_bai
    File SV_vcf
    File smalls_vcf
    String normal_name
    String tumour_name
  }

  parameter_meta {
    tumour_bam: "Input tumor file (bam)."
    tumour_bai: "Input tumor file index (bai)."
    tumour_name: "Name of tumour sample"
    normal_bam: "Input normal file (bam)."
    normal_bai: "Input normal file index (bai)."
    normal_name: "Name of normal sample"
    SV_vcf: "somatic Structural Variant File (.vcf) from gridss."
    smalls_vcf: "somatic small (SNV+indel) Variant File (.vcf) [tested with mutect2]."
  }

  call amber {
    input:
      tumour_bam = tumour_bam,
      tumour_bai = tumour_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = normal_name,
      tumour_name = tumour_name
  }

  call cobalt {
    input:
      tumour_bam = tumour_bam,
      tumour_bai = tumour_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = normal_name,
      tumour_name = tumour_name
  }

  call filterVCF as filterSV {
    input: 
      vcf = SV_vcf,
      tumour_name = tumour_name
  }

  call filterVCF as filterSMALL {
    input: 
      vcf = smalls_vcf,
      tumour_name = tumour_name
  }

  call purple {
    input:
      normal_name = normal_name,
      tumour_name = tumour_name,
      amber_directory = amber.output_directory,
      cobalt_directory = cobalt.output_directory,
      SV_vcf = filterSV.filtered_vcf,
      smalls_vcf = filterSMALL.filtered_vcf
  }

  meta {
    author: "Felix Beaudry"
    email: "fbeaudry@oicr.on.ca"
    description: "performs purity and ploidy estimation"
    dependencies: [
    {
      name: "PURPLE",
      url: "https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md"
    }
    ]
  }

  output {
    File purple_directory = "~{tumour_name}.purple.zip"
  }
}

task amber {
  input {
    String tumour_name
    File tumour_bam
    File tumour_bai
    String normal_name
    File normal_bam
    File normal_bai
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

    mkdir ~{tumour_name}.amber  

    ~{amberScript} \
      -reference ~{normal_name} -reference_bam ~{normal_bam} \
      -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
      -output_dir ~{tumour_name}.amber/ \
      -loci ~{PON} \
      -ref_genome_version ~{genomeVersion}

    zip -r ~{tumour_name}.amber.zip ~{tumour_name}.amber/

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File output_directory = "~{tumour_name}.amber.zip"
  }
}

task cobalt {
  input {
    String tumour_name
    File tumour_bam
    File tumour_bai
    String normal_name
    File normal_bam
    File normal_bai
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

    mkdir ~{tumour_name}.cobalt 

    ~{colbaltScript} \
      -reference ~{normal_name} -reference_bam ~{normal_bam} \
      -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
      -output_dir ~{tumour_name}.cobalt/ \
      -gc_profile ~{gcProfile} \
      -pcf_gamma ~{gamma}

    zip -r ~{tumour_name}.cobalt.zip ~{tumour_name}.cobalt/

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File output_directory = "~{tumour_name}.cobalt.zip"
  }
}

task filterVCF {
  
  input {
    String tumour_name
    File vcf
    String modules = "bcftools"
    String bcftoolsScript = "$BCFTOOLS_ROOT/bin/bcftools view"
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }
    command <<<
    set -euo pipefail

    ~{bcftoolsScript} -f 'PASS' ~{vcf}  >~{tumour_name}.PASS.vcf

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File filtered_vcf = "~{tumour_name}.PASS.vcf"
  }
}

task purple {
  input {
    String normal_name
    String tumour_name
    File amber_directory
    File cobalt_directory
    File SV_vcf
    File smalls_vcf
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String refFasta = "$HMFTOOLS_DATA_ROOT/hg38_random.fa"
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

    unzip ~{amber_directory} 
    unzip ~{cobalt_directory} 
    mkdir ~{tumour_name}.purple 

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

      zip -r ~{tumour_name}.purple.zip ~{tumour_name}.purple/
  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File purple_directory = "~{tumour_name}.purple.zip"
  }
}
