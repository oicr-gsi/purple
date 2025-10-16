version 1.0

struct GenomeResources {
    String PON
    String modules
    String gatkModules
    String gcProfile
    String ensemblDir
    String refFasta
    String refFai
    String pon_sgl_file
    String pon_sv_file
    String known_hotspot_file
    String repeat_mask_file
    String knownfusion
}

workflow purple {
  input {
    File tumour_bam
    File tumour_bai
    File normal_bam
    File normal_bai
    String genomeVersion = "38"
    Boolean doSV = true
    Boolean doSMALL = true
  }

  parameter_meta {
    tumour_bam: "Input tumor file (bam)"
    tumour_bai: "Input tumor file index (bai)"
    normal_bam: "Input normal file (bam)"
    normal_bai: "Input normal file index (bai)"
    genomeVersion: "Genome Version, only 38 supported"
    doSV: "include somatic structural variant calls, true/false"
    doSMALL: "include somatic small (SNV+indel) calls, true/false"
  }

Map[String,GenomeResources] resources = {
  "38": {
    "modules": "hmftools/1.1 hg38/p12 hmftools-data/53138",
    "gatkModules": "hg38-gridss-index/1.0 gatk/4.1.6.0",
    "refFasta": "$HG38_ROOT/hg38_random.fa",
    "refFai": "$HG38_GRIDSS_INDEX_ROOT/hg38_random.fa.fai",
    "PON" : "$HMFTOOLS_DATA_ROOT/copy_number/GermlineHetPon.38.vcf.gz",
    "ensemblDir": "$HMFTOOLS_DATA_ROOT/ensembl_data",
    "gcProfile": "$HMFTOOLS_DATA_ROOT/copy_number/GC_profile.1000bp.38.cnp",
    "pon_sgl_file": "$HMFTOOLS_DATA_ROOT/sv/sgl_pon.38.bed.gz",
    "pon_sv_file": "$HMFTOOLS_DATA_ROOT/sv/sv_pon.38.bedpe.gz",
    "known_hotspot_file": "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe",
    "repeat_mask_file": "$HMFTOOLS_DATA_ROOT/sv/repeat_mask_data.38.fa.gz",
    "knownfusion": "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe"
  }
}

  call extractName as extractTumorName {
    input:
    refFasta = resources [ genomeVersion ].refFasta,
    refFai = resources [ genomeVersion ].refFai,
    modules = resources [ genomeVersion ].gatkModules,
    inputBam = tumour_bam,
    inputBai = tumour_bai
  }

  call extractName as extractNormalName {
    input:
    refFasta = resources [ genomeVersion ].refFasta,
    refFai = resources [ genomeVersion ].refFai,
    modules = resources [ genomeVersion ].gatkModules,
    inputBam = normal_bam,
    inputBai = normal_bai
  }

  call amber {
    input:
      tumour_bam = tumour_bam,
      tumour_bai = tumour_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = extractNormalName.input_name,
      tumour_name = extractTumorName.input_name,
      genomeVersion = genomeVersion,
      modules = resources [ genomeVersion ].modules,
      PON = resources [ genomeVersion ].PON
  }

  call cobalt {
    input:
      tumour_bam = tumour_bam,
      tumour_bai = tumour_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = extractNormalName.input_name,
      tumour_name = extractTumorName.input_name,
      modules = resources [ genomeVersion ].modules,
      gcProfile = resources [ genomeVersion ].gcProfile
  }

  if(doSV) {
    call filterSV {
      input: 
        normal_name = extractNormalName.input_name,
        tumour_name = extractTumorName.input_name,
        genomeVersion = genomeVersion,
        refFasta = resources [ genomeVersion ].refFasta,
        pon_sgl_file = resources [ genomeVersion ].pon_sgl_file,
        pon_sv_file = resources [ genomeVersion ].pon_sv_file,
        known_hotspot_file = resources [ genomeVersion ].known_hotspot_file,
        repeat_mask_file = resources [ genomeVersion ].repeat_mask_file,
        modules = resources [ genomeVersion ].modules
    }
  }

  if(doSMALL) {
    call filterSMALL {
      input: 
        normal_name = extractNormalName.input_name,
        tumour_name = extractTumorName.input_name
    }
  }

  call runPURPLE {
    input:
      normal_name = extractNormalName.input_name,
      tumour_name = extractTumorName.input_name,
      amber_directory = amber.output_directory,
      cobalt_directory = cobalt.output_directory,
      SV_vcf = filterSV.filtered_vcf,
      smalls_vcf = filterSMALL.filtered_vcf,
      genomeVersion = genomeVersion,
      modules = resources [ genomeVersion ].modules,
      gcProfile = resources [ genomeVersion ].gcProfile,
      ensemblDir = resources [ genomeVersion ].ensemblDir,
      refFasta = resources [ genomeVersion ].refFasta
  }

  call expandAlternates {
  }

  scatter (alternate in expandAlternates.alternates_ploidies) {
      call runPURPLE as runPURPLEAlternates {
        input:
          solution_name = alternate[0],
          min_ploidy = alternate[0],
          max_ploidy = alternate[1],
          normal_name = extractNormalName.input_name,
          tumour_name = extractTumorName.input_name,
          amber_directory = amber.output_directory,
          cobalt_directory = cobalt.output_directory,
          SV_vcf = filterSV.filtered_vcf,
          smalls_vcf = filterSMALL.filtered_vcf,
          genomeVersion = genomeVersion,
          modules = resources [ genomeVersion ].modules,
          gcProfile = resources [ genomeVersion ].gcProfile,
          ensemblDir = resources [ genomeVersion ].ensemblDir,
          refFasta = resources [ genomeVersion ].refFasta
      }
  }

  call group_alternates {
    input:
      tumour_name = extractTumorName.input_name,
      alternate_solutions = select_all(runPURPLEAlternates.purple_directory)
  }

  if(doSV) {
    call LINX{
      input:
          tumour_name = extractTumorName.input_name,
          ensemblDir = resources [ genomeVersion ].ensemblDir,
          genomeVersion = genomeVersion, 
          fusions_file = resources [ genomeVersion ].knownfusion,
          purple_dir = runPURPLE.purple_directory,
          modules = resources [ genomeVersion ].modules
    }
  }

  meta {
    author: "Felix Beaudry, Lawrence Heisler"
    email: "fbeaudry@oicr.on.ca, lheisler@oicr.on.ca"
    description: "performs purity and ploidy estimation"
    dependencies: [
    {
      name: "PURPLE",
      url: "https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md"
    },
    {
      name: "AMBER",
      url: "https://github.com/hartwigmedical/hmftools/blob/master/amber/README.md"
    },
    {
      name: "COBALT",
      url: "https://github.com/hartwigmedical/hmftools/blob/master/cobalt/README.md"
    },
    {
      name: "Linx",
      url: "https://github.com/hartwigmedical/hmftools/blob/master/linx/README.md"
    },
    {
      name: "bcftools",
      url: "https://www.htslib.org/doc/1.9/bcftools.html"
    }
    ]
    output_meta: {
    purple_directory: {
        description: "Zipped results from PURPLE",
        vidarr_label: "purple_directory"
    },
    purple_qc: {
        description: "QC results from PURPLE",
        vidarr_label: "purple_qc"
    },
    purple_purity: {
        description: "tab seperated Purity estimate from PURPLE",
        vidarr_label: "purple_purity"
    },
    purple_purity_range: {
        description: "tab seperated range of Purity estimate from PURPLE",
        vidarr_label: "purple_purity_range"
    },
    purple_segments: {
        description: "tab seperated segments estimated by PURPLE",
        vidarr_label: "purple_segments"
    },
    purple_cnv: {
        description: "tab seperated somatic copy number variants from PURPLE",
        vidarr_label: "purple_cnv"
    },
    purple_cnv_gene: {
        description: "tab seperated somatic gene-level copy number variants from PURPLE",
        vidarr_label: "purple_cnv_gene"
    },
    purple_SV_index: {
        description: "Structural Variant .vcf index edited by PURPLE",
        vidarr_label: "purple_SV_index"
    },
    purple_SV: {
        description: "Structural Variant .vcf edited by PURPLE",
        vidarr_label: "purple_SV"
    },
    purple_SMALL_index: {
        description: "SNV+IN/DEL .vcf index edited by PURPLE",
        vidarr_label: "purple_SMALL_index"
    },
    purple_SMALL: {
        description: "SNV+IN/DEL .vcf edited by PURPLE",
        vidarr_label: "purple_SMALL"
    },
    purple_alternate_directory: {
        description: "Directory for alternate solution files",
        vidarr_label: "purple_alternate_directory"
    }
}
  }

  output {
    File? purple_directory = runPURPLE.purple_directory
    File purple_alternate_directory = group_alternates.purple_alternate_directory
    File purple_qc = runPURPLE.purple_qc
    File purple_purity = runPURPLE.purple_purity
    File purple_purity_range = runPURPLE.purple_purity_range
    File purple_segments = runPURPLE.purple_segments
    File purple_cnv = runPURPLE.purple_cnv
    File purple_cnv_gene = runPURPLE.purple_cnv_gene
    File? purple_SV_index = runPURPLE.purple_SV_index
    File? purple_SV = runPURPLE.purple_SV
    File? purple_SMALL_index = runPURPLE.purple_SMALL_index
    File? purple_SMALL = runPURPLE.purple_SMALL
  }
}


# =========================================
# Job to extract names from input bam files
# =========================================
task extractName {
  input {
    String modules
    String refFasta 
    String refFai 
    File inputBam
    File inputBai
    Int memory = 4
    Int timeout = 4
  }

  parameter_meta {
    inputBam: "input .bam file"
    inputBai: "input .bai file"
    refFasta: "Reference FASTA file"
    refFai: "Reference fai index"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    if [ -f "~{inputBam}" ]; then
      gatk --java-options "-Xmx1g" GetSampleName -R ~{refFasta} -I ~{inputBam} -O input_name.txt -encode
    fi

    cat input_name.txt
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      input_name: "name of the input"
    }
  }

  output {
    String input_name = read_string(stdout()) 
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
    String PON
    String genomeVersion
    Int min_mapping_quality = 30
    Int min_base_quality = 25
    String modules
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    tumour_bam: "Tumour bam"
    tumour_bai: "Matching bai for Tumour bam"
    normal_name: "Name for Normal sample"
    normal_bam: "Normal bam"
    normal_bai: "Matching bai for Normal bam"
    amberScript: "location of AMBER script"
    PON: "Panel of Normal (PON) file, generated for AMBER"
    genomeVersion: "genome version (37 or 38)"
    min_mapping_quality: "Minimum mapping quality for an alignment to be used"
    min_base_quality: "Minimum quality for a base to be considered"
		modules: "Required environment modules"
		memory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

  command <<<
    set -euo pipefail

    mkdir ~{tumour_name}.amber  

    ~{amberScript} \
      -reference ~{normal_name} -reference_bam ~{normal_bam} \
      -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
      -output_dir ~{tumour_name}.amber/ \
      -loci ~{PON} \
      -ref_genome_version ~{genomeVersion} \
      -min_mapping_quality ~{min_mapping_quality} \
      -min_base_quality ~{min_base_quality} 

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

  meta {
		output_meta: {
			output_directory: "Zipped AMBER results directory"
		}
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
    String gcProfile
    String gamma = 300
    Int min_mapping_quality = 30
    String modules
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    tumour_bam: "Tumour bam"
    tumour_bai: "Matching bai for Tumour bam"
    normal_name: "Name for Normal sample"
    normal_bam: "Normal bam"
    normal_bai: "Matching bai for Normal bam"
    colbaltScript: "location of COBALT script"
    gcProfile: "GC profile, generated for COBALT"
    gamma: "gamma (penalty) value for segmenting"
    min_mapping_quality: "Minimum mapping quality for an alignment to be used"
		modules: "Required environment modules"
		memory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

  command <<<
    set -euo pipefail

    mkdir ~{tumour_name}.cobalt 

    ~{colbaltScript} \
      -reference ~{normal_name} -reference_bam ~{normal_bam} \
      -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
      -output_dir ~{tumour_name}.cobalt/ \
      -gc_profile ~{gcProfile} \
      -pcf_gamma ~{gamma} \
      -min_quality ~{min_mapping_quality}

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

  meta {
		output_meta: {
			output_directory: "Zipped COBALT results directory"
		}
	}

}

task filterSV {
  
  input {
    String normal_name
    String tumour_name
    File? vcf
    String gripssScript = "java -Xmx80G -jar $HMFTOOLS_ROOT/gripss.jar"
    String refFasta
    String genomeVersion
    String pon_sgl_file
    String pon_sv_file
    String known_hotspot_file
    String repeat_mask_file
    Int hard_min_tumor_qual = 500
    String filter_sgls = "-filter_sgls"
    String modules
    Int memory = 80
    Int threads = 1
    Int timeout = 100
  }

  parameter_meta {
    normal_name:  "Name for normal sample"
    tumour_name: "Name for Tumour sample"
    vcf: "VCF file for filtering"
    gripssScript: "location and java call for gripss jar"
    refFasta: "reference fasta"
    genomeVersion: "version of the genome"
    pon_sgl_file: "panel of normals, single breakends"
    pon_sv_file: "panel of normals, germline structural variants"
    known_hotspot_file: "known/common hotspots for structural variants"
    repeat_mask_file: "repeating masking information"
    hard_min_tumor_qual: "Any variant with QUAL less than x is filtered "
    filter_sgls: "include filtering of single breakends"
		modules: "Required environment modules"
		memory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}
  
  command <<<
    set -euo pipefail

    mkdir gripss

    ~{gripssScript} \
        -vcf ~{vcf}  \
        -sample ~{tumour_name} -reference ~{normal_name} \
        -ref_genome_version ~{genomeVersion} \
        -ref_genome ~{refFasta} \
        -pon_sgl_file ~{pon_sgl_file} \
        -pon_sv_file ~{pon_sv_file} \
        -known_hotspot_file ~{known_hotspot_file} \
        -repeat_mask_file ~{repeat_mask_file} \
        -output_dir gripss/ \
        -hard_min_tumor_qual ~{hard_min_tumor_qual} \
        ~{filter_sgls}

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File soft_filtered_vcf = "gripss/~{tumour_name}.gripss.vcf.gz"
    File filtered_vcf = "gripss/~{tumour_name}.gripss.filtered.vcf.gz"
  }

  meta {
		output_meta: {
			filtered_vcf: "high confidence structural variant VCF",
      soft_filtered_vcf: "structural variant VCF after first filtering"
		}
	}

}

task filterSMALL {
  
  input {
    String tumour_name
    String normal_name
    File? vcf
    File? vcf_index
    String bcftoolsScript = "$BCFTOOLS_ROOT/bin/bcftools"
    String genome = "$HG38_ROOT/hg38_random.fa"
    String regions = "chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr2,chr20,chr21,chr22,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chrX"
    String difficultRegions = "--targets-file $HG38_DAC_EXCLUSION_ROOT/hg38-dac-exclusion.v2.bed"
    String tumorVAF = "0.01"
    String modules = "bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0"
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  parameter_meta {
    normal_name:  "Name for normal sample"
    tumour_name: "Name for Tumour sample"
    vcf: "VCF file for filtering"
    vcf_index: "index of VCF file for filtering"
    bcftoolsScript: "location for bcftools"
    genome: "reference fasta"
    regions: "regions/chromosomes to include"
    difficultRegions: "regions to exclude because they are difficult"
    tumorVAF: "minimum variant allele frequency for tumour calls to pass filter"
		modules: "Required environment modules"
		memory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    echo ~{normal_name} >samples.txt
    echo ~{tumour_name} >>samples.txt

     ~{bcftoolsScript} view -f "PASS" -S samples.txt -r ~{regions} ~{difficultRegions} ~{vcf} |\
     ~{bcftoolsScript} norm --multiallelics - --fasta-ref ~{genome} |\
     ~{bcftoolsScript} filter -i "(FORMAT/AD[1:1])/(FORMAT/AD[1:0]+FORMAT/AD[1:1]) >= ~{tumorVAF}"  > ~{tumour_name}.PASS.vcf

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

  meta {
    output_meta: {
      filtered_vcf: "Filtered SNV+in/del VCF"
    }
  }

}

task runPURPLE {
  input {
    String normal_name
    String tumour_name
    String solution_name = "Primary"
    String outfilePrefix = tumour_name + ".sol" + solution_name
    File amber_directory
    File cobalt_directory
    File? SV_vcf
    File? smalls_vcf
    String ensemblDir
    String refFasta
    String genomeVersion
    String gcProfile 
    Int min_diploid_tumor_ratio_count = 60
    String purpleScript = "java -Xmx8G -jar $HMFTOOLS_ROOT/purple.jar"
    String? min_ploidy
    String? max_ploidy
    String? min_purity
    String? max_purity
    String? ploidy_penalty_factor
    String? ploidy_penalty_standard_deviation
    String modules
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    normal_name: "Name for Normal sample"
    solution_name: "Name of solution"
    outfilePrefix: "Prefix of output file"
    amber_directory: "zipped output from AMBER"
    cobalt_directory: "zipped output from COBALT"
    SV_vcf: "filtered structural variant (SV) vcf"
    smalls_vcf: "filtered SNV and indel (smalls) vcf"
    ensemblDir: "Directory of Ensembl data for PURPLE"
    refFasta: "fasta of reference genome"
    gcProfile: "GC profile, generated for COBALT"
    min_diploid_tumor_ratio_count: "smooth over contiguous segments which are fewer than this number of depth windows long and which have no SV support on either side and which are bounded on both sides by copy number regions which could be smoothed together using our normal smoothing rules."
    genomeVersion: "genome version for AMBER, default set to V38"
    purpleScript: "location of PURPLE script"
    min_ploidy: "minimum ploidy"
    max_ploidy: "max ploidy"
    min_purity: "mininimum purity"
    max_purity: "max purity"
    ploidy_penalty_factor: "multiplies aggregate event penalty by this factor"
    ploidy_penalty_standard_deviation: "not entirely sure what this does"
    modules: "Required environment modules"
		memory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

  command <<<
    set -euo pipefail

    unzip ~{amber_directory} 
    unzip ~{cobalt_directory} 
    mkdir ~{outfilePrefix}.purple 

    ~{purpleScript} \
      -ref_genome_version ~{genomeVersion} \
      -ref_genome ~{refFasta}  \
      -gc_profile ~{gcProfile} \
      -ensembl_data_dir ~{ensemblDir}  \
      -reference ~{normal_name} -tumor ~{tumour_name}  \
      -amber ~{tumour_name}.amber -cobalt ~{tumour_name}.cobalt \
      ~{"-ploidy_penalty_factor" + ploidy_penalty_factor} \
      ~{"-ploidy_penalty_standard_deviation" + ploidy_penalty_standard_deviation} \
      ~{"-somatic_sv_vcf " + SV_vcf} \
      ~{"-somatic_vcf " + smalls_vcf} \
      ~{"-min_ploidy " + min_ploidy} \
      ~{"-max_ploidy " + max_ploidy} \
      ~{"-min_purity " + min_purity} \
      ~{"-max_purity " + max_purity} \
      -no_charts \
      -min_diploid_tumor_ratio_count ~{min_diploid_tumor_ratio_count} \
      -output_dir ~{outfilePrefix}.purple 

    zip -r ~{outfilePrefix}.purple.zip ~{outfilePrefix}.purple/

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File? purple_directory = "~{outfilePrefix}.purple.zip"
    File purple_qc = "~{outfilePrefix}.purple/~{tumour_name}.purple.qc"
    File purple_purity = "~{outfilePrefix}.purple/~{tumour_name}.purple.purity.tsv"
    File purple_purity_range = "~{outfilePrefix}.purple/~{tumour_name}.purple.purity.range.tsv"
    File purple_segments = "~{outfilePrefix}.purple/~{tumour_name}.purple.segment.tsv"
    File purple_cnv = "~{outfilePrefix}.purple/~{tumour_name}.purple.cnv.somatic.tsv"
    File purple_cnv_gene = "~{outfilePrefix}.purple/~{tumour_name}.purple.cnv.gene.tsv"
    File? purple_SV_index = "~{outfilePrefix}.purple/~{tumour_name}.purple.sv.vcf.gz.tbi"
    File? purple_SV = "~{outfilePrefix}.purple/~{tumour_name}.purple.sv.vcf.gz"
    File? purple_SMALL_index = "~{outfilePrefix}.purple/~{tumour_name}.purple.somatic.vcf.gz.tbi"
    File? purple_SMALL = "~{outfilePrefix}.purple/~{tumour_name}.purple.somatic.vcf.gz"
  }

  meta {
		output_meta: {
			purple_directory: "Zipped Output PURPLE directory",
      purple_qc: "QC results from PURPLE",
      purple_purity: "tab seperated Purity estimate from PURPLE",
      purple_purity_range: "tab seperated range of Purity estimate from PURPLE",
      purple_segments: "tab seperated segments estimated by PURPLE",
      purple_cnv: "tab seperated somatic copy number variants from PURPLE",
      purple_cnv_gene: "tab seperated somatic gene-level copy number variants from PURPLE",
      purple_SV_index: "Structural Variant .vcf index edited by PURPLE",
      purple_SV: "Structural Variant .vcf edited by PURPLE",
      purple_SMALL_index: "SNV+IN/DEL .vcf index edited by PURPLE",
      purple_SMALL: "SNV+IN/DEL .vcf edited by PURPLE"
		}
	}
}

task expandAlternates {
  input {
    Int min_alternate_ploidy = 1
    Int max_alternate_ploidy = 8
    Int alternate_ploidy_step = 1
    Int jobMemory = 1
		Int timeout = 1
    String modules = "python/3.10.6"
  }
  parameter_meta {
    min_alternate_ploidy : "Minimum of alternative ploidy"
    max_alternate_ploidy : "Maximum of alternative ploidy"
    alternate_ploidy_step: "NUmber of steps"
    modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		timeout: "Hours before task timeout"
	}

  command <<<
    python<<CODE

    for ploidy in range(~{min_alternate_ploidy},~{max_alternate_ploidy},~{alternate_ploidy_step}):
      ploidy_plus_one = ploidy+~{alternate_ploidy_step}
      print(str(ploidy)+"\t"+str(ploidy_plus_one)+"\n")

    CODE
	>>>

  runtime {
		memory:  "~{jobMemory} GB"
		timeout: "~{timeout}"
    modules: "~{modules}"
	}

	output {
		Array[Array[String]] alternates_ploidies = read_tsv(stdout())
	}
}


task group_alternates {
  input {
    Array[File] alternate_solutions 
    String tumour_name
    Int jobMemory = 1
    Int threads = 1
    Int timeout = 2
  }

  parameter_meta {
    alternate_solutions: "Name for alternate_solutions"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}
  command <<<
    set -euo pipefail

    mkdir ~{tumour_name}.purple_alternates

    cp ~{sep=' ' alternate_solutions} ~{tumour_name}.purple_alternates/

    zip -r ~{tumour_name}.purple_alternates.zip ~{tumour_name}.purple_alternates/


  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File purple_alternate_directory = "~{tumour_name}.purple_alternates.zip"
  }
}


task LINX {
  input {
    File? purple_dir
    String tumour_name
    String ensemblDir
    String genomeVersion
    String fusions_file
    String linxScript = "java -Xmx32G -cp $HMFTOOLS_ROOT/linx.jar com.hartwig.hmftools.linx.LinxApplication"
    String modules
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    ensemblDir: "Directory of Ensembl data for PURPLE"
    genomeVersion: "genome version for AMBER, default set to V38"
    linxScript: "location of LINX script"
    modules: "Required environment modules"
		memory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

  command <<<
    set -euo pipefail

    unzip ~{purple_dir} 
    mkdir ~{tumour_name}.linx 
     
    ~{linxScript} \
      -sample ~{tumour_name} \
      -ref_genome_version ~{genomeVersion} \
      -ensembl_data_dir ~{ensemblDir}  \
      -check_fusions \
      -known_fusion_file ~{fusions_file} \
      -purple_dir ~{tumour_name}.solPrimary.purple \
      -output_dir ~{tumour_name}.linx 

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File linx_fusions = "~{tumour_name}.linx/~{tumour_name}.linx.fusion.tsv"
    File linx_svs = "~{tumour_name}.linx/~{tumour_name}.linx.svs.tsv"
    File drivers_svs = "~{tumour_name}.linx/~{tumour_name}.linx.drivers.tsv"
  }

  meta {
      output_meta: {
      linx_fusions: "tab seperated LINX fusions",
      linx_svs: "tab seperated LINX Structural Variants",
      drivers_svs: "tab seperated LINX Driver Structural Variants"
      }
  }
}
