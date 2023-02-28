version 1.0

struct GenomeResources {
    String PON
    String amberModules
    String gcProfile
    String cobaltModules
    String ensemblDir
    String refFasta
    String gcProfile
    String runPURPLEModules
}

workflow purple {
  input {
    File tumour_bam
    File tumour_bai
    File normal_bam
    File normal_bai
    String normal_name
    String tumour_name
    String genomeVersion = "V38"
    Boolean doSV = true
    Boolean doSMALL = true
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

Map[String,GenomeResources] resources = {
  "V38": {
    "amberModules": "hmftools/1.1 hg38/p12 hmftools-data/hg38",
    "cobaltModules": "hmftools/1.1 hg38/p12 hmftools-data/hg38",
    "runPURPLEModules": "hmftools/1.1 hg38/p12 hmftools-data/hg38",
    "PON" : "$HMFTOOLS_DATA_ROOT/GermlineHetPon.38.vcf.gz",
    "gcProfile": "$HMFTOOLS_DATA_ROOT/GC_profile.1000bp.38.cnp",
    "ensemblDir": "$HMFTOOLS_DATA_ROOT/ensembl",
    "refFasta": "$HMFTOOLS_DATA_ROOT/hg38_random.fa",
    "gcProfile": "$HMFTOOLS_DATA_ROOT/GC_profile.1000bp.38.cnp"
  }
}

  call amber {
    input:
      tumour_bam = tumour_bam,
      tumour_bai = tumour_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = normal_name,
      tumour_name = tumour_name,
      genomeVersion = genomeVersion,
      modules = resources [ genomeVersion ].amberModules,
      PON = resources [ genomeVersion ].PON
  }

  call cobalt {
    input:
      tumour_bam = tumour_bam,
      tumour_bai = tumour_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = normal_name,
      tumour_name = tumour_name,
      modules = resources [ genomeVersion ].cobaltModules,
      gcProfile = resources [ genomeVersion ].gcProfile
  }

  if(doSV) {
    call filterSV {
      input: 
        tumour_name = tumour_name
    }
  }

  if(doSMALL) {
    call filterSMALL {
      input: 
        tumour_name = tumour_name
    }
  }

  call runPURPLE {
    input:
      normal_name = normal_name,
      tumour_name = tumour_name,
      amber_directory = amber.output_directory,
      cobalt_directory = cobalt.output_directory,
      SV_vcf = filterSV.filtered_vcf,
      smalls_vcf = filterSMALL.filtered_vcf,
      genomeVersion = genomeVersion,
      modules = resources [ genomeVersion ].runPURPLEModules,
      gcProfile = resources [ genomeVersion ].gcProfile,
      ensemblDir = resources [ genomeVersion ].ensemblDir,
      refFasta = resources [ genomeVersion ].refFasta,
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
    output_meta: {
			purple_directory: "zipped directory containing all PURPLE output"
		}
  }

  output {
    File purple_qc = runPURPLE.purple_qc
    File purple_purity = runPURPLE.purple_purity
    File purple_purity_range = runPURPLE.purple_purity_range
    File purple_segments = runPURPLE.purple_segments
    File purple_cnv = runPURPLE.purple_cnv
    File? purple_SV_index = runPURPLE.purple_SV_index
    File? purple_SV = runPURPLE.purple_SV
    File? purple_SMALL_index = runPURPLE.purple_SMALL_index
    File? purple_SMALL = runPURPLE.purple_SMALL

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
    genomeVersion: "genome version for AMBER, default set to V38"
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
    String gamma = 100
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

  meta {
		output_meta: {
			output_directory: "Zipped COBALT results directory"
		}
	}

}

task filterSV {
  
  input {
    String tumour_name
    File? vcf
    String bcftoolsScript = "$BCFTOOLS_ROOT/bin/bcftools view"
    String modules = "bcftools"
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    vcf: "VCF file for filtering"
    bcftoolsScript: "location of BCFTOOLS script, including view command"
		modules: "Required environment modules"
		memory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
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

  meta {
		output_meta: {
			filtered_vcf: "Filtered VCF"
		}
	}

}

task filterSMALL {
  
  input {
    String tumour_name
    File? vcf
    String bcftoolsScript = "$BCFTOOLS_ROOT/bin/bcftools view"
    String modules = "bcftools"
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    vcf: "VCF file for filtering"
    bcftoolsScript: "location of BCFTOOLS script, including view command"
		modules: "Required environment modules"
		memory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
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

  meta {
		output_meta: {
			filtered_vcf: "Filtered VCF"
		}
	}

}

task runPURPLE {
  input {
    String normal_name
    String tumour_name
    File amber_directory
    File cobalt_directory
    File? SV_vcf
    File? smalls_vcf
    String ensemblDir
    String refFasta
    String genomeVersion
    String gcProfile 
    String purpleScript = "java -Xmx8G -jar $HMFTOOLS_ROOT/purple.jar"
    String modules
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    normal_name: "Name for Normal sample"
    amber_directory: "zipped output from AMBER"
    cobalt_directory: "zipped output from COBALT"
    SV_vcf: "filtered structural variant (SV) vcf"
    smalls_vcf: "filtered SNV and indel (smalls) vcf"
    ensemblDir: "Directory of Ensembl data for PURPLE"
    refFasta: "fasta of reference genome"
    gcProfile: "GC profile, generated for COBALT"
    genomeVersion: "genome version for AMBER, default set to V38"
    purpleScript: "location of PURPLE script"
    modules: "Required environment modules"
		memory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

  String SV_vcf_arg = if(defined(SV_vcf))then
                                  "-structural_vcf ~{SV_vcf}"
                                 else
                                  ""

  String smalls_vcf_arg = if(defined(smalls_vcf))then
                                  "-somatic_vcf ~{smalls_vcf}"
                                 else
                                  ""

  command <<<
    set -euo pipefail

    unzip ~{amber_directory} 
    unzip ~{cobalt_directory} 
    mkdir purple_output 

    ~{purpleScript} \
      -ref_genome_version ~{genomeVersion} \
      -ref_genome ~{refFasta}  \
      -gc_profile ~{gcProfile} \
      -ensembl_data_dir ~{ensemblDir}  \
      -reference ~{normal_name} -tumor ~{tumour_name}  \
      -amber ~{tumour_name}.amber -cobalt ~{tumour_name}.cobalt \
      ~{SV_vcf_arg} ~{smalls_vcf} \
      -output_dir purple_output/

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File purple_qc = "purple_output/~{tumour_name}.purple.qc"
    File purple_purity = "purple_output/~{tumour_name}.purple.purity.tsv"
    File purple_purity_range = "purple_output/~{tumour_name}.purple.purity.range.tsv"
    File purple_segments = "purple_output/~{tumour_name}.purple.segment.tsv"
    File purple_cnv = "purple_output/~{tumour_name}.purple.cnv.somatic.tsv"
    File? purple_SV_index = "purple_output/~{tumour_name}.purple.sv.vcf.gz.tbi"
    File? purple_SV = "purple_output/~{tumour_name}.purple.sv.vcf.gz"
    File? purple_SMALL_index = "purple_output/~{tumour_name}.purple.somatic.vcf.gz.tbi"
    File? purple_SMALL = "purple_output/~{tumour_name}.purple.somatic.vcf.gz"
  }

  meta {
		output_meta: {
			purple_directory: "Zipped Output PURPLE directory"
		}
	}
}
