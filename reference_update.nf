#!/usr/bin/env nextflow
nextflow.enable.dsl=1

VERSION = "0.1.1"

process create_orthodb {
  afterScript 'rm *.fa'

  input:
    val orthodb_path from "${params.orthodb_path}"

  output:
    path "orthodb.fasta" into orthodb_fasta

  """
  touch orthodb.fasta

  url_regex='^(https?|ftp|file)://[-A-Za-z0-9\\+&@#/%?=~_|!:,.;]*[-A-Za-z0-9\\+&@#/%=~_|]\\.[-A-Za-z0-9\\+&@#/%?=~_|!:,.;]*[-A-Za-z0-9\\+&@#/%=~_|]\$'

  if [[ ${orthodb_path} =~ \$url_regex ]]; then wget ${orthodb_path}; else ln -s ${orthodb_path}; fi

  if [[ `ls *.fa.gz 2> /dev/null` ]]; then gunzip *.fa.gz ; cat *.fa > orthodb.fasta; fi
  """
}

if (!params.from_local){
  process get_reference_species {
      conda 'environment.yml'

      output:
        stdout species_ch

      """
      IFS=\$'\\n'; ${params.EUPATHWS_SCRIPTS_PATH}/get_reference_species -l ${params.veupathdb_username}:${params.veupathdb_password} ${params.veupathdb_domain}
      """
  }

  species_ch.splitText().map{it -> it.trim()}.set { org }
  process get_organism {
      conda 'environment.yml'
      errorStrategy 'retry'
      maxRetries 5
      maxForks 8
      publishDir "${params.REFERENCE_PATH}", mode: 'copy'

      input:
        val org
        path "orthodb.fasta" from orthodb_fasta

      output:
        path 'clean_gff/*.gff3' into org_gff3
        path '*_Proteins.fasta' into org_prots
        path '*_Genome.fasta' into org_fasta
        path '*.gaf'
        path '*_chr.json' into org_chrs
        path '*_meta.json' into org_meta
        stdout org_ch

      """
      ${params.EUPATHWS_SCRIPTS_PATH}/get_organism -l ${params.veupathdb_username}:${params.veupathdb_password} ${params.veupathdb_domain} \"${org}\"
      rename "s/-| /_/g"  *
      rename "s/[\\)\\(]//g" *
      mkdir clean_gff
      for x in *.gff3 ; do gt gff3 -sort -retainids -tidy \$x > clean_gff/\$x & done
      ls *.gff3 | sed 's/.gff3//g'
      """
  }

  org_prots.concat(orthodb_fasta).collectFile(name: "all_annotated_proteins.fasta", storeDir: "${params.REFERENCE_PATH}")

} else {
  species_ch = Channel.fromPath( String.format( "%s/*", params.from_local ) )
  process get_all_local_organisms {
    publishDir "${params.REFERENCE_PATH}", mode: 'copy'

    input:
      path "*" from species_ch.collect()
      path "orthodb.fasta" from orthodb_fasta

    output:
      path 'clean_gff/*.gff3' into org_gff3
      path '*_Proteins.fasta', includeInputs: true
      path '*_Genome.fasta', includeInputs: true into org_fasta
      path '*.gaf', includeInputs: true
      path '*_chr.json', includeInputs: true into org_chrs
      path '*_meta.json', includeInputs: true into org_meta
      path 'all_annotated_proteins.fasta'
      stdout org_ch

    """
    rename "s/-| /_/g"  *
    rename "s/[\\)\\(]//g" *
    mkdir clean_gff
    for x in *.gff3 ; do gt gff3 -sort -retainids -tidy \$x > clean_gff/\$x & done
    cat *_Proteins.fasta orthodb.fasta > all_annotated_proteins.fasta
    ls *.gff3 | sed 's/.gff3//g'
    """
  }
}

org_fasta.into{
  org_fasta_augustus
  org_fasta_validation
}

org_ch.into{
  augustus_org_ch
  validation_org_ch
}

if (params.do_augustus) {
  augustus_org_ch.splitText().map{it -> it.trim()}.set { org }
  process train_augustus {
      input:
        val x from org
        path "*" from org_gff3
        path "*" from org_fasta_augustus

      output:
        path "${x}.gff3" into gff3_out

      """
      AUGUSTUS_CONFIG_PATH=${params.AUGUSTUS_CONFIG_PATH} \
      ${params.AUGUSTUS_SCRIPTS_PATH}/new_species.pl --species=${x} --ignore
      ${params.AUGUSTUS_SCRIPTS_PATH}/gff2gbSmallDNA.pl ${x}.gff3 ${x}_Genome.fasta 1000 sampled.gb
      ${params.AUGUSTUS_BIN_PATH}/etraining --species=${x} sampled.gb
      """
  }
} else {
  org_gff3.set { gff3_out }
}

process prepare_references {
    input:
      path "*" from gff3_out.collect()
      path "*" from org_chrs.collect()
      path "*" from org_meta.collect()

    output:
      file "references-in-*.json" into references_in
      path "groups.txt" into groups

    """
    for x in *.gff3; do python3 ${baseDir}/bin/parse_chromosomes.py \$x; done > ChromosomeFile.txt
    ls *.gff3 | awk -F'[_]' '{print \$1}' | sort | uniq -c | awk '{if (\$1 > 0) print \$2}' > groups.txt
    for x in `cat groups.txt`; do ls \$x*.gff3 | perl ${baseDir}/bin/generateReferences-i.json.pl ${params.REFERENCE_PATH} \$x ${params.AUGUSTUS_CONFIG_PATH} ${VERSION} > references-in-\$x.json; done
    """
}

if (!params.do_all_vs_all){
  groups.splitText().map{it -> it.trim()}.set { group }
  process update_references{
      debug true
      publishDir "${params.REFERENCE_PATH}/Ref_${x}", mode: 'copy',
        saveAs: { filename ->
          if(filename.startsWith("${x}") & filename.endsWith(".fasta")) {
            null
          } else {
            filename
          }
        }

      input:
        val x from group
        file "*" from references_in

      output:
        path "*", type: "dir"
        path "references.json" into refs
        tuple val(x), path("${x}*.fasta") into prots

      """
      cp references-in-${x}.json references-in.json
      ${params.COMPANION_BIN_PATH}/update_references.lua
      for y in ${x}_*/proteins.fasta
      do
        species="\$(dirname \$y)"
        cp \$y \$species.fa
        ${params.COMPANION_BIN_PATH}/adjust_fasta_header.pl \$species \$species.fa 1
      done
      """
  }
} else {
  process update_all_references{
    debug true
    publishDir "${params.REFERENCE_PATH}/Reference", mode: 'copy'

    input:
      file "*" from references_in

    output:
      path "*", type: "dir"
      path "references.json" into refs
      tuple val(null), path("*.fasta") into prots

    """
    echo "{\\"groups\\":{},\\"species\\":{}}" >> references.tmp
    for x in references-in-*.json; do cp \$x references-in.json; ${params.COMPANION_BIN_PATH}/update_references.lua; jq -s '.[0] * .[1]' references.tmp references.json | sponge references.tmp; done
    for y in */proteins.fasta; do cp \$y "\$(dirname \$y).fasta"; done
    mv references.tmp references.json
    """
  }
}

if (params.do_orthofinder) {
  process run_orthofinder {
    cpus params.diamond_threads
    debug true
    publishDir params.do_all_vs_all ? "${params.REFERENCE_PATH}/Reference/_all" : "${params.REFERENCE_PATH}/Ref_${grp}/_all", mode: 'copy'

    input:
      tuple val(grp), path("0.input_faa/*") from prots

    output:
      path "all_orthomcl.out"

    """
    orgs=`ls 0.input_faa | wc -l`
    if [[ \$orgs -gt 1 ]]
    then
      ${params.ORTHOFINDER_PATH}/orthofinder.py -f 0.input_faa/ -o results -t ${task.cpus}

      # filter out clusters with single gene.
      awk 'BEGIN { FS="[ ]" }; { if (\$3) print \$0 }' results/*/Orthogroups/Orthogroups.txt > all_orthomcl.out
    else
      touch all_orthomcl.out
    fi
    """
  }
} else {
  process empty_orthogroups {
    publishDir params.do_all_vs_all ? "${params.REFERENCE_PATH}/Reference/_all" : "${params.REFERENCE_PATH}/Ref_${grp}/_all", mode: 'copy'

    input:
      tuple val(grp), path from prots

    output:
      path "all_orthomcl.out"

    """
    touch all_orthomcl.out
    """
  }
}

if (params.validate_refs) {
  if (params.validate_refs == "sampled") {
    process sample_fasta_for_validation {
      errorStrategy 'ignore'

      input:
        path "*" from org_fasta_validation

      output:
        path '*.sampled.fa' into org_fasta_sampled

      """
      in=\$(basename *_Genome.fasta)
      out=\$in.sampled.fa
      seqkit seq -m 10000 -M 500000 \$in > 1
      reformat.sh in=1 out=\$out reads=2 minlength=10000 overwrite=true

      if [ ! -s \$out ]
      then
        exit 1
      fi
      """
    }
  } else if (params.validate_refs == "full") {
    process full_fasta_for_validation {
      input:
        path "*" from org_fasta_validation

      output:
        path '*.sampled.fa' into org_fasta_sampled

      """
      in=\$(basename *_Genome.fasta)
      out=\$in.sampled.fa
      cp \$in \$out
      """
    }
  }

  process prepare_validation {
    cache false
    beforeScript """rm -f ${ref_d}/failed_refs.txt"""

    if (params.do_all_vs_all) {
      ref_dir = "${params.REFERENCE_PATH}/Reference"
    } else {
      ref_dir = "${params.REFERENCE_PATH}"
    }

    input:
      val ref_d from ref_dir
      path "references.json" from refs.collect()

    output:
      val ref_d into ref_dir

    """
    touch ${ref_d}/failed_refs.txt
    touch ${ref_d}/good_refs.txt
    """
  }

  ref_dir.into{
    ref_dir_val
    ref_dir_filt
  }

  validation_org_ch.splitText().map{it -> it.trim()}.set { org }
  process run_validation {
    errorStrategy 'ignore'
    afterScript """
    if [[ \${nxf_main_ret:=0} != 0 ]]
    then
      echo ${ref_species} >> ${ref_d}/failed_refs.txt
    else
      echo ${ref_species} >> ${ref_d}/good_refs.txt
    fi
    """

    input:
      path "*" from org_fasta_sampled.collect()
      val ref_species from org
      val ref_d from ref_dir_val

    output:
      path ".exitcode" into exitcode

    """
    if [ ! -d ${ref_d}/${ref_species} ]
    then
      ref_dir="\$(dirname ${ref_d}/*/${ref_species})"
      echo \$ref_dir
    else
      ref_dir=${ref_d}
    fi

    nextflow run ${params.COMPANION_SCRIPT_PATH} -profile docker -c ${baseDir}/params_companion.config -with-trace -resume \
     --inseq ${ref_species}*.sampled.fa \
     --ref_dir \$ref_dir \
     --ref_species ${ref_species} \
    """
  }

  process filter_failed_refs {
    input:
      path "*" from exitcode.collect()
      val ref_d from ref_dir_filt

    """
    for x in ${ref_d}/failed_refs.txt; do
      refs=`grep -Ril ${params.REFERENCE_PATH}/*/references.json -e '\$x'`
      jq '.species.\$x += { "validated": false }' \$refs
    done
    """
  }
}
