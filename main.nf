#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
         ELEMENTO - MELNICK LABS ||  ATAC-SEQ   BEST   PRACTICE
========================================================================================
#### Authors
Ashley Stephen Doane <asd2007@med.cornell.edu>

*/
//params.index = '/home/asd2007/Scripts/nf/fripflow/sindex.tsv'



version = 1.1

// SET PARAMS

params.name = 'ecadd4547'
params.index = 'indexpooled.csv'
   //params.index = 'indexpooled.tsv'
       //params.index = 'Sample_Ly7_pooled_100k'
       //params.index = "$baseDir/indexecad.csv"
//params.index = "$baseDir/indexTest.csv"
       //params.index = 'sampleIndexjc.csv'
       // Configurable variables
       params.name = false
       params.project = false
       params.genome = 'hg38'
       params.genomes = []
       params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
       params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
       params.blacklist = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
              //genome = file(params.genome)
              //index = file(params.index)
       params.chrsz = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes"
              //params.bwa_index = "/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
              //params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
              //params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
       params.notrim = false
       params.saveReference = true
       params.saveTrimmed = false
       params.saveAlignedIntermediates = false
       params.broad = false
       params.outdir = "$baseDir/resultsNB"
       params.email = 'ashley.doane@gmail.com'
       params.chromsizes = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes"
       params.lncaprefpeak = "$baseDir/data/lncapPeak.narrowPeak" 
       params.bcellrefpeak = "$baseDir/data/gcb.tn5.broadPeak"
       params.nucleoatac = true

/*
params.DNASE_BED="/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
params.BLACK="/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
params.PICARDCONF="/athena/elementolab/scratch/asd2007/reference/hg38/picardmetrics.conf"
params.REF="/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.REFbt2="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.bwt2_idx="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.RG="hg38"
params.SPEC="hs"
params.REFGen="/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.seq="/athena/elementolab/scratch/asd2007/reference/hg38/seq"
params.gensz="hs"
params.bwt2_idx="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.REF_FASTA="/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.species_browser='hg38'
params.TSS_ENRICH="/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_gencode_tss_unique.bed.gz"
params.DNASE='/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz'
params.PROM='/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz'
params.ENH='athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz'
params.REG2MAP='/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_dnase_avg_fseq_signal_formatted.txt.gz'
params.REG2MAP_BED="/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_celltype_compare_subsample.bed.gz"
params.ROADMAP_META="/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt"
*/

if( params.bwa_index ){
    Channel
        .fromPath(params.bwa_index)
        .ifEmpty { exit 1, "BWA index not found: ${params.bwa_index}" }
    .into { bwa_index; refgenome }
} else if ( params.fasta ){
    fasta = file(params.fasta)
        if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
                              } else {
    exit 1, "No reference genome specified!"
}



custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}

// Header log info
log.info "=================================================================="
log.info " Elemento-Melnick-Labs-ATACseq: ATAC-Seq Best Practice v${version}"
log.info "==================================================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
       //summary['Reads']        = params.reads
summary['Genome']       = params.genome
if(params.bwa_index)  summary['BWA Index'] = params.bwa_index
else if(params.fasta) summary['Fasta Ref'] = params.fasta
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Save Reference'] = params.saveReference
summary['Save Trimmed']   = params.saveTrimmed
summary['Save Intermeds'] = params.saveAlignedIntermediates
if(params.email) summary['E-mail Address'] = params.email
if(workflow.commitId) summary['Pipeline Commit']= workflow.commitId
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "===================================="

index = file(params.index)

       //bwaref = file(params.ref)
       results_path = file(params.outdir)

blacklist=file(params.blacklist)

/*
 * PREPROCESSING - Build BWA index
 */

////// Check input parameters //////

if (!params.index) {
  exit 1, "Please specify the input table file"
}


sizefactors = Channel.from(1)



       /*
 *atacs = Channel
 *       .from(index.readLines())
 *       .map { line ->
 *       def list = line.split()
 *              def bed = file(list[0])
 *              def peaks = file(list[1])
 *              def sname = list[2]
 *              def dprefix = file(list[3])
 *       println bed
 *       println peaks
 *       [ sname, bed, peaks, dprefix ]
 *}
 */



bcellrefpeaks = file(params.bcellrefpeak)
lncaprefpeaks = file(params.lncaprefpeak)
chrsz = file(params.chrsz)
dnase = file(params.DNASE)
tss_enrich = file(params.TSS_ENRICH)
prom = file(params.PROM)
enh = file(params.ENH)
reg2map = file(params.REG2MAP)
roadmap_meta = file(params.ROADMAP_META)

  Channel
      .from(index.readLines())
      .map { line ->
             def list = line.split(',')
             def Sample = list[0]
             def samplebam = file(list[1])
             def sizef = list[2]
             [ Sample, samplebam, sizef ]
}
    .set{ bam } 










/**
process processbam {
    tag "$Sample"
    publishDir "$results_path/$Sample/$Sample", mode: 'copy'

    executor 'sge'
    clusterOptions '-l h_vmem=4G -pe smp 8 -l h_rt=16:00:00 -l athena=true'
    scratch true
    // cpus 8


    input:
    set Sample, file(nbam) from newbam
    file(BLACK) from blacklist

    output:
    set Sample, file("${Sample}.sorted.bam") into sortedbam
    set Sample, file("${Sample}.sorted.bam") into sortedbamqc
    set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into finalbam
    set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into finalbamforqc
    set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into bamforsignal
    set Sample, file("${Sample}*.pbc.qc") into pbcqc
    set Sample, file("${Sample}*.dup.qc") into dupqc
    file("*nsort.fixmate.bam") into fixmatebam
    file("*window500.hist_data") into hist_data
    file("*window500.hist_graph.pdf") into fragsizes
    set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbam
        // set Sample, file("${Sample}.sorted.nodup.noM.black.bam"), file("${Sample}.sorted.nodup.noM.black.bam.bai") into bamforsignal
    set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbamforqc

    script:
    """
    #!/bin/bash -l
    set -o pipefail
    spack load jdk
    processAlignment.nf.sh ${nbam} ${BLACK} 8
    """
}
*/

//hist_data.subscribe { println "Received: " + file(hist_data)}

//fragsizes.subscribe { println "Received: " + file(fragsizes)}




process nosup {
     tag "$Sample"
     publishDir  "$results_path/$Sample/$Sample", mode: 'copy'
     cpus 8

     input:
     set Sample, file(finalbam), val(sz) from bam

     output:
     set Sample, file('*.nosup.bam') into bamnosup
     set Sample, file('*.nosup.bam') into bamfornuc
     set Sample, file('*.nosup.bam'), val(sz) into bamforsignal


    script:
    """
    nosup.bam.sh ${finalbam} ${task.cpus}
    """
            }

process bam2bed {
    tag "$Sample"
    publishDir  "$results_path/$Sample/$Sample", mode: 'copy'

    cpus 8

    input:
    set Sample, file(finalbamnosup) from bamnosup

    output:
    set Sample, file("${Sample}.nodup.tn5.tagAlign.gz") into finalbedqc
    set Sample, file("${Sample}.nodup.tn5.tagAlign.gz"), file(finalbamnosup) into finalbed
    set Sample, file("${Sample}.nodup.bedpe.gz") into finalbedpe


    script:
    """
    #!/bin/bash -l
    sambamba sort --memory-limit 86GB -n -t ${task.cpus} --out ${Sample}.nsorted.nodup.noM.bam ${finalbamnosup}
    samtools fixmate -@8 ${Sample}.nsorted.nodup.noM.bam ${Sample}.nsorted.fixmate.nodup.noM.bam
    convertBAMtoBED.sh ${Sample}.nsorted.fixmate.nodup.noM.bam
    cp ${Sample}.nsorted.fixmate.nodup.noM.tn5.tagAlign.gz  ${Sample}.nodup.tn5.tagAlign.gz
    cp ${Sample}.nsorted.fixmate.nodup.noM.bedpe.gz ${Sample}.nodup.bedpe.gz
    """
        }




process callpeaks {
    tag "$Sample"
    publishDir  "$results_path/$Sample/$Sample", mode: 'copy'
        // conda '/athena/elementolab/scratch/asd2007/projectshg38/datasets/atacseq/bcellpooled_atacflow/idp2.yml'
        //conda "/home/asd2007/intel/intelpython3/envs/idp2"
        //conda 'bioconda:macs2 python=2.7'
        // conda '/home/asd2007/intel/intelpython3/envs/idp2'

    input:
    set Sample, file(rbed), file(bamfile) from finalbed

    output:
    set Sample, file("${Sample}.tn5.narrowPeak.gz") into narrowpeak
    set Sample, file("${Sample}.narrowPeak") into narrowpeaks
    set Sample, file("${Sample}.tag.narrow_summits.bed") into summits
    set Sample, file("${Sample}.tn5.broadPeak.gz") into broadpeak
    set Sample, file("${Sample}.tn5.broadPeak.gz"), file("${Sample}.sorted.nodup.noM.black.nosup.bam") into broadPeaknuc

    script:
    """
    #!/bin/bash -l
    source activate atacFlow
    macs2 callpeak -t ${rbed} -f BED -n ${Sample}.tag.narrow -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-3

    cp ${Sample}.tag.narrow_peaks.narrowPeak ${Sample}.narrowPeak

    /home/asd2007/ATACseq/narrowpeak.py ${Sample}.tag.narrow_peaks.narrowPeak ${Sample}.tn5.narrowPeak 

    macs2 callpeak -t  ${rbed} -f BED -n ${Sample}.tag.broad -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --broad --broad-cutoff 0.1

    /home/asd2007/ATACseq/broadpeak.py  ${Sample}.tag.broad_peaks.broadPeak ${Sample}.tn5.broadPeak 
    """
        }


process signalTrack {
    tag "$Sample"

    publishDir "$results_path/$Sample/$Sample", mode: 'copy'

        //executor 'sge'
        //penv 'smp'
        //memory { 3.GB * task.attempt }
        //clusterOptions '-l h_rt=48:00:00 -l athena=true'
        //scratch true
    cpus 12
        //errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }

    input:
    set Sample, file(sbam), val(sz) from bamforsignal
        //val sz from sizefactors


    output:
        set Sample, file("${Sample}.deseqsf.bw") into insertionTrackbw

    script:
    """
    #!/bin/bash
    . /home/asd2007/.spackloads.sh
    ##source activate deep
    getbamcov.sh ${sbam} ${Sample} ${sz} ${task.cpus}
    """
    }


finalbedpe.mix(broadpeak)
    .groupTuple(sort: true)
    .set{ fripin }

/*
process frip {
    tag "$Sample"

    publishDir "$results_path/$Sample/frip", mode: 'copy'

    input:
    set Sample, file(file_list) from fripin
        //set Sample, file(peaks) from broadpeak
    file(lncapref) from lncaprefpeaks
    file(bcellref) from bcellrefpeaks

    output:
    set Sample, file("${Sample}.frip.txt") into frips
    set Sample, file("${Sample}.lncapref.frip.txt") into frips2
    set Sample, file("${Sample}.bcellref.frip.txt") into frips3

    script:
    """
    getFripQC.py --bed ${Sample}.nodup.bedpe.gz --peaks ${Sample}.tn5.broadPeak.gz --out ${Sample}.frip.txt

    getFripQC.py --bed ${Sample}.nodup.bedpe.gz --peaks ${lncapref} --out ${Sample}.lncapref.frip.txt

    getFripQC.py --bed ${Sample}.nodup.bedpe.gz --peaks ${bcellref} --out ${Sample}.bcellref.frip.txt
    """
}
*/




//group broadPeak and bam files
//bamfornuc.mix(broadPeaknuc)
//    .groupTuple(sort: true)
//    .view()
//    .set{ grouped_peak_bam_ch }



if (params.nucleoatac == true) {

    process nucleoatac {
        tag "$Sample"
        publishDir "$results_path/$Sample/nucleoAtac", mode: 'copy'
            //conda 'intel::cython=0.28.2=py27_intel_4 numpy=1.14.3=py27_intel_0'
        // conda '/home/asd2007/intel/intelpython3/envs/nuc'
        //conda '/athena/elementolab/scratch/asd2007/projectshg38/datasets/atacseq/bcellpooled_atacflow/nuc.yml'
        //conda '$baseDir/nucnextflow.yml'
        cpus 24
        //executor 'sge'
        //penv 'smp'
        //clusterOptions '-l h_vmem=4G -l h_rt=96:00:00 -l athena=true -R y'
        //scratch true

        input:
        file chrsz from chrsz
        file index from refgenome.first()
        set Sample, file(peaks), file(bam) from broadPeaknuc
            //    set Sample, file(filelist) from grouped_peak_bam_ch

        output:
        file '*.{bed.gz,txt,eps,gz,bw,tbi}' into nucleoAtac_results

        script:
            //sampleId = replicateId.replaceAll(/[12]$/,'')
            //prefix = bam[0].toString() - ~/(\.sorted)?(\.nodup)?(\.noM)?(\.black)?(\.bam)?$/
        """
        #!/bin/bash -l
        ##spack load samtools
        ##spack load bedtools2
        ##myspack
        source activate nuc
        samtools index ${Sample}.sorted.nodup.noM.black.nosup.bam
        gunzip -cf ${Sample}.tn5.broadPeak.gz > ${Sample}.tn5.broadPeak

        bedtools slop -i ${Sample}.tn5.broadPeak -g ${chrsz} -b 1000 > ${Sample}.slop1k.bed

        run_nucleoatac.nf.sh ${Sample}.sorted.nodup.noM.black.nosup.bam ${Sample}.slop1k.bed ${Sample} ${task.cpus}

        """
            }
}



/*
process nucleoatac {
    tag "$Sample"

    publishDir "$results_path/$Sample/nucleoAtac", mode: 'copy'

    executor 'sge'
    clusterOptions '-l h_vmem=4G -pe smp 24-32 -l h_rt=96:00:00 -l athena=true -R y'
    scratch true

    input:
    file chrsz from chrsz
    file index from refgenome.first()
    set Sample, file(broadPeak), file(bam) from grouped_peak_bam_ch 

    output:
    file '*.{bed.gz,txt,eps,gz,bw,tbi}' into nucleoAtac_results


    script:
        //sampleId = replicateId.replaceAll(/[12]$/,'')
        //prefix = bam[0].toString() - ~/(\.sorted)?(\.nodup)?(\.noM)?(\.black)?(\.bam)?$/
    """
    source activate idp2
    spack load samtools
    samtools index ${bam}
    gunzip ${Sample}.tn5.broadPeak.gz

    bedtools slop -i ${Sample}.tn5.broadPeak -g ${chrsz} -b 1000 > ${Sample}.slop1k.bed

    nucleoatac run --bed ${Sample}.slop1k.bed \\
        --bam ${bam} \\
        --fasta ${index}/genome.fa \\
        --out ${Sample} \\
        --write_all \\
        --cores \${NSLOTS}
    """

}

*/





workflow.onComplete {
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
    def subject = 'pipeline execution'
    def recipient = 'ashley.doane@gmail.com'

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}




/*
 * Completion e-mail notification
 *
 *workflow.onComplete {
 *
 *    // Set up the e-mail variables
 *    def subject = "[OElab-ATACseq] Successful: $workflow.runName"
 *    if(!workflow.success){
 *      subject = "[OElab-ATACseq] FAILED: $workflow.runName"
 *    }
 *    def email_fields = [:]
 *    email_fields['version'] = version
 *    email_fields['runName'] = custom_runName ?: workflow.runName
 *    email_fields['success'] = workflow.success
 *    email_fields['dateComplete'] = workflow.complete
 *    email_fields['duration'] = workflow.duration
 *    email_fields['exitStatus'] = workflow.exitStatus
 *    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
 *    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
 *    email_fields['commandLine'] = workflow.commandLine
 *    email_fields['projectDir'] = workflow.projectDir
 *    email_fields['summary'] = summary
 *    email_fields['summary']['Date Started'] = workflow.start
 *    email_fields['summary']['Date Completed'] = workflow.complete
 *    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
 *    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
 *    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp
 *    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
 *    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
 *    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
 *    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
 *    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
 *    if(workflow.container) email_fields['summary']['Singularity image'] = workflow.container
 *
 *    // Render the TXT template
 *    def engine = new groovy.text.GStringTemplateEngine()
 *    def tf = new File("$baseDir/assets/email_template.txt")
 *    def txt_template = engine.createTemplate(tf).make(email_fields)
 *    def email_txt = txt_template.toString()
 *
 *    // Render the HTML template
 *    def hf = new File("$baseDir/assets/email_template.html")
 *    def html_template = engine.createTemplate(hf).make(email_fields)
 *    def email_html = html_template.toString()
 *
 *    // Render the sendmail template
 *    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
 *    def sf = new File("$baseDir/assets/sendmail_template.txt")
 *    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
 *    def sendmail_html = sendmail_template.toString()
 *
 *    // Send the HTML e-mail
 *    if (params.email) {
 *        try {
 *          // Try to send HTML e-mail using sendmail
 *          [ 'sendmail', '-t' ].execute() << sendmail_html
 *          log.debug "[NGI-ChIPseq] Sent summary e-mail using sendmail"
 *        } catch (all) {
 *          // Catch failures and try with plaintext
 *          [ 'mail', '-s', subject, params.email ].execute() << email_txt
 *          log.debug "[OElab-ATACseq] Sendmail failed, failing back to sending summary e-mail using mail"
 *        }
 *        log.info "[OElab-ATACseq] Sent summary e-mail to $params.email"
 *    }
 *
 *    // Write summary e-mail HTML to a file
 *    def output_d = new File( "${params.outdir}/Documentation/" )
 *    if( !output_d.exists() ) {
 *      output_d.mkdirs()
 *    }
 *    def output_hf = new File( output_d, "pipeline_report.html" )
 *    output_hf.withWriter { w -> w << email_html }
 *    def output_tf = new File( output_d, "pipeline_report.txt" )
 *    output_tf.withWriter { w -> w << email_txt }
 *
 *    log.info "[OElab-ATACseq] Pipeline Complete"
 *}
 */

