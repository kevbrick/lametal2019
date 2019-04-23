#!/usr/bin/env nextflow

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "LAM ET AL. 2019 PIPELINE (Version 1.0)                                    "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run pipe_LamEtAl.groovy \\"
  log.info " -with-trace -with-timeline -with-report"
  log.info " "
  log.info "HELP: nextflow run pipe_LamEtAl.groovy --help"
  log.info " "
  log.info "================================================================================================================="
  log.info "Required Arguments:"
  log.info " --projectdir    Lam et al. root dir (default = .)"
  log.info " "
  log.info "Output Arguments"
  log.info "  --outdir       Output dir "
  log.info " "
  log.info "================================================================================================================"
  exit 1
  }

// Define arg defaults

//number of threads
params.threads        = 6
params.mem            = "16G"
params.debugmode      = ""

//output and tempory directories
params.projectdir       = "."
params.tmpdir           = "/lscratch/\$SLURM_JOBID"
params.codedir          = "${params.projectdir}/accessoryFiles/scripts/"
params.codeRdir         = "${params.projectdir}/accessoryFiles/scripts/R"
params.datadir          = "${params.projectdir}/accessoryFiles/data/"
params.genomedir        = "${params.projectdir}/accessoryFiles/genomeFiles/"
params.timecoursedir    = "${params.projectdir}/accessoryFiles/data/timeCourse/"
params.allHMdir         = "${params.projectdir}/accessoryFiles/data/allHistoneMods/"
params.rnadir           = "${params.projectdir}/accessoryFiles/RNASeq/"
params.blacklist        = "${params.projectdir}/accessoryFiles/data/blacklist/mm10_hotspot_blackList.bed"
params.table1Init       = "${params.projectdir}/accessoryFiles/data/table1/table1Data.tab"

params.outdir      		  = ""
params.outdirPeaks 		  = "${params.outdir}/peakCalls"
params.outdirAnnot 	 	  = "${params.outdir}/annotation"
params.outdirKmetBG 	  = "${params.outdir}/coverage"
params.outdirCoverage 	= "${params.outdir}/coverage"
params.outdirSlices     = "${params.outdir}/slices"
params.outdirFig4Slices = "${params.outdir}/slices/fig4"
params.outdirTables 	  = "${params.outdir}/tables"
params.outdirRTables    = "${params.outdir}/Rtables"
params.outdirImages 	  = "${params.outdir}/figs"
params.outbam        	  = "${params.outdir}/bam"

params.outdir_tmp  	    = "/tmp"

params.timecoursebams        = "${params.timecoursedir}/*KmetStat.bam"
params.timecourseH3K4me3bams = "${params.timecoursedir}/*R1.H3K4me3*.bam*"
params.allHMBAMs             = "${params.allHMdir}/[aH]*_SCP3pos_H1Tneg.bam"
params.allHMBAMsI            = "${params.allHMdir}/*put_SCP3pos_H1Tneg.bam"
params.allHMBAMsCTRL         = "${params.allHMdir}/I*_SCP3pos_H1Tneg.bam"
params.everyHMBAM            = "${params.allHMdir}/*_SCP3pos_H1Tneg.bam"
params.everyHMBAMI           = "${params.allHMdir}/*_SCP3pos_H1Tneg.bam.bai"

params.getRNASeqData    = true

def mm10FA              = "${params.genomedir}/mm10_genome.fa"
def mm10IDX             = "${params.genomedir}/mm10_genome.fa.fai"
def mm10KMetFA          = "${params.genomedir}/mm10_KmetStat_genome.fa"
def mm10KMetIDX         = "${params.genomedir}/mm10_KmetStat_genome.fa.fai"

//output and tempory directories
log.info "===================================================================="
log.info "Lam et al. PIPELINE : "
log.info "===================================================================="
log.info " "
log.info "- pipeline args ----------------------------------------------------"
log.info "Project dir         : ${params.projectdir}"
log.info "threads             : ${params.threads}"
log.info "mem                 : ${params.mem}"
log.info " "
log.info "- annotation data --------------------------------------------------"
log.info "code dir            : ${params.codedir}"
log.info "data dir            : ${params.datadir}"
log.info "genome dir          : ${params.genomedir}"
log.info "timecourse BAM dir  : ${params.timecoursedir}"
log.info "HistMod BAM dir     : ${params.allHMdir}"
log.info "RNA-Seq dir         : ${params.rnadir}"
log.info " "
log.info "- output directories -----------------------------------------------"
log.info "outdir              : ${params.outdir}"
log.info "outdir [Peaks]      : ${params.outdirPeaks}"
log.info "outdir [Annotation] : ${params.outdirAnnot}"
log.info "outdir [Kmet BG]    : ${params.outdirKmetBG}"
log.info "outdir [BG and BW]  : ${params.outdirCoverage}"
log.info "outdir [Images]     : ${params.outdirImages}"
log.info " "
log.info "--------------------------------------------------------------------"

//Create timeCourse input channels
Channel
     .fromPath(params.timecoursebams)
     .ifEmpty { exit 1, "TimeCourse BAM files not found or mis-named" }
     .map { sample -> tuple(getTCType(sample), getTCRep(sample), getTCName(sample), file(getIDX(sample)),  sample) }
     .into {bamTC; bamTCalt}

 def getTCType( file ) {
   def nm="${file.name}"
   def nRet = nm.replaceFirst(/.R.+bam/,"")

   //println("type = $nRet")
   return nRet
 }

 def getTCRep( file ) {
   def nm="${file.name}"
   def nRet = nm.replaceFirst(".{2}.(R[0-9]).+bam","\$1")

   //println("rep = $nRet")
   return nRet
 }

 def getTCName( file ) {
   def nm="${file.name}"
   def nRet = nm.replaceFirst("(.{2}.R[0-9].(H3K9Ac|H3K4me3|abK4me|input)).+bam","\$1")

   println("name = $nRet")
   return nRet
 }

 def getIDX( file ) {
    def nm="${file}"
    return nm.replaceFirst(/.bam/,".bam.bai")
 }

//Create timeCourse H3K4me3 bam input channels
Channel
    .fromPath(params.timecourseH3K4me3bams)
    .ifEmpty { exit 1, "TimeCourse BAM files not found or mis-named" }
    .set {bamH3K4m3}

// >> Done with input channels

//// BUILD ANNOTATION FILES

process getAnnotationFiles {

    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:200 '
    echo true
    cpus 1
    memory "8G"

    module 'macs/2.1.2'
    module 'bedtools/2.27.1'
    module 'R/3.5.2'
    module 'meme/5.0.1'
    module 'ucsc/373'

 	  time { 2.hour }
    errorStrategy { 'retry' }
    maxRetries 1

    tag{sampleID}

    publishDir params.outdirAnnot, mode: 'copy', overwrite: true

    input:

    output:
    file 'B6_maleHS.1bp.bedgraph'       into hotspotBG
    file 'B6_maleHS.500bp.bedgraph'     into hotspotBG500
    file 'B6_maleHS.oneMotif.500bp.bed' into hotspotOneMotif500a, hotspotOneMotif500b, hotspotOneMotif500c, hotspotOneMotif500d, hotspotOneMotif500e
    file 'B6_maleHS.1Kb.bed'            into hotspotBED
    file 'gencodeTSS.1bp.noMerge.bed'   into tssNoMergeBED1bp
    file 'gencodeTSS.1Kb.noMerge.bed'   into tssNoMergeBED1kb_a, tssNoMergeBED1kb_b
    file 'gencodeTSS.1KbDets.bed'       into tssBEDDets
    file 'gencodeTSS.1Kb.bed'           into tssBEDa  , tssBEDb  , tssBEDc  , tssBEDd  , tssBEDe  , tssBEDf  , tssBEDg  , tssBEDh, tssBEDi, tssBEDj, tssBEDk, tssBEDl, tssBEDm, tssBEDn
    file 'gencodeTES.1Kb.bed'           into gencodeTESBED
    file 'gencodeGene.bed'              into gencodeGeneBED
    file 'HS_and_TSS.1KbDets.bed'       into hstssBEDDets
    file 'HS_and_TSS.1Kb.bed'           into hstssBEDa, hstssBEDb, hstssBEDc, hstssBEDd, hstssBEDe, hstssBEDf, hstssBEDg, hstssBEDh, hstssBEDi, hstssBEDj, hstssBEDk
    file 'gencode.vM20.annotation.gtf'  into gencodeGTFa, gencodeGTFb, gencodeGTFc, gencodeGTFd
    file 'gencodeTSS.3Kb.noMerge.bed'   into tss3Kb
    file 'B6_maleHS.3Kb.bed'            into hs3Kb, hs3Kbb, hs3Kbc
    file 'B6_Pr9KO_maleHS.3Kb.bed'      into hsPrKO3Kb, hsPrKO3Kbb
    file 'refseqTSS.bed'                into refseqTSSBED,refseqTSSBEDa,refseqTSSBEDb
    file 'refseqTES.bed'                into refseqTESBED
    file 'refseqGene.bed'               into refseqGeneBED
    file 'refSeqEnhancers.bed'          into refseqEnhancersBED
    file 'repmaskerLINE.bed'            into rmLINEBED

 	  script:
 	  """
 	  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2664nnn/GSM2664275/suppl/GSM2664275_Testis_SSDS_T1.DSBhotspots.bedgraph.gz
 	  gunzip -c GSM2664275_Testis_SSDS_T1.DSBhotspots.bedgraph.gz |cut -f1-3,6 |grep -P \'^chr[0-9]+\' >B6_maleHS.bedgraph

   	bedtools slop -l -0.5 -r -0.5 -pct -i B6_maleHS.bedgraph -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >B6_maleHS.1bp.bedgraph
    cut -f1-3 B6_maleHS.1bp.bedgraph                                                                                           >B6_maleHS.1bp.bed

   	bedtools slop -l 250  -r 250       -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6_maleHS.500bp.bedgraph
   	bedtools slop -l 500  -r 500       -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6_maleHS.1Kb.bedgraph
    bedtools slop -l 500  -r 500       -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX} |cut -f1-3 >B6_maleHS.1Kb.bed
    bedtools slop -l 1500 -r 1500      -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX} |cut -f1-3 >B6_maleHS.3Kb.bed
   	cat B6_maleHS.1Kb.bedgraph |perl -lane \'print join("\\t",@F[0..3],"HS","+")\'                 >B6_maleHS.1Kb.forMerge.bed

    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2664nnn/GSM2664291/suppl/GSM2664291_Testis_SSDS_Prdm9ko.DSBhotspots.bedgraph.gz
    gunzip -c GSM2664291_Testis_SSDS_Prdm9ko.DSBhotspots.bedgraph.gz |cut -f1-3,6 |grep -P \'^chr[0-9]+\' >B6_Pr9KO_maleHS.bedgraph

    bedtools slop -l -0.5 -r -0.5 -pct -i B6_Pr9KO_maleHS.bedgraph     -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >B6_Pr9KO_maleHS.1bp.bedgraph
    cut -f1-3 B6_Pr9KO_maleHS.1bp.bedgraph                                                                                               >B6_Pr9KO_maleHS.1bp.bed

    bedtools slop -l 250  -r 250       -i B6_Pr9KO_maleHS.1bp.bedgraph -g ${mm10IDX}            >B6_Pr9KO_maleHS.500bp.bedgraph
    bedtools slop -l 500  -r 500       -i B6_Pr9KO_maleHS.1bp.bedgraph -g ${mm10IDX}            >B6_Pr9KO_maleHS.1Kb.bedgraph
    bedtools slop -l 500  -r 500       -i B6_Pr9KO_maleHS.1bp.bedgraph -g ${mm10IDX} |cut -f1-3 >B6_Pr9KO_maleHS.1Kb.bed
    bedtools slop -l 1500 -r 1500      -i B6_Pr9KO_maleHS.1bp.bedgraph -g ${mm10IDX} |cut -f1-3 >B6_Pr9KO_maleHS.3Kb.bed

    perl -lane \'\$nm=join("_",@F[0..2]); print join("\\t",@F[0..2],\$nm,\$nm,\$F[3])' B6_maleHS.500bp.bedgraph >B6HS500forFIMO.bed
    bedtools getfasta -fi ${mm10FA} -bed B6HS500forFIMO.bed -name -fo B6_maleHS.500bp.fa

    fimo --max-stored-scores 1000000 --thresh 1e-3 --o fimo ${params.datadir}/PRDM9motif/PRBS_B6.MEMEv4.pwm B6_maleHS.500bp.fa

    perl ${params.codedir}/getHotspotsWithSingleMotif.pl --fimo ./fimo/fimo.tsv --w 250 --out B6_maleHS.oneMotif.500bp.bed

    ##GENCODE
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf.gz
    ##TSS
    perl ${params.codedir}/gencodeGTFtoTSS.pl gencode.vM20.annotation.gtf.gz |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >gencodeTSS.1bp.noMerge.bed

    bedtools slop -l 500 -r 500   -i gencodeTSS.1bp.noMerge.bed -g ${mm10IDX} >gencodeTSS.1Kb.noMerge.bed
    bedtools slop -l 1500 -r 1500 -i gencodeTSS.1bp.noMerge.bed -g ${mm10IDX} >gencodeTSS.3Kb.noMerge.bed
    mergeBed -i gencodeTSS.1Kb.noMerge.bed -c 4,5,6 -o distinct,distinct,first |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX}            >gencodeTSS.1KbDets.bed
    mergeBed -i gencodeTSS.1Kb.noMerge.bed -c 4,5,6 -o distinct,distinct,first |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} |cut -f1-3 >gencodeTSS.1Kb.bed

    cat gencodeTSS.1KbDets.bed |perl -lane \'print join("\\t",@F[0..2],0,@F[4..5])\'                           >gencodeTSS.1Kb.forMerge.bed

    ##TES
    perl ${params.codedir}/gencodeGTFtoTES.pl gencode.vM20.annotation.gtf.gz |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} |grep -P "\\s+" >gencodeTES.1bp.noMerge.bed

    bedtools slop -l 500 -r 500   -i gencodeTES.1bp.noMerge.bed -g ${mm10IDX} >gencodeTES.1Kb.noMerge.bed
    mergeBed -i gencodeTES.1Kb.noMerge.bed -c 4,5,6 -o distinct,distinct,first |\
                                            ${params.codedir}/sortBEDByFAI.pl - \
                                            ${mm10KMetIDX} >gencodeTES.1Kb.bed

    ##GENES
    perl ${params.codedir}/gencodeGTFtoCDS.pl gencode.vM20.annotation.gtf.gz |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >gencodeGene.noMerge.bed
    mergeBed -s -i gencodeGene.noMerge.bed -c 4,5,6 \
             -o distinct,distinct,first |${params.codedir}/sortBEDByFAI.pl - \
            ${mm10KMetIDX} |grep -P "\\s+" >gencodeGene.bed

    ##HS and TSS
    cat B6_maleHS.1Kb.forMerge.bed gencodeTSS.1Kb.forMerge.bed |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX}            >HS_and_TSS.1KbDets.bed
    cat B6_maleHS.1Kb.forMerge.bed gencodeTSS.1Kb.forMerge.bed |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} |cut -f1-3 >HS_and_TSS.1Kb.bed

    gunzip gencode.vM20.annotation.gtf.gz

    ##REFSEQ
    wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ncbiRefSeqCurated.txt.gz
    perl ${params.codedir}/parseRefSeq.pl ncbiRefSeqCurated.txt.gz TSS  |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >refseqTSS.bed
    perl ${params.codedir}/parseRefSeq.pl ncbiRefSeqCurated.txt.gz TES  |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >refseqTES.bed
    perl ${params.codedir}/parseRefSeq.pl ncbiRefSeqCurated.txt.gz Gene |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >refseqGene.bed

    ##ENHANCERS
    bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/mm10/ncbiRefSeq/refSeqFuncElems.bb stdout |cut -f1-3 |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >refSeqEnhancers.bed

    ##LINES
    wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
    zcat rmsk.txt.gz | perl -lane 'print join("\\t",@F[5..7],\$F[10],\$F[11],\$F[9]) if (\$F[1] > 20000 && \$_ =~ /LINE/)' |grep -P "\\s+" |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >repmaskerLINE.bed

 	  """
  }

process getSpo11OligoData {

  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:200 '
  echo true
  cpus 1
  memory "8G"

  module 'bedtools/2.27.1'

  time { 2.hour }
  errorStrategy { 'retry' }
  maxRetries 2

  tag{sampleID}

  publishDir params.outdirAnnot, mode: 'copy', overwrite: true

  input:

  output:
	file 'B6_spo11Oligo.bedgraph'           into spo11BGa, spo11BGb

 	script:
 	"""
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2247nnn/GSM2247727/suppl/GSM2247727_B6_unique_as_unique_RPM.txt.gz
 	gunzip -c GSM2247727_B6_unique_as_unique_RPM.txt.gz |perl -lane \'print join("\\t",\$F[0],\$F[1]-1,\$F[1],\$F[2])\' |grep -P \'^chr[0-9]+\' |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >B6_spo11Oligo.bedgraph
 	"""
  }

if (params.getRNASeqData){
    process getRNASeqData {

    	scratch '/lscratch/$SLURM_JOBID'
    	clusterOptions ' --gres=lscratch:200 '
    	echo true
    	cpus 12
    	memory "16G"

    	module 'kallisto/0.45.0'
    	module 'sratoolkit/2.9.2'

      	time { 2.hour }
    	errorStrategy { 'retry' }
    	maxRetries 2

      	input:
      	file(gtf) from gencodeGTFa

    	output:
    	file('*fastq.gz') 					into fastqRNASeq
    	file('gencode.vM20.transcripts.fa') into transcriptFA
    	file('gencode.vM20.kallistoIDX')    into kallistoIDX

     	script:
     	"""
     	fastq-dump --stdout --gzip SRX3928871 >LE.fastq.gz
     	fastq-dump --stdout --gzip SRX3928868 >ZY.fastq.gz
     	fastq-dump --stdout --gzip SRX3928869 >PA.fastq.gz
     	fastq-dump --stdout --gzip SRX3928866 >DI.fastq.gz

    	perl ${params.codedir}/catGTFgeneNames.pl <${gtf} >catNames.gtf
     	${params.codedir}/gtf_to_fasta catNames.gtf ${mm10FA} gencode.vM20.transcripts.fa
     	perl -pi -e \'s/^\\>\\d+\\s+(\\S+)\\s+chr.+/>\$1/' gencode.vM20.transcripts.fa

    	kallisto index -i gencode.vM20.kallistoIDX --make-unique gencode.vM20.transcripts.fa

 	"""
  }
 }else{
    Channel
         .fromPath("${params.rnadir}/fastq/*gz")
         .ifEmpty { exit 1, "RNA-Seq fastq files not found [${params.rnadir}/fastq/*gz] " }
         .set {fastqRNASeq}

    Channel
         .fromPath("${params.rnadir}/gencode/gencode.vM20.transcripts.fa")
         .ifEmpty { exit 1, "Gencode fasta file not found " }
         .set {transcriptFA}

    Channel
         .fromPath("${params.rnadir}/gencode/gencode.vM20.kallistoIDX")
         .ifEmpty { exit 1, "Kallisto index file not found " }
         .set {kallistoIDX}
 }

process runKallistoForRNASeq {

    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:200 '
    echo true
    cpus 6
    memory "16G"

    module 'kallisto/0.45.0'
    module 'sratoolkit/2.9.2'

  	time { 2.hour }
    errorStrategy { 'retry' }
    maxRetries 1

    publishDir params.outdirAnnot, mode: 'copy', overwrite: true

  	input:
  	//file(gtf)  from gencodeGTFc
  	file(idx)  from kallistoIDX
  	file(fa)   from transcriptFA
  	file(tss)  from tssNoMergeBED1kb_a
    file(fq)   from fastqRNASeq.collect()

    output:
    file('spermatogenesisRNASeq.kallisto.tsv') into kallistoExpr

 	script:
 	"""
    kallisto quant --index=${idx} --output-dir=LE --threads=12 --single --single-overhang --fragment-length=250 \
                   --sd=100 --bootstrap-samples=6 --seed=42 --plaintext LE.fastq.gz

    kallisto quant --index=${idx} --output-dir=ZY --threads=12 --single --single-overhang --fragment-length=250 \
                   --sd=100 --bootstrap-samples=6 --seed=42 --plaintext ZY.fastq.gz

    kallisto quant --index=${idx} --output-dir=PA --threads=12 --single --single-overhang --fragment-length=250 \
                   --sd=100 --bootstrap-samples=6 --seed=42 --plaintext PA.fastq.gz

    kallisto quant --index=${idx} --output-dir=DI --threads=12 --single --single-overhang --fragment-length=250 \
                   --sd=100 --bootstrap-samples=6 --seed=42 --plaintext DI.fastq.gz

    paste  LE/abundance.tsv ZY/abundance.tsv PA/abundance.tsv DI/abundance.tsv |cut -f1,5,10,15,20 |grep -v target_id >allTSV.tab
    perl ${params.codedir}/alignExpressionWithTSSs.pl --t ${tss} --e allTSV.tab |sort -k1,1 -k2n,2n -k3n,3n >spermatogenesisRNASeq.kallisto.tsv

 	"""
  }

//// HERE IS THE PROCESSING SECTION FOR EACH KMETSTAT BAM
process cleanBAMKMet {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:400 --partition=norm'
    echo true
    cpus 12
    memory '24g'

    time { 16.hour }
    errorStrategy { 'retry' }
    maxRetries 2

    module 'samtools/1.9'
    module 'picard/2.17.11'
    module 'bedtools/2.27.1'
    module 'bamtools/2.5.1'

    tag { bamOut }

    input:
    set val(type), val(rep), val(bamOut), file(idx), file(bam) from bamTC

    output:
  	set val(bamOut), file("*cleaned.bam"), file("*cleaned.bam.bai") into cleanBAM
  	set val(bamOut), file("*COclean.bam"), file("*bamtype")         into coSortBAM

    script:
    """
    ### Necessary to allow for K9Ac Rep2 SR data !!!
    if [ `samtools view -hb ${bam} chr1:40000000-45000000 |samtools view -c -f 1 - 2>/dev/null` == '0' ];
    then
       bamtools filter -in ${bam} -isFailedQC false -isMapped true -isPrimaryAlignment true \
                                  -forceCompression -out SROK.bam

       java -jar \$PICARDJAR SortSam I=SROK.bam O=COsorted.bam SO=queryname VALIDATION_STRINGENCY=LENIENT
       java -jar \$PICARDJAR CleanSam I=COsorted.bam O=${bamOut}.COclean.bam VALIDATION_STRINGENCY=LENIENT
       echo "SR" >SR.bamtype
    else
       bamtools filter -in ${bam} -isFailedQC false -isMapped true -isMateMapped true \
                                  -isProperPair true -isPaired true -isPrimaryAlignment true \
                                  -forceCompression -out pairedAndOK.bam

       java -jar \$PICARDJAR SortSam I=pairedAndOK.bam O=COsorted.bam SO=queryname VALIDATION_STRINGENCY=LENIENT
       java -jar \$PICARDJAR FixMateInformation I=COsorted.bam O=fixed_mate.bam ADD_MATE_CIGAR=true VALIDATION_STRINGENCY=LENIENT
       java -jar \$PICARDJAR CleanSam I=fixed_mate.bam O=${bamOut}.COclean.bam VALIDATION_STRINGENCY=LENIENT
       echo "PE" >PE.bamtype
    fi

    java -jar \$PICARDJAR SortSam I=${bamOut}.COclean.bam O=${bamOut}.initial.cleaned.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT
    samtools index ${bamOut}.initial.cleaned.bam
    """
    }

process makeWinFilesKMet {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:10'
    echo true
    cpus 2
    memory '4g'

    time { 1.hour }
    errorStrategy { 'retry' }
    maxRetries 2

    module 'bedtools/2.27.1'

    input:

    output:
    file "win*bed" into winFilesa,winFilesb,winFilesc,winFilesd,winFilese,winFilesf

    script:
    """
    bedtools makewindows -g ${mm10KMetIDX} -w 25          | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 25   && \$_ !~ /Kmet/)\' >win25.bed
    bedtools makewindows -g ${mm10KMetIDX} -w 147         | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 147  && \$_ !~ /Kmet/)\' >win147.bed
    bedtools makewindows -g ${mm10KMetIDX} -w 200         | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 200  && \$_ !~ /Kmet/)\' >win200.bed
    bedtools makewindows -g ${mm10KMetIDX} -w 500         | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500  && \$_ !~ /Kmet/)\' >win500.bed
    bedtools makewindows -g ${mm10KMetIDX} -w 1000        | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 1000 && \$_ !~ /Kmet/)\' >win1000.bed
    bedtools makewindows -g ${mm10KMetIDX} -w 1000 -s 147 | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 1000 && \$_ !~ /Kmet/)\' >win1ks147.bed
    """
    }

process makeFragmentsBEDKMet {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:200'
    echo true
    cpus 2
    memory '12g'

    time { 3.hour }
    errorStrategy { 'retry' }
    maxRetries 2

    module 'bedtools/2.27.1'

  tag { bamOut }

    //publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    set val(bamOut), file(CObam), file(bamType) from coSortBAM

    output:
    set val(bamOut), file("*cleaned.bed"), file("*bedtype")     into cleanBEDa,cleanBEDb

    script:
    """
       if [ -e SR.bamtype ]; then
          bedtools bamtobed        -i ${CObam} >cleaned.tmp.bed
          cp SR.bamtype SR.bedtype
       else
          bedtools bamtobed -bedpe -i ${CObam} | perl -lane \'print join("\\t",\$F[0],\$F[1],\$F[5],@F[6..8]) if (\$F[0] eq \$F[3])\' >cleaned.tmp.bed
          cp PE.bamtype PE.bedtype
       fi
       perl ${params.codedir}/sortBEDByFAI.pl cleaned.tmp.bed ${mm10KMetIDX} | perl -lane \'print \$_ if (\$F[4] >= 20)\' >${bamOut}.initFragments.cleaned.bed
    """
    }

process monoNucBGKMet {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:400 --partition=norm'
    echo true
    cpus 2
    memory '12g'

    time { 8.hour }
    errorStrategy { 'retry' }
    maxRetries 2

    module 'bedtools/2.27.1'

    tag { bamOut }

    publishDir params.outdirKmetBG, mode: 'copy', overwrite: true

    input:
    set val(bamOut), file(bed), file(srpeType) from cleanBEDa
    file winFiles from winFilesa.collect()

    output:
    set val(bamOut), file("*monoNucleosomes*bedgraph") into monoBG
    set val(bamOut), file("*monoNucleosomes.bed")      into monoBED
    file("*coverStats.tab")                            optional true into coverStats
    file("*monoNucleosomes.ws25bp.bedgraph")           optional true into monoBG25all
    file("*H3K4*monoNucleosomes.ws25bp.bedgraph")      optional true into bgH3K4
    file("*abK4*monoNucleosomes.ws25bp.bedgraph")      optional true into bgabK4
    file("*H3K9*monoNucleosomes.ws25bp.bedgraph")      optional true into bgH3K9
    file("*input*monoNucleosomes.ws25bp.bedgraph")     optional true into bgInput
    file("*K*monoNucleosomes.bed")                     optional true into monoBEDallChIP
    file("*H3K4*monoNucleosomes.bed")                  optional true into monoBEDH3K4
    file("*abK4*monoNucleosomes.bed")                  optional true into monoBEDabK4
    file("*H3K9*monoNucleosomes.bed")                  optional true into monoBEDH3K9
    file("*input*monoNucleosomes.bed")                 optional true into monoBEDInput

    script:
    """
    if [ -e SR.bedtype ]; then
      cp ${bed} ${bamOut}.monoNucleosomes.bed
    else
      perl -lane \'\$sz = (\$F[2]-\$F[1]); print \$_ if (\$sz > 100 && \$sz < 200)\' ${bed} >${bamOut}.monoNucleosomes.bed
    fi

    intersectBed -a win25.bed     -b ${bamOut}.monoNucleosomes.bed -c -sorted |bedtools slop -l -1 -r -0     -g ${mm10KMetIDX} -i - >${bamOut}.monoNucleosomes.ws25bp.bedgraph
    intersectBed -a win147.bed    -b ${bamOut}.monoNucleosomes.bed -c -sorted |bedtools slop -l -1 -r -0     -g ${mm10KMetIDX} -i - >${bamOut}.monoNucleosomes.ws147bp.bedgraph
    intersectBed -a win200.bed    -b ${bamOut}.monoNucleosomes.bed -c -sorted |bedtools slop -l -1 -r -0     -g ${mm10KMetIDX} -i - >${bamOut}.monoNucleosomes.ws200bp.bedgraph
    intersectBed -a win500.bed    -b ${bamOut}.monoNucleosomes.bed -c -sorted |bedtools slop -l -1 -r -0     -g ${mm10KMetIDX} -i - >${bamOut}.monoNucleosomes.ws500bp.bedgraph
    intersectBed -a win1000.bed   -b ${bamOut}.monoNucleosomes.bed -c -sorted |bedtools slop -l -1 -r -0     -g ${mm10KMetIDX} -i - >${bamOut}.monoNucleosomes.ws1000bp.bedgraph
    intersectBed -a win1ks147.bed -b ${bamOut}.monoNucleosomes.bed -c -sorted |bedtools slop -l -427 -r -427 -g ${mm10KMetIDX} -i - >${bamOut}.monoNucleosomes.w1000s147bp.bedgraph

    cut -f1 ${bamOut}.monoNucleosomes.bed |uniq -c |perl -lane \'\$c += \$F[0]; print join("\\t",\$F[1],\$F[0],\$c)\' >${bamOut}.monoNucleosomes.coverStats.tab
    """
    }

process makeBAMsKMet {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:400 --partition=norm'
    echo true
    cpus 2
    memory '12g'

    time { 8.hour }
    errorStrategy { 'retry' }
    maxRetries 2

    module 'bedtools/2.27.1'
    module 'samtools/1.9'
    module 'picard/2.17.11'

    tag { bamOut }

    publishDir params.outdirKmetBG, mode: 'copy', overwrite: true

    input:
    set val(bamOut), file(bam), file(bai) from cleanBAM

    output:
    set val(bamOut), file("*ucleo*bam"), file("*ucleo*bam.bai"), file("*bamStats.tab")     into finalbam
  	file("*.monoNucleosomes.bam") 									   into justMonoBAM
  	file("*K*monoNucleosomes.bam")                     optional true into (monoBAMChIP, monoBAMChIP_b)
    file("*H3K4*monoNucleosomes.bam")                  optional true into monoBAMH3K4
    file("*abK4*monoNucleosomes.bam")                  optional true into monoBAMabK4
    file("*H3K9*monoNucleosomes.bam")                  optional true into monoBAMH3K9
    file("*input*monoNucleosomes.bam")                 optional true into monoBAMInput

    script:
    """
    if [ `samtools view -hb ${bam} chr1:40000000-45000000 |samtools view -c -f 1 - 2>/dev/null` == '0' ];
    then
       cp ${bam} ${bamOut}.monoNucleosomes.bam
    else
       samtools view -h ${bam} |perl -lane \'print \$_ if (\$_ =~ /^@/ || (abs(\$F[8]) >= 100 && abs(\$F[8]) <= 200))\' |samtools view -Shb - >${bamOut}.monoNucleosomes.bam
    fi

    samtools index ${bamOut}.monoNucleosomes.bam
    samtools idxstats ${bamOut}.monoNucleosomes.bam >${bamOut}.monoNucleosomes.bamStats.tab
    """
    }

process correctKMetBedgraphs {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:40'
    echo true
    cpus 1

    //memory '16g'
    memory { 12.GB * task.attempt }

    errorStrategy { 'retry' }
    maxRetries 2

    time { 3.hour }

    tag  { bgT }

    module 'ucsc/373'
    module 'bedtools/2.27.1'

    publishDir params.outdirCoverage, mode: 'copy', overwrite: true

    input:
    file(bgT)     from bgH3K4
    file(allBGI)  from bgInput.collect()
    file(allStat) from coverStats.collect()

    output:
    file '*.monoCorrected.ws25bp.bedgraph' into correctedH3K4BGa, correctedH3K4BGb
    file '*.monoCorrected.ws25bp.bigwig'   into correctedH3K4BWa, correctedH3K4BWb

  	script:
    def chipStat  = bgT.name.replaceFirst(".monoNucleosomes.ws25bp.bedgraph"       ,".monoNucleosomes.coverStats.tab")
    def outBG     = bgT.name.replaceFirst(".monoNucleosomes.ws25bp.bedgraph"       ,".monoCorrected.ws25bp.bedgraph")
    def outBW     = bgT.name.replaceFirst(".monoNucleosomes.ws25bp.bedgraph"       ,".monoCorrected.ws25bp.bigwig")
  	def inputBG   = bgT.name.replaceFirst("H3K4me3.monoNucleosomes.ws25bp.bedgraph","input.monoNucleosomes.ws25bp.bedgraph")
  	def inputStat = bgT.name.replaceFirst("H3K4me3.monoNucleosomes.ws25bp.bedgraph","input.monoNucleosomes.coverStats.tab")
  	"""

    perl ${params.codedir}/correctBGwithKmetStatPanel.pl --t ${bgT} --c ${inputBG} --sT ${chipStat} --sC ${inputStat} --out initCorrected.bedgraph --fai ${mm10KMetIDX}

    ${params.codedir}/sortBEDByFAI.pl initCorrected.bedgraph ${mm10KMetIDX} >${outBG}
    sort -k1,1 -k2n,2n -k3n,3n ${outBG} >forbigwig.bg

    bedGraphToBigWig forbigwig.bg ${mm10KMetIDX} ${outBW}
    """
  }

process callPeaksFromTimeCourseData {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:200 '
    echo true
    cpus 1
    memory "16G"

    module 'macs/2.1.2'
    module 'bedtools/2.27.1'
    module 'samtools/1.9'
    module 'R/3.5.2'

  	time { 2.hour }
    errorStrategy { 'retry' }
    maxRetries 2

    tag{chipBAM}

    publishDir params.outdirPeaks, mode: 'copy', overwrite: true

    input:
    file (chipBAM)   from monoBAMChIP
    file (inputBAMs) from monoBAMInput.collect()

    output:
    	file '*peaks.count'  into peakCountsTC
    	file '*peaks.xls'    into peaksXLSTC
    	file '*peaks.bed'    into peaksBEDTCa, peaksBEDTCb, peaksBEDTCc
    	file '*.model.pdf'   into peakModelsPDFTC
    	file '*model.r'      into peakModelsTC

    script:
    	def outName  = chipBAM.name.replaceFirst(".monoNucleosomes.bam","")
    	def outStage = chipBAM.name.replaceFirst("^(.{2}).R\\d..+.monoNucleosomes.bam","\$1")
    	def outRep   = chipBAM.name.replaceFirst("^.{2}.(R\\d)..+.monoNucleosomes.bam","\$1")
    	def outType  = chipBAM.name.replaceFirst("^.{2}.R\\d.(.+).monoNucleosomes.bam","\$1")
    	//def inputBAM = chipBAM.name.replaceFirst("(\\w\\.\\w)\\.\\w\\.monoNucleosomes.bam","\$1.input.monoNucleosomes.bam")
      def inputBAM = "${outStage}.${outRep}.input.monoNucleosomes.bam"
    	def outBP    = "${outName}_peaks.broadPeak"
    	def outPeaks = "${outName}.peaks.bed"
    	def outCount = "${outName}.peaks.count"
    	def outModel = "${outName}_model.r"
    	def outPDF   = "${outName}.model.pdf"
    	"""
    	macs2 callpeak -q 0.1 -n ${outName} -g mm --broad -t ${chipBAM} -c ${inputBAM}

    	cut -f1-3 ${outBP} |grep -P \'chr\\d+\' |perl -lane \'print join("\\t",@F,"\'${outStage}\'","\'${outRep}\'","\'${outType}\'")\' >tmpPeaks.bed

      intersectBed -a tmpPeaks.bed -b ${params.blacklist} -v >${outPeaks}

      npeaks=`cat ${outPeaks} |wc -l`
    	echo -e "$outName\t${outStage}\t${outRep}\t${outType}\t\$npeaks" >${outCount}

      R --silent --vanilla <${outModel} >${outPDF}
    	"""
  }

process mergePeaksFromTimeCourseData {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:100 '
    echo true
    cpus 1
    memory "4G"

    module 'bedtools/2.27.1'

    time { 0.25.hour }
    errorStrategy { 'retry' }
    maxRetries 2

    publishDir params.outdirPeaks, mode: 'copy', overwrite: true

    input:
    file peaksK4me3  from peaksBEDTCa.collect()

    output:
    file '*K4me3*peaks.bedgraph' into k4me3PeaksBG
    file '*K4me3*peaks.bed'      into k4me3Peaksa, k4me3Peaksb, k4me3Peaksc, k4me3Peaksd, k4me3Peakse, k4me3Peaksf, k4me3Peaksg, k4me3Peaksh, k4me3Peaksi, k4me3Peaksk, k4me3Peaksl
    file '*abK4*peaks.bedgraph'  into abk4PeaksBG
    file '*abK4*peaks.bed'       into abk4Peaksa, abk4Peaksb, abk4Peaksc, abk4Peaksd, abk4Peakse, abk4Peaksf, abk4Peaksg, abk4Peaksh, abk4Peaksi, abk4Peaksj, abk4Peaksk, abk4Peaksl
    file '*K9Ac*peaks.bedgraph'  into k9acPeaksBG
    file '*K9Ac*peaks.bed'       into k9acPeaksa, k9acPeaksb, k9acPeaksc, k9acPeaksd, k9acPeakse, k9acPeaksf, k9acPeaksg, k9acPeaksh, k9acPeaksi, k9acPeaksj, k9acPeaksk, k9acPeaksl
    script:
    """
    sort -k1,1 -k2n,2n -k3n,3n *H3K4me3.peaks.bed                           >init.sorted.bed
    bedtools slop -l -0.5 -r -0.5 -pct -i init.sorted.bed -g ${mm10KMetIDX} >init.1bp.bed
    bedtools slop -l 500 -r 500 -i init.1bp.bed           -g ${mm10KMetIDX} >init.1k.tmp
    sort -k1,1 -k2n,2n -k3n,3n init.1k.tmp                                  >init.1k.bed
    mergeBed -i init.1k.bed -d -750 -c 3 -o count                           >init.merge.bed
    ${params.codedir}/sortBEDByFAI.pl init.merge.bed ${mm10KMetIDX}         >timecourseH3K4me3peaks.bedgraph
    cut -f 1-3 timecourseH3K4me3peaks.bedgraph                              >timecourseH3K4me3peaks.bed

    sort -k1,1 -k2n,2n -k3n,3n *H3K9Ac.peaks.bed                            >init.sorted.bed
    bedtools slop -l -0.5 -r -0.5 -pct -i init.sorted.bed -g ${mm10KMetIDX} >init.1bp.bed
    bedtools slop -l 500 -r 500 -i init.1bp.bed           -g ${mm10KMetIDX} >init.1k.tmp
    sort -k1,1 -k2n,2n -k3n,3n init.1k.tmp                                  >init.1k.bed
    mergeBed -i init.1k.bed -d -750 -c 3 -o count                           >init.merge.bed
    ${params.codedir}/sortBEDByFAI.pl init.merge.bed ${mm10KMetIDX}         >timecourseH3K9Acpeaks.bedgraph
    cut -f 1-3 timecourseH3K9Acpeaks.bedgraph                               >timecourseH3K9Acpeaks.bed

    sort -k1,1 -k2n,2n -k3n,3n *abK4*.peaks.bed                             >init.sorted.bed
    bedtools slop -l -0.5 -r -0.5 -pct -i init.sorted.bed -g ${mm10KMetIDX} >init.1bp.bed
    bedtools slop -l 500 -r 500 -i init.1bp.bed           -g ${mm10KMetIDX} >init.1k.tmp
    sort -k1,1 -k2n,2n -k3n,3n init.1k.tmp                                  >init.1k.bed
    mergeBed -i init.1k.bed -d -750 -c 3 -o count                           >init.merge.bed
    ${params.codedir}/sortBEDByFAI.pl init.merge.bed ${mm10KMetIDX}         >timecourseabK4mepeaks.bedgraph
    cut -f 1-3 timecourseabK4mepeaks.bedgraph                               >timecourseabK4mepeaks.bed
    """
  }

process makeTable1Data {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:100 '
  echo true
  cpus 1
  memory "12G"

  module 'bedtools/2.27.1'
  module 'samtools/1.9'

  time { 1.hour * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 2

  tag { bam }
  //publishDir params.outdirPeaks, mode: 'copy', overwrite: true

  input:
  file (bam)        from monoBAMChIP_b
  file (hs)         from hs3Kbc
  file (peaksK4me3) from peaksBEDTCc.collect()
  file (k4me3Peaks) from k4me3Peaksl.collect()
  file (k9acPeaks)  from k9acPeaksl.collect()
  file (abK4Peaks)  from abk4Peaksl.collect()

  output:
  file '*table1Data.tab' optional true into table1Data

  script:
  def stage = bam.name.replaceFirst("^(.{2}).R\\d..+.monoNucleosomes.bam","\$1")
  def rep   = bam.name.replaceFirst("^.{2}.(R\\d)..+.monoNucleosomes.bam","\$1")
  def type  = bam.name.replaceFirst("^.{2}.R\\d.(.+).monoNucleosomes.bam","\$1")
  def allPeaks = "timecourse${type}peaks.bed"
  def myPeaks  = "${stage}.${rep}.${type}.peaks.bed"
  """
  if [[ ${bam} =~ H3K9Ac ]]; then
    echo "Skipping ..."
  else
    ln -s `readlink -f ${bam}`.bai ${bam}.bai
    perl ${params.codedir}/sortBEDByFAI.pl ${myPeaks}  ${mm10IDX} >pkMe.bed
    perl ${params.codedir}/sortBEDByFAI.pl ${allPeaks} ${mm10IDX} >pkAll.bed
    perl ${params.codedir}/sortBEDByFAI.pl ${hs} ${mm10IDX} >hs.bed

    myHScount=`cat pkMe.bed |wc -l`
    allHScount=`cat pkAll.bed |wc -l`
    dsbHScount=`cat hs.bed |wc -l`

    nTot=`samtools idxstats ${bam} |cut -f3 | paste -sd+ | bc`

    nDSBHS=`intersectBed -a ${bam} -b hs.bed -sorted -bed |wc -l`
    dsbSPOT=`echo -e "\$nDSBHS\\t\$nTot" |perl -lane 'print sprintf("%4.2f",(\$F[0]/\$F[1]*100))'`

    nMeHS=`intersectBed -a ${bam} -b pkMe.bed -sorted -bed |wc -l`
    meSPOT=`echo -e "\$nMeHS\\t\$nTot" |perl -lane 'print sprintf("%4.2f",(\$F[0]/\$F[1]*100))'`

    nAllHS=`intersectBed -a ${bam} -b pkAll.bed -sorted -bed |wc -l`
    allSPOT=`echo -e "\$nAllHS\\t\$nTot" |perl -lane 'print sprintf("%4.2f",(\$F[0]/\$F[1]*100))'`

    grep -P '${stage}\\s+${rep}\\s+${type}' ${params.table1Init} \
       |perl -lane 'print join("\\t",@F,"'\$nTot'","'\$nDSBHS'","'\$nMeHS'","'\$nAllHS'","'\$dsbSPOT'","'\$meSPOT'","'\$allSPOT'","'\$myHScount'","'\$allHScount'")' \
       >${stage}.${rep}.${type}.table1Data.tab
  fi
  """
  }

process makeTable1 {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:100 '
  echo true
  cpus 1
  memory "2G"

  module 'bedtools/2.27.1'

  time { 0.25.hour }
  errorStrategy { 'retry' }
  maxRetries 2

  publishDir params.outdirTables, mode: 'copy', overwrite: true

  input:
  file (tab1Data)  from table1Data.collect()

  output:
  file 'table1.tab' into table1Out

  script:

  """
  echo -e "stage\\trep\\ttype\\tab\\tpurity\\tDNA\\ttotReads\\tHotspotreads\\tmyPkreads\\tallPkreads\\tSPOThs\\tSPOTme\\tSPOTall\\tmyPeaks\\tallPeaks" >table1.tab
  cat *table1Data.tab >>table1.tab
  """
  }

process overlapKM {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '8g'

  time { 1.hour }

  module 'bedtools/2.27.1'

  tag {kBG}
  //publishDir params.outdirTables, mode: 'copy', overwrite: true

  input:
  file(tss)   from tssBEDb
  file(hstss) from hstssBEDb
  file(me3pk) from k4me3Peaksa
  file(abpk)  from abk4Peaksa
  file(k9pk)  from k9acPeaksa
  file(kBG)   from correctedH3K4BGa

  output:
  file('*.TSS.ol')    into tssOLKM
  file('*.HSTSS.ol')  into hstssOLKM
  file('*.K4ME3.ol')  into k4me3OLKM
  file('*.ABK4ME.ol') into abk4OLKM
  file('*.K9AC.ol')   into k9acOLKM

	script:
	def bgName = kBG.name.replaceFirst(".monoCorrected.ws25bp.bedgraph",".KM")
	"""
	echo ${bgName} >header.txt

  mapBed -a ${tss}   -b ${kBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${bgName}.TSS.tmp
  mapBed -a ${hstss} -b ${kBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${bgName}.HSTSS.tmp
  mapBed -a ${me3pk} -b ${kBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${bgName}.K4ME3.tmp
  mapBed -a ${abpk}  -b ${kBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${bgName}.ABK4ME.tmp
  mapBed -a ${k9pk}  -b ${kBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${bgName}.K9AC.tmp

  cat header.txt ${bgName}.TSS.tmp    >${bgName}.TSS.ol
  cat header.txt ${bgName}.HSTSS.tmp  >${bgName}.HSTSS.ol
  cat header.txt ${bgName}.K4ME3.tmp  >${bgName}.K4ME3.ol
  cat header.txt ${bgName}.ABK4ME.tmp >${bgName}.ABK4ME.ol
  cat header.txt ${bgName}.K9AC.tmp   >${bgName}.K9AC.ol
  """
  }

//// NCIS AND STABLE PROMOTER CORRECTIONS
process runNCIScorrection {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:40'
    echo true
    cpus 1
    memory '16g'

    time { 3.hour }

    tag  { chipBED }

    module 'R'
    module 'ucsc/373'

    publishDir params.outdirCoverage, mode: 'copy', overwrite: true

    input:
    file(chipBED) from monoBEDallChIP
    file(allBG)   from monoBG25all.collect()
    file(allIn)   from monoBEDInput.collect()

    output:
    file '*.ncisCorrected.ws25bp.bedgraph' into ncisBGa, ncisBGb, ncisBGc
    file '*.ncisCorrected.ws25bp.bigwig'   into ncisBWa, ncisBWb, ncisBWc

  	script:

    def outStage = chipBED.name.replaceFirst("^(.{2}).R\\d..+.monoNucleosomes.bed","\$1")
    def outRep   = chipBED.name.replaceFirst("^.{2}.(R\\d)..+.monoNucleosomes.bed","\$1")
    def outType  = chipBED.name.replaceFirst("^.{2}.R\\d.(.+).monoNucleosomes.bed","\$1")
    def inputBED = "${outStage}.${outRep}.input.monoNucleosomes.bed"
    def inputBG  = "${outStage}.${outRep}.input.monoNucleosomes.ws25bp.bedgraph"

    def chipBG    = chipBED.name.replaceFirst(".monoNucleosomes.bed"       ,".monoNucleosomes.ws25bp.bedgraph")
    def outBG     = chipBED.name.replaceFirst(".monoNucleosomes.bed"       ,".ncisCorrected.ws25bp.bedgraph")
    def outBW     = chipBED.name.replaceFirst(".monoNucleosomes.bed"       ,".ncisCorrected.ws25bp.bigwig")
    def outNCIS   = chipBED.name.replaceFirst(".monoNucleosomes.bed"       ,".ncis.tab")

    //def inputBG   = chipBED.name.replaceFirst("\\w\\.\\w\\..+K.+.monoNucleosomes.bed" ,".input.monoNucleosomes.ws25bp.bedgraph")
    //def inputBED  = chipBED.name.replaceFirst("\\w\\.\\w\\..+K.+.monoNucleosomes.bed" ,".input.monoNucleosomes.bed")
    //def outName  = chipBAM.name.replaceFirst(".monoNucleosomes.bam","")

  	"""
  	grep -w chr1 ${chipBED} >chip.bed
  	grep -w chr1 ${inputBED} >input.bed

  	perl ${params.codedir}/runNCIS.pl           --t chip.bed --c input.bed --lib ${params.codedir}/NCIS --out ${outNCIS}
  	ncis=`cut -f1 ${outNCIS}`

    perl ${params.codedir}/correctBGwithNCIS.pl --t ${chipBG} --c ${inputBG} --ncis \$ncis --out ${outBG}

    grep -P \'^chr[\\d]+\\s\' ${outBG}  |sort -k1,1 -k2n,2n -k3n,3n -T "." >bwSort.bedgraph

    bedGraphToBigWig bwSort.bedgraph ${mm10KMetIDX} ${outBW}
    """
  }

process overlapNCISraw {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '8g'

  time { 1.hour }

  module 'bedtools/2.27.1'
  tag {kBG}
  //publishDir params.outdirTables, mode: 'copy', overwrite: true

  input:
  file(tss)  from tssBEDa
  file(m3pk) from k4me3Peaksb
  file(abpk) from abk4Peaksb
  file(k9pk) from k9acPeaksb
  file(kBG)  from ncisBGa

  output:
  file('*.TSS.olNCISraw')   into tssOLNCISRAW
  file('*.K4ME3.olNCISraw') into k4me3OLNCISRAW
  file('*.ABK4ME.olNCISraw') into abk4OLNCISRAW
  file('*.K9AC.olNCISraw') into k9acOLNCISRAW

  script:
  def bgName = kBG.name.replaceFirst(".ncisCorrected.ws25bp.bedgraph",".NCIS")
  """
  echo ${bgName} >header.txt

  mapBed -a ${tss}  -b ${kBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${bgName}.TSS.tmp
  mapBed -a ${m3pk} -b ${kBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${bgName}.K4ME3.tmp
  mapBed -a ${abpk} -b ${kBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${bgName}.ABK4ME.tmp
  mapBed -a ${k9pk} -b ${kBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${bgName}.K9AC.tmp

  cat header.txt ${bgName}.TSS.tmp    >${bgName}.TSS.olNCISraw
  cat header.txt ${bgName}.K4ME3.tmp  >${bgName}.K4ME3.olNCISraw
  cat header.txt ${bgName}.ABK4ME.tmp >${bgName}.ABK4ME.olNCISraw
  cat header.txt ${bgName}.K9AC.tmp   >${bgName}.K9AC.olNCISraw
  """
  }

process getStablePromoters {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '8g'

  time { 2.hour }

  module 'R/3.5.2'
  module 'bedtools/2.27.1'

  publishDir params.outdirTables, mode: 'copy', overwrite: true, pattern: "*.bed"

  input:
  file(tss)   from tssBEDc
  file(hstss) from hstssBEDh
  file(m3pk)  from k4me3Peaksc
  file(abpk)  from abk4Peaksc
  file(k9pk)  from k9acPeaksc
  file(nOL)   from tssOLNCISRAW.collect()

  output:
  file("stableTSSs.abK4me.bed")  into stableTSSsabK4me
  file("stableTSSs.H3K4me3.bed") into stableTSSsH3K4me3
  file("stableTSSs.H3K9Ac.bed")  into stableTSSsH3K9Ac
  file("tss*stable.ol")          into tssOLstable
  file("hstss*stable.ol")        into hstssOLstable
  file("k4me3*stable.ol")        into k4me3OLstable
  file("abk4*stable.ol")         into abk4OLstable
  file("k9ac*stable.ol")         into k9acOLstable

  script:
  //paste ${tss} *H3K4me3.NCIS* |grep -v NCIS |grep -P \'^chr[\\d]+\\s\' >h3k4me3.txt
  //paste ${tss} *H3K9Ac.NCIS*  |grep -v NCIS |grep -P \'^chr[\\d]+\\s\' >h3k9ac.txt
  //paste ${tss} *abK4me.NCIS*  |grep -v NCIS |grep -P \'^chr[\\d]+\\s\' >abk4me.txt
  """
  echo -e "cs\tfrom\tend" >tss.bed
  cat ${tss} >>tss.bed
  paste tss.bed *H3K4me3.NCIS* |grep -v NCIS |intersectBed -a - -b ${m3pk} |uniq -f 3 >h3k4me3.txt
  paste tss.bed *H3K9Ac.NCIS*  |grep -v NCIS |intersectBed -a - -b ${k9pk} |uniq -f 3 >h3k9ac.txt
  paste tss.bed *abK4me.NCIS*  |grep -v NCIS |intersectBed -a - -b ${abpk} |uniq -f 3 >abk4me.txt

  R -q --no-save --no-restore --slave -f ${params.codedir}/genStableTSS.R --args 1.8 h3k4me3.txt stableTSSs.H3K4me3.bed
  R -q --no-save --no-restore --slave -f ${params.codedir}/genStableTSS.R --args 1.8 h3k9ac.txt  stableTSSs.H3K9Ac.bed
  R -q --no-save --no-restore --slave -f ${params.codedir}/genStableTSS.R --args 1.8 abk4me.txt  stableTSSs.abK4me.bed

  echo "stableH3K4me3" >tss.H3K4me3.stable.ol
  echo "stableH3K4me3" >hstss.H3K4me3.stable.ol
  echo "stableH3K4me3" >k4me3pk.H3K4me3.stable.ol
  echo "stableH3K4me3" >abk4pk.H3K4me3.stable.ol
  echo "stableH3K4me3" >k9acpk.H3K4me3.stable.ol

  echo "stableH3K9Ac" >tss.H3K9Ac.stable.ol
  echo "stableH3K9Ac" >hstss.H3K9Ac.stable.ol
  echo "stableH3K9Ac" >k4me3pk.H3K9Ac.stable.ol
  echo "stableH3K9Ac" >abk4pk.H3K9Ac.stable.ol
  echo "stableH3K9Ac" >k9acpk.H3K9Ac.stable.ol

  echo "stableabK4me" >tss.abK4me.stable.ol
  echo "stableabK4me" >hstss.abK4me.stable.ol
  echo "stableabK4me" >k4me3pk.abK4me.stable.ol
  echo "stableabK4me" >abk4pk.abK4me.stable.ol
  echo "stableabK4me" >k9acpk.abK4me.stable.ol

  mapBed -a ${tss}   -b stableTSSs.H3K4me3.bed -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>tss.H3K4me3.stable.ol
  mapBed -a ${hstss} -b stableTSSs.H3K4me3.bed -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>hstss.H3K4me3.stable.ol
  mapBed -a ${m3pk}  -b stableTSSs.H3K4me3.bed -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>k4me3pk.H3K4me3.stable.ol
  mapBed -a ${abpk}  -b stableTSSs.H3K4me3.bed -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>abk4pk.H3K4me3.stable.ol
  mapBed -a ${k9pk}  -b stableTSSs.H3K4me3.bed -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>k9acpk.H3K4me3.stable.ol

  mapBed -a ${tss}   -b stableTSSs.H3K9Ac.bed  -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>tss.H3K9Ac.stable.ol
  mapBed -a ${hstss} -b stableTSSs.H3K9Ac.bed  -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>hstss.H3K9Ac.stable.ol
  mapBed -a ${m3pk}  -b stableTSSs.H3K9Ac.bed  -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>k4me3pk.H3K9Ac.stable.ol
  mapBed -a ${abpk}  -b stableTSSs.H3K9Ac.bed  -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>abk4pk.H3K9Ac.stable.ol
  mapBed -a ${k9pk}  -b stableTSSs.H3K9Ac.bed  -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>k9acpk.H3K9Ac.stable.ol

  mapBed -a ${tss}   -b stableTSSs.abK4me.bed  -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>tss.abK4me.stable.ol
  mapBed -a ${hstss} -b stableTSSs.abK4me.bed  -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>hstss.abK4me.stable.ol
  mapBed -a ${m3pk}  -b stableTSSs.abK4me.bed  -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>k4me3pk.abK4me.stable.ol
  mapBed -a ${abpk}  -b stableTSSs.abK4me.bed  -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>abk4pk.abK4me.stable.ol
  mapBed -a ${k9pk}  -b stableTSSs.abK4me.bed  -c 3 -o count -g ${mm10KMetIDX} |cut -f4 |perl -lane \'print ((\$F[0]>0)?"TRUE":"FALSE")\' >>k9acpk.abK4me.stable.ol
  """
  }

process correctNCISwithStablePromoters {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1

  memory { 8.GB * task.attempt }
  time   { 2.hour * task.attempt }

  errorStrategy { 'retry' }
  maxRetries 3

  tag  { bg }

  module 'bedtools/2.27.1'
  module 'ucsc/373'

  publishDir params.outdirCoverage, mode: 'copy', overwrite: true, pattern: "*.ncisCorrectedByStablePromoters.ws25bp.b*"

  input:
  file(stableH3K4me3) from stableTSSsH3K4me3
  file(stableH3K9Ac)  from stableTSSsH3K9Ac
  file(stableabK4me)  from stableTSSsabK4me
  file(tss)           from tssBEDd
  file(hstss)         from hstssBEDd
  file(m3pk)          from k4me3Peaksd
  file(abpk)          from abk4Peaksd
  file(k9pk)          from k9acPeaksd
  file(bg)            from ncisBGb

  output:
  file("*.ncisCorrectedByStablePromoters.ws25bp.bedgraph")           optional true into (stableTSSscorrectedBGa, stableTSSscorrectedBGb, stableTSSscorrectedBGc)
  file("*.ncisCorrectedByStablePromoters.ws25bp.bigwig")             optional true into (stableTSSscorrectedBWa, stableTSSscorrectedBWb, stableTSSscorrectedBWc)
  file("LE*.H3K4me3.ncisCorrectedByStablePromoters.ws25bp.bedgraph") optional true into (stableTSSscorrectedLEBG_K4Me3)
  file("LE*.H3K9Ac.ncisCorrectedByStablePromoters.ws25bp.bedgraph")  optional true into (stableTSSscorrectedLEBG_K9Ac)
  file("LE*.abK4me.ncisCorrectedByStablePromoters.ws25bp.bedgraph")  optional true into (stableTSSscorrectedLEBG_abK4)
  file('*.NCISTSS.TSS.ol')                                                         into tssOLNCISTSS
  file('*.NCISTSS.HSTSS.ol')                                                       into hstssOLNCISTSS
  file('*.NCISTSS.K4ME3.ol')                                                       into k4me3OLNCISTSS
  file('*.NCISTSS.ABK4ME.ol')                                                      into abk4OLNCISTSS
  file('*.NCISTSS.K9AC.ol')                                                        into k9acOLNCISTSS

  script:
  def outBG  = bg.name.replaceFirst(".ncisCorrected.ws25bp.bedgraph" ,".ncisCorrectedByStablePromoters.ws25bp.bedgraph")
  def outBW  = bg.name.replaceFirst(".ncisCorrected.ws25bp.bedgraph" ,".ncisCorrectedByStablePromoters.ws25bp.bigwig")
  def olName = bg.name.replaceFirst(".ncisCorrected.ws25bp.bedgraph",".NCISTSS")

  """
  if [[ ${bg} =~ H3K9Ac ]]; then
    perl ${params.codedir}/correctBGbyStablePromoters.pl --bg ${bg} --tss ${stableH3K9Ac} --out ${outBG}
  fi

  if [[ ${bg} =~ H3K4me3 ]]; then
    perl ${params.codedir}/correctBGbyStablePromoters.pl --bg ${bg} --tss ${stableH3K4me3} --out ${outBG}
  fi

  if [[ ${bg} =~ abK4me ]]; then
    perl ${params.codedir}/correctBGbyStablePromoters.pl --bg ${bg} --tss ${stableabK4me} --out ${outBG}
  fi

  echo ${olName} >header.txt

  mapBed -a ${tss}   -b ${outBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${olName}.TSS.tmp
  mapBed -a ${hstss} -b ${outBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${olName}.HSTSS.tmp
  mapBed -a ${m3pk}  -b ${outBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${olName}.K4ME3.tmp
  mapBed -a ${abpk}  -b ${outBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${olName}.ABK4ME.tmp
  mapBed -a ${k9pk}  -b ${outBG} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${olName}.K9AC.tmp

  cat header.txt ${olName}.TSS.tmp    >${olName}.TSS.ol
  cat header.txt ${olName}.HSTSS.tmp  >${olName}.HSTSS.ol
  cat header.txt ${olName}.K4ME3.tmp  >${olName}.K4ME3.ol
  cat header.txt ${olName}.ABK4ME.tmp >${olName}.ABK4ME.ol
  cat header.txt ${olName}.K9AC.tmp   >${olName}.K9AC.ol

  sort -k1,1 -k2n,2n -k3n,3n ${outBG} |grep -P \'^chr[\\dXY]+\' >tmp.bg
  bedGraphToBigWig tmp.bg ${mm10KMetIDX} ${outBW}
  """
  }

// THIS SIMPLY COPIES ZYGOTENE H3K4me3 A FILE CALLED H3K4me3_SCP3pos_H1Tneg.bam
// SIMPLIFIES LATER ANALYSES
process getEarlyH3K4me3Bam {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:400 --partition=norm'
  echo true
  cpus 4
  memory '12g'

  time { 2.hour }
  errorStrategy { 'retry' }
  maxRetries 3

  module 'samtools/1.9'
  module 'picard/2.17.11'

  publishDir params.outbam, mode: 'copy', overwrite: false

  input:
  	file(bam) from bamH3K4m3.collect()

  output:
    file('H3K4me3_SCP3pos_H1Tneg.bam')     into earlyK4m3Bam
    file('H3K4me3_SCP3pos_H1Tneg.bam.bai') into earlyK4m3Bai
    val "OK" into k4me3ok_a,k4me3ok_b,k4me3ok_c,k4me3ok_d

  script:
    """
    zyBam=`ls ZY.R1.H3K4me3.*.bam`
    cp \$zyBam H3K4me3_SCP3pos_H1Tneg.bam
    samtools index H3K4me3_SCP3pos_H1Tneg.bam
    """
 }

// MAKE CHANNELS FOR ALL HISTONE MODIFICATIONS AND INPUTS
Channel
     .fromPath(params.allHMBAMs)
     .ifEmpty { exit 1, "HM BAM files not found or mis-named" }
     .into {allHistModBAMs; allHistModBAMs_b; allHistModBAMs_c}

Channel
    .fromPath(params.allHMBAMsI)
    .ifEmpty { exit 1, "HM (input) BAM files not found or mis-named" }
    .into {allHistModBAMsI; allHistModBAMsIa}

allHistoneModBAMswithK4me3merge = earlyK4m3Bam.concat(allHistModBAMs)

Channel
     .fromPath(params.everyHMBAM)
     .ifEmpty { exit 1, "HM BAM files not found or mis-named" }
     .into {everyHistModBAMs; everyHistModBAMs_a}

allHistoneModBAMswithK4me3merge_a = earlyK4m3Bam.concat(everyHistModBAMs)
allHistoneModBAMswithK4me3merge_b = earlyK4m3Bam.concat(everyHistModBAMs_a)

process runNCIScorrection_allHMs {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '16g'

  time { 3.hour }

  tag  { bam }

  module 'bedtools/2.27.1'
  module 'R/3.5.2'
  module 'ucsc/373'

  publishDir params.outdirCoverage, mode: 'copy', overwrite: true, pattern: "*ws25bp.b*"

  input:
  file(bam)       from allHistoneModBAMswithK4me3merge
  file(tss)       from tssBEDe
  file(hstss)     from hstssBEDe
  file(m3pk)      from k4me3Peakse
  file(abpk)      from abk4Peakse
  file(k9pk)      from k9acPeakse
  file(bamm3I)    from earlyK4m3Bai.collect()
  file(bamI)      from allHistModBAMsI.collect()
  file(wins)      from winFilesd.collect()

  output:
  file '*.ws25bp.bedgraph' into allHMBG
  file '*.ws25bp.bigwig'   into allHMBWa, allHMBWb, allHMBWc, allHMBWd, allHMBWe, allHMBWf
  file '*.TSS.ol'          into tssOLHM
  file '*.HSTSS.ol'        into hstssOLHM
  file '*.K4ME3.ol'        into k4me3OLHM
  file '*.ABK4ME.ol'       into abk4OLHM
  file '*.K9AC.ol'         into k9acOLHM

  script:
  def name   = bam.name.replaceFirst("_SCP3pos_H1Tneg.bam" ,"")
  def bed    = bam.name.replaceFirst(".bam" ,".bed")
  def bedI   = bam.name.replaceFirst(".bam" ,".input.bed")
  def ncis   = bam.name.replaceFirst(".bam" ,".NCIS.tab")
  def bg     = bam.name.replaceFirst(".bam" ,".NCIS.ws25bp.bedgraph")
  def bw     = bam.name.replaceFirst(".bam" ,".NCIS.ws25bp.bigwig")
  def bgI    = bam.name.replaceFirst(".bam" ,".input.NCIS.ws25bp.bedgraph")
  def bwI    = bam.name.replaceFirst(".bam" ,".input.NCIS.ws25bp.bigwig")
  """
  bedtools bamtobed -i ${bam}  >${bed}
  bedtools bamtobed -i ${bamI} >${bedI}

  grep -w chr1 ${bed}  >chip.bed
  grep -w chr1 ${bedI} >input.bed

  perl ${params.codedir}/runNCIS.pl --t chip.bed --c input.bed --lib ${params.codedir}/NCIS --out ${ncis}
  ncis=`cut -f1 ${ncis}`

  mapBed -a win25.bed -b ${bed}  -c 3 -o count -g ${mm10KMetIDX} >${bg}.tmp
  mapBed -a win25.bed -b ${bedI} -c 3 -o count -g ${mm10KMetIDX} >${bgI}

  perl ${params.codedir}/correctBGwithNCIS.pl --t ${bg}.tmp --c ${bgI} --ncis \$ncis --out ${bg}

  grep -P \'^chr[\\d]+\\s\' ${bg} |sort -k1,1 -k2n,2n -k3n,3n -T "." >bwSort.bedgraph

  bedGraphToBigWig bwSort.bedgraph ${mm10KMetIDX} ${bw}

  echo ${name} >header.txt

  mapBed -a ${tss}   -b ${bg} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${name}.TSS.tmp
  mapBed -a ${hstss} -b ${bg} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${name}.HSTSS.tmp
  mapBed -a ${m3pk}  -b ${bg} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${name}.K4ME3.tmp
  mapBed -a ${abpk}  -b ${bg} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${name}.ABK4ME.tmp
  mapBed -a ${k9pk}  -b ${bg} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${name}.K9AC.tmp

  cat header.txt ${name}.TSS.tmp    >${name}.TSS.ol
  cat header.txt ${name}.HSTSS.tmp  >${name}.HSTSS.ol
  cat header.txt ${name}.K4ME3.tmp  >${name}.K4ME3.ol
  cat header.txt ${name}.ABK4ME.tmp >${name}.ABK4ME.ol
  cat header.txt ${name}.K9AC.tmp   >${name}.K9AC.ol
  """
  }

Channel
    .fromPath(params.allHMBAMsCTRL)
    .ifEmpty { exit 1, "HM (IgG/Input) BAM files not found or mis-named" }
    .set {allHistModBAMsCTRL}

// BUILD OVERLAP FILES FOR R TABLES
process ctrlBAMstoBW {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '16g'

  time { 3.hour }

  tag  { bam }

  module 'bedtools/2.27.1'
  module 'R/3.5.2'
  module 'ucsc/373'

  publishDir params.outdirCoverage, mode: 'copy', overwrite: true, pattern: "*ws25bp.b*"

  input:
  file(bam)       from allHistModBAMsCTRL
  file(tss)       from tssBEDk
  file(hstss)     from hstssBEDk
  file(m3pk)      from k4me3Peaksk
  file(abpk)      from abk4Peaksk
  file(k9pk)      from k9acPeaksk
  file(wins)      from winFilese.collect()

  output:
  file '*.ws25bp.bedgraph' into allHMiBG
  file '*.ws25bp.bigwig'   into allHMiBW, allHMiBWb
  file '*.TSS.ol'          into tssOLi
  file '*.HSTSS.ol'        into hstssOLi
  file '*.K4ME3.ol'        into k4me3OLi
  file '*.ABK4ME.ol'       into abk4OLi
  file '*.K9AC.ol'         into k9acOLi

  script:
  def name   = bam.name.replaceFirst("_SCP3pos_H1Tneg.bam" ,"")
  def bed    = bam.name.replaceFirst(".bam" ,".bed")
  def bg     = bam.name.replaceFirst(".bam" ,".ws25bp.bedgraph")
  def bw     = bam.name.replaceFirst(".bam" ,".ws25bp.bigwig")
  """
  bedtools bamtobed -i ${bam}  >${bed}

  mapBed -a win25.bed -b ${bed} -c 3 -o count -g ${mm10KMetIDX} >${bg}

  grep -P \'^chr[\\d]+\\s\' ${bg} |sort -k1,1 -k2n,2n -k3n,3n -T "." >bwSort.bedgraph

  bedGraphToBigWig bwSort.bedgraph ${mm10KMetIDX} ${bw}

  echo ${name} >header.txt

  mapBed -a ${tss}   -b ${bg} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${name}.TSS.tmp
  mapBed -a ${hstss} -b ${bg} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${name}.HSTSS.tmp
  mapBed -a ${m3pk}  -b ${bg} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${name}.K4ME3.tmp
  mapBed -a ${abpk}  -b ${bg} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${name}.ABK4ME.tmp
  mapBed -a ${k9pk}  -b ${bg} -c 4 -o sum -g ${mm10KMetIDX} |cut -f4 >${name}.K9AC.tmp

  cat header.txt ${name}.TSS.tmp    >${name}.TSS.ol
  cat header.txt ${name}.HSTSS.tmp  >${name}.HSTSS.ol
  cat header.txt ${name}.K4ME3.tmp  >${name}.K4ME3.ol
  cat header.txt ${name}.ABK4ME.tmp >${name}.ABK4ME.ol
  cat header.txt ${name}.K9AC.tmp   >${name}.K9AC.ol
  """
  }

process overlapAnnotation {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '8g'

  time { 1.hour }

  module 'bedtools/2.27.1'

  input:
  file(tss)       from tssBEDf
  file(hstss)     from hstssBEDf
  file(m3pk)      from k4me3Peaksf
  file(abpk)      from abk4Peaksf
  file(k9pk)      from k9acPeaksf
  file(hsU500)    from hotspotOneMotif500a
  file(expr)      from kallistoExpr
  file(spo11)     from spo11BGa
  file(hs3Kb)     from hs3Kb
  file(hsPrKO3Kb) from hsPrKO3Kb
  file(tss3Kb)    from tss3Kb

  output:
  file('*.TSS.ol')   into tssOLAnnot
  file('*.HSTSS.ol') into hstssOLAnnot
  file('*.K4ME3.ol') into k4me3OLAnnot
  file('*.ABK4ME.ol') into abk4OLAnnot
  file('*.K9AC.ol') into k9acOLAnnot

  script:
  """
  ## COMMENT
  echo "HS" >hotspots.TSS.ol
  echo "HS" >hotspots.HSTSS.ol
  echo "HS" >hotspots.K4ME3.ol
  echo "HS" >hotspots.ABK4ME.ol
  echo "HS" >hotspots.K9AC.ol

  echo "HSPrdm9KO" >hotspotsPrKO.TSS.ol
  echo "HSPrdm9KO" >hotspotsPrKO.HSTSS.ol
  echo "HSPrdm9KO" >hotspotsPrKO.K4ME3.ol
  echo "HSPrdm9KO" >hotspotsPrKO.ABK4ME.ol
  echo "HSPrdm9KO" >hotspotsPrKO.K9AC.ol

  echo "HSone" >hotspotsOneMotif.TSS.ol
  echo "HSone" >hotspotsOneMotif.HSTSS.ol
  echo "HSone" >hotspotsOneMotif.K4ME3.ol
  echo "HSone" >hotspotsOneMotif.ABK4ME.ol
  echo "HSone" >hotspotsOneMotif.K9AC.ol

  echo "TSS" >tssALL.TSS.ol
  echo "TSS" >tssALL.HSTSS.ol
  echo "TSS" >tssALL.K4ME3.ol
  echo "TSS" >tssALL.ABK4ME.ol
  echo "TSS" >tssALL.K9AC.ol

  echo -e "name\\tnumExp\\teLE\\teZY\\tePA\\teDI" >expr.TSS.ol
  echo -e "name\\tnumExp\\teLE\\teZY\\tePA\\teDI" >expr.HSTSS.ol
  echo -e "name\\tnumExp\\teLE\\teZY\\tePA\\teDI" >expr.K4ME3.ol
  echo -e "name\\tnumExp\\teLE\\teZY\\tePA\\teDI" >expr.ABK4ME.ol
  echo -e "name\\tnumExp\\teLE\\teZY\\tePA\\teDI" >expr.K9AC.ol

  echo "spo11" >spo11.TSS.ol
  echo "spo11" >spo11.HSTSS.ol
  echo "spo11" >spo11.K4ME3.ol
  echo "spo11" >spo11.ABK4ME.ol
  echo "spo11" >spo11.K9AC.ol

  ${params.codedir}/sortBEDByFAI.pl ${hsU500} ${mm10KMetIDX} >singleMotifHS.bed
  grep -vP "LE\\tZY" ${expr} |${params.codedir}/sortBEDByFAI.pl - ${mm10KMetIDX} >expr.bed

  mapBed -a ${tss}   -b ${hs3Kb}          -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspots.TSS.ol
  mapBed -a ${tss}   -b singleMotifHS.bed -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspotsOneMotif.TSS.ol
  mapBed -a ${tss}   -b ${hsPrKO3Kb}      -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspotsPrKO.TSS.ol
  mapBed -a ${tss}   -b ${tss3Kb}         -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>tssALL.TSS.ol

  mapBed -a ${tss}   -b expr.bed          -c 4,4,5,6,7,8 -o distinct,count,sum,sum,sum,sum -g ${mm10KMetIDX}|perl -pi -e \'s/\\t\\./\\t0/g\' |cut -f4-9 >>expr.TSS.ol
  mapBed -a ${tss}   -b ${spo11}          -c 4 -o sum -g ${mm10KMetIDX} |perl -pi -e \'s/\\t\\./\\t0/g\' |cut -f4 >>spo11.TSS.ol

  mapBed -a ${hstss} -b ${hs3Kb}          -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspots.HSTSS.ol
  mapBed -a ${hstss} -b singleMotifHS.bed -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspotsOneMotif.HSTSS.ol
  mapBed -a ${hstss} -b ${hsPrKO3Kb}      -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspotsPrKO.HSTSS.ol
  mapBed -a ${hstss} -b ${tss3Kb}         -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>tssALL.HSTSS.ol

  mapBed -a ${hstss} -b expr.bed          -c 4,4,5,6,7,8 -o distinct,count,sum,sum,sum,sum -g ${mm10KMetIDX}|perl -pi -e \'s/\\t\\./\\t0/g\' |cut -f4-9 >>expr.HSTSS.ol
  mapBed -a ${hstss} -b ${spo11}          -c 4 -o sum -g ${mm10KMetIDX} |perl -pi -e \'s/\\t\\./\\t0/g\' |cut -f4 >>spo11.HSTSS.ol

  mapBed -a ${m3pk}  -b ${hs3Kb}          -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspots.K4ME3.ol
  mapBed -a ${m3pk}  -b singleMotifHS.bed -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspotsOneMotif.K4ME3.ol
  mapBed -a ${m3pk}  -b ${hsPrKO3Kb}      -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspotsPrKO.K4ME3.ol
  mapBed -a ${m3pk}  -b ${tss3Kb}         -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>tssALL.K4ME3.ol

  mapBed -a ${m3pk}  -b expr.bed          -c 4,4,5,6,7,8 -o distinct,count,sum,sum,sum,sum -g ${mm10KMetIDX}|perl -pi -e \'s/\\t\\./\\t0/g\' |cut -f4-9 >>expr.K4ME3.ol
  mapBed -a ${m3pk}  -b ${spo11}          -c 4 -o sum -g ${mm10KMetIDX} |perl -pi -e \'s/\\t\\./\\t0/g\' |cut -f4 >>spo11.K4ME3.ol

  mapBed -a ${abpk}  -b ${hs3Kb}          -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspots.ABK4ME.ol
  mapBed -a ${abpk}  -b singleMotifHS.bed -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspotsOneMotif.ABK4ME.ol
  mapBed -a ${abpk}  -b ${hsPrKO3Kb}      -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspotsPrKO.ABK4ME.ol
  mapBed -a ${abpk}  -b ${tss3Kb}         -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>tssALL.ABK4ME.ol

  mapBed -a ${abpk}  -b expr.bed          -c 4,4,5,6,7,8 -o distinct,count,sum,sum,sum,sum -g ${mm10KMetIDX}|perl -pi -e \'s/\\t\\./\\t0/g\' |cut -f4-9 >>expr.ABK4ME.ol
  mapBed -a ${abpk}  -b ${spo11}          -c 4 -o sum -g ${mm10KMetIDX} |perl -pi -e \'s/\\t\\./\\t0/g\' |cut -f4 >>spo11.ABK4ME.ol

  mapBed -a ${k9pk}  -b ${hs3Kb}          -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspots.K9AC.ol
  mapBed -a ${k9pk}  -b singleMotifHS.bed -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspotsOneMotif.K9AC.ol
  mapBed -a ${k9pk}  -b ${hsPrKO3Kb}      -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>hotspotsPrKO.K9AC.ol
  mapBed -a ${k9pk}  -b ${tss3Kb}         -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>tssALL.K9AC.ol

  mapBed -a ${k9pk}  -b expr.bed          -c 4,4,5,6,7,8 -o distinct,count,sum,sum,sum,sum -g ${mm10KMetIDX}|perl -pi -e \'s/\\t\\./\\t0/g\' |cut -f4-9 >>expr.K9AC.ol
  mapBed -a ${k9pk}  -b ${spo11}          -c 4 -o sum -g ${mm10KMetIDX} |perl -pi -e \'s/\\t\\./\\t0/g\' |cut -f4 >>spo11.K9AC.ol
  """
  }

process overlapH3K4me3Peaks {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '8g'

  time { 1.hour }

  module 'bedtools/2.27.1'

  input:
  file (k4me3Peaks) from peaksBEDTCb
  file(tss)       from tssBEDh
  file(hstss)     from hstssBEDa
  file(m3pk)      from k4me3Peaksg
  file(abpk)      from abk4Peaksg
  file(k9pk)      from k9acPeaksg

  output:
  file('*.TSS.ol')    into tssOLPeaks
  file('*.HSTSS.ol')  into hstssOLPeaks
  file('*.K4ME3.ol')  into k4me3OLPeaks
  file('*.ABK4ME.ol') into abk4OLPeaks
  file('*.K9AC.ol')   into k9acOLPeaks

  script:
  def name   = k4me3Peaks.name.replaceFirst(".bed" ,"")
  """
  ## COMMENT
  echo "${name}" >${name}.TSS.ol
  echo "${name}" >${name}.HSTSS.ol
  echo "${name}" >${name}.K4ME3.ol
  echo "${name}" >${name}.ABK4ME.ol
  echo "${name}" >${name}.K9AC.ol

  ${params.codedir}/sortBEDByFAI.pl ${k4me3Peaks} ${mm10KMetIDX} >peaks.bed

  mapBed -a ${tss}   -b peaks.bed         -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>${name}.TSS.ol
  mapBed -a ${hstss} -b peaks.bed         -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>${name}.HSTSS.ol
  mapBed -a ${m3pk}  -b peaks.bed         -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>${name}.K4ME3.ol
  mapBed -a ${abpk}  -b peaks.bed         -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>${name}.ABK4ME.ol
  mapBed -a ${k9pk}  -b peaks.bed         -c 3 -o count -g ${mm10KMetIDX}|cut -f4 |perl -lane \'print ((\$F[0] eq 0)?"FALSE":"TRUE")\' >>${name}.K9AC.ol
  """
  }

// MAKE OUTPUT TABLES AND SLICE DATA
Channel.from( ["chr1", 176790000, 177130000, "chr1_177Mb"],
              ["chr1",  21875000,  22125000, "chr1_22Mb"],
              ["chr8", 111517000, 111676000, "chr8_115Mb"],
              ["chr9", 115050000, 115350000, "chr9_115Mb"],
              ["chr10", 92520000,  92730000, "chr10_92Mb"])
       .set {slices2use}

process makeSlices {
    scratch '/lscratch/$SLURM_JOBID'
    clusterOptions ' --gres=lscratch:40'
    echo true
    cpus 1
    memory '8g'

    time { 1.hour }

    module 'bedtools/2.27.1'
    module 'ucsc/373'
    module 'R/3.5.2'

    publishDir params.outdirSlices, mode: 'copy', overwrite: true

    input:
    set val(nCS), val(nFrom), val(nTo), val(sName) from slices2use
    //file(tss)                                      from tssNoMergeBED1kb_b
    file(tss)                                      from refseqTSSBEDb
    file(hs)                                       from hotspotBED
    file(win4bg)                                   from winFilesb.collect()
    file(bgKMet)                                   from correctedH3K4BGb.collect()
    file(bgNCIS)                                   from stableTSSscorrectedBGa.collect()

    output:
    file("*DZ.tab")        into sliceDZ
    file("*chr*bedgraph")  into (sliceBG, sliceBGa)
    file("HS*bed")         into sliceHS
    file("TSS*bed")        into sliceTSS

    script:
    """
    echo -e "${nCS}\\t${nFrom}\\t${nTo}" >region.bed

    intersectBed -a ${tss}  -b region.bed -wa -u |cut -f1-3 >TSS.${sName}.bed
    intersectBed -a ${hs}   -b region.bed -wa -u |cut -f1-3 >HS.${sName}.bed

    intersectBed -a win1ks147.bed -b region.bed -wa -u >wins.${sName}.bed

    for stage in 'LE.R1.H3K4me3' 'ZY.R1.H3K4me3' 'EP.R1.H3K4me3' 'LP.R1.H3K4me3' 'DI.R1.H3K4me3'; do
        inBG=\$stage".monoCorrected.ws25bp.bedgraph"
        outBG=\$stage"_KM.${sName}.bedgraph"
        outDZ=\$stage"_KM.${sName}.DZ.tab"

        mapBed -a wins.${sName}.bed -b \$inBG -c 4 -o sum -g ${mm10KMetIDX} >\$outBG

		    perl -lane \'\$min = int((\$F[1]+\$F[2])/2); print join("\\t",\$F[0],\$min,\$F[3])\' \$outBG >\$outDZ
    done

    for stage in 'LE.R1.H3K4me3' 'ZY.R1.H3K4me3' 'EP.R1.H3K4me3' 'LP.R1.H3K4me3' 'DI.R1.H3K4me3' \
                 'LE.R1.H3K9Ac' 'ZY.R1.H3K9Ac' 'EP.R1.H3K9Ac' 'LP.R1.H3K9Ac' 'DI.R1.H3K9Ac'; do
        inBG=\$stage".ncisCorrectedByStablePromoters.ws25bp.bedgraph"
        outBG=\$stage"_ncis.${sName}.bedgraph"
        outDZ=\$stage"_ncis.${sName}.DZ.tab"

        mapBed -a wins.${sName}.bed -b \$inBG -c 4 -o sum -g ${mm10KMetIDX} >\$outBG

        perl -lane \'\$min = int((\$F[1]+\$F[2])/2); print join("\\t",\$F[0],\$min,\$F[3])\' \$outBG >\$outDZ
    done
    """
    }

process makeRtablesTSS {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '8g'

  time { 1.hour }

  publishDir params.outdirRTables, mode: 'copy', overwrite: true

  input:
  file(tss) from tssBEDg
  file(tOL) from tssOLKM.collect()
  file(nOL) from tssOLNCISTSS.collect()
  file(aOL) from tssOLAnnot.collect()
  file(sOL) from tssOLstable.collect()
  file(hOL) from tssOLHM.collect()
  file(pOL) from tssOLPeaks.collect()

  output:
  file("tssTable.forR.tab") into (tssTableForR, tssTableForRa, tssTableForRb)

  script:
  """
  echo -e "cs\\tfrom\\tto" >tssTable.tab

  cat ${tss}  >>tssTable.tab

  paste tssTable.tab *TSS.ol |grep -P "^(cs|chr[\\d]+)\\s" >tssTable.forR.tab
  """
  }

process makeRtablesHSTSS {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '8g'

  time { 1.hour }

  publishDir params.outdirRTables, mode: 'copy', overwrite: true

  input:
  file(hstss) from hstssBEDDets
  file(tOL)   from hstssOLKM.collect()
  file(nOL)   from hstssOLNCISTSS.collect()
  file(aOL)   from hstssOLAnnot.collect()
  file(sOL)   from hstssOLstable.collect()
  file(hOL)   from hstssOLHM.collect()
  file(pOL)   from hstssOLPeaks.collect()

  output:
  file('hstssTable.forR.tab') into (hstssTableForR, hstssTableForRa, hstssTableForRb)

	script:
	"""
	echo -e "cs\\tfrom\\tto\\tHSstrength\\ttype\\tstrand" >hstssTable.tab
  cat ${hstss}  >>hstssTable.tab

  paste hstssTable.tab *HSTSS.ol |grep -P "^(cs|chr[\\d]+)\\s" >hstssTable.forR.tab
  """
  }

process makeRtablesK4ME3 {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '8g'

  time { 1.hour }

  publishDir params.outdirRTables, mode: 'copy', overwrite: true

  input:
  file(m3pk) from k4me3PeaksBG
  file(tOL)  from k4me3OLKM.collect()
  file(rOL)  from k4me3OLNCISRAW.collect()
  file(nOL)  from k4me3OLNCISTSS.collect()
  file(aOL)  from k4me3OLAnnot.collect()
  file(sOL)  from k4me3OLstable.collect()
  file(hOL)  from k4me3OLHM.collect()
  file(pOL)  from k4me3OLPeaks.collect()

  output:
  file 'k4me3Table.forR.tab' into (k4me3TableForR, k4me3TableForRa, k4me3TableForRb)

  script:
  """

  echo -e "cs\\tfrom\\tto\\tnpeaks" >k4me3Table.tab
  cat ${m3pk}  >>k4me3Table.tab

  paste k4me3Table.tab *K4ME3.ol |grep -P "^(cs|chr[\\d]+)\\s" >k4me3Table.forR.tab
  """
  }

process makeRtablesABK4ME {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '8g'

  time { 1.hour }

  publishDir params.outdirRTables, mode: 'copy', overwrite: true

  input:
  file(pk)   from abk4PeaksBG
  file(tOL)  from abk4OLKM.collect()
  file(rOL)  from abk4OLNCISRAW.collect()
  file(nOL)  from abk4OLNCISTSS.collect()
  file(aOL)  from abk4OLAnnot.collect()
  file(sOL)  from abk4OLstable.collect()
  file(hOL)  from abk4OLHM.collect()
  file(pOL)  from abk4OLPeaks.collect()

  output:
  file 'abk4meTable.forR.tab' into abk4meTableForR

  script:
  """

  echo -e "cs\\tfrom\\tto\\tnpeaks" >abk4meTable.tab
  cat ${pk}  >>abk4meTable.tab

  paste abk4meTable.tab *ABK4ME.ol |grep -P "^(cs|chr[\\d]+)\\s" >abk4meTable.forR.tab
  """
  }

process makeRtablesK9AC {
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '8g'

  time { 1.hour }

  publishDir params.outdirRTables, mode: 'copy', overwrite: true

  input:
  file(pk)  from k9acPeaksBG
  file(tOL) from k9acOLKM.collect()
  file(rOL) from k9acOLNCISRAW.collect()
  file(nOL) from k9acOLNCISTSS.collect()
  file(aOL) from k9acOLAnnot.collect()
  file(sOL) from k9acOLstable.collect()
  file(hOL) from k9acOLHM.collect()
  file(pOL) from k9acOLPeaks.collect()

  output:
  file 'k9acTable.forR.tab' into k9acTableForR

  script:
  """
  echo -e "cs\\tfrom\\tto\\tnpeaks" >k9acTable.tab
  cat ${pk}  >>k9acTable.tab

  paste k9acTable.tab *K9AC.ol |grep -P "^(cs|chr[\\d]+)\\s" >k9acTable.forR.tab
  """
  }

// MAKE ALL HEATMAPS
process makeHeatmapsForH3K9Ac{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 8
  memory '12g'

  time { 1.hour }

  module 'deeptools/3.0.1'
  module 'bedtools/2.27.1'

  publishDir params.outdirImages, mode: 'copy', overwrite: true

  input:
  file (lepBW) from stableTSSscorrectedLEBG_K9Ac // ONLY used to assure we run on each replicate
  file bw      from stableTSSscorrectedBWa.collect()
  file bwO     from ncisBWa.collect()
  file (hs)    from hotspotOneMotif500b
  file (tss)   from tssBEDi

  output:
  file "*matrix.gz"  into hmOutMatrixK9Ac
  file "*png"        into hmOutPNGK9Ac
  file "*svg"        into hmOutSVGK9Ac

  script:
  def nRA = lepBW.name.replaceFirst(".H3K9Ac.ncisCorrectedByStablePromoters.ws25bp.bedgraph", "")
  def rep = nRA.replaceFirst("LE.", "")
  """
  echo "OK ... "
  intersectBed -a ${hs} -b ${tss} -v -wa >hs.oneMotif.bed

  computeMatrix reference-point --referencePoint center -R LE ZY EP LP DI -a 3000 -b 3000 \
                                -R hs.oneMotif.bed \
                                -S LE.${rep}.H3K9Ac.ncisCorrected.ws25bp.bigwig \
                                   ZY.${rep}.H3K9Ac.ncisCorrected.ws25bp.bigwig \
                                   EP.${rep}.H3K9Ac.ncisCorrected.ws25bp.bigwig \
                                   LP.${rep}.H3K9Ac.ncisCorrected.ws25bp.bigwig \
                                   DI.${rep}.H3K9Ac.ncisCorrected.ws25bp.bigwig \
                                   -p max -o H3K9Ac.${rep}.ncisCorrectedRaw.HS.matrix.gz \
                                   --binSize 25

  plotHeatmap -m H3K9Ac.${rep}.ncisCorrectedRaw.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean H3K4me3" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat png \
              -o H3K9Ac.${rep}.ncisCorrectedRaw.HS.png

  plotHeatmap -m H3K9Ac.${rep}.ncisCorrectedRaw.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean H3K4me3" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat svg \
              -o H3K9Ac.${rep}.ncisCorrectedRaw.HS.svg

  computeMatrix reference-point --referencePoint center -R LE ZY EP LP DI -a 3000 -b 3000 \
                  -R hs.oneMotif.bed \
                  -S LE.${rep}.H3K9Ac.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                     ZY.${rep}.H3K9Ac.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                     EP.${rep}.H3K9Ac.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                     LP.${rep}.H3K9Ac.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                     DI.${rep}.H3K9Ac.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                     -p max -o H3K9Ac.${rep}.ncisCorrectedByStablePromoters.HS.matrix.gz \
                     --binSize 25

  plotHeatmap -m H3K9Ac.${rep}.ncisCorrectedByStablePromoters.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean H3K9Ac" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat png \
              -o H3K9Ac.${rep}.ncisCorrectedByStablePromoters.HS.png

  plotHeatmap -m H3K9Ac.${rep}.ncisCorrectedByStablePromoters.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean H3K9Ac" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat svg \
              -o H3K9Ac.${rep}.ncisCorrectedByStablePromoters.HS.svg

  """
  }

process makeHeatmapsForH3K4me3{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 12
  memory '16g'

  time { 1.hour }

  module 'deeptools/3.0.1'
  module 'bedtools/2.27.1'

  publishDir params.outdirImages, mode: 'copy', overwrite: true

  input:
  file (lepBW) from stableTSSscorrectedLEBG_K4Me3 // ONLY used to assure we run on each replicate
  file bwN     from stableTSSscorrectedBWb.collect()
  file bwM     from correctedH3K4BWa.collect()
  file bwO     from ncisBWb.collect()
  file (hs)    from hotspotOneMotif500c
  file (tss)   from tssBEDj

  output:
  file "*matrix.gz"  into hmOutMatrixK4me3
  file "*png"     into hmOutPNGK4me3
  file "*svg"     into hmOutSVGK4me3

  script:
  def nRA = lepBW.name.replaceFirst(".H3K4me3.ncisCorrectedByStablePromoters.ws25bp.bedgraph", "")
  def rep = nRA.replaceFirst("LE.", "")
  """
  echo "OK"
  intersectBed -a ${hs} -b ${tss} -v -wa >hs.oneMotif.bed

  computeMatrix reference-point --referencePoint center -R LE ZY EP LP DI -a 3000 -b 3000 \
                                -R hs.oneMotif.bed \
                                -S LE.${rep}.H3K4me3.ncisCorrected.ws25bp.bigwig \
                                   ZY.${rep}.H3K4me3.ncisCorrected.ws25bp.bigwig \
                                   EP.${rep}.H3K4me3.ncisCorrected.ws25bp.bigwig \
                                   LP.${rep}.H3K4me3.ncisCorrected.ws25bp.bigwig \
                                   DI.${rep}.H3K4me3.ncisCorrected.ws25bp.bigwig \
                                   -p max -o H3K4me3.${rep}.ncisCorrectedRaw.HS.matrix.gz \
                                   --binSize 25

  plotHeatmap -m H3K4me3.${rep}.ncisCorrectedRaw.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean H3K4me3" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat png \
              -o H3K4me3.${rep}.ncisCorrectedRaw.HS.png

  plotHeatmap -m H3K4me3.${rep}.ncisCorrectedRaw.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean H3K4me3" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat svg \
              -o H3K4me3.${rep}.ncisCorrectedRaw.HS.svg

  computeMatrix reference-point --referencePoint center -R LE ZY EP LP DI -a 3000 -b 3000 \
                                -R hs.oneMotif.bed \
                                -S LE.${rep}.H3K4me3.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                                   ZY.${rep}.H3K4me3.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                                   EP.${rep}.H3K4me3.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                                   LP.${rep}.H3K4me3.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                                   DI.${rep}.H3K4me3.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                                   -p max -o H3K4me3.${rep}.ncisCorrectedByStablePromoters.HS.matrix.gz \
                                   --binSize 25

  plotHeatmap -m H3K4me3.${rep}.ncisCorrectedByStablePromoters.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean H3K4me3" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat png \
              -o H3K4me3.${rep}.ncisCorrectedByStablePromoters.HS.png

  plotHeatmap -m H3K4me3.${rep}.ncisCorrectedByStablePromoters.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean H3K4me3" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat svg \
              -o H3K4me3.${rep}.ncisCorrectedByStablePromoters.HS.svg

  computeMatrix reference-point --referencePoint center -R LE ZY EP LP DI -a 3000 -b 3000 \
                                -R hs.oneMotif.bed  \
                                -S LE.${rep}.H3K4me3.monoCorrected.ws25bp.bigwig \
                                   ZY.${rep}.H3K4me3.monoCorrected.ws25bp.bigwig \
                                   EP.${rep}.H3K4me3.monoCorrected.ws25bp.bigwig \
                                   LP.${rep}.H3K4me3.monoCorrected.ws25bp.bigwig \
                                   DI.${rep}.H3K4me3.monoCorrected.ws25bp.bigwig \
                                   -p max -o H3K4me3.${rep}.monoCorrected.HS.matrix.gz \
                                   --binSize 25

  plotHeatmap -m H3K4me3.${rep}.monoCorrected.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean H3K4me3" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat png \
              -o H3K4me3.${rep}.monoCorrected.HS.png

  plotHeatmap -m H3K4me3.${rep}.monoCorrected.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean H3K4me3" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat svg \
              -o H3K4me3.${rep}.monoCorrected.HS.svg

  """
  }

process makeHeatmapsForAbK4me{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 8
  memory '12g'

  time { 1.hour }

  module 'deeptools/3.0.1'
  module 'bedtools/2.27.1'

  publishDir params.outdirImages, mode: 'copy', overwrite: true

  input:
  file (lepBW) from stableTSSscorrectedLEBG_abK4 // ONLY used to assure we run on each replicate
  file bw      from stableTSSscorrectedBWc.collect()
  file bwO     from ncisBWc.collect()
  file (hs)    from hotspotOneMotif500b
  file (tss)   from tssBEDi

  output:
  file "*matrix.gz"  into hmOutMatrixabK4
  file "*png"     into hmOutPNGabK4
  file "*svg"     into hmOutSVGabK4

  script:
  def nRA = lepBW.name.replaceFirst(".abK4me.ncisCorrectedByStablePromoters.ws25bp.bedgraph", "")
  def rep = nRA.replaceFirst("LE.", "")
  """
  echo "OK ... "
  intersectBed -a ${hs} -b ${tss} -v -wa >hs.oneMotif.bed

  computeMatrix reference-point --referencePoint center -R LE ZY EP LP DI -a 3000 -b 3000 \
                                -R hs.oneMotif.bed \
                                -S LE.${rep}.abK4me.ncisCorrected.ws25bp.bigwig \
                                   ZY.${rep}.abK4me.ncisCorrected.ws25bp.bigwig \
                                   EP.${rep}.abK4me.ncisCorrected.ws25bp.bigwig \
                                   LP.${rep}.abK4me.ncisCorrected.ws25bp.bigwig \
                                   DI.${rep}.abK4me.ncisCorrected.ws25bp.bigwig \
                                   -p max -o abK4me.${rep}.ncisCorrectedRaw.HS.matrix.gz \
                                   --binSize 25

  plotHeatmap -m abK4me.${rep}.ncisCorrectedRaw.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean abK4me" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat png \
              -o abK4me.${rep}.ncisCorrectedRaw.HS.png

  plotHeatmap -m abK4me.${rep}.ncisCorrectedRaw.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean abK4me" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat svg \
              -o abK4me.${rep}.ncisCorrectedRaw.HS.svg

  computeMatrix reference-point --referencePoint center -R LE ZY EP LP DI -a 3000 -b 3000 \
                  -R hs.oneMotif.bed \
                  -S LE.${rep}.abK4me.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                     ZY.${rep}.abK4me.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                     EP.${rep}.abK4me.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                     LP.${rep}.abK4me.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                     DI.${rep}.abK4me.ncisCorrectedByStablePromoters.ws25bp.bigwig \
                     -p max -o abK4me.${rep}.ncisCorrectedByStablePromoters.HS.matrix.gz \
                     --binSize 25

  plotHeatmap -m abK4me.${rep}.ncisCorrectedByStablePromoters.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean abK4me" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat png \
              -o abK4me.${rep}.ncisCorrectedByStablePromoters.HS.png

  plotHeatmap -m abK4me.${rep}.ncisCorrectedByStablePromoters.HS.matrix.gz --samplesLabel LE ZY EP LP DI \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean abK4me" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat svg \
              -o abK4me.${rep}.ncisCorrectedByStablePromoters.HS.svg

  """
  }

Channel
     .fromPath(params.everyHMBAM)
     .ifEmpty { exit 1, "HM BAM files not found or mis-named" }
     .into {everyHMBAM; everyHMBAMb}

process makeRPKMhistoneBigWigs{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 16
  memory '16g'

  time { 4.hour * task.attempt }
  tag  {bam}

  module 'deeptools/3.0.1'
  module 'bedtools/2.27.1'
  module 'ucsc/373'

  publishDir params.outdirCoverage, mode: 'copy', overwrite: true

  input:
  //file (bam) from everyHMBAM
  file (bam) from allHistoneModBAMswithK4me3merge_a

  output:
  file "*RPKM.bigwig" into rpkmBWHMa, rpkmBWHMb

  script:
  def bw = bam.name.replaceFirst("bam","150bp.RPKM.bigwig")
  """
  samtools index ${bam}
  bamCoverage --bam ${bam} --outFileName ${bw} --outFileFormat bigwig --binSize 150 -p max --ignoreForNormalization chrX chrM chrY --normalizeUsing RPKM
  """
  }

Channel.from(["abK4me","H3K4me3" ,"H3K9ac"  ,"H3K36me3","H3K4me2","H3K4me1",
              "H3K27ac" ,"H4ac5"   ,"H4K20me3","H4K12ac","H4K8ac",
              "H3K79me1","H3K27me3","H3K79me3","H3K4ac" ,"H3K27me1",
              "H3K9me3" ,"H3K9me2" ,"H3"      ,"Input"  ,"IgG"])
       .into {hm4supplementOrder; hm4supplementOrder2; hm4supplementOrder3}

process makeIntervalsForHistoneModsSupp{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 8
  memory '12g'

  time { 0.5.hour }
  tag  { "ALLHM" }

  module 'deeptools/3.0.1'
  module 'bedtools/2.27.1'
  module 'ucsc/373'

  //publishDir params.outdirAnnot, mode: 'copy', overwrite: true, pattern: "*bed"

  input:
  file (hs)   from hotspotOneMotif500e
  file (tss)  from refseqTSSBED
  file (tes)  from refseqTESBED
  file (gene) from refseqGeneBED
  file (line) from rmLINEBED
  file (enh)  from refseqEnhancersBED

  output:
  file("hm*bed") into suppHMIntervals

  script:
  """
  shuf ${hs}   |head -n 10000 |sort -k1,1 -k2n,2n >hmHS.bed
  shuf ${tss}  |head -n 10000 |sort -k1,1 -k2n,2n >hmTSS.bed
  shuf ${tes}  |head -n 10000 |sort -k1,1 -k2n,2n >hmTES.bed
  shuf ${gene} |head -n 10000 |sort -k1,1 -k2n,2n >hmCDS.bed
  shuf ${line} |head -n 10000 |sort -k1,1 -k2n,2n >hmLINE.bed
  shuf ${enh}  |head -n 10000 |sort -k1,1 -k2n,2n >hmENH.bed
  """
  }

process makeHeatmapsForHistoneModsSuppFig{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 16
  memory '16g'

  time { 4.hour }
  tag  { ibed }

  module 'deeptools/3.0.1'
  module 'bedtools/2.27.1'
  module 'ucsc/373'

  publishDir params.outdirImages, mode: 'copy', overwrite: true, pattern: "*png"
  publishDir params.outdirTables, mode: 'copy', overwrite: true, pattern: "*matrix.gz"

  input:
  each file (ibed) from suppHMIntervals
  val  (hm)        from hm4supplementOrder2.collect()
  file (bwall)     from rpkmBWHMa.collect()

  output:
  file("*matrix.gz") into matrixSuppHM
  file("*png")       into heatmapSuppHM

  script:
  def nm       = ibed.name.replaceFirst(".bed","")
  def matrix   = ibed.name.replaceFirst(".bed",".matrix.gz")
  def pngout   = ibed.name.replaceFirst(".bed",".heatmap.png")
  def hmNames  = hm.join(" ")
  def bw       = hm.join("_SCP3pos_H1Tneg.150bp.RPKM.bigwig ")
  """
  if [[ "${nm}" == "hmLINE" || "${nm}" == "hmCDS" ]] ; then
    computeMatrix scale-regions   -o ${matrix} --regionsFileName ${ibed} -S ${bw}_SCP3pos_H1Tneg.150bp.RPKM.bigwig -a 5000 -b 5000 -p max  -bs 50 -m 10000 --skipZeros --missingDataAsZero
  else
    if [ "${nm}" == "hmHS" ] ; then
      computeMatrix reference-point -o ${matrix} --regionsFileName ${ibed} -S ${bw}_SCP3pos_H1Tneg.150bp.RPKM.bigwig --referencePoint center -a 5000 -b 5000 -p max
    else
      computeMatrix reference-point -o ${matrix} --regionsFileName ${ibed} -S ${bw}_SCP3pos_H1Tneg.150bp.RPKM.bigwig --referencePoint center -a 5000 -b 5000 -p max
    fi
  fi

  plotHeatmap -m ${matrix} --colorMap hot --outFileName ${pngout} --samplesLabel ${hmNames} --regionsLabel ${nm} --heatmapHeight 12 --xAxisLabel distance --dpi 300 --plotFileFormat png

  """
  }

Channel.from( ["set1","abK4me" ,"H3K9ac"  ,"H3K36me3","H3K4me2","H3K4me1"],
              ["set2","H3K27ac" ,"H4ac5"   ,"H4K20me3","H4K12ac","H4K8ac"],
              ["set3","H3K79me1","H3K27me3","H3K79me3","H3K4ac" ,"H3K27me1"],
              ["set4","H3K9me3" ,"H3K9me2" ,"H3"      ,"Input"  ,"IgG"])
       .set {hm4heatmaps}

process makeHeatmapsForHistoneMods{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 16
  memory '16g'

  time { 4.hour }
  tag  { nSet }

  module 'deeptools/3.0.1'
  module 'bedtools/2.27.1'

  publishDir params.outdirImages, mode: 'copy', overwrite: true

  input:
  set val(nSet), val(hm1), val(hm2), val(hm3), val(hm4), val(hm5) from hm4heatmaps
  file (bwAll)                                                    from rpkmBWHMb.collect()
  file (hsBED)                                                    from hotspotOneMotif500d

  output:
  file "*matrix.gz"  into hmOutMatrixHMx5
  file "*png"     into hmOutPNGHMx5
  file "*svg"     into hmOutSVGHMx5

  script:
  def bw1 = "${hm1}_SCP3pos_H1Tneg.150bp.RPKM.bigwig"
  def bw2 = "${hm2}_SCP3pos_H1Tneg.150bp.RPKM.bigwig"
  def bw3 = "${hm3}_SCP3pos_H1Tneg.150bp.RPKM.bigwig"
  def bw4 = "${hm4}_SCP3pos_H1Tneg.150bp.RPKM.bigwig"
  def bw5 = "${hm5}_SCP3pos_H1Tneg.150bp.RPKM.bigwig"
  """
  computeMatrix reference-point --referencePoint center -a 2500 -b 2500 \
                                --regionsFileName ${hsBED} \
                                -S ${bw1} \
                                   ${bw2} \
                                   ${bw3} \
                                   ${bw4} \
                                   ${bw5} \
                                -p max -o ${nSet}.ncisCorrectedRaw.HS.matrix.gz \
                                --binSize 25

  plotHeatmap -m ${nSet}.ncisCorrectedRaw.HS.matrix.gz --samplesLabel ${hm1} ${hm2} ${hm3} ${hm4} ${hm5} \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean signal" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat png --heatmapHeight 12 --yMax 9.7 --zMax 9.7  \
              -o histoneMods.${nSet}.HS.png

  plotHeatmap -m ${nSet}.ncisCorrectedRaw.HS.matrix.gz --samplesLabel ${hm1} ${hm2} ${hm3} ${hm4} ${hm5} \
              --colorMap RdGy_r --refPointLabel HS --yAxisLabel "Mean signal" --xAxisLabel 'distance (Kb)' \
              --regionsLabel "DSB hotspots" --plotFileFormat svg  --heatmapHeight 12 --yMax 9.7 --zMax 9.7 \
              -o histoneMods.${nSet}.HS.svg
  """
  }

// DRAW FIGURES
Channel.from( ["slice1","chr17:35100000-35140000"],
              ["slice2","chr4:139350000-139383000"],
              ["slice3","chr17:31842000-32125000"],
              ["slice4","chr17:31269000-31548000"],
              ["slice5","chr5:129994114-130055807"],
              ["slice6","chr8:111517000-111676000"],
              ["slice7","chr1:86042000-86236000"],
              ["slice8","chr10:92526000-92727000"],
              ["slice9","chr9:83788000-83941000"],)
       .set {fig4Slices}

process makeSliceDataForFigure4{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 2
  memory '8g'

  time { 0.5.hour }
  tag  { slice }

  module 'bedtools/2.27.1'
  module 'samtools/1.9'

  publishDir params.outdirFig4Slices, mode: 'copy', overwrite: true

  input:
  set  val(slice), val(coords) from fig4Slices
  file (tss)                   from refseqTSSBEDa
  file (hs)                    from hs3Kbb
  file (hsKO)                  from hsPrKO3Kbb
  file (bams)                  from allHistoneModBAMswithK4me3merge_b.collect()
  file (winz)                  from winFilesf.collect()

  output:
  file "slice*bed"     into sliceAnnotfig4
  file "*ALL.bedgraph" into sliceBGfig4a, sliceBGfig4b, sliceBGfig4c

  script:
  """
  echo $coords |perl -pi -e 's/[:\\-]/\t/g' >coord.bed

  intersectBed -a ${tss}   -b coord.bed -wa -u               |cut -f1-3 |sort -k1,1 -k2n,2n >${slice}.TSSgencode.bed
  intersectBed -a ${hs}    -b coord.bed -wa -u               |cut -f1-3 |sort -k1,1 -k2n,2n >${slice}.HSB6.bed
  intersectBed -a ${hsKO}  -b coord.bed -wa -u               |cut -f1-3 |sort -k1,1 -k2n,2n >${slice}.HSKO.bed

  intersectBed -a win1ks147.bed -b coord.bed -wa >coordWindows.bed

  for bam in *.bam; do
    nm=\${bam/_SCP3pos_H1Tneg.bam/}
    sliceBam="${slice}.\$nm.bam"
    sliceDZ="${slice}.\$nm.DZtpm.bedgraph"

    bamLink=`readlink -f \$bam`
    ln -s \$bamLink.bai \$bam.bai

    corrFactor=`samtools idxstats \$bam |perl -lane '\$c += \$F[2]; print \$c/1000000' |tail -n1`

    samtools view -hb \$bam $coords >\$sliceBam
    samtools index \$sliceBam

    intersectBed -a coordWindows.bed -b \$sliceBam -c |perl -lane \'print join("\\t",\$F[0],int((\$F[1]+\$F[2])/2),\$F[3]/1000/'\$corrFactor',"'\$nm'")' >\$sliceDZ
  done

  cat *DZtpm.bedgraph >${slice}.ALL.bedgraph
  """
  }

process drawFigure3{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 8
  memory '24g'

  time { 8.hour }
  tag  { "Drawing Figure 3" }

  module 'R/3.5.2'

  publishDir params.outdirImages, mode: 'copy', overwrite: true

  input:
  file(peakTab)  from k4me3TableForRa
  file(hstssTab) from hstssTableForRa
  file(tssTab)   from tssTableForRa
  file(s3Done)   from sliceBGa.collect()

  output:
  file "*png"    into figure3PNGs
  file "*pdf"    into figure3PDFs

  script:
  """
  ## Add globals to .Renviron (so we can see them inside R!)
  echo KBPIPEOUTDIR="${params.outdir}/"      >>~/.Renviron
  echo KBPIPEWORKDIR="${params.projectdir}/" >>~/.Renviron

  R --no-save --no-site-file --no-init-file --no-restore --silent --slave <${params.codeRdir}/plotFigure3.R ||true
  """
  }

process drawKmetStatFig{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 2
  memory '8g'

  time { 1.hour }

  module 'R/3.5.2'

  publishDir params.outdirImages, mode: 'copy', overwrite: true, pattern: '*png'
  publishDir params.outdirImages, mode: 'copy', overwrite: true, pattern: '*pdf'
  publishDir params.outdirTables, mode: 'copy', overwrite: true, pattern: '*tab'

  input:

  output:
  file "KmetStat_Enrichment.tab" into dataKmetStatEnrichment
  file "*png"                    into figureKmetStatPNG
  file "*png"                    into figureKmetStatPDF

  script:

  """
  ## Add globals to .Renviron (so we can see them inside R!)
  echo KBPIPEOUTDIR="${params.outdir}/"      >>~/.Renviron
  echo KBPIPEWORKDIR="${params.projectdir}/" >>~/.Renviron

  echo -e "stage\\trep\\tab\\tsi\\tn" >KmetStat_Enrichment.tab

  for b in ${params.timecoursebams};
    do samtools idxstats \$b |perl -lane 'print join("\t","'\$b'",\$_)' |grep KmetStat_ |perl -pi -e 's/^.+(..)\\.(R[12])\\.(\\S+?)\\..+bam\\s+KmetStat_(\\S+)\\s+\\d+\\s+(\\d+)\\s+\\d+/\$1\\t\$2\\t\$3\\t\$4\\t\$5/' 2>e.e >>KmetStat_Enrichment.tab
  done

  R --no-save --no-site-file --no-init-file --no-restore --silent --slave <${params.codeRdir}/sortingPaper_plotKmetStatEnrichment.R  ||true
  """
  }

process drawFigure4{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 8
  memory '12g'

  time { 4.hour }
  tag  { "Drawing Figure 4" }

  module 'R/3.5.2'

  publishDir params.outdirImages, mode: 'copy', overwrite: true

  input:
  file(peakTab)  from k4me3TableForRb
  file(hstssTab) from hstssTableForRb
  file(tssTab)   from tssTableForRb
  file(hMat)     from matrixSuppHM.collect()
  file(s4Done)   from sliceBGfig4b.collect()

  output:
  file "*png"     into figure4PNGs
  file "*pdf"     into figure4PDFs

  script:

  """
  ## Add globals to .Renviron ; (so we can see them inside R!)
  ##
  echo KBPIPEOUTDIR="${params.outdir}/"      >>~/.Renviron
  echo KBPIPEWORKDIR="${params.projectdir}/" >>~/.Renviron

  R -q --no-save --no-site-file --no-init-file --no-restore --silent --slave <${params.codeRdir}/plotFigure4.R ||true
  """
  }

process drawFigureS5{
  scratch '/lscratch/$SLURM_JOBID'
  clusterOptions ' --gres=lscratch:40'
  echo true
  cpus 1
  memory '2g'

  time { 0.1.hour }
  tag  { "Drawing Figure S5" }

  publishDir params.outdirImages, mode: 'copy', overwrite: true

  input:
  file(initPNG)  from hmOutPNGHMx5.collect()

  output:
  file "Figure_S5*png"    into figureS5PNGs

  script:

  """
  convert -append histoneMods.set1.HS.png histoneMods.set2.HS.png Figure_S5a.png
  convert -append histoneMods.set3.HS.png histoneMods.set4.HS.png Figure_S5b.png
  """
  }

// //// UNCOMMENT THIS PROCESS IF YOU WANT TO MODIFY PIPELINE AND ASSURE IT DOES  NOT COMPLETE !!
// //// VERY USEFUL FOR TWEAKING THINGS & DEBUGGING
// process failNoGood {
//     scratch '/lscratch/$SLURM_JOBID'
//     clusterOptions ' --gres=lscratch:4'
//     echo true
//     cpus 1
//     memory '1g'
//
//     time { 1.hour }
//
//     input:
//     file failta  from allCounts
//   	file failtd  from tssTableForR
//   	file failte  from hstssTableForR
//     file failte2 from k4me3TableForR
//     file failtg  from figure3PNGs
//     file failth  from figure4PNGs
//     file failti  from figureS5PNGs
//
//
//     script:
//     """
//     ## Clean up .Renviron
//     grep -vP "KBPIPE(WORK|OUT)DIR" ~/.Renviron  >~/.Renviron
//
//     bedtools ls proon
//     """
//   }
