configfile: "config/config.yaml"
print(config)
configfile: config['samplesList']
configfile: config['pairList']

include: "ichorCNA.snakefile"
include: "getAlleleCounts.snakefile"
include: "GDCDownload.snakefile"
import os.path

CHRS = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
CLUST = {1:[1], 2:[1,2], 3:[1,2,3], 4:[1,2,3,4], 5:[1,2,3,4,5], 6:[1,2,3,4,5,6], 7:[1,2,3,4,5,6,7], 8:[1,2,3,4,5,6,7,8], 9:[1,2,3,4,5,6,7,8,9], 10:[1,2,3,4,5,6,7,8,9,10]}
PLOIDY = {2:[2], 3:[2,3], 4:[2,3,4]}


# rule all:
# 	input: 
# 		expand("results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
# 		#expand("results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
# 		expand("results/titan/hmm/{tumor}/optimalClusterSolution.txt",tumor=config["pairings"].keys())
rule all:
	input: 
		expand("results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config['run'], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		#expand("results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		expand("results/titan/hmm/{tumor}/optimalClusterSolution.txt",tumor=config["run"])			
rule makeOutDir:
	output:
		"results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/"
	shell:
		"mkdir -p {output}"
		
rule runTitanCNA:
	input:
		alleleCounts="results/titan/tumCounts/{tumor}.tumCounts.txt",
		corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt"		
	output:
		outRoot="results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}/",
		titan="results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
		param="results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
		segTxt="results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt",
		seg="results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.seg"
	params:
		titanRscript=config["TitanCNA_rscript"],
		libdir=config["TitanCNA_libdir"],
		numCores=config["TitanCNA_numCores"],
		normal=config["TitanCNA_normalInit"],
		chrs=config["TitanCNA_chrs"],
		estimatePloidy=config["TitanCNA_estimatePloidy"],
		estimateClonality=config["TitanCNA_estimateClonality"],
		estimateNormal=config["TitanCNA_estimateNormal"],
		centromere=config["TitanCNA_centromere"],
		alphaK=config["TitanCNA_alphaK"],
		alphaKHigh=config["TitanCNA_alphaKHigh"],
		#alphaR=config["TitanCNA_alphaR"],
		#alleleModel=config["TitanCNA_alleleModel"],
		cytobandFile=config["TitanCNA_cytobandFile"],
		txnExpLen=config["TitanCNA_txnExpLen"],
		plotYlim=config["TitanCNA_plotYlim"],
		genome=config["TitanCNA_genomeStyle"]
	log:
		"logs/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.log"
	shell:
		"Rscript {params.titanRscript} --genomeStyle {params.genome} --hetFile {input.alleleCounts} --cnFile {input.corrDepth} --outFile {output.titan} --outSeg {output.segTxt} --outParam {output.param} --outIGV {output.seg} --outPlotDir {output.outRoot} --libdir {params.libdir} --id {wildcards.tumor} --numClusters {wildcards.clustNum} --numCores {params.numCores} --normal_0 {params.normal} --ploidy_0 {wildcards.ploidy} --chrs \"{params.chrs}\" --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateClonality {params.estimateClonality}  --centromere {params.centromere} --alphaK {params.alphaK} --alphaKHigh {params.alphaKHigh} --txnExpLen {params.txnExpLen} --plotYlim \"{params.plotYlim}\" --cytobandFile {params.cytobandFile} > {log} 2> {log}"
	
#--alleleModel {params.alleleModel} --alphaR {params.alphaR}
	
				
rule selectSolution:
	input:
		ploidyDirs= lambda wildcards:expand("results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]],tumor=wildcards.tumor),
		resultFiles=lambda wildcards: expand("results/titan/hmm/{tumor}/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]],tumor=wildcards.tumor)
	output:
		"results/titan/hmm/{tumor}/optimalClusterSolution.txt"
	params:
		solutionRscript=config["TitanCNA_selectSolutionRscript"],
		threshold=config["TitanCNA_solutionThreshold"]
	log:
		"logs/titan/{tumor}_selectSolution.log"
	shell:
		"""
		if [ -d results/titan/hmm/{wildcards.tumor}/titanCNA_ploidy3/ ]; then
			ploidyRun3=results/titan/hmm/{wildcards.tumor}/titanCNA_ploidy3/
		else
			ploidyRun3=NULL
		fi
		if [ -d results/titan/hmm/{wildcards.tumor}/titanCNA_ploidy4/ ]; then
			ploidyRun4=results/titan/hmm/{wildcards.tumor}/titanCNA_ploidy4/
		else
			ploidyRun4=NULL
		fi
		Rscript {params.solutionRscript} --ploidyRun2 {input.ploidyDirs[0]} --ploidyRun3 $ploidyRun3 --ploidyRun4 $ploidyRun4 --threshold {params.threshold} --outFile {output} > {log} 2> {log}
		"""
#	run:
#		if "results/titan/hmm/{tumor}/titanCNA_ploidy3" in input:
#			ploidyRun3 = input[1]
#		else:
#			ploidyRun3 = "NULL"
#		if "results/titan/hmm/{tumor}/titanCNA_ploidy4" in input:
#			ploidyRun4 = input[2]
#		else:
#			ploidyRun4 = "NULL"	
#		os.system("Rscript params.solutionRscript --ploidyRun2 input[0] --ploidyRun3 ploidyRun3 --ploidyRun4 ploidyRun4 --threshold params.threshold --outFile output")
	
		
