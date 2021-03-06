configfile: "config/config.yaml"
configfile: config['samplesList']
configfile: config['pairList']

def loadBarcode2BaitKit(meta):
	dic = dict()
	with open(meta,'r') as fi:
		for line in fi:
			fs = line.strip().split('\t')
			bc = fs[6]
			kit = fs[19]
			dic[bc] = kit
	return dic
config['barcode2kit'] = loadBarcode2BaitKit(config['bamMeta'])

rule correctDepth:
  input: 
  	expand("results/ichorCNA/{tumor}/{tumor}.cna.seg", tumor=config["pairings"]),
  	expand("results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt", tumor=config["pairings"]),
  	expand("results/readDepth/{samples}.bin{binSize}.wig", samples=config["samples"], binSize=str(config["binSize"]))

rule read_counter:
	input:
		done = 'data/gdc/{samples}.done',
	output:
		"results/readDepth/{samples}.bin{binSize}.wig"		
	params:
		readCounter=config["readCounterScript"],
		binSize=config["binSize"],
		qual="20",
		chrs=config["chrs"],
		bam=lambda wildcards: config["samples"][wildcards.samples]
	resources:
		mem=4
	log:
		"logs/readDepth/{samples}.bin{binSize}.log"
	shell:
		"{params.readCounter} {params.bam} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule ichorCNA:
	input:
		tum="results/readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig",
		norm=lambda wildcards: "results/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
	output:
		corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
		#param="results/ichorCNA/{tumor}/{tumor}.params.txt",
		cna="results/ichorCNA/{tumor}/{tumor}.cna.seg",
		#segTxt="results/ichorCNA/{tumor}/{tumor}.seg.txt",
		#seg="results/ichorCNA/{tumor}/{tumor}.seg",
		#rdata="results/ichorCNA/{tumor}/{tumor}.RData",
		outDir="results/ichorCNA/{tumor}/",
	params:
		rscript=config["ichorCNA_rscript"],
		id="{tumor}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
		gcwig=config["ichorCNA_gcWig"],
		mapwig=config["ichorCNA_mapWig"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		centromere=config["ichorCNA_centromere"],
		exons=lambda wildcards: config["ichorCNA_exons"][config['barcode2kit'][wildcards.tumor]],
		txnE=config["ichorCNA_txnE"],
		txnStrength=config["ichorCNA_txnStrength"],
		fracReadsChrYMale="0.001",
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"]
	resources:
		mem=4
	log:
		"logs/ichorCNA/{tumor}.log"	
	shell:
		"Rscript {params.rscript} --id {params.id} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --NORMWIG {input.norm} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --fracReadsInChrYForMale {params.fracReadsChrYMale} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {output.outDir} > {log} 2> {log}"
