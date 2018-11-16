configfile: "config/config.yaml"
configfile: config['samplesList']
configfile: config['pairList']

CHRS = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]

rule tumCounts:
	input: 
		expand("results/titan/tumCounts/{tumor}/{tumor}.tumCounts.{chr}.txt", tumor=config["pairings"], chr=CHRS),		
		expand("results/titan/tumCounts/{tumor}.tumCounts.txt", tumor=config["pairings"])
	
rule getHETsites:
	input:
		done = lambda wildcards: 'data/gdc/{s}.done'.format(s=config["pairings"][wildcards.tumor])
		
	output:
		"results/titan/hetPosns/{tumor}/{tumor}.{chr}.vcf"
	params:
		refFasta=config["refFasta"],
		snpDB=config["snpDB"],
		bam=lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
	log:
		"logs/titan/hetPosns/{tumor}/{tumor}.{chr}.log"
	shell:
		"samtools mpileup -uv -I -f {params.refFasta} -r {wildcards.chr} -l {params.snpDB} {params.bam} | bcftools call -v -c - | grep -e '0/1' -e '#' > {output} 2> {log}"


rule getAlleleCountsByChr:
	input:
		done = 'data/gdc/{tumor}.done',
		hetSites="results/titan/hetPosns/{tumor}/{tumor}.{chr}.vcf",
	output:
		"results/titan/tumCounts/{tumor}/{tumor}.tumCounts.{chr}.txt"
	params:
		countScript=config["pyCountScript"],
		#pyEnv=config["pyEnv"],
		#refFasta=config["refFasta"],
		mapQ=config["map_quality"],
		baseQ=config["base_quality"],
		vcfQ=config["vcf_quality"],
		tumBam=lambda wildcards: config["samples"][wildcards.tumor]
	log:
		"logs/titan/tumCounts/{tumor}/{tumor}.{chr}.log"
	shell:
		"python {params.countScript} {wildcards.chr} {input.hetSites} {params.tumBam} {params.baseQ} {params.mapQ} {params.vcfQ} > {output} 2> {log}"

rule catAlleleCountFiles:
	input:
		expand("results/titan/tumCounts/{{tumor}}/{{tumor}}.tumCounts.{chr}.txt", chr=CHRS)
	output:
		"results/titan/tumCounts/{tumor}.tumCounts.txt"
	log:
		"logs/titan/tumCounts/{tumor}/{tumor}.cat.log"
	shell:
		"cat {input} | grep -v Chr > {output} 2> {log}"






