configfile: "config/config.yaml"
configfile: config['samplesList']
configfile: config['pairList']

def loadBamManifest(meta):
    dic = dict()
    with open(meta,'r') as fi:
        for line in fi:
            fs = line.strip().split()
            bc = fs[6]
            uuid = fs[0]
            dic[bc] = uuid
    return dic
config['barcode2uuid'] = loadBamManifest(config['bamMeta'])

rule download_bam:
    output:
        'data/gdc/{sample}.done'
    params:
         uuid = lambda wildcards: config['barcode2uuid'][wildcards.sample]
    threads:
        4
    shell:
    # It will timeout after 15m maximum run for gdc-client
    # if more than 15m, stop at gdc-client and will not touch
        "timeout 5h gdc-client download -t {config[token]} "
        "-d data/ -n {threads} {params.uuid} "
        "&& touch {output} && bash fixindex.sh data/{params.uuid}"
