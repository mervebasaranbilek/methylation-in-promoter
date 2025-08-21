import os
grouphome = os.environ.get("GROUPHOME", "")
samples_file = os.path.join(grouphome, "Valafar Lab/TSS_Project/Known_TSS/samples.csv")
SAMPLES = [
    i for i in open(samples_file).read().splitlines() if len(i) > 0
]

rule all:
    input:
        expand("results/{sample}_methylation_in_tss.csv", sample=SAMPLES)

rule methylation_overlap:
    input:
        tss="Valafar Lab/TSS_Project/Known_TSS/{sample}.csv",
        meth="Valafar Lab/TSS_Project/Methylation/{sample}_methylation.csv"
    output:
        "results/{sample}_methylation_in_tss.csv"
    shell:
        """
        python3 scripts/tss_methyl.py \
          -t {input.tss} \
          -m {input.meth} \
          -o {output}
        """
