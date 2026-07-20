# bamCompare runs per sample against that sample's assigned control (control_for).
# Only samples that have a control can produce a .peaks.bw, enforced by the
# wildcard constraint; when no sample has a control the constraint never matches
# and the rule is effectively inert.
rule bamcompare:
    input:
        sample_bam  = OUT + "/bam/{sample}.bam",
        sample_bai  = OUT + "/bam/{sample}.bam.bai",
        control_bam = lambda wc: OUT + f"/bam/{control_for(wc)}.bam",
        control_bai = lambda wc: OUT + f"/bam/{control_for(wc)}.bam.bai",
    output:
        OUT + "/bigWig/{sample}.peaks.bw"
    wildcard_constraints:
        sample = CONTROL_TREATMENT_CONSTRAINT,
    params:
        extra                = config["bamcompare"].get("extra", ""),
        bin_size             = config["bamcompare"]["bin_size"],
        operation            = config["bamcompare"]["operation"],
        scale_factors_method = config["bamcompare"]["scale_factors_method"],
        n                    = config["bamcompare"]["n"],
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["bamcompare"]["mem_mb"],
        runtime         = config["resources"]["bamcompare"]["runtime"],
    log:
        OUT + "/logs/bamcompare/{sample}.log"
    shell:
        """
        bamCompare \
          -b1 {input.sample_bam} -b2 {input.control_bam} \
          -o {output} \
          --binSize {params.bin_size} --operation {params.operation} \
          --scaleFactorsMethod {params.scale_factors_method} -n {params.n} \
          -p {threads} {params.extra} 2>{log}
        """
