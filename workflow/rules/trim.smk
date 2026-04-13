rule trim_se:
    wildcard_constraints:
        sample = "|".join(sorted(SE_SAMPLES)) if SE_SAMPLES else "(?!)",
    input:
        r1  = lambda wc: get_r1(wc.sample),
    output:
        r1           = OUT + "/trimmed/{sample}.R1.fastq.gz",
        total_frags  = OUT + "/stats/{sample}.total_frags.txt",
        subsampled_frags = OUT + "/stats/{sample}.subsampled_frags.txt",
    params:
        max_frags = config["bbduk"].get("max_frags") or "None",
        adapters  = config["bbduk"].get("adapters") or "adapters",
        mem_gb    = config["bbduk"].get("mem_gb", 8),
        extra     = config["bbduk"].get("extra", ""),
        k         = config["bbduk"].get("k", 21),
        mink      = config["bbduk"].get("mink", 11),
        ktrim     = config["bbduk"].get("ktrim", "r"),
        qtrim     = config["bbduk"].get("qtrim", "r"),
        trimq     = config["bbduk"].get("trimq", 6),
        maq       = config["bbduk"].get("maq", 10),
        ow        = config["bbduk"].get("ow", "t"),
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["trim_align"]["mem_mb"],
        runtime         = config["resources"]["trim_align"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        subsample = OUT + "/logs/bbduk/{sample}.subsample.log",
        trim      = OUT + "/logs/bbduk/{sample}.trim.log",
    shell:
        """
        R1_FILES=({input.r1})

        TOTAL_LINES=$(
          for fq in "${{R1_FILES[@]}}"; do
            if gzip -t "$fq" 2>/dev/null; then gzip -dc "$fq"; else cat "$fq"; fi
          done | wc -l
        )
        TOTAL_FRAGS=$((TOTAL_LINES / 4))
        echo "$TOTAL_FRAGS" > {output.total_frags}

        MAXFRAGS={params.max_frags}
        if [ "$MAXFRAGS" = "None" ] || [ -z "$MAXFRAGS" ]; then
          SUB_FRAC=1
        else
          SUB_FRAC=$(echo "$MAXFRAGS / $TOTAL_FRAGS" | bc -l)
          if (( $(echo "$SUB_FRAC > 1" | bc -l) )); then SUB_FRAC=1; fi
        fi

        (
          for fq in "${{R1_FILES[@]}}"; do
            if gzip -t "$fq" 2>/dev/null; then gzip -dc "$fq"; else cat "$fq"; fi
          done
        ) \
        | reformat.sh \
          in=stdin.fq out=stdout.fq \
          int=f \
          -Xmx{params.mem_gb}g \
          samplerate=$SUB_FRAC \
          sampleseed=9 \
          2>{log.subsample} \
        | bbduk.sh \
          -Xmx{params.mem_gb}g \
          threads={threads} \
          int=f \
          in=stdin.fq out={output.r1} \
          ref={params.adapters} \
          k={params.k} mink={params.mink} ktrim={params.ktrim} qtrim={params.qtrim} \
          trimq={params.trimq} maq={params.maq} ow={params.ow} \
          {params.extra} \
          2>{log.trim}

        grep "^Output:" {log.subsample} \
          | awk '{{print int($2)}}' > {output.subsampled_frags}
        """


rule trim_pe:
    wildcard_constraints:
        sample = "|".join(sorted(PE_SAMPLES)) if PE_SAMPLES else "(?!)",
    input:
        r1 = lambda wc: get_r1(wc.sample),
        r2 = lambda wc: config["samples"][wc.sample]["r2"],
    output:
        r1               = OUT + "/trimmed/{sample}.R1.fastq.gz",
        r2               = OUT + "/trimmed/{sample}.R2.fastq.gz",
        total_frags      = OUT + "/stats/{sample}.total_frags.txt",
        subsampled_frags = OUT + "/stats/{sample}.subsampled_frags.txt",
    params:
        max_frags = config["bbduk"].get("max_frags") or "None",
        adapters  = config["bbduk"].get("adapters") or "adapters",
        mem_gb    = config["bbduk"].get("mem_gb", 8),
        extra     = config["bbduk"].get("extra", ""),
        k         = config["bbduk"].get("k", 21),
        mink      = config["bbduk"].get("mink", 11),
        ktrim     = config["bbduk"].get("ktrim", "r"),
        qtrim     = config["bbduk"].get("qtrim", "r"),
        trimq     = config["bbduk"].get("trimq", 6),
        maq       = config["bbduk"].get("maq", 10),
        ow        = config["bbduk"].get("ow", "t"),
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["trim_align"]["mem_mb"],
        runtime         = config["resources"]["trim_align"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        subsample = OUT + "/logs/bbduk/{sample}.subsample.log",
        trim      = OUT + "/logs/bbduk/{sample}.trim.log",
    shell:
        """
        R1_FILES=({input.r1})
        R1="${{R1_FILES[0]}}"
        R2={input.r2}

        TOTAL_FRAGS=$(
          if gzip -t "$R1" 2>/dev/null; then gzip -dc "$R1"; else cat "$R1"; fi | wc -l
        )
        TOTAL_FRAGS=$((TOTAL_FRAGS / 4))
        echo "$TOTAL_FRAGS" > {output.total_frags}

        MAXFRAGS={params.max_frags}
        if [ "$MAXFRAGS" = "None" ] || [ -z "$MAXFRAGS" ]; then
          SUB_FRAC=1
        else
          SUB_FRAC=$(echo "$MAXFRAGS / $TOTAL_FRAGS" | bc -l)
          if (( $(echo "$SUB_FRAC > 1" | bc -l) )); then SUB_FRAC=1; fi
        fi

        reformat.sh \
          in1=$R1 in2=$R2 \
          out=stdout.fq \
          -Xmx{params.mem_gb}g \
          int=t \
          samplerate=$SUB_FRAC \
          sampleseed=9 \
          2>{log.subsample} \
        | bbduk.sh \
          -Xmx{params.mem_gb}g \
          threads={threads} \
          int=t \
          in=stdin.fq out1={output.r1} out2={output.r2} \
          ref={params.adapters} \
          k={params.k} mink={params.mink} ktrim={params.ktrim} tbo tpe qtrim={params.qtrim} trimq={params.trimq} maq={params.maq} ow={params.ow} \
          {params.extra} \
          2>{log.trim}

        grep "^Output:" {log.subsample} \
          | awk '{{print int($2/2)}}' > {output.subsampled_frags}
        """
