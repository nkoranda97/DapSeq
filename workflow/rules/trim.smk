rule trim_se:
    wildcard_constraints:
        sample = "|".join(sorted(SE_SAMPLES)) if SE_SAMPLES else "(?!)",
    input:
        r1  = lambda wc: get_r1(wc.sample),
    output:
        r1           = OUT + "/trimmed/{sample}.R1.fastq.gz",
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
        # bbmap 39.81's parallel FASTQ reader (FastqStreamer) has a thread-safety
        # bug that deadlocks at high thread counts (>=12 observed) with spurious
        # "missing plus" assertions on valid reads. Cap trimming at 8 threads;
        # alignment keeps config["threads"].
        min(config["threads"], 8)
    resources:
        mem_mb          = config["resources"]["trim_align"]["mem_mb"],
        runtime         = config["resources"]["trim_align"]["runtime"],
    log:
        subsample = OUT + "/logs/bbduk/{sample}.subsample.log",
        trim      = OUT + "/logs/bbduk/{sample}.trim.log",
    shell:
        """
        set -euo pipefail
        R1_FILES=({input.r1})

        TOTAL_LINES=$(
          for fq in "${{R1_FILES[@]}}"; do
            if gzip -t "$fq" 2>/dev/null; then gzip -dc "$fq"; else cat "$fq"; fi
          done | wc -l
        )
        TOTAL_FRAGS=$((TOTAL_LINES / 4))

        MAXFRAGS={params.max_frags}
        if [ "$MAXFRAGS" = "None" ] || [ -z "$MAXFRAGS" ]; then
          printf "Input:\t%s reads\nOutput:\t%s reads (100.00%%)\n" "$TOTAL_FRAGS" "$TOTAL_FRAGS" > {log.subsample}

          (
            for fq in "${{R1_FILES[@]}}"; do
              if gzip -t "$fq" 2>/dev/null; then gzip -dc "$fq"; else cat "$fq"; fi
            done
          ) \
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
        else
          SUB_FRAC=$(echo "$MAXFRAGS / $TOTAL_FRAGS" | bc -l)
          if (( $(echo "$SUB_FRAC > 1" | bc -l) )); then SUB_FRAC=1; fi

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
        fi

        """


rule trim_pe:
    wildcard_constraints:
        sample = "|".join(sorted(PE_SAMPLES)) if PE_SAMPLES else "(?!)",
    input:
        r1 = lambda wc: get_r1(wc.sample),
        r2 = lambda wc: get_r2(wc.sample),
    output:
        r1               = OUT + "/trimmed/{sample}.R1.fastq.gz",
        r2               = OUT + "/trimmed/{sample}.R2.fastq.gz",
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
        join_dir  = OUT + "/temp/{sample}_lanes",
    threads:
        # bbmap 39.81's parallel FASTQ reader (FastqStreamer) has a thread-safety
        # bug that deadlocks at high thread counts (>=12 observed) with spurious
        # "missing plus" assertions on valid reads. Cap trimming at 8 threads;
        # alignment keeps config["threads"].
        min(config["threads"], 8)
    resources:
        mem_mb          = config["resources"]["trim_align"]["mem_mb"],
        runtime         = config["resources"]["trim_align"]["runtime"],
    log:
        subsample = OUT + "/logs/bbduk/{sample}.subsample.log",
        trim      = OUT + "/logs/bbduk/{sample}.trim.log",
    shell:
        """
        set -euo pipefail
        shopt -s inherit_errexit  # a failed cat/gzip inside $(join_lanes ...) must abort the job
        R1_FILES=({input.r1})
        R2_FILES=({input.r2})
        JOIN_DIR="{params.join_dir}"
        trap 'rm -rf "$JOIN_DIR"' EXIT

        # Join a mate's lanes into one file so bbduk sees a single in1/in2.
        # One lane is used as-is. Same-compression lanes are byte-concatenated
        # (concatenated gzip members are a valid gzip stream); mixed lanes are
        # decompressed. The extension matches the content, which BBTools reads.
        join_lanes() {{
          local out=$1; shift
          if [ $# -eq 1 ]; then echo "$1"; return; fi
          local n_gz=0 fq
          for fq in "$@"; do
            # gzip magic bytes 1f 8b; avoids decompressing each lane just to detect it
            if [ "$(head -c2 "$fq" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]; then n_gz=$((n_gz + 1)); fi
          done
          mkdir -p "$JOIN_DIR"
          if [ "$n_gz" -eq "$#" ]; then
            cat "$@" > "$out.fastq.gz"; echo "$out.fastq.gz"
          elif [ "$n_gz" -eq 0 ]; then
            cat "$@" > "$out.fastq"; echo "$out.fastq"
          else
            for fq in "$@"; do
              if gzip -t "$fq" 2>/dev/null; then gzip -dc "$fq"; else cat "$fq"; fi
            done > "$out.fastq"
            echo "$out.fastq"
          fi
        }}
        R1=$(join_lanes "$JOIN_DIR/R1" "${{R1_FILES[@]}}")
        R2=$(join_lanes "$JOIN_DIR/R2" "${{R2_FILES[@]}}")

        TOTAL_FRAGS=$(
          if gzip -t "$R1" 2>/dev/null; then gzip -dc "$R1"; else cat "$R1"; fi | wc -l
        )
        TOTAL_FRAGS=$((TOTAL_FRAGS / 4))

        MAXFRAGS={params.max_frags}
        if [ "$MAXFRAGS" = "None" ] || [ -z "$MAXFRAGS" ]; then
          TOTAL_READS=$((TOTAL_FRAGS * 2))
          printf "Input:\t%s reads\nOutput:\t%s reads (100.00%%)\n" "$TOTAL_READS" "$TOTAL_READS" > {log.subsample}

          bbduk.sh \
            -Xmx{params.mem_gb}g \
            threads={threads} \
            int=t \
            in1="$R1" in2="$R2" out1={output.r1} out2={output.r2} \
            ref={params.adapters} \
            k={params.k} mink={params.mink} ktrim={params.ktrim} tbo tpe qtrim={params.qtrim} trimq={params.trimq} maq={params.maq} ow={params.ow} \
            {params.extra} \
            2>{log.trim}
        else
          SUB_FRAC=$(echo "$MAXFRAGS / $TOTAL_FRAGS" | bc -l)
          if (( $(echo "$SUB_FRAC > 1" | bc -l) )); then SUB_FRAC=1; fi

          reformat.sh \
            in1="$R1" in2="$R2" \
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
        fi

        """
