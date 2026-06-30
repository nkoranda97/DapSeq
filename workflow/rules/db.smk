rule update_db:
    input:
        report = OUT + "/stats/report.csv",
    output:
        flag = OUT + "/stats/db_updated.flag",
    params:
        db_path          = os.path.join(os.path.dirname(workflow.basedir), "pipeline_db.db"),
        output_dir       = OUT,
        samples_cfg      = config["samples"],
        control          = config.get("control"),
        genome_ref       = config["genome_ref"],
        genome_size      = config["genome_size"],
        threads          = config["threads"],
        mapq             = config["samtools"]["mapq"],
        max_frags        = config["bbduk"].get("max_frags"),
        macs3_format     = config["macs3"]["format"],
        macs3_foldch_levels  = config["macs3"]["foldch_levels"],
        macs3_meme_foldch_level = MEME_FOLD_IDX,
        meme_nmotifs     = config["meme"]["nmotifs"],
        meme_minw        = config["meme"]["minw"],
        meme_maxw        = config["meme"]["maxw"],
        meme_maxpeaks    = config["meme"]["maxpeaks"],
        fimo_thresh      = config["fimo"]["thresh"],
        treatment_samples = TREATMENT_SAMPLES,
        author           = config.get("author", ""),
        gene_annotation  = config.get("gene_annotation") or "",
    resources:
        mem_mb          = config["resources"]["update_db"]["mem_mb"],
        runtime         = config["resources"]["update_db"]["runtime"],
    log:
        OUT + "/logs/update_db.log"
    script:
        "../scripts/update_db.py"
