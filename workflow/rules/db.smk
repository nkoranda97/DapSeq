rule update_db:
    input:
        qc_summary = OUT + "/stats/qc_summary.tsv",
    output:
        flag = OUT + "/stats/db_updated.flag",
    params:
        db_path          = os.path.join(os.path.dirname(workflow.basedir), "pipeline_db.csv"),
        output_dir       = OUT,
        samples_cfg      = config["samples"],
        input_control    = config.get("input_control"),
        genome_ref       = config["genome_ref"],
        genome_size      = config["genome_size"],
        threads          = config["threads"],
        mapq             = config["samtools"]["mapq"],
        max_frags        = config["bbduk"].get("max_frags"),
        macs3_format     = config["macs3"]["format"],
        macs3_min_foldch = config["macs3"]["min_foldch"],
        meme_nmotifs     = config["meme"]["nmotifs"],
        meme_minw        = config["meme"]["minw"],
        meme_maxw        = config["meme"]["maxw"],
        meme_maxpeaks    = config["meme"]["maxpeaks"],
        fimo_thresh      = config["fimo"]["thresh"],
        treatment_samples = TREATMENT_SAMPLES,
    resources:
        mem_mb          = config["resources"]["update_db"]["mem_mb"],
        runtime         = config["resources"]["update_db"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/update_db.log"
    script:
        "../scripts/update_db.py"
