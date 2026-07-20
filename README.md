# DAP-seq Pipeline

Snakemake implementation of the JGI DAP-seq analysis pipeline. Runs on SLURM clusters or local workstations via Apptainer.

## Setup

Install [Pixi](https://pixi.sh):

```sh
curl -fsSL https://pixi.sh/install.sh | sh
```

## Running

### HPC (SLURM cluster)

```sh
module load apptainer
pixi shell --frozen # activates the environment
snakemake --profile profiles/slurm --configfile /path/to/your/config.yaml
```

Set `slurm_partition` and `slurm_account` in the slurm profile

### Local workstation

Build the Apptainer SIF once (requires the pixi environment):

```sh
pixi run apptainer build apptainer_build/dapseq.sif apptainer_build/dapseq.def
```

Then run using the local profile (no SLURM, 4 parallel jobs by default):

```sh
pixi run snakemake --profile profiles/local --configfile /path/to/your/config.yaml
```

`slurm_partition` and `slurm_account` can be omitted or left as `null` in your config when using the local profile. Override parallelism with `--jobs N`.

### Shared installation (multi-user)

The pipeline can be installed once on a shared filesystem and used by multiple users simultaneously. Each user:

1. Copies `config.yaml` from the repo root to their own location and fills in their samples, paths, and parameters.
2. Runs from the shared repo root, passing their own config with `--configfile`:


Concurrent runs are safe as long as each user sets `output_dir` to a unique location. Snakemake's locks are keyed per output directory, so simultaneous runs do not interfere with each other.

The shared database (`pipeline_db.db` in the repo root) accumulates results from all users' runs. Set `db_path` in your config to redirect to a private database if needed.

## Configuration

Copy `config.yaml` from the repo root and fill in the required fields. Defaults for all tool parameters are in `config/config.yaml`.

### Required fields

```yaml
author: "Your Name"  # optional

samples:
  sample_name:
    r1: /path/to/sample.R1.fastq.gz
    r2: /path/to/sample.R2.fastq.gz  # null for single-end
    control: input_DNA               # name of another sample; null for no background
    experiment_date: "2026-07-14"    # optional, per-sample (free-form)
    gdna_batch: "batch-07"           # optional, per-sample (free-form)
  input_DNA:
    r1: /path/to/control.R1.fastq.gz
    r2: null
    control: null                    # a control declares no control of its own

output_dir: /path/to/output/
genome_ref: /path/to/genome.fa
genome_size: "3000000000"
gene_annotation: /path/to/annotation.gtf  # null to skip HOMER annotation

slurm_partition: "caslake"
slurm_account: "pi-yourlab"
```

#### Running multiple experiments in one run

Control is assigned **per sample** via each sample's `control:` key, so you can run
several independent experiments in a single pipeline invocation — just point each
sample at its own control:

```yaml
samples:
  TF_A:    { r1: ..., control: input_A }   # experiment A
  TF_B:    { r1: ..., control: input_B }   # experiment B
  input_A: { r1: ..., control: null }
  input_B: { r1: ..., control: null }
```

Each treatment is peak-called against its assigned control; each control is
peak-called against itself and shown as a flagged report row. There is exactly one
control per sample (the `control:` value is a single sample name, not a list). The
previous multi-control **merge** behavior (a top-level `control:` list merged into
one `merged_control`) is no longer supported.

### Notable options

| Option | Default | Notes |
|---|---|---|
| `samples.<name>.control` | `null` | Per-sample: name of another sample to use as this sample's MACS3 `-c` background and bamCompare `-b2`. Exactly one control per sample; `null` peak-calls the sample with no background. A sample named as a control by others is peak-called **twice**: once with no input control (QC self peak-call), and once **against itself** (`-t control -c control`, expected to yield few or no peaks). The against-itself run flows through the full treatment pipeline (filtering, MEME/FIMO motifs, stats) and appears as a flagged row in the report. Assign different controls to different samples to run several experiments in one pipeline invocation. **The old top-level `control:` field has been removed** — a config that still sets it fails fast with a migration message |
| `samples.<name>.experiment_date` | `null` | Optional per-sample date (free-form); shown as a report column and stored per-sample in the results database |
| `samples.<name>.gdna_batch` | `null` | Optional per-sample gDNA batch label (free-form); shown as a report column and stored per-sample in the results database |
| `aligner` | `bowtie2` | `bowtie2` or `bwa_mem2` |
| `bbduk.k` | `21` | Kmer length for adapter/contaminant trimming (program default: `31`) |
| `bbduk.mink` | `11` | Enables short-kmer matching at read tips for adapter trimming (program default: `0`, disabled) |
| `bbduk.ktrim` | `r` | Trim direction for kmer-matched adapters (program default: `f`, no trimming) |
| `bbduk.qtrim` | `r` | Trim direction for quality trimming (program default: `f`, no trimming) |
| `bbduk.maq` | `10` | Minimum average read quality after trimming (program default: `0`, no filter) |
| `samtools.mapq` | `30` | MAPQ filter applied after alignment (program default: `0`, no filter) |
| `bamcoverage.normalize_using` | `RPGC` | Coverage normalization method (program default: `None`) |
| `bamcoverage.bin_size` | `1` | bigwig bin size, in bp (program default: `50`) |
| `bamcoverage.extend_reads` | `true` | Extend reads to fragment size; PE only (program default: `f`) |
| `bamcoverage.max_fragment_length` | `600` | Max fragment length for PE read/pair inclusion (program default: `0`, no limit) |
| `bamcoverage.ignore_duplicates` | `true` | Count duplicate reads only once (program default: `f`) |
| `chrom_filter` | `[]` | Chromosomes excluded from peaks before MEME (e.g. `[chrEBV]`) |
| `complexity_filter.enabled` | `false` | Drop low-complexity peak FASTA sequences before MEME |
| `complexity_filter.min_entropy` | `3.0` | Minimum 3-mer Shannon entropy in bits when `complexity_filter.enabled: true`; lower is more permissive, max possible is `6.0` |
| `macs3.min_foldch` | `2.0` | Peak fold-change filter |
| `macs3.format` | `BAMPE` | Set to `BAM` for single-end data |
| `meme.nmotifs` | `2` | Number of motifs to search for (program default: `1`) |
| `meme.maxw` | `32` | Maximum motif width (program default: `50`) |
| `meme.mod` | `anr` | Motif site distribution model (program default: `zoops`) |
| `meme.markov_order` | `0` | Background Markov model order built from input sequences; try `2` to suppress composition-driven repeats |
| `meme.maxsize` | `10000000` | Max total size (bp) of input sequences searched for motifs (`searchsize`; program default: sampling above `100000`) |
| `meme.summit_extend` | `50` | bp around the peak summit used for motif search in "summits" mode |
| `meme.base_colors` | unset | Optional hex color overrides for sequence logos |
| `fimo.thresh` | `1.0e-5` | p-value threshold for FIMO motif scanning |

### Migration note

The old `tandem_filter` block (`enabled`, `k`, `k_max`) has been renamed and replaced by `complexity_filter` (`enabled`, `min_entropy`). Existing user config files that still use `tandem_filter` must be updated; the old key is ignored and the filter will stay off.

### Extra arguments

Any tool accepts an `extra` field for additional flags:

```yaml
macs3:
  extra: "--qvalue 0.01"
```

---

## Default Pipeline Parameters

>Only biologically relevent parameters are shown. Please refere to original documentation for additional parameters

### [bbduk](https://bbmap.org/tools/bbduk)

- `qin=auto` — Input quality offest: `33`, `64`, or `auto` — type: `enum` — program default: `auto`
- `k=21` — Kmer length used for finding contaminants. Contaminants shorter than k will not be found — type: `int` — program default: `31`
- `rcomp=t` — Look for reverse-complements of kmers in addition to forward kmers — type: `bool` — program default: `t`
- `maskmiddle=t` — Tread the middle base of a kmer as a wildcard to increase sensetivity in the presence of an error — type: `bool`|`enum`|`int` (1 for odd-length kmers | 2 for even-length kmers | n to mask n many bp's) — program default: `t`
- `minkmerhits=1` — Reads need at leasy this many matching kmers to be considered as matching the reference — type: `int` — program default: `1`
- `minkmerfraction=0.0` — Reads needs at least this fraction of its total kmers to hit a ref, in order to be considered a match. If this and minkmerhits are set, the greater is used — type: `float` — program default: `0.0`
- `mincovfraction=0.0` — Reads needs at least this fraction of its total bases to be covered by ref kmers to be considered a match. If specified, mcf overrides mkh and mkf — type: `float` — program default: `0.0`
- `hammingdistance=0` — Maximum Hamming distance for ref kmers (subs only). Memory use scales as (3×K)^hdist: E. coli with k=31, hdist=1 uses ~15GB vs 140MB at hdist=0 — type: `int` — program default: `0`
- `qhdist=0` — Hamming distance for query kmers; impacts speed, not memory. Alternative to hdist for memory-constrained systems - mutates read kmers instead of reference kmers — type: `int` — program default: `0`
- `editdistance=0` — Maximum edit distance from ref kmers (subs and indels). Memory use scales as (8×K)^edist. Use qhdist instead of edist/hdist if insufficient memory available — type: `int` — program default: `0`
- `hammingdistance2=0` — Sets hdist for short kmers, when using mink — type: `int` — program default: `0`
- `qhdist2=0` — Sets qhdist for short kmers, when using mink — type: `int` — program default: `0`
- `editdistance2=0` — Sets edist for short kmers, when using mink — type: `int` — program default: `0`
- `forbidn=f` — Forbids matching of read kmers containing N. By default, these will match a reference 'A' if hdist>0 or edist>0, to increase sensitivity — type: `bool` — program default: `f`
- `removeifeitherbad=t` — Paired reads get sent to 'outmatch' if either is match (or either is trimmed shorter than minlen). Set to false to require both — type: `bool` — program default: `t`
- `trimfailures=f` — Instead of discarding failed reads, trim them to 1bp. This makes the statistics a bit odd — type: `bool` — program default: `f`
- `findbestmatch=f` — If multiple matches, associate read with sequence sharing most kmers. Reduces speed — type: `bool` — program default: `f`
- `skipr1=f` — Don't do kmer-based operations on read 1 — type: `bool` — program default: `f`
- `skipr2=f` — Don't do kmer-based operations on read 2 — type: `bool` — program default: `f`
- `ecco=f` — For overlapping paired reads only. Performs error-correction with BBMerge prior to kmer operations — type: `bool` — program default: `f`
- `recalibrate=f` — Recalibrate quality scores. Requires calibration matrices generated by CalcTrueQuality — type: `bool` — program default: `f`
- `sam=<file,file>` — If recalibration is desired, and matrices have not already been generated, BBDuk will create them from the sam file — type: `string` — program default: unset
- `amino=f` — Run in amino acid mode. Some features have not been tested, but kmer-matching works fine. Maximum k is 12 — type: `bool` — program default: `f`
- `ktrim=r` — Trim reads to remove bases matching reference kmers, plus all bases to the left or right. Values: f (don't trim), r (trim to the right), l (trim to the left) — type: `enum` — program default: `f`
- `ktrimtips=0` — Set this to a positive number to perform ktrim on both ends, examining only the outermost X bases — type: `int` — program default: `0`
- `kmask=` — Replace bases matching ref kmers with another symbol. Allows any non-whitespace character, and processes short kmers on both ends if mink is set. 'kmask=lc' will convert masked bases to lowercase — type: `string` — program default: unset
- `maskfullycovered=f` — Only mask bases that are fully covered by kmers — type: `bool` — program default: `f`
- `ksplit=f` — For single-ended reads only. Reads will be split into pairs around the kmer. If the kmer is at the end of the read, it will be trimmed instead. Singletons will go to out, and pairs will go to outm. Do not use ksplit with other operations such as quality-trimming or filtering — type: `bool` — program default: `f`
- `mink=11` — Look for shorter kmers at read tips down to this length, when k-trimming or masking. 0 means disabled. Essential for adapter trimming when adapter remnants are shorter than k. Enabling this will disable maskmiddle — type: `int` — program default: `0`
- `qtrim=r` — Trim read ends to remove bases with quality below trimq. Performed AFTER looking for kmers. Uses Phred algorithm for more accurate trimming than naive approaches. Values: rl (trim both ends), f (neither end), r (right end only), l (left end only), w (sliding window) — type: `enum` — program default: `f`
- `trimq=6` — Regions with average quality BELOW this will be trimmed, if qtrim is set to something other than f. Can be a floating-point number like 7.3 — type: `float` — program default: `6`
- `quantize` — Bin quality scores to reduce file size. quantize=2 will eliminate all odd quality scores, while quantize=0,10,37 will only allow quality scores of 0, 10, or 37 — type: `int` — program default: unset
- `trimclip=f` — Trim soft-clipped bases from sam files — type: `bool` — program default: `f`
- `minlength=10` — Reads shorter than this after trimming will be discarded. Pairs will be discarded if both are shorter — type: `int` — program default: `10`
- `mlf=0` — Reads shorter than this fraction of original length after trimming will be discarded — type: `float` — program default: `0`
- `maxlength=` — Reads longer than this after trimming will be discarded — type: `int` — program default: unset
- `minavgquality=0` — Reads with average quality (after trimming) below this will be discarded — type: `float` — program default: `0`
- `maqb=10` — If positive, calculate maq from this many initial bases — type: `int` — program default: `0`
- `minbasequality=0` — Reads with any base below this quality (after trimming) will be discarded — type: `int` — program default: `0`
- `maxns=-1` — If non-negative, reads with more Ns than this (after trimming) will be discarded — type: `int` — program default: `-1`
- `mcb=0` — Discard reads without at least this many consecutive called bases — type: `int` — program default: `0`
- `ottm=f` — Output reads trimmed to shorter than minlength to outm rather than discarding — type: `bool` — program default: `f`
- `tp=0` — Trim this much extra around matching kmers — type: `int` — program default: `0`
- `tbo=f` — Trim adapters based on where paired reads overlap. Uses BBMerge internally and does not require known adapter sequences — type: `bool` — program default: `f`
- `strictoverlap=t` — Adjust sensitivity for trimbyoverlap mode — type: `bool` — program default: `t`
- `minoverlap=14` — Require this many bases of overlap for detection — type: `int` — program default: `14`
- `mininsert=40` — Require insert size of at least this for overlap. Should be reduced to 16 for small RNA sequencing — type: `int` — program default: `40`
- `tpe=f` — When kmer right-trimming, trim both reads to the minimum length of either. Useful when adapter kmer detected in only one read of a pair — type: `bool` — program default: `f`
- `forcetrimleft=0` — If positive, trim bases to the left of this position (exclusive, 0-based) — type: `int` — program default: `0`
- `forcetrimright=0` — If positive, trim bases to the right of this position (exclusive, 0-based) — type: `int` — program default: `0`
- `forcetrimright2=0` — If positive, trim this many bases on the right end — type: `int` — program default: `0`
- `forcetrimmod=0` — If positive, right-trim length to be equal to zero, modulo this number. Use ftm=5 to remove inaccurate extra bases from 151bp reads (converts to 150bp), common in Illumina runs — type: `int` — program default: `0`
- `restrictleft=0` — If positive, only look for kmer matches in the leftmost X bases — type: `int` — program default: `0`
- `restrictright=0` — If positive, only look for kmer matches in the rightmost X bases — type: `int` — program default: `0`
- `mingc=0` — Discard reads with GC content below this — type: `float` — program default: `0`
- `maxgc=1` — Discard reads with GC content above this — type: `float` — program default: `1`
- `gcpairs=t` — Use average GC of paired reads. Also affects gchist — type: `bool` — program default: `t`
- `tossjunk=f` — Discard reads with invalid characters as bases — type: `bool` — program default: `f`
- `swift=f` — Trim Swift sequences: Trailing C/T/N R1, leading G/A/N R2 — type: `bool` — program default: `f`
- `entropy=-1` — Set between 0 and 1 to filter reads with entropy below that value. Higher is more stringent — type: `float` — program default: `-1`
- `entropywindow=50` — Calculate entropy using a sliding window of this length — type: `int` — program default: `50`
- `entropyk=5` — Calculate entropy using kmers of this length — type: `int` — program default: `5`
- `minbasefrequency=0` — Discard reads with a minimum base frequency below this — type: `float` — program default: `0`
- `entropytrim=f` — Values: f: (false) Do not entropy-trim. r: (right) Trim low entropy on the right end only. l: (left) Trim low entropy on the left end only. rl: (both) Trim low entropy on both ends — type: `enum` — program default: `f`
- `entropymask=f` — Values: f: (filter) Discard low-entropy sequences. t: (true) Mask low-entropy parts of sequences with N. lc: Change low-entropy parts of sequences to lowercase — type: `enum` — program default: `f`
- `entropymark=f` — Mark each base with its entropy value. This is on a scale of 0-41 and is reported as quality scores, so the output should be fastq or fasta+qual — type: `bool` — program default: `f`
- `cardinality=f` — Count unique kmers using the LogLog algorithm. Provides accurate cardinality estimation using minimal memory, with no upper bound on kmer length — type: `bool` — program default: `f`
- `cardinalityout=f` — Count unique kmers in output reads — type: `bool` — program default: `f`
- `loglogk=31` — Use this kmer length for counting — type: `int` — program default: `31`
- `loglogbuckets=2048` — Use this many buckets for counting — type: `int` — program default: `2048`
- `khist=<file>` — Kmer frequency histogram; plots number of kmers versus kmer depth. This is approximate — type: `string` — program default: unset
- `khistout=<file>` — Kmer frequency histogram for output reads — type: `string` — program default: unset

---

### [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

- `N=0` — Number of mismatches allowed in a seed alignment during multiseed alignment (0 or 1). Higher values increase sensitivity but slow alignment — type: `int` — program default: `0`
- `L=22` — Length of seed substrings used in multiseed alignment. Smaller values are slower but more sensitive — type: `int` — program default: `22`
- `i=S,1,1.15` — Function governing the interval between extracted seed substrings, as a function of read length (e.g. S,1,2.5 = 1 + 2.5×sqrt(read length)) — type: `func` — program default: `S,1,1.15`
- `n-ceil=L,0,0.15` — Function governing the maximum number of ambiguous characters (e.g. Ns) allowed in a read, as a function of read length. Reads exceeding this are filtered out — type: `func` — program default: `L,0,0.15`
- `dpad=15` — Number of columns to pad dynamic programming problems on either side, to allow for gaps — type: `int` — program default: `15`
- `gbar=4` — Disallow gaps within this many positions of the read ends — type: `int` — program default: `4`
- `ignore-quals=f` — Treat all quality values as the highest possible when calculating mismatch penalties. Always on for -f/-r/-c modes — type: `bool` — program default: `f`
- `nofw=f` — Disable alignment to the forward (Watson) reference strand — type: `bool` — program default: `f`
- `norc=f` — Disable alignment to the reverse-complement (Crick) reference strand — type: `bool` — program default: `f`
- `no-1mm-upfront=f` — Skip Bowtie2's default attempt at an exact or 1-mismatch end-to-end alignment before the multiseed heuristic, for predictable behavior with N/L. Slower — type: `bool` — program default: `f`
- `end-to-end=t` — Require the entire read to align without soft-clipping either end. Mutually exclusive with local — type: `bool` — program default: `t`
- `local=f` — Allow soft-clipping of read ends to maximize alignment score. Mutually exclusive with end-to-end — type: `bool` — program default: `f`
- `ma=2` — Match bonus added to the alignment score per matching position. Only used in local mode — type: `int` — program default: `2`
- `mp=6,2` — Maximum and minimum mismatch penalties, scaled by base quality unless ignore-quals is set — type: `int,int` — program default: `6,2`
- `np=1` — Penalty for positions where the read, reference, or both contain an ambiguous character (e.g. N) — type: `int` — program default: `1`
- `rdg=5,3` — Read gap open and extend penalties. A gap of length N costs int1 + N*int2 — type: `int,int` — program default: `5,3`
- `rfg=5,3` — Reference gap open and extend penalties. A gap of length N costs int1 + N*int2 — type: `int,int` — program default: `5,3`
- `score-min=L,-0.6,-0.6` — Function governing the minimum alignment score for an alignment to be considered valid, as a function of read length — type: `func` — program default: `L,-0.6,-0.6`

---

### [bwa-mem2](https://bio-bwa.sourceforge.net/bwa.shtml)

>bwa-mem2 has identical parameters and output but is optimized for parallelization via SIMD. Link is for bwa-mem1

- `k=19` — Minimum seed length. Matches shorter than this will be missed. Alignment speed is usually insensitive to this value unless it significantly deviates from 20 — type: `int` — program default: `19`
- `w=100` — Band width. Gaps longer than this will not be found, though the maximum gap length is also affected by the scoring matrix and hit length — type: `int` — program default: `100`
- `d=100` — Off-diagonal X-dropoff (Z-dropoff). Stop extension when the difference between the best and current extension score exceeds |i-j|*A+d — type: `int` — program default: `100`
- `r=1.5` — Trigger re-seeding for a MEM longer than minSeedLen*r. Larger values yield fewer seeds: faster but less accurate — type: `float` — program default: `1.5`
- `c=10000` — Discard a MEM if it has more than this many occurrences in the genome — type: `int` — program default: `10000`
- `P=f` — In paired-end mode, perform SW to rescue missing hits only, without trying to find hits that fit a proper pair — type: `bool` — program default: `f`
- `A=1` — Matching score — type: `int` — program default: `1`
- `B=4` — Mismatch penalty. Approximate sequence error rate is .75 * exp[-log(4) * B/A] — type: `int` — program default: `4`
- `O=6` — Gap open penalty — type: `int` — program default: `6`
- `E=1` — Gap extension penalty. A gap of length k costs O + k*E — type: `int` — program default: `1`
- `L=5` — Clipping penalty. If the best score reaching the end of the query exceeds the best SW score minus this penalty, clipping is not applied (the SAM AS tag still reports the best SW score) — type: `int` — program default: `5`
- `U=9` — Penalty for an unpaired read pair. Used to compare scoreRead1+scoreRead2-U against the paired score to decide whether to force pairing — type: `int` — program default: `9`
- `p=f` — Assume the first input query file is interleaved paired-end FASTA/Q — type: `bool` — program default: `f`
- `R=` — Complete read group header line, e.g. `@RG\tID:foo\tSM:bar`. The read group ID is attached to every output read — type: `string` — program default: unset
- `T=30` — Don't output alignments with a score lower than this. Only affects output — type: `int` — program default: `30`
- `a=f` — Output all found alignments for single-end or unpaired paired-end reads, flagged as secondary alignments — type: `bool` — program default: `f`
- `C=f` — Append FASTA/Q comment to SAM output, to transfer read metadata (e.g. barcode). Comments must conform to the SAM spec (e.g. BC:Z:CGTAC) or output will be malformed — type: `bool` — program default: `f`
- `H=f` — Use hard clipping in the SAM output. May dramatically reduce output redundancy when mapping long contigs or BACs — type: `bool` — program default: `f`
- `M=f` — Mark shorter split hits as secondary, for Picard compatibility — type: `bool` — program default: `f`
- `v=3` — Verbosity level: 0 disables stderr output, 1 errors only, 2 warnings and errors, 3 all normal messages, 4+ for debugging (output is no longer SAM at this level) — type: `int` — program default: `3`

---

### [SAMtools view](https://www.htslib.org/doc/samtools-view.html)

- `mapq=30 (q=30)` — Skip alignments with MAPQ smaller than this — type: `int` — program default: `0`
- `F=1804` — Do not output alignments with any of these bits set in the FLAG field. Can be specified in hex (`0x...`), octal (`0...`), decimal, or as a comma-separated list of flag names — type: `int`|`string` — program default: `0`

---

### [bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)

- `scaleFactor=1.0` — The computed scaling factor (or 1, if not applicable) is multiplied by this — type: `float` — program default: `1.0`
- `MNase=f` — Determine nucleosome positions from MNase-seq data, counting only the 3 central nucleotides of each fragment. Only fragments 130-200bp are considered unless overridden by minFragmentLength/maxFragmentLength. Requires paired-end data; a bin size of 1 is recommended — type: `bool` — program default: `f`
- `Offset=` — Use this offset (or range of offsets) inside each read as the signal, e.g. for RiboSeq/GROseq. Negative values offset from the read end; 0 is not permitted — type: `int`|`int,int` — program default: unset
- `filterRNAstrand=` — Select RNA-seq reads originating from genes on the given strand (`forward` or `reverse`), assuming a standard dUTP library prep. Consider samExcludeFlag instead for other strand-filtering contexts — type: `enum` — program default: unset
- `binSize=1` — Size of the bins, in bases, for the output bigwig/bedgraph file — type: `int` — program default: `50`
- `region=` — Limit the operation to a region of the genome, in `chr:start:end` format. Useful for testing parameters quickly — type: `string` — program default: unset
- `blackListFileName=` — A BED or GTF file of regions to exclude from all analyses, by rejecting genomic chunks that overlap an entry. Adjust effectiveGenomeSize accordingly if used — type: `string` — program default: unset
- `numberOfProcessors=1` — Number of processors to use. Can be `max` or `max/2` — type: `int`|`string` — program default: `1`
- `verbose=f` — Show processing messages — type: `bool` — program default: `f`
- `effectiveGenomeSize=` — The portion of the genome that is mappable, excluding large stretches of Ns and adjusting for excluded repetitive regions. See deepTools' effectiveGenomeSize table for values — type: `int` — program default: unset
- `normalizeUsing=RPGC` — Normalization method for reads per bin: RPKM, CPM, BPM, RPGC, or None. RPGC requires effectiveGenomeSize. Each read is considered independently; use samFlagInclude/samFlagExclude to count only one mate of a pair — type: `enum` — program default: `None`
- `exactScaling=f` — Process all reads to compute exact scaling factors instead of sampling. More accurate but significantly slower — type: `bool` — program default: `f`
- `ignoreForNormalization=` — Space-delimited list of chromosome names to exclude when computing normalization, e.g. for samples with unequal coverage across chromosomes — type: `string` — program default: unset
- `skipNonCoveredRegions=f` — Skip regions with no overlapping reads instead of treating them as zero — type: `bool` — program default: `f`
- `smoothLength=0` — Window size, larger than binSize, over which to average reads for each bin (including neighbors). Values smaller than binSize are ignored — type: `int` — program default: `0`
- `extendReads=t` — Extend reads to fragment size. For single-end data a fragment length must be given (or estimated from data); for paired-end data, mates are extended to match the fragment defined by both ends. Not recommended for spliced-read data like RNA-seq — type: `bool`|`int` — program default: `f`
- `ignoreDuplicates=t` — Count reads with the same orientation, start position (and mate position, if paired) only once — type: `bool` — program default: `f`
- `minMappingQuality=` — If set, only reads with at least this mapping quality are considered — type: `int` — program default: unset
- `centerReads=f` — Center reads with respect to fragment length, producing a sharper signal around enriched regions — type: `bool` — program default: `f`
- `samFlagInclude=` — Include only reads matching this SAM flag, e.g. `64` for first-mate reads only — type: `int` — program default: unset
- `samFlagExclude=` — Exclude reads matching this SAM flag, e.g. `16` to exclude reverse-strand reads — type: `int` — program default: unset
- `minFragmentLength=0` — Minimum fragment length required for read/pair inclusion. Useful for filtering mono-/di-nucleosome fragments in ATAC-seq — type: `int` — program default: `0`
- `maxFragmentLength=600` — Maximum fragment length required for read/pair inclusion — type: `int` — program default: `0`

---

### [MACS3 callpeak](https://macs3-project.github.io/MACS/docs/callpeak.html)

- `gsize=hs` — Mappable/effective genome size (smaller than the raw genome size due to repeats). Precompiled shortcuts: `hs` (~2.9e9, human), `mm` (~2.65e9, mouse), `ce` (~1e8, C. elegans), `dm` (~1.4e8, fly). Check deepTools for other assemblies, or estimate by removing Ns and simple repeats from the genome — type: `string`|`int` — program default: `hs`
- `tsize=` — Sequencing tag size. If unset, MACS3 determines it from the first 10 reads of the treatment file — type: `int` — program default: unset (auto-detected)
- `qvalue=0.05` — q-value (FDR) cutoff for significant regions, computed via Benjamini-Hochberg. Try `0.01` for broad marks — type: `float` — program default: `0.05`
- `pvalue=` — p-value cutoff. If set, overrides qvalue — type: `float` — program default: unset
- `min-length=` — Minimum length of a called peak. If unset, defaults to the predicted fragment size d. Setting it shorter than the fragment size may have no effect — type: `int` — program default: unset (= fragment size)
- `max-gap=` — Maximum gap between nearby regions to merge into one peak. If unset, defaults to the detected read length. With `--broad`, weaker/broad peaks are merged using 4x this value — type: `int` — program default: unset (= read length)
- `nolambda=f` — Use the background lambda as the local lambda everywhere, ignoring local bias at candidate peaks. Recommended when calling peaks without a control — type: `bool` — program default: `f`
- `slocal=1000` — Size of the small local region (bp) used to compute local lambda — type: `int` — program default: `1000`
- `llocal=10000` — Size of the large local region (bp) used to compute local lambda, capturing long-range bias like open chromatin domains — type: `int` — program default: `10000`
- `nomodel=f` — Skip building the shifting model. Combine with extsize and shift to control fragment extension manually — type: `bool` — program default: `f`
- `extsize=` — When nomodel is set, extend reads in the 5'->3' direction to this fixed fragment size — type: `int` — program default: unset
- `shift=0` — Shift read positions by this many bp before extension (when nomodel is set). Negative shifts move 5' ends 3'->5'. Must be 0 for paired-end (BAMPE/BEDPE) data. For DNase-seq-style cutting-site enrichment use `shift=-100 extsize=200`; for nucleosome centers use `shift=37 extsize=73` — type: `int` — program default: `0`
- `keep-dup=1` — How to handle duplicate reads at the exact same coordinate and strand: an integer caps the count kept, `auto` calculates the max via a binomial p-value cutoff of 1e-5, `all` keeps every read. Ignored for FRAG format — type: `int`|`enum` — program default: `1`
- `broad=f` — Call broad peaks (UCSC gappedPeak), nesting weaker/broad peaks (cutoff via broad-cutoff) within stronger/narrow peaks (cutoff via -p/-q). The merge gap for broad peaks is 4x max-gap — type: `bool` — program default: `f`
- `broad-cutoff=0.1` — Cutoff for broad regions when broad is set. A p-value cutoff if pvalue is set, otherwise a q-value cutoff — type: `float` — program default: `0.1`
- `scale-to=small` — When `large`, linearly scale the smaller dataset up to the larger dataset's depth; when `small` (default), scale the larger dataset down. Scaling up is more prone to false positives — type: `enum` — program default: `small`
- `call-summits=f` — Reanalyze the signal shape within each peak to deconvolve overlapping subpeaks, useful for detecting adjacent binding events. Subpeaks share the parent peak's boundaries but get distinct summits/scores — type: `bool` — program default: `f`
- `min_foldch=2.0` — Minimum fold-change (narrowPeak column 7) for a peak to be kept in `*_peaks_filt.narrowPeak`. This is a pipeline-level post-filter applied after callpeak, not a MACS3 option — type: `float` — program default: `2.0`

---

### [MEME Motif](https://meme-suite.org/meme/doc/meme.html)

- `objfun=classic` — Objective function used to choose motif width/sites, score significance, and pick the best motif. Options: `classic` (original E-value approach), `de`/`se` (enrichment vs. control sequences), `cd`/`ce` (enrichment near sequence centers, for ChIP/CLIP-seq), `nc` (numerically-correct classic, very slow) — type: `enum` — program default: `classic`
- `test=mhg` — Statistical test for motif enrichment when objfun is `de` or `se`: multiple hypergeometric (`mhg`), binomial (`mbn`), or rank-sum (`mrs`) — type: `enum` — program default: `mhg`
- `use_llr=f` — Use the log-likelihood ratio method (instead of the E-value method) to evaluate EM starting points. Only affects the `classic` objective function — type: `bool` — program default: `f`
- `neg=` — FASTA file of control sequences (same alphabet as primary sequences), used by enrichment-based objective functions. If unset and required, MEME shuffles the primary sequences to create controls — type: `string` — program default: unset
- `shuf=2` — Kmer size (1-6) whose frequencies are preserved when shuffling primary sequences to create control sequences. Avoid 1 with masked (e.g. RepeatMasker) sequences — type: `int` — program default: `2`
- `hsfrac=0.5` — Fraction of primary sequences held out for estimating motif width, site count, and significance, for objective functions that require hold-out sets — type: `float` — program default: `0.5`
- `cefrac=0.25` — Fraction of each sequence's length used as the central region for the Central Enrichment (and Central Distance) objective functions — type: `float` — program default: `0.25`
- `searchsize=` — Max total size (letters) of primary sequences searched for motifs; larger inputs are randomly subsampled. `0` disables sampling (uses everything, very slow). Doesn't affect reported motif frequencies/sites — type: `int` — program default: sampling occurs above `100000` letters
- `norand=f` — Don't randomize sequence order before applying searchsize sampling. Use if sequences are already ordered by confidence — type: `bool` — program default: `f`
- `csites=1000` — Max sites used for E-value computation with the `classic` objective function; larger sets are subsampled (oops/zoops) or capped (anr) — type: `int` — program default: `1000`
- `seed=0` — Random seed used for shuffling and sampling sequences — type: `int` — program default: `0`
- `alph=` — Custom alphabet definition file (may include wildcard symbols `*`/`-`). Incompatible with `dna`/`rna`/`protein` — type: `string` — program default: unset (protein alphabet assumed)
- `dna=f` — Use the standard DNA alphabet, treating ambiguous characters as unknown. Incompatible with `alph`/`rna`/`protein` — type: `bool` — program default: `f`
- `rna=f` — Use the standard RNA alphabet. Incompatible with `alph`/`dna`/`protein` — type: `bool` — program default: `f`
- `protein=f` — Use the standard protein alphabet. Incompatible with `alph`/`dna`/`rna` — type: `bool` — program default: `f`
- `revcomp=f` — Search both the given strand and its reverse complement (DNA/RNA only) — type: `bool` — program default: `f`
- `pal=f` — Only search for palindromes, averaging letter frequencies across complementary motif columns during EM — type: `bool` — program default: `f`
- `mod=anr` — Distribution of motif sites per sequence: `oops` (exactly one per sequence, fastest), `zoops` (zero or one, ~2x slower), `anr` (any number of repeats, ~10x slower; `tcm` is a synonym) — type: `enum` — program default: `zoops`
- `nmotifs=2` — Number of motifs to find before stopping — type: `int` — program default: `1`
- `evt=` — Stop searching once the most recently found motif has an E-value above this (or after 1000 motifs), unless nmotifs is also set — type: `float` — program default: unset
- `time=` — Stop searching for additional motifs once at least one has been found and continuing is estimated to exceed this many seconds — type: `int` — program default: unset
- `nsites=` — Force motifs to have exactly this many sites (shorthand for setting minsites and maxsites equal). Ignored under the `oops` model — type: `int` — program default: unset
- `minsites=2` — Minimum number of sites a motif may have. Ignored under the `oops` model — type: `int` — program default: `2`
- `maxsites=` — Maximum number of sites a motif may have. Ignored under `oops`; defaults to the number of sequences under `zoops`, or 5x that under `anr` — type: `int` — program default: unset
- `wnsites=0.8` — Strength (0-1) of the bias toward motifs with exactly the expected number of sites (per nsites/minsites/maxsites). Closer to 1 is a stronger bias — type: `float` — program default: `0.8`
- `w=` — Search only for motifs of exactly this width, instead of the minw-maxw range — type: `int` — program default: unset
- `minw=8` — Minimum motif width to search for — type: `int` — program default: `8`
- `maxw=32` — Maximum motif width to search for — type: `int` — program default: `50`
- `allw=f` — Find EM starting points for every width from minw to maxw (instead of a geometric progression). Slow for large width ranges — type: `bool` — program default: `f`
- `nomatrim=f` — Don't trim motif width using multiple alignments. Only affects the `classic` objective function, which trims by default — type: `bool` — program default: `f`
- `wg=11` — Gap opening cost for the alignments used to trim motif width — type: `int` — program default: `11`
- `ws=1` — Gap extension cost for the alignments used to trim motif width — type: `int` — program default: `1`
- `noendgaps=f` — Don't penalize end gaps in the alignments used to trim motif width — type: `bool` — program default: `f`
- `bfile=` — Markov background model file, used as the EM null model and for log-likelihood/significance/PSSM calculations. If unset, a model is estimated from the primary sequences — type: `string` — program default: unset (estimated from input)
- `markov_order=0` — Maximum order of the Markov background model to use/create. If bfile is unset, MEME builds this order of model from the primary sequences with a 0.1 pseudocount — type: `int` — program default: `0`
- `psp=` — Position-specific priors file biasing where motif sites are expected within sequences — type: `string` — program default: unset
- `maxiter=50` — Maximum EM iterations per starting point (or until convergence, see `distance`) — type: `int` — program default: `50`
- `distance=0.001` — EM convergence threshold: stop iterating when the change (Euclidean distance) between successive frequency matrices drops below this — type: `float` — program default: `0.001`
- `prior=` — Type of prior on model parameters: `dirichlet` (simple, default for DNA), `dmix`/`megap`/`mega` (Dirichlet mixtures, default for protein depending on mod), or `addone` — type: `enum` — program default: unset (chosen automatically based on alphabet/mod)
- `b=` — Strength of the prior on model parameters (`0` = use the prior's intrinsic strength for `dmix`). Default is `0.01` for `dirichlet`, `0` for `dmix` — type: `float` — program default: unset
- `plib=` — Dirichlet mixtures prior library file. Default depends on alphabet (DNA/protein) — type: `string` — program default: unset
- `spfuzz=` — Fuzziness of the sequence-to-theta mapping; meaning depends on spmap (uniform prior weight, or PAM mutation distance) — type: `float` — program default: unset (depends on spmap)
- `spmap=` — Mapping function for estimating theta: `uni` (uniform prior, default for DNA) or `pam` (PAM matrices, default for protein) — type: `enum` — program default: unset (chosen automatically based on alphabet)
- `cons=` — (May be repeated.) Seed a starting point from this consensus sequence instead of sampling, suppressing sampling for that many motifs. Pad short DNA/RNA consensus sequences with `N`s to width 6 — type: `string` — program default: unset
- `maxpeaks=100` — Maximum number of peaks written to the input FASTA, keeping those with the highest fold-change. This is a pipeline-level option for `narrow_peak_to_fasta.py`, not a MEME option — type: `int` — program default: `100`
- `summit_extend=50` — Number of bp on each side of the peak summit to include in the input FASTA for the "summits" mode (use `all` to use the full peak instead). This is a pipeline-level option for `narrow_peak_to_fasta.py`, not a MEME option — type: `int`|`enum` — program default: `50`
- `complexity_filter.enabled=false` — Enable low-complexity filtering of peak FASTA sequences before MEME. This is a pipeline-level option for `narrow_peak_to_fasta.py`, not a MEME option — type: `bool` — program default: `false`
- `complexity_filter.min_entropy=3.0` — Minimum 3-mer Shannon entropy in bits for retained peak FASTA sequences when `complexity_filter.enabled` is true. Lower values keep more repetitive sequence; higher values are more aggressive. Max possible is `6.0`; calibrate against your own data — type: `float` — program default: `3.0`
