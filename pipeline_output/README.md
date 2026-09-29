# Pipeline output files

These are the files written by the bioinformatic pipeline of Witsø et al. (2025), PLOS One
20(9): e0330754 (<https://doi.org/10.1371/journal.pone.0330754>), for the 36 plastisphere
metagenomes plus one negative control. You build every analysis table from them. The full
description of the methods is in the paper (`docs/`), this file explains the files and their
formats.

Sample naming: `R1`-`R24` river (Lier river, two locations), `V1`-`V6` raw wastewater,
`V7`-`V12` treated wastewater (VEAS treatment plant), `Blank` = extraction and library
negative control. Every file name starts with the sequencing run ID
(`<sample>_<library>_<flowcell>_<lane>`, e.g. `R1_EKDN230028429-1A_HKTC2DSX7_L3`). Sample V8
was sequenced twice: `V8_..._HHGM2DSX7_L1` (4.9 million read pairs) and `V8_..._HJH72DSX7_L2`
(47 million). `metadata/sample_metadata.csv` links samples, runs and the study design.

The pipeline, in order:

| step | tool | output folder |
|---|---|---|
| Sequencing: Illumina NovaSeq 6000, 2 x 150 bp, adapter and quality trimmed by the provider | | |
| Read quality control | FastQC 0.12.1 + MultiQC | `01_read_qc` |
| Read mapping to resistance genes | KMA 1.4.9 vs CARD protein homolog models (3.2.8, see below) | `02_kma_card`, `03_card_database` |
| Taxonomic profiling of reads | Kraken2 2.1.3 (standard DB 2023-10-09) + Bracken 2.9 | `04_kraken2_bracken` |
| Assembly (MEGAHIT 1.2.9, contigs >= 1000 bp), BLAST of contigs vs CARD, flanking-region filter, Kraken2 on ARG contigs | | `05_contigs_blast_card` |
| Binning (MetaBAT2 2.15) and bin quality (CheckM 1.2.2) | | `06_mags_checkm` |
| Taxonomy of high-quality MAGs | GTDB-Tk 2.3.2, GTDB release 214 | `07_mags_gtdbtk` |
| Gene prediction and ARG screening of high-quality MAGs | Prokka + AMRFinderPlus; BLAST vs CARD | `08_mags_amrfinder` |
| One complete example MAG | | `09_example_mag` |

Not included because of size: raw reads (NCBI BioProject PRJNA1219086), Kraken2 report
files, assemblies, bin FASTA files, full Prokka annotations, KMA alignment and variant
files.

## 01_read_qc

- `multiqc_general_stats.txt`: one row per read file (`<run>_1`, `<run>_2` are the two
  mates); `total_sequences` = number of reads in that file = read pairs of the run,
  `percent_duplicates`, `percent_gc`, `avg_sequence_length`.
- `multiqc_fastqc.txt`: the FastQC module pass/warn/fail flags per file.
- `multiqc_report.html` and `plots/`: the MultiQC report; open in a browser.

## 02_kma_card

KMA was run as `kma -ipe <run>_1.fq.gz <run>_2.fq.gz -mem_mode -ef -1t1 -cge -nf -vcf -t_db nt_protein_homolog`
(index built from CARD `nucleotide_fasta_protein_homolog_model.fasta`). Three files per run:

- `<run>.res`: one row per reference gene ("template") with reads. Columns: `#Template`,
  `Score`, `Expected`, `Template_length` (bp), `Template_Identity` (% of the template
  identical to the consensus of the mapped reads), `Template_Coverage` (% of the template
  covered by reads; > 100 with insertions), `Query_Identity` (% identity within the aligned
  part), `Query_Coverage`, `Depth` (mean read depth), `q_value`, `p_value`. The template
  name contains the ARO accession: `gb|M14039.1|+|412-898|ARO:3002835|lnuA [Staphylococcus haemolyticus]`.
- `<run>.mapstat`: six comment lines (`## method`, `## version`, `## database`,
  `## fragmentCount` = fragments processed, `## date`, `## command`), then a header
  (`# refSequence ...`) and one row per template: `readCount`, `fragmentCount`, `mapScoreSum`,
  `refCoveredPositions`, `refConsensusSum`, `bpTotal`, `depthVariance`,
  `nucHighDepthVariance`, `depthMax`, `snpSum`, `insertSum`, `deletionSum`, `readCountAln`,
  `fragmentCountAln`.
- `<run>.fsa`: the consensus sequence of every template in FASTA format.

Full definitions: `docs/KMAspecification.pdf`.

## 03_card_database

From the CARD 3.2.8 release (<https://card.mcmaster.ca>, Alcock et al. 2023, NAR 51:D690).

- `aro_index.tsv`: one row per reference model. `ARO Accession` is the key; `CARD Short
  Name` is the gene name used in figures; `AMR Gene Family`, `Drug Class` (several classes
  separated by `;`), `Resistance Mechanism`.
- `aro_categories_index.tsv`, `shortname_antibiotics.tsv`, `shortname_pathogens.tsv`,
  `CARD-Download-README.txt` (licence: free for academic use).

## 04_kraken2_bracken

- `<run>.species.bracken` and `<run>.genus.bracken`: Bracken output, tab-separated:
  `name`, `taxonomy_id`, `taxonomy_lvl` (`S` or `G`), `kraken_assigned_reads` (reads
  Kraken2 placed exactly on the taxon), `added_reads` (reads Bracken re-distributed from
  higher ranks), `new_est_reads` (the final estimate), `fraction_total_reads`. Human reads
  removed. Kraken2 counts a read pair as one sequence. Only the deep run of V8; no Blank.
- `bact_count.txt`: run and number of read pairs assigned to the domain Bacteria (from the
  Kraken2/Bracken report). This is the denominator of the RPKM normalisation in the paper.

The per-sample files were re-created from the merged Bracken table of the original analysis
because the original per-sample files were no longer available; the content is the Bracken
output.

## 05_contigs_blast_card

- `blast_output/<run>_blast_output.txt`: BLASTn of all contigs against CARD (outfmt 6 with
  a header line): `qseqid` (contig), `sseqid` (CARD template), `salltitles`, `length`
  (alignment), `pident`, `evalue`, `qlen`, `slen` (template length), `qstart`, `qend`,
  `sstart`, `send`.
- `blast_output/filtered_<run>_blast_output.txt`: hits with identity >= 80 % and alignment
  covering >= 80 % of the template (`length_slen_ratio`).
- `blast_output/flanking_<run>_blast_output.txt`: the filtered hits with the length of
  contig sequence on each side of the gene (`flanking_left`, `flanking_right`). Files exist
  only for runs with such hits.
- `contig_names_amr_flank.tsv`, `combined_flanking_names.tsv`, `contig_and_file.txt`: helper
  lists of contigs with >= 1500 bp flanks and the run they came from.
- `contigs_flank.kraken_out.txt`: Kraken2 classification of those contigs (standard Kraken2
  output: `C`/`U`, contig, `name (taxid N)`, length, k-mer assignments).
- `contigs_flank.kraken_report.txt`: the Kraken2 report for the same contigs.
- `contigs_with_taxonomy.txt`: contig and classification, two columns.

## 06_mags_checkm

- `checkm_output.csv`: CheckM `lineage_wf` summary for all 3282 bins, `;`-separated
  (there is a trailing `;`, which gives an empty last column): `Bin_Id`, `Marker lineage`,
  `# genomes`, `# markers`, `# marker sets`, `0` `1` `2` `3` `4` `5+` (marker gene
  counts), `Completeness`, `Contamination`. Bin IDs are `<run>.<bin number>`.
- `bin_stats.analyze.tsv`: bin ID and a Python dictionary string with genome statistics
  (`'Genome size'`, `'# contigs'`, `'N50 (contigs)'`, `'GC'`, `'Coding density'`,
  `'# predicted genes'`, ...).

## 07_mags_gtdbtk

- `gtdbtk.bac120.summary.tsv` (531 bacterial MAGs) and `gtdbtk.ar53.summary.tsv` (6
  archaeal MAGs): `user_genome`, `classification` (`d__...;p__...;c__...;o__...;f__...;g__...;s__...`,
  empty rank = not resolved), the closest reference genome and its ANI, `classification_method`,
  `red_value`, warnings.
- `gtdbtk.json`, `gtdbtk.log`: the command and run log.

## 08_mags_amrfinder

- `<bin>_amrfinder_protein.txt`: AMRFinderPlus result on the Prokka proteins of each
  high-quality MAG (537 files, including 4 bins of the shallow V8 run). Standard AMRFinderPlus
  columns; a file with only the header line means no gene was found. `Element type` is
  `AMR` for resistance genes; `Class`/`Subclass` give the drug class.
- `prokka_summary/<bin>.txt`: Prokka's summary (contigs, bases, CDS, rRNA, tRNA).
- `blast_card/<bin>_blast_output.txt`: BLASTn of the MAG's contigs against CARD (same
  columns as in 05, without header line).

## 09_example_mag

`V7_EKDN230028459-1A_HKTTNDSX7_L4.11` (treated WW; GTDB: *Cloacibacterium*; CheckM: 96.5 %
complete, 2.9 % contamination; one class A beta-lactamase): contigs (`.fa`), predicted
proteins (`.faa`), Prokka annotation table (`.tsv`: locus tag, type, length, gene, EC
number, product, in contig order), Prokka summary (`.txt`) and its AMRFinderPlus result.
