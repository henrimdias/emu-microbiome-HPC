
# High-fidelity microbiome profiling in soil and plants using a reproducible EMU workflow based on HPC environment

Accurately profiling soil and root-associated bacterial communities generates vast volumes of sequencing data that can overwhelm standard desktop computers. This protocol harnesses powerful high-performance computing clusters to process millions of full-length ribosomal ribonucleic acid gene reads, ensuring rapid quality filtering, organelle removal, and exact sequence matching. By providing user-friendly command-line scripts and automated workflows, it lowers the barrier to handling massive datasets, minimizes errors in data processing, and enables transparent, reproducible analyses. Researchers can make confident decisions based on robust community profiles of soil and plant-associated bacteria, accelerating insights into soil health, plant growth promotion, and sustainable agricultural practices.

### How it works

Before launching into the four protocols outlined above, users should take time to organize their computational environment, define key parameters, and map out the overall workflow. The diagram in Figure 1 provides a high‐level flowchart of the EMU pipeline from raw Nanopore FASTQ files through quality control, organelle filtering, exact‐matching taxonomic assignment, and downstream ecological analyses. Below we summarize the main considerations and options that will ensure a smooth, reproducible execution of the protocol.

![alt text](Figure_1_workflow.png)

#### Project directory layout (recommended)

```text
/$USER/scratch/emu_pipeline/
├── raw_data/
├── nanoplot_reports/
├── nanofilt_reports/
├── nanoplot_after_nanofilt_data/
├── qc_summary_reports/
├── filtered_data/
├── kraken2_db/
├── no_organelle_filtered_data/
├── emu_db/
├── emu_results/
├── downstream_analysis/
└── scripts/
```

### Full protocol

All steps, parameters, and troubleshooting notes are described in the published protocol.
If you use this protocol for your analysis, please cite it using the following [doi: xxxxxxx](https://www.biorxiv.org/content/10.1101/2024.11.12.623179v1)


