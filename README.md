# Hackathon LLM evals



## How to run skills with an LLM

[Docs](https://code.claude.com/docs/en/skills#discovery-from-parent-and-nested-directories)
_Claude Code loads project skills from .claude/skills/ in the directory where you start it and in every parent directory up to the repository root.[...]_



## Examples: 
---


### 🧬 ATAC-seq to BAM files pipeline

**Goal**: Provided two ATAC-seq files (or givin a path where these live) write a script to generate bam files mapped to the mouse genome.

**Prompt**:
```
Use available skills you have access to whenever possible. 
task: I would like to process some ATAC-Seq data (that lives in the sandbox folder in this project), can you please write me a python script to generate bam files mapped to the mouse genome Mus_musculus.GRCm39.116.gtf.gz. Note: use a pixi environment to use python 3.12 and write also a README.md to see the steps and the compute requirements before starting with this task. 
```

**Skills Used**: exploratory-data-analysis, database-lookup, genomic-coordinates, genomic-intelligence, bulk-rnaseq

---

