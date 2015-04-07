mnegri-ov170
============

This repository holds scripts and raw data used by the publication:
    
L. Beltrame, M. Di Marino, R. Fruscio, E. Calura, B. Chapman, L. Clivio, F. Sina, C. Mele, P. Iatropoulos, T. Grassi, V. Fotia, C. Romualdi, P. Martini, M. Noris, L. Paracchini, Ilaria Craparotta, M. Petrillo, R. Milani, P. Perego, A. Ravaggi, A. Zambelli, E. Ronchetti, M. D'Incalci, and S. Marchini
**Profiling cancer gene mutations in longitudinal epithelial ovarian cancer biopsies by targeted next-generation sequencing: a retrospective study**. *Annals of Oncology*, 2015.

and it is divided in several subdirectories:

- ``config``: Includes configuration data for bcbio-nextgen
- ``data``: Includes results (VCFs) and other data files
- ``programs``: Includes the programs mentioned in the paper (Supplementary Section 2)

### Programs

``dump_counts`` and ``extract_shared_mutations`` require Python 2.7, while the other programs require Python 3.4.

Other dependencies:

- ``cruzdb``
- ``matplotlib``
- ``sqlalchemy``
- ``pandas``
- ``gemini``
- ``sarge``
- ``fastcluster``
- ``numpy``
- ``scipy``
- ``pathlib`` (for Python 2.7 programs)

All programs have command-line help. 

### Data files

The data files included are as follows:

- ``all_somatic_mutations.vcf.gz``: bgzipped VCF file with all the found somatic mutations
- ``all_somatic_mutations.vcf.gz.tbi``: tabix index for the above file
- ``Ft-S_mutation_counts.txt``: Mutation counts for Ft-S samples
- ``Sd-S_mutation_counts.txt``: Mutation counts for Sd-S samples
- ``phenotype_color_table.txt``: Table to use with ``mutation_heatmap`` to plot OS, PFS, and other phenotypes
- ``pathway_genes_and_colors.txt``: Table to use with ``mutation_heatmap`` to plot the genes associated with pathways
