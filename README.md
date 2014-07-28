mnegri-ov170
============

Scripts and raw data used by the XXX publication.

This repository has several subdirectories:

- config: Includes configuration data for bcbio-nextgen
- data: Includes results (VCFs) and other data files
- programs: Includes the programs mentioned in the paper (Supplementary Section 2)

### Programs

``dump_counts`` and ``extract_shared_mutations`` require Python 2, while the other programs require Python 3.4.

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

### Data files

The data files included are as follows:

- ``all_somatic_mutations.vcf.gz``: bgzipped VCF file with all the found somatic mutations
- ``all_somatic_mutations.vcf.gz.tbi``: tabix index for the above file
- ``Ft-S_mutation_counts.txt``: Mutation counts for Ft-S samples
- ``Sd-S_mutation_counts.txt``: Mutation counts for Sd-S samples
- ``phenotype_color_table.txt``: Table to use with ``mutation_heatmap`` to plot OS, PFS, and other phenotypes
- ``pathway_genes_and_colors.txt``: Table to use with ``mutation_heatmap`` to plot the genes associated with pathways

