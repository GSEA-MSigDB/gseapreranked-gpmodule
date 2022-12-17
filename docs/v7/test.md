# GSEAPreranked (v7.4.x)

Runs the gene set enrichment analysis against a user-supplied ranked
list of genes.

**Author:** Chet Birger, David Eby; Broad Institute

**Contact:**

[See the GSEA forum](https://groups.google.com/forum/#!forum/gsea-help)
for GSEA questions.

[Contact the GenePattern
team](http://software.broadinstitute.org/cancer/software/genepattern/contact)
for GenePattern issues.

**GSEA Version:** 4.3._x_

## Introduction

GSEAPreranked runs **Gene Set Enrichment Analysis** (GSEA) against a
user-supplied, ranked list of genes. It determines whether *a priori*
defined sets of genes show statistically significant enrichment at
either end of the ranking. A statistically significant enrichment
indicates that the biological activity (e.g., biomolecular pathway)
characterized by the gene set is correlated with the user-supplied
ranking.

## Details

**Gene Set Enrichment Analysis** (GSEA) is a powerful analytical method
for interpreting gene expression data. It evaluates cumulative changes
in the expression of groups of multiple genes defined based on prior
biological knowledge.

The GSEAPreranked module can be used to conduct gene set enrichment
analysis on data that do not conform to the typical GSEA scenario. For
example, it can be used when the ranking metric choices provided by the
GSEA module are not appropriate for the data, or when a ranked list of
genomic features deviates from traditional microarray expression data
(e.g., GWAS results, ChIP-Seq, RNA-Seq, etc.).

The user provides GSEAPreranked with a pre-ranked gene list. Paired
with each gene in the list is the numeric ranking statistic, which
GSEAPreranked uses to rank order genes in descending order.
GSEAPreranked calculates an enrichment score for each gene set. A gene
set’s enrichment score reflects how often members of that gene set occur
at the top or bottom of the ranked data set (for example, in expression
data, in either the most highly expressed genes or the most
underexpressed genes).

### The ranked list must not contain duplicate ranking values.

Duplicate ranking values may lead to arbitrary ordering of genes and to
erroneous results. Therefore, it is important to make sure that the
ranked list contains no duplicate ranking values.

### Permutation test

In GSEAPreranked, permutations are always done by gene set. In standard
GSEA, you can choose to set the parameter *Permutation type* to
*phenotype* (the default) or *gene set*, but GSEAPreranked does not
provide this
option.

### Understand and keep in mind how GSEAPreranked computes enrichment scores.

The GSEA PNAS 2005 paper introduced a method where a running sum
statistic is incremented by the absolute value of the ranking metric
when a gene belongs to the set. This method has proven to be efficient
and facilitates intuitive interpretation of ranking metrics that reflect
correlation of gene expression with phenotype. In the case of
GSEAPreranked, you should make sure that this weighted scoring scheme
applies to your choice of ranking statistic. If in doubt, we recommend
using a more conservative scoring approach by setting *scoring
scheme* parameter to *classic;* however, the scoring scheme parameter’s
default value is *weighted*, the default value employed by the GSEA
module. Please refer to the GSEA PNAS 2005 paper for further details.

## References

Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA,
Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP. Gene set
enrichment analysis: A knowledge-based approach for interpreting
genome-wide expression profiles. *PNAS*. 2005;102(43);15545-15550.
([link](http://www.pnas.org/content/102/43/15545.full.pdf.html))

Mootha VK, Lindgren CM, Eriksson K-F, Subramanian A, Sihag S, Lehar J,
Puigserver P, Carlsson E, Ridderstrale M, Laurila E, Houstis N, Daly MJ,
Patterson N, Mesivor JP, Golub TR, Tamayo P, Spiegelman B, Lander ES,
Hirschhorn JN, Altshuler D, Groop LC. PGC-1-α responsive genes involved
in oxidative phosphorylation are coordinately downregulated in human
diabetes. *Nat Genet*. 2003;34:267-273.
([link](http://www.nature.com/ng/journal/v34/n3/full/ng1180.html))

GSEA User Guide:
<http://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html>

GSEA website: <http://www.gsea-msigdb.org/>

This version of the module is based on the GSEA v4.1.x code base. See
the [Release
Notes](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/GSEA_v4.1.x_Release_Notes)
for new features and other notable changes.

## Parameters

**NOTE**: Certain parameters are considered to be "advanced"; that is,
they control details of the GSEAPreranked algorithm that are typically
not changed. You should not override the default values unless you are
conversant with the algorithm. These parameters are marked "Advanced"
in the parameter descriptions.

- **ranked list**<span style="color:red;">*</span>
    - This is a file in [RNK](http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29) format that contains the rank ordered gene (or feature) list.
- **gene sets database**<span style="color:red;">*</span>
    - This parameter's drop-down allows you to select gene sets from the [Molecular Signatures Database (MSigDB)](http://www.gsea-msigdb.org/gsea/msigdb/index.jsp) on the GSEA website. This drop-down provides access to only the most current version of MSigDB. You can also upload your own gene set file(s) in [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29), [GMX](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29), or [GRP](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GRP:_Gene_set_file_format_.28.2A.grp.29) format. 
    - If you want to use files from an earlier version of MSigDB you will need to download them from the archived releases on the [website](http://www.gsea-msigdb.org/gsea/downloads.jsp).
- **number of permutations**<span style="color:red;">*</span>
    - Specifies the number of permutations to perform in assessing the statistical significance of the enrichment score. It is best to start with a small number, such as 10, in order to check that your analysis will complete successfully (e.g., ensuring you have gene sets that satisfy the minimum and maximum size requirements). After the analysis completes successfully, run it again with a full set of permutations. The recommended number of permutations is 1000. 
    - Default: 1000
- **collapse dataset**<span style="color:red;">*</span>
    - Select whether to collapse each probe set in the expression dataset into a single vector for the gene, which gets identified by its gene symbol. It is also possible to remap symbols from one namespace to another without collapsing (an error will occur if multiple source genes map to a single destination gene). 
    - When using the *Collapse* or *Remap_Only* mode with an annotated CHIP (such as those from MSigDB), the resulting reports will also be annotated.
    - *No_Collapse* will use the dataset as-is, with its native feature identifiers. When you select this option, the chip annotation file (*chip platform* parameter) is ignored and you must specify a gene set file (*gene sets database file* parameter) that identify genes using the same feature (gene or probe) identifiers as is used in your expression dataset. 
    - Default: *Remap_Only*
- **chip platform**
    - This drop-down allows you to specify the chip annotation file, which lists each probe on a chip and its matching HUGO gene symbol, used for the expression array. This parameter is required if *collapse dataset* is set to true. The chip files listed here are from the [GSEA website](http://www.gsea-msigdb.org/gsea/downloads.jsp). If you used a file not listed here, you will need to provide it (in [CHIP](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CHIP:_Chip_file_format_.28.2A.chip.29) format) using 'Upload your own file'. 
    - Please see the [MSigDB 7.0 Release Notes](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/MSigDB_v7.0_Release_Notes) for information about symbol remapping.
- **scoring scheme**<span style="color:red;">*</span>
    - The enrichment statistic. This parameter affects the running-sum statistic used for the enrichment analysis, controlling the value of p used in the enrichment score calculation. Options are:
      - classic Kolmorogorov-Smirnov: p=0 
      - weighted (default): p=1; a running sum statistic that is incremented by the absolute value of the ranking metric when a gene belongs to the set (see the [2005 PNAS paper](http://www.gsea-msigdb.org/gsea/doc/subramanian_tamayo_gsea_pnas.pdf) for details)
      - weighted_p2: p=2 
      - weighted_p1.5: p=1.5
- **max gene set size**<span style="color:red;">*</span>
    - After filtering from the gene sets any gene not in the expression dataset, gene sets larger than this are excluded from the analysis. 
    - Default: 500 
- **min gene set size**<span style="color:red;">*</span>
    - After filtering from the gene sets any gene not in the expression dataset, gene sets smaller than this are excluded from the analysis. 
    - Default: 15
- **collapsing mode for probe sets with more than one match**<span style="color:red;">*</span>
    - Collapsing mode for sets of multiple probes for a single gene. Used only when the *collapse dataset* parameter is set to *Collapse*. Select the expression values to use for the single probe that will represent all probe sets for the gene. For custom ranking metrics, be very cautious when selecting any of these modes to be sure it is compatible with your metric. Options are:
      - Abs_max_of_probes (default): For each sample, use the expression value for the probe set with the maximum **absolute value**.  Note that each value retains its original sign but is chosen based on absolute value. In other words, the largest magnitude value is used. While this method is useful with computational-based input datasets it is generally **not recommended** for use with quantification-based expression measures such as counts or microarray fluorescence.
      - Max_probe: For each sample, use the maximum expression value for the probe set. That is, if there are three probes that map to a single gene, the expression value that will represent the collapsed probe set will be the maximum expression value from those three probes.
      - Median_of_probes: For each sample, use the median expression value for the probe set. 
      - Mean_of_probes: For each sample, use the mean expression value for the probe set. 
      - Sum_of_probes: For each sample, sum all the expression values of the probe set.
- **normalization mode**<span style="color:red;">*</span>
    - Method used to normalize the enrichment scores across analyzed gene sets. Options are:
      - meandiv (default): GSEA normalizes the enrichment scores as described in [Normalized Enrichment Score (NES)](http://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Normalized_Enrichment_Score) in the GSEA User Guide. 
      - None: GSEA does not normalize the enrichment scores.
- **omit features with no symbol match**<span style="color:red;">*</span>
    - Used only when *collapse dataset* is set to *Collapse*. By default (*true*), the new dataset excludes probes/genes that have no gene symbols. Set to *false* to have the new dataset contain all probes/genes that were in the original dataset.
- **make detailed gene set report**<span style="color:red;">*</span>
    - Create detailed gene set report (heat map, mountain plot, etc.) for each enriched gene set. 
    - Default: true
- **num top sets**<span style="color:red;">*</span>
    - GSEAPreranked generates summary plots and detailed analysis results for the top x genes in each phenotype, where x is 20 by default. The top genes are those with the largest normalized enrichment scores. 
    - Default: 20
- **random seed**<span style="color:red;">*</span>
    - Seed used to generate a random number for phenotype and gene_set permutations. Timestamp is the default. Using a specific integer-valued seed generates consistent results, which is useful when testing software.
- **output file name**<span style="color:red;">*</span> 
    - Name of the output file. The name cannot include spaces. 
    - Default: <expression.dataset_basename>.zip
- **create svgs**<span style="color:red;">*</span>
    - Whether to create SVG images (compressed) along with PNGs. Saving PNGs requires **a lot of storage**; therefore, this parameter is set to false by default.
- **selected gene sets**
    - Semicolon-separated list of gene sets from the provided gene sets database files (GMT/GMX/GRP). If you are using multiple files then you *must* prefix each selected gene set with its file name followed by '#' (like "my_file1.gmt#selected_gene_set1,my_file2.gmt#selected_gene_set2"). With a single file only the names are necessary. Leave this blank to select all gene sets.
- **alt delim**
    - Optional alternate delimiter character for gene set names instead of comma for use with selected.gene.sets. If used, a semicolon is recommended.

<span style="color:red;">*</span> = required

## Input Files

1. *ranked
list:* [RNK](http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29) file

This file contains the rank ordered gene (or feature) list.

2. *gene sets database
file:* [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29), [GMX](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29),
or [GRP](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GRP:_Gene_set_file_format_.28.2A.grp.29) file

Gene set files, either your own or from the listed MSigDB files.

3. *chip platform:* an
optional [CHIP](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CHIP:_Chip_file_format_.28.2A.chip.29)
file may be provided if you do not select a *chip platform* from the
drop-down

## Output Files

1. Enrichment Report archive: ZIP

ZIP file containing the result files. For more information on
interpreting these results, see [Interpreting GSEA
Results](http://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Interpreting_GSEA_Results)
in the GSEA User Guide. Note that in prior versions the ZIP bundle was
created as the only output file. This behavior has been changed to give
direct access to the results without the need for a download.

2. Enrichment Report: HTML and PNG images

The GSEA Enrichment Report. As above, see the GSEA User Guide for more
info.

3. Optional SVG images (compressed)

Identical to the PNGs in the Enrichment Report, but in SVG format for
higher resolution. These are GZ compressed to reduce space usage; they
can be decompressed using 'gunzip' on Mac or Linux and 7-Zip on Windows

## Platform Dependencies

**Task Type:**  
Gene List Selection

**CPU Type:**  
any

**Operating System:**  
any

**Language:**  
Java

## Version Comments

- **7.4.0** (2022-10-2): Updated to Human MSigDB v2022.1.Hs and Mouse MSigDB 2022.1.Mm. 
- **7.3.7** (2022-9-26): Fixed a manifest typo. Updated to Human MSigDB v2022.1.Hs. Direct support for Mouse MSigDB 2022.1.Mm is not yet available.  
- **7.3.6** (2022-9-15): Updated to Human MSigDB v2022.1.Hs. Direct support for Mouse MSigDB 2022.1.Mm is not yet available. 
- **7.3.5** (2022-3-22): Removed Log4J entirely from the code base. Fixed weighted_p1.5 computation. Added min dataset size warnings. 
- **7.3.4** (2022-1-20): Updated Log4J to 2.17.1. 
- **7.3.3** (2022-1-19): Updated to MSigDB v7.5.1. 
- **7.3.2** (2022-1-12): Updated to MSigDB v7.5. 
- **7.3.1** (2021-12-23): Updated with the GSEA Desktop 4.2.1 code base. Updated to Log4J 2.17.0.
- **7.3.0** (2021-12-17): Updated with the GSEA Desktop 4.2.0 code base with numerous bug fixes. Adds the Abs_max_of_probes collapse mode. Fixes some issues handling datasets with missing values. Improved warnings and logging. Changed the FDR q-value scale on the NES vs Significance plot. Fixed bugs in weighted_p1.5 scoring.
- **7.2.4** (2021-4-22): Fixed minor typo.
- **7.2.3** (2021-4-2): Updated to MSigDB v7.4.
- **7.2.2** (2021-3-22): Updated to MSigDB v7.3.
- **7.2.1** (2020-10-27): Fixed a bug in the Collapse Sum mode. 
- **7.2.0** (2020-9-23): Updated to MSigDB v7.2. Updated to use dedicated Docker container.
- **7.1.0** (2020-7-30): Updated to use the GSEA v4.1.0 code base.
- **7.0.4** (2020-4-2): Updated to use the GSEA v4.0.3 code base. Updated to give access to MSigDB v7.1.
- **7.0.3** (2019-10-24): Updated to use the GSEA v4.0.2 code base. Updated to give access to MSigDB v7.0. OpenJDK 11 port. Java code moved into the GSEA Desktop code base.
- **6.0.12** (2019-10-10): Updated to use the GSEA v3.0 open-source code base. Updated to give access to MSigDB v6.2. Unified the Gene Set DB selector parameters and better downloading of MSigDB files. Added selected.gene.sets, alt.delim and create.svgs parameters. Better temp file clean-up and other internal code improvements.
- **5** (2017-05-18): Updated to give access to MSigDB v6.0. 
- **4** (2016-02-04): Updated to give access to MSigDB v5.1.
- **3** (2015-12-04): Updating the GSEA jar to deal with an issue with FTP access. Fixes an issue for GP@IU.
- **2** (2015-06-16): Updated for MSigDB v5.0 and hallmark gene sets support.
- **1** (2013-06-17): Initial Release.

Copyright © 2003-2022 Broad Institute, Inc., Massachusetts Institute of
Technology, and Regents of the University of California. All rights
reserved.

