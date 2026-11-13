# Blind cross-tissue single-cell analysis of one publicly available HCA-donor
A reanalysis of five tissues from a single Tabula Sapiens donor (TSP14), run as a blind benchmark: the published cell-type annotations were sealed before any analysis, the analysis plan and the label mapping were committed before any hypothesis was tested, and the reference labels were unsealed only at the final step.

### 1. The biological problem
A T cell in the gut is not the same as a T cell in fat, even when both are memory CD4 T cells by lineage. The same is true of macrophages, endothelial cells and fibroblasts. Cells are shaped by the tissue they live in. However, how that happens, and how much of it is a shared "tissue environment" versus something each cell type does on its own, is not settled.

Two molecular conversations plausibly drive it: the signals cells send each other, and how antigen is presented locally. This project asks three questions about them in one individual, where genetics, age and medical history are held constant:

1. Do memory T cells group by the kind of tissue they inhabit? Specifically, do cells from barrier tissues (large and small intestine, continuously exposed to gut contents) resemble each other more than they resemble cells from non-barrier tissues (fat, skeletal muscle)?
2. Which cells produce the secreted signals in a tissue, and does the same cell type secrete differently depending on where it lives?
3. Which cells carry the antigen-presentation machinery, and does that vary by tissue?
Running alongside these is a methodological question that turned out to matter more than expected: do the standard thresholds and corrections in a single-cell pipeline survive contact with this data?

Running alongside these is a methodological question that turned out to matter more than expected: do the standard thresholds and corrections in a single-cell pipeline survive contact with this data?

Most single-cell reanalyses apply thresholds inherited from tutorials and validate the result against the labels that shipped with the data. Both habits hide errors. Here:

1. Every threshold is derived from the data and justified in the notebook that sets it. Where no natural cutoff exists, the value is reported per gene rather than applied as a filter.
2. Annotation is blind. The published labels were written to a separate file in the first notebook, and notebooks 2–6 assert they are absent from the object. A mapping from analysis labels to acceptable reference labels was committed before unsealing, so the agreement score could not be tuned afterwards.
3. Statistics match the design. One donor means cells are pseudoreplicates. All differential expression is pseudobulk with 10x channels as replicates; comparisons with a single channel on one side are reported as effect sizes with no p-values.

### 2. The data

| | |
|---|---|
| Source | Tabula Sapiens (Science 376, eabl4896, 2022), donor TSP14 |
| Tissues | Blood, fat, skeletal muscle, large intestine, small intestine |
| Sampling sites | Fat: mesenteric + subcutaneous. Muscle: abdomen + diaphragm (two lanes). Gut: proximal + distal for both intestines. Blood: six libraries, two of them 5′ chemistry |
| Channels | 15 10x libraries |
| Input | 47,923 cells × 43,981 genes, raw counts plus a decontX-corrected layer |
| After QC | 44,847 cells × 37,628 genes |
| Annotated | 44,319 cells → 73 cell types in 34 groups; 528 excluded as probable cell pairs |

Cells per tissue at input: blood 9,962, fat 7,670, large intestine 11,439, muscle 8,394,
small intestine 10,458.

**Used in the analysis:** donor, sex, tissue, sampling site, 10x channel, chemistry.
**Withheld until the final notebook:** all cell-type annotations, the authors' QC flags,
and their per-cell statistics; fifteen metadata columns in total.

No data are stored in this repository.

<img width="1347" height="976" alt="umap_tissue" src="https://github.com/user-attachments/assets/2ed7ec6f-1ee3-4c0c-8ea4-0b450a99e4ee" />

### 3. The pipeline
#### 3.1 Annotation and integrity 
**What it does.** Reads count matrices and metadata exported from R, transposes them into the cells × genes layout, assembles a single annotated object, and separates the published annotations into a sealed file that nothing downstream reads.

**Tools.** scanpy, anndata, scipy (Matrix Market I/O).

**Results.** Three checks that the data is what it claims to be, all passed:
1. The counts are integers, so this  is raw data and not something already normalised.
2. Every row is matched to its own metadata. Recomputing library size and detected genes per cell and correlating them against the supplied statistics gives a correlation of exactly 1.0000, but only against the decontX layer. That revealed something not stated in the original metadata: the authors' per-cell statistics were computed on corrected counts, not raw ones.
3. The donor is male in all five tissues. XIST is detected in at most 1.2% of cells and the Y-linked genes are present everywhere, which also rules out a sample swap between tissues.



#### 3.2 Quality control from the data, not from convention
**What it does.** Computes per-cell quality metrics, inspects their distributions per tissue and per channel, sets every threshold from what those distributions show, scores dissociation stress, detects doublets, and removes cells only after a provisional clustering makes within-cluster comparisons possible.

**Tools.** scanpy, Scrublet, scipy (median absolute deviation), the Leiden algorithm, UMAP.

**Results**
1. The matrix was already pre-filtered at roughly 2,500 UMIs per cell. Every tissue has almost the same floor (1st percentile 2,621–3,188; minima 2,524–2,655), which a natural distribution would not produce. An additional count threshold would therefore do nothing.
2. Red blood cells were found rather than assumed. Blood splits into nucleated cells (≤17% haemoglobin) and a haemoglobin-dominated population (40–97%), with a nearly empty valley between them. The threshold was set in that valley at 30%, flagging 1,639 cells; the choice is robust, since 20%, 30% and 40% differ by only ten cells. Those cells then separated into two groups by sequencing chemistry, not biology. The same cell type measured two ways.
3. No fixed mitochondrial cutoff was applied, and this was the most consequential decision. 43% of large-intestine and 25% of small-intestine cells exceed 20% mitochondrial content. They are not damaged: they have more detected genes than unflagged cells (1,605 vs 1,533 in colon), two to three times more counts, and not one of them has fewer than 1,000 non-mitochondrial counts. Marker genes later identified them as intestinal epithelium. A conventional cutoff would have deleted nearly half of one tissue.
4. Scrublet's automatic thresholds were rejected. The simulated-score distributions were not bimodal, and in four channels the cutoff fell beyond nearly every simulated doublet, calling essentially none. A single threshold of 0.25 was set from the distributions, flagging 1.44% of cells.
5. Dissociation stress is kept as a covariate rather than removed, because it is confounded with tissue: muscle 1.75 > fat 1.09 > large intestine 0.78 > small intestine 0.63 > blood 0.47, with tissue explaining 36% of the score's variance.
6. 3,076 cells removed (6.42%). Blood loses 17.7%, almost entirely red cells; the solid tissues lose 3.0–4.1% each, so their relative composition is preserved.

<img width="4136" height="490" alt="blood_low_complexity_bands" src="https://github.com/user-attachments/assets/50a8aba1-934e-40f0-abb2-301926a64cfa" />


#### 3.3 Clustering and blind annotation 
**What it does.** Rebuilds the embedding from scratch on the cleaned cells, decides which count layer to cluster on, corrects one technical effect in one place, assigns lineages, re-clusters each lineage separately, annotates every cluster from marker genes, and cross-checks the immune annotation against an external classifier.

**Tools.** scanpy, the Leiden algorithm, UMAP, Harmony (chemistry only), CellTypist (Immune_All_Low), scikit-learn (nearest-neighbour reassignment).

**Results**
1. Raw counts were kept over the supplied decontX layer. Clustering on corrected counts dissolved an artefactual cluster, but it changed the global structure more than chance does (agreement of 0.73 against a seed-to-seed baseline of 0.89–0.91), including in tissues where ambient RNA was not a problem. decontX was also run by the original authors using their own cluster assignments, so clustering on it risked importing the annotation under validation.
2. One batch correction, deliberately scoped. Blood myeloid cells were splitting by library chemistry. Because both chemistries were run on the same blood samples, this is purely technical and not confounded with tissue, so Harmony was applied within the immune lineage only. No cross-tissue integration was performed — tissue is the signal here, not a batch.
3. Clustering is hierarchical, because one resolution cannot serve everything. Immune subtypes differ subtly while epithelial and muscle cells differ enormously. Six lineages were assigned from canonical markers, then each was re-clustered with its own variable genes.
4. Three populations turned out not to be what they first appeared. About 460 immune cells had been pulled into the epithelial lineage by the epithelial RNA floating in the gut samples; they were recovered by matching each one to its nearest immune neighbours in a gene space that excludes those contaminating genes. Nearly 400 further cells carry both an endothelial and a pericyte program in the same droplet — pericytes still attached to the capillary they wrap around — which the doublet detector missed because such pairs are common enough to look like a genuine population. A smaller group of probable immune–epithelial pairs was excluded on the same logic.
5. The external classifier agreed on most populations and was useful where it disagreed. It correctly flagged a regulatory T-cell population that had been missed, which marker genes then confirmed and which was split out — recorded as prompted by the classifier. It also called every neutrophil state a classical monocyte with high confidence; that call was rejected on evidence, since those cells carry neutrophil-specific genes in 83–99% of cases and essentially no MHC class II, while monocytes carry it in 58%.
6. The outcome is 73 cell types in 34 groups, every one with its marker evidence recorded.

<img width="498" height="520" alt="umap_final_cell_groups" src="https://github.com/user-attachments/assets/04e2a90e-9161-4a13-bbf5-e5aafa8ee9e9" />

<img width="1742" height="1182" alt="dotplot__lineage_markers" src="https://github.com/user-attachments/assets/c2fc2f61-f13a-48db-a07d-42addbcd5463" />


**UMAP for immune cells**

<img width="1134" height="1048" alt="umap_immune_harmony" src="https://github.com/user-attachments/assets/848ad239-2893-473d-9c82-3c1ea13fab37" />


#### 3.4 Feasibility analysis
**What it does.** Counts cells per cell type per channel and fixes the gene universe for the secretome analysis. It produces no results.

**Tools.** scanpy, pandas, the UniProt REST interface.

**Results** The decision, made before any test, of which comparisons have replication on both sides and which can only be descriptive. Memory CD4 T cells clear the minimum in two channels in all four solid tissues. Fat versus muscle is testable with channel replication for endothelial cells, macrophages and CD4 T cells, but not for fibroblasts or pericytes. Eight cell groups are present in both fat depots. And 1,799 of 2,144 annotated human secreted proteins are present in the data.

#### 3.5 T-cell adaptation and MHC survey 
**What it does.** Tests whether memory T cells group by tissue class, finds the genes behind it, controls for the alternative explanations, and surveys antigen-presentation genes across all cell types.

**Tools.** scanpy, PyDESeq2, scipy (hierarchical clustering), numpy (permutation testing).

**Results**
1. Blood was dropped as a baseline from the data, not from preference. Only 10 memory CD4 T cells were captured in blood; that compartment consists almost entirely of naive and regulatory T cells. A blood-versus-tissue contrast would have compared naive with memory cells and "confirmed" residency for a trivial reason, since several residency genes are also naive-cell genes.
2. Two analyses had to be rebuilt. Selecting variable genes the usual way gave unstable results — the two muscle samples correlated −0.023 at full depth but +0.270 at matched depth, opposite conclusions from nearly the same cells because selection was dominated by the much larger gut samples. And correlations tracked cell number almost perfectly, so depth was masquerading as biology. Fixing the gene set to genes detected in every sample, and subsampling all samples to the size of the smallest, resolved both.
3. Tissue class structures the transcriptome. The samples separate into a barrier block and a non-barrier block, with a within-class minus between-class correlation gap of 0.079–0.108. A permutation test gives p = 0.028 in every variant — the floor of this design, since eight samples admit only 35 possible four-against-four splits. The barrier split is the best of all 35. It survives depth matching, gene-set choice, ambient filtering, and restriction to the single memory subset present in all four tissues.
4. Anatomical proximity does not make mesenteric fat gut-like. Mesenteric fat sits against the intestine, so its T cells were predicted to resemble gut T cells more than subcutaneous ones do. They do not: mean correlation with gut is 0.852 versus 0.847, and the direction flips depending on the analysis. A clean null, reproduced three ways.
5. The first differential expression result was entirely ambient RNA and the top genes were epithelial and stromal markers, not T-cell biology.
6. After filtering, the class difference is a difference in residency depth. Fat and muscle T cells retain KLF2 and S1PR1, the genes a cell needs in order to leave a tissue. Gut T cells have switched those off and carry the homing receptors instead — CCR9, the small-intestine-homing receptor whose ligand is made by small-intestine epithelium, and GPR15, the colon-tropic one. 
7. Residency is graded rather than binary: small intestine > large intestine > fat > muscle. The ordering holds only when CD69 is excluded, because CD69 is both a residency marker and a gene induced by dissociation, and fat and muscle were the most stressed tissues. Comparing only cells matched for stress level, the ordering is identical in every one of six stress bins and the gap widens: stress was masking the difference, not creating it.
8. These are tissue cells, not blood cells trapped in vessels. Fat and muscle are highly vascularised and were not perfused, so circulating T cells would carry exactly this signature. But they have lost SELL and CCR7, which circulating cells retain, while expressing CXCR6 and ITGAE at gut-like levels.
9. The finding replicates independently in CD8 T cells, where stress happens to be equal across the three usable tissues, so no correction is needed at all. Two residency flavours emerge: intraepithelial cells in gut, and a GZMK⁺ population in fat.
10. Residency genes are core, homing genes are tissue-specific. Testing each tissue against the others, CD69, ITGAE, CXCR6 and PDCD1 are invariant across all four. What varies is the egress module and the homing receptors.
11. Antigen presentation is dominated by cell type, modulated by tissue. Class II scores span a wide range between cell types but at most a narrow range within a cell type across tissues. The notable result is that vascular endothelium in fat and muscle expresses the complete class II module above the resident macrophages beside it — the transcription factor CIITA included. That low macrophage signal independently matches the LYVE1⁺ perivascular subtype the blind annotation had assigned from unrelated markers.

<img width="1016" height="567" alt="cd4_sample_dendrogram" src="https://github.com/user-attachments/assets/ae3523f1-f387-46cf-b103-b7627ac9b5b3" />


#### 3.6 Secretome analysis
**What it does.** Takes every annotated secreted protein as the search space, maps each one to a producing cell lineage, compares the same cell types between fat and muscle, and compares the two fat depots.

**Tools.** scanpy, PyDESeq2, UniProt secreted-protein annotation.

**Results**
1. Roughly half of the expressed secreted genes can be attributed to a source lineage in each tissue, with positive controls landing correctly: complement from myeloid cells, von Willebrand factor from endothelium, granzymes from lymphocytes, decorin from stroma.
2. The first specificity metric was wrong and had to be replaced. Comparing a gene's top cell type against the runner-up penalised anything shared by two related populations: lymphocyte genes scored as non-specific purely because both CD4 and CD8 T cells were present. Comparing against the best cell type from a different lineage fixed the controls and raised the attributable set by half.
3. Two of my own hypotheses were tested and rejected. Specificity does not fall with abundance, and a lineage's secretome does not become harder to attribute the more that lineage dominates the sample. What actually governs attributability is how distinctive a lineage's secretory program is: immune lineages score highest because chemokines and granzymes are made by nothing else, while stromal, mural and vascular cells all produce overlapping matrix proteins.
4. Tissue adaptation is cell-type-specific, not tissue-wide. Of all the secreted genes differing between fat and muscle in endothelial cells, macrophages and CD4 T cells, only two are shared by all three — and both are heat-shock genes. None is shared in the fat direction. There is no common tissue-environment signature; the only shared component is technical, and the analysis caught it.
5. Stress affects only one direction. Stress and digestion genes make up a sizeable fraction of every muscle-up set and none at all of the three fat-up sets. An initial stress control was circular by construction — the gene sets were defined as differing between tissues, and stress also differs between tissues — and recomputing it within each tissue separately broke the circularity and reversed the apparent result for two of the three cell types.
6. Fat macrophages carry a coherent lipid-handling secretome containing no stress-flagged genes, matching the lipid-associated subtype annotated blind two notebooks earlier.
7. The two fat depots give the opposite answer to the two tissues. They differ by a set of developmental positional genes shared across nearly every cell type — subcutaneous and mesenteric markers that fell out of a comparison using nothing but the depot label parsed from a file name. The program is confined to the mesenchymal compartment and is never highest in an immune cell type, so the signal in immune cells is spill from the mesenchymal producers.
8. Cell-type-specific depot differences are not supported. The number of apparently depot-specific genes scales almost perfectly inversely with how many cells that type has, so it tracks sampling noise rather than biology

**Secretome: ambient contribution per gene, per library**
<img width="1615" height="490" alt="secretome_margin" src="https://github.com/user-attachments/assets/b5e7d6a8-132a-4545-b31d-0c3d658f7c2f" />

**Mesenchymal correlation**
<img width="2128" height="1397" alt="h3_mesenchymal_correlation" src="https://github.com/user-attachments/assets/734e46e9-55d9-40fc-bfbb-e95d074138a7" />


#### 3.7 Validation
**What it does.** Opens the withheld annotations for the first time and compares them to the blind labels using the metric and mapping fixed in advance.

**Tools.** pandas, scikit-learn (adjusted Rand index, normalised mutual information), seaborn.

**Results** 91.2% of evaluable cells carry a reference label inside the allowed set for their analysis label. Averaged across cell groups rather than cells the figure is lower, where most cells agree, while several small groups disagree systematically, and both are reported. The annotation step adds real information: the annotated cell groups agree with the reference substantially better than the raw clusters they were built from.

### 4. The ambient RNA problem
1. Unfiltered, the first cross-tissue T-cell contrast returned epithelial and stromal genes. Not T-cell biology, but the dominant ambient transcripts of each tissue. Every droplet receives a little of the free-floating RNA from cells that broke open during dissociation, and each tissue's soup has its own composition, so a cross-tissue comparison can silently become a comparison of which cell type dominates each sample.
2. The supplied decontX correction is insufficient here. It removes only 0.5–2.7% of counts per sample, and the worst epithelial contaminant still showed a fold change of roughly 2⁻⁸ in a T-cell contrast after correction.
3. A per-library "foreign margin" was built instead: a gene's level in the target cells minus its level in the other cells that shared the same sequencing channel. Threshold sensitivity showed a smooth continuum with no natural cutoff, so the margin is reported per gene rather than applied as a hard filter, and the interpretable set is the intersection of both defenses: 748 of 1,072 differentially expressed genes. The two approaches disagree on 161 genes, so neither is sufficient alone.
4. Ambient RNA from a cell type that is present in the data looks concentrated in that cell type; ambient from a cell type that was never captured looks uniform across everything. In these samples the uniform genes are immunoglobulins, from plasma cells, and haemoglobin, from red cells: neither of which forms a profile in fat. That distinction is what separates contamination from genes that are genuinely expressed by many cell types.

### 5. New relative to Tabula Sapiens paper
TSP14 is part of the published dataset. Where this analysis agrees with theirs, that is method concordance on the same cells, not independent replication.

**Recovered independently (their findings, reached blind and by a different route):** CCR9 and GPR15 as the region-specific gut T-cell trafficking receptors; endothelial cells clustering by tissue with fat and muscle most distinct; macrophages shared across tissues with subtle differences including epiregulin; the intestinal proliferative hierarchy; and IgA plasma cells concentrated in the gut.

**Absent from the paper, and new here:**
1. Adipose depot identity. The paper treats fat as one tissue. Separating mesenteric from subcutaneous recovers a developmental positional program, confined to the mesenchymal compartment, in two depots of the same person and refutes the specific prediction that gut-adjacent fat hosts more gut-like T cells.
2. A residency gradient across four tissues. The paper compares T cells between the two intestines only. Extending to fat and muscle shows residency is graded rather than binary, that fat and muscle T cells retain the egress machinery gut cells silence, and that the same pattern appears independently in CD8 T cells.
3. The secretome as a category. Which lineage produces which secreted protein, and the finding that tissue adaptation is cell-type-specific with no shared tissue signature.
4. MHC class II distribution across cell types, including vascular endothelium out-expressing the macrophages beside it.
5. The ambient RNA findings in section 4, which bear on how such atlases are processed rather than on this donor specifically.

### 6. Validation of labels
Where the marker evidence favours the blind call. The reference labels several hundred BEST4⁺ colonic cells as Paneth cells. The genes that define Paneth cells, the defensins, are undetectable across every colonic sub-cluster in this data, while those cells carry a coherent BEST4 program; Paneth cells are normally restricted to the small intestine, though they can appear metaplastically in proximal colon. Separately, a small group of genuine Paneth cells was split out of a small-intestine goblet cluster on exactly those defensin genes; the reference calls them goblet cells. Both disagreements are reported with their evidence, and neither label was changed after unsealing.

Where the difference is granularity, not error. Smooth muscle cells and pericytes map to the reference's broader fibroblast category; stem, transit-amplifying and enteroendocrine cells map to enterocyte labels. Neither vocabulary is wrong.

Where the analysis had already flagged the ambiguity. Platelets were labelled "platelet, neutrophil-associated" here because the platelet markers co-occurred with neutrophil genes; the reference calls them neutrophils. Both describe the same observation.

One disagreement is self-inflicted and was left in place. "Paneth cell" was deliberately not accepted as a reference equivalent for BEST4⁺ cells when the mapping was written. Allowing it would have raised the headline figure. The mapping was committed before unsealing and was not revised, which is the only reason the number means anything.

<img width="2573" height="1954" alt="validation_confusion_compact" src="https://github.com/user-attachments/assets/afd226f5-6a98-4959-b7c1-6d66a293a658" />


### 7. Conclusions 
1. Memory T cells are organised by the kind of tissue they inhabit, not merely by tissue identity. Barrier and non-barrier tissues separate at the statistical floor this design permits, and the genes behind the split are about tissue retention and homing.
2. Residency is a gradient, not a switch. Fat and muscle T cells have entered and adapted to their tissues. They have lost the lymph-node homing receptors and gained residency integrins and, yet, they retain the machinery for leaving, which gut T cells have shut down. Two independent T-cell lineages show this.
3. Tissues and depots differ in different ways. Different tissues impose different programs on each resident cell type, with essentially nothing shared between cell types. Two depots of the same tissue differ instead by one program shared across cell types and confined to a single compartment. The same pipeline produced both answers.
4. Anatomical proximity does not predict immunological similarity. A specific, pre-registered prediction about mesenteric fat was tested and refuted three ways.
5. Vascular endothelium is a stronger antigen-presentation surface than the macrophages beside it, in both fat and muscle. The result that also independently confirmed a macrophage subtype assigned blind from unrelated markers.
6. Published corrections need testing, not trusting. The ambient-RNA correction supplied with this dataset removes a small fraction of counts and leaves contamination large enough to dominate a cross-tissue comparison. A simple comparison against the cells sharing each library caught what it left behind.
7. The blind design worked. Annotation performed without the published labels agreed with expert annotation for 91.2% of cells, and the disagreements that remain are interpretable: two of them cases where the marker evidence favours the blind call.

### 8. Limitations
**Design**
1. One donor. Replicates are sampling sites and lanes within one person, so every result describes this donor's tissues and none of it generalises to humans.
2. Tissue is confounded with tissue processing. Each tissue was dissociated differently, so part of any difference is handling. Stress was scored, excluded from clustering, flagged in results and controlled by stress-matching — which bounds the problem but does not remove it.
3. Agreement with the source paper is method concordance on the same cells, since TSP14 is in that dataset and one of its analyses uses this exact donor.

**Statistics**
1. Eight pseudobulk samples give a permutation floor of p = 0.028, so the main structural test is at its ceiling by construction.
2. Only three cell types had channel replication on both sides of the fat–muscle comparison, and the muscle replicates include two lanes of one diaphragm sample, which measures lane noise rather than handling variability.
3. One library per adipose depot, so depot and library effects cannot be separated by replication; the internal control across shared cell types is a substitute, not a solution.
4. Muscle is the weakest arm: few memory CD4 T cells, the highest stress score, and the only tissue where residency correlates positively with stress.
5. Stress-matched comparisons use only the upper half of the stress range, where all four tissues overlap, so they say nothing about low-stress cells

**Scope**
1. The cells that define each tissue's secretome are absent. Adipocytes are floated off with the stromal vascular fraction and myofibres are too large for droplets, so the secretome described here is the stromal-vascular one. mRNA also indicates secretory potential, not secreted protein.
2. Composition is capture composition. Gut macrophages, dendritic cells and blood lymphocytes are all systematically under-represented by this protocol.
3. The published annotations were inspected once while assembling the data: category names, not cell-level labels. So blinding was procedural rather than absolute. This is disclosed rather than claimed otherwise.
4. Not yet replicated in an external multi-donor atlas

### 9. References
**Data**
The Tabula Sapiens Consortium. The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans. Science 376, eabl4896 (2022). doi:10.1126/science.abl4896

UniProt Consortium. UniProt: the Universal Protein Knowledgebase. Secreted-protein annotation retrieved via keyword KW-0964 (reviewed human entries), accessed 2026-10-11. https://www.uniprot.org

**Tools**
1. Wolf, F.A., Angerer, P. & Theis, F.J. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology 19, 15 (2018). doi:10.1186/s13059-017-1382-0
2. Virshup, I. et al. anndata: Annotated data. (AnnData / scverse).
3. Wolock, S.L., Lopez, R. & Klein, A.M. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell Systems 8, 281–291 (2019). doi:10.1016/j.cels.2018.11.005
4. Traag, V.A., Waltman, L. & van Eck, N.J. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports 9, 5233 (2019). doi:10.1038/s41598-019-41695-z
5. McInnes, L., Healy, J. & Melville, J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv:1802.03426 (2018).
6. Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nature Methods 16, 1289–1296 (2019). doi:10.1038/s41592-019-0619-0
7. Domínguez Conde, C. et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science 376, eabl5197 (2022). doi:10.1126/science.abl5197 (CellTypist)
8. Muzellec, B., Teleńczuk, M., Cabeli, V. & Andreux, M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics 39, btad547 (2023). doi:10.1093/bioinformatics/btad547
9. Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014). doi:10.1186/s13059-014-0550-8
10. Pedregosa, F. et al. Scikit-learn: Machine Learning in Python. Journal of Machine Learning Research 12, 2825–2830 (2011).
11. Yang, S. et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biology 21, 57 (2020). doi:10.1186/s13059-020-1950-6 (Applied to the data by the original authors; evaluated and found insufficient here — see section 4.)
12. IPD-IMGT/HLA Database. Allele statistics. https://www.ebi.ac.uk/ipd/imgt/hla/about/statistics

**Gene sets**
1. van den Brink, S.C. et al. Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations. Nature Methods 14, 935–936 (2017). doi:10.1038/nmeth.4437 (dissociation-stress gene set, human orthologs)
2. Kumar, B.V. et al. Human Tissue-Resident Memory T Cells Are Defined by Core Transcriptional and Functional Signatures. Cell Reports 20, 2921–2934 (2017). doi:10.1016/j.celrep.2017.08.078 (source for the core tissue-residency markers; the gene set used here was assembled from published TRM markers rather than copied verbatim from a supplementary table)
3. Marker genes used for annotation are canonical and are recorded per cluster with their evidence alongside the labels.

**Data access.** Tabula Sapiens is public. No raw or processed data are redistributed here; see the Tabula Sapiens portal for access terms.
