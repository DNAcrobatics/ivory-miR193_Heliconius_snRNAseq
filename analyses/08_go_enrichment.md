# Gene mapping, GO enrichment and plotting

**Paper:** Figure 7C/F.

**Status:** Plotting supplied; enrichment provenance needed.

Recover gene-to-FlyBase mapping, submitted gene lists and enrichment settings.

## Inputs

- Final DE tables, `HomologGeneNames.txt`, `fbgn_fbtr_fbpp_expanded_fb_2024_01.tsv`.
- UpGenes/DownGenes and late equivalents, including FBpp submission lists.
- Four DirectGoEnrichment_PANGEA CSV exports named in the plotting code.

## Outputs

- GO enrichment result tables and early/late PDF plots.

## Software

R and ggplot2 for plotting; enrichment service/database version and submission settings need recovery.

## Execution and interpretation

These blocks plot previously downloaded results; they do not execute the enrichment test. Recover the universe, submitted identifiers, term library, correction and export settings. Early filtering mixes Benjamini–Yekutieli and Benjamini–Hochberg, while later filters use BY. Reconcile this with the paper’s permutation/FDR description before choosing final code.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Code

### Early upregulated GO plot

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 1352–1370](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L1352-L1370).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
UpGenePANGEA <- read.csv('UpGenes_DirectGoEnrichment_PANGEA.csv',h=T)

UpGenePANGEA<- (subset(UpGenePANGEA, Benjamini...Yekutieli < 0.05))

a <- ggplot(UpGenePANGEA, aes(x = Count.Overlap.Gene, y = reorder(Gene.Set.Name, Fold.Enrichment), 
               size = Fold.Enrichment, color = Benjamini...Yekutieli)) +
  geom_point() +
#  facet_grid(Cluster~., scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  labs(title = "GO Enrichment by Cluster",
       x = "Count",
       y = "DirectGO Term",
       color = "Adjusted p-value",
       size = "Fold Enrichment") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))
```

### Early downregulated GO plot

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 1373–1396](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L1373-L1396).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
DownGenePANGEA <- read.csv('DownGenes_DirectGoEnrichment_PANGEA.csv',h=T)

DownGenePANGEA<- (subset(DownGenePANGEA, Benjamini...Hochberg < 0.05))

b <- ggplot(DownGenePANGEA, aes(x = Count.Overlap.Gene, y = reorder(Gene.Set.Name, Fold.Enrichment), 
               size = Fold.Enrichment, color = Benjamini...Hochberg)) +
  geom_point() +
#  facet_grid(Cluster~., scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  labs(title = "GO Enrichment by Cluster",
       x = "Count",
       y = "DirectGO Term",
       color = "Adjusted p-value",
       size = "Fold Enrichment") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))


plot_grid(a,b,ncol=1)
dev.copy(pdf, 'GoTermsOut.pdf',width=7,height=7)
dev.off()
```


[Return to analysis index](../README.md#analysis-index).
