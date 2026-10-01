# Cell-cycle scoring and phase plots

**Paper:** Figure 5E/F.

**Status:** Scoring code and gene mapping supplied.

Score S and G2/M phases with the supplied cross-species gene mapping.

## Inputs

- Annotated `Hm.sops.sct` and early/late objects `test` and `test1`.
- [CellCycleGenes.txt](../archive/original/ivory_proj_transfer_aug2025/CellCycleGenes.txt).

## Outputs

- S.Score, G2M.Score and Phase metadata; early/late phase proportion plots.

## Software

R, Seurat and ggplot2, with the scoring assay and package version still to confirm.

## Execution and interpretation

The mapping contains 66 rows: 29 S and 37 G2/M. Retain its original NO_HIT entries and duplicate mappings in the archive. Confirm how the final scoring run handled absent and repeated features; the snippet passes the HmID column directly. The phase plot fragments occur before scoring in the pasted session, so their dependencies must be restored before making a runnable workflow.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Code

### Recovered scoring

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 254–261](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L254-L261).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
### cell cycle scoring
### import list of canonical cell cycle markers
ccgenes <- read.table('CellCycleGenes.txt',h=T)
### make objects for each class
s.genes <- ccgenes$HmID[ccgenes$stage == 'S']
g2m.genes <- ccgenes$HmID[ccgenes$stage == 'G2/M']
### assign each cell a cell cycle score
Hm.sops.sct <- CellCycleScoring(Hm.sops.sct, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
```

### Recovered late phase plot and PDF export

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 234–252](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L234-L252).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
PhasePlotClusterLate <- table(test1$mut[test1$celltypes =='scale'], test1@meta.data$Phase[test1$celltypes =='scale' ])
hsize <- 1
okabe5 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",'#D55E00')
PhaseClustersL <- ggplot(as.data.frame(PhasePlotClusterLate), aes(fill=Var2, x=hsize, y=Freq))+
           geom_bar(position='fill',stat='identity',)+
           coord_polar("y")+
           xlim(c(0.2,hsize+.5))+
           scale_fill_manual(values=okabe5[4:6])+
           theme_void()+
           theme(axis.text=element_text(colour='black'))+
           facet_grid(~Var1)


PhaseClustersE
dev.copy(pdf, 'PhaseClustersEarly.pdf',width=7,height=3)
dev.off()
PhaseClustersL
dev.copy(pdf, 'PhaseClustersLate.pdf',width=7,height=3)
dev.off()
```


[Return to analysis index](../README.md#analysis-index).
