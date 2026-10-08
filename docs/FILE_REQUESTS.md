# Remaining source and release requirements

## Already available

The supplied `ivory_reanalysis(1).r` contains lineage labeling and downstream analysis code. Saved-object exports provide subset metadata, assays, identities, contrast candidates and complete DE tables. Supplementary tables S1–S8 and paper documents have been located; selected DE statistics match the saved results at publication precision. Ribosomal/mitochondrial lists and the reduced-ivory counting GTF are received.

The Pegasus miR intake supplies the mature sequence, compressed target UTRs, transcript list, command, raw predictions, gene summary and predicted-gene set. Clean parsing reproduces all 1,265 gene summaries; comparative annotation reproduces the Table S8 hit counts and conserved-candidate rules. These files do not need to be sent again.

## Needed for complete upstream reproduction

1. **Final object construction and execution records:** the code/settings and successful logs for the actual `Hm_merged`/combined object used in the paper, including corrected sample metadata, integration parameters, assay preparation and original package versions. The March integration tail is present but its upstream setup is commented; do not resend the downstream lineage script.
2. **Final sample/run mapping and contrast confirmation:** run accessions, exclusions and final sample sheet; confirmation of the exact executed DE cohorts/assays where saved subset labels differ from parent labels. Saved metadata is available, so only missing mappings or branch confirmation are needed.
3. **Reference and BAM provenance:** reference FASTA source/checksum, counting-GTF generation and STAR-index commands, successful alignment logs, actual CB-tagged BAM-producing command and final locus inputs. The supplied GTF and gene lists are already available.
4. **miR upstream preparation:** code that selected longest transcripts and generated 3′ UTR FASTA; exact source annotation/FASTA identifiers and checksums. Original parsing code is useful provenance but is no longer a blocker for reconstructing the supplied summaries.
5. **Comparative DE provenance:** original Bicyclus/Junonia input tables, source/version and filtering/orthology rules that supplied the Tian/Fandino annotations. S8 provides the resulting flags and conserved-candidate rule, not the primary analyses.
6. **Enrichment provenance:** original PANGEA/DirectGO or modPhEA job settings, background universe, annotation version and unfiltered results where absent. Existing plot inputs do not establish the generating statistical test.
7. **Release details:** code license, contributor credits, final paper citation and archived release DOI. Confirm final figure assemblies from the already located documents rather than uploading them again.

Nanopore processing was performed by a collaborator and is outside the user's code-intake scope. Credit and link that contribution when available; no reconstruction is requested from the user.

## Scientific descriptions requiring reconciliation

- Supplied RNAhybrid command/output cutoff: −25 kcal/mol; manuscript: −20 kcal/mol. Supplied job logs show activation errors and do not establish the successful environment/version.
- Nucleus-level `FindMarkers` code versus pseudobulk wording in manuscript descriptions.
- Historical integration settings versus the settings reported in Methods; distinct branches remain separate.
- Parent ivory labels versus saved contrast-subset assignments; the early WT/mutant source call uses all of `test` despite its object name.
- S8 KrT95D has no energy value although its archived prediction is −26.4 kcal/mol; originals remain unchanged.

These are specific provenance or editorial gaps, not requests to replace already supplied code or results.
