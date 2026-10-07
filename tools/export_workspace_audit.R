# Export a compact audit from an ALREADY LOADED R session, avoiding another
# full copy of a large workspace in RAM. This does not recompute any analysis.
# Usage in the existing session:
# source("tools/export_workspace_audit.R")
# export_workspace_audit(outdir = "ivory_workspace_audit")
export_workspace_audit <- function(outdir, envir = .GlobalEnv, max_table_mb = 20) {
  if (!is.environment(envir)) stop("envir must be an environment")
  if (file.exists(outdir)) stop("Output directory already exists; choose a new path")
  if (!is.numeric(max_table_mb) || length(max_table_mb) != 1L || !is.finite(max_table_mb) || max_table_mb <= 0) stop("Invalid table size limit")
  dir.create(outdir, recursive = TRUE)
  objects <- ls(envir, all.names = TRUE)
  inventory <- list(); exports <- list(); assays <- list(); bindings <- list()
  record <- function(object, kind, file, status, detail = "") {
    exports[[length(exports) + 1L]] <<- data.frame(object, kind, file, status, detail, stringsAsFactors = FALSE)
  }
  # Explicit row names retain gene IDs / full sample-prefixed cell IDs.
  table_export <- function(x, object, kind, file) {
    if (!is.data.frame(x)) return(invisible(NULL))
    if (as.numeric(object.size(x)) > max_table_mb * 1024^2) {
      record(object, kind, file, "skipped_size", paste(nrow(x), "rows; raise max_table_mb if needed")); return(invisible(NULL))
    }
    tryCatch({
      id_column <- ".saved_rowname"
      while (id_column %in% colnames(x)) id_column <- paste0(id_column, "_")
      ids <- data.frame(rownames(x), stringsAsFactors = FALSE); names(ids) <- id_column
      con <- gzfile(file.path(outdir, file), "wt")
      tryCatch(write.table(cbind(ids, x), con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA"), finally = close(con))
      record(object, kind, file, "exported", paste(nrow(x), "rows; row identifier column", id_column))
    }, error = function(e) record(object, kind, file, "error", conditionMessage(e)))
  }
  for (i in seq_along(objects)) {
    name <- objects[i]
    # Active bindings can run code when read; inventory them without invoking.
    if (bindingIsActive(name, envir)) {
      bindings[[length(bindings) + 1L]] <- data.frame(object = name, status = "active_binding_skipped")
      next
    }
    x <- get(name, envir = envir, inherits = FALSE)
    id <- sprintf("object_%04d", i)
    dims <- tryCatch(dim(x), error = function(e) NULL)
    inventory[[length(inventory) + 1L]] <- data.frame(index = i, object = name, class = paste(class(x), collapse = ";"), dimensions = paste(dims, collapse = "x"), bytes = as.numeric(object.size(x)), stringsAsFactors = FALSE)
    if (is.data.frame(x)) table_export(x, name, "table", paste0(id, ".tsv.gz"))
    if (inherits(x, "Seurat") && isS4(x)) {
      sn <- slotNames(x)
      if ("meta.data" %in% sn) {
        md <- slot(x, "meta.data")
        if (is.data.frame(md) && "active.ident" %in% sn) {
          identities <- slot(x, "active.ident")
          column <- ".saved_active_ident"
          while (column %in% colnames(md)) column <- paste0(column, "_")
          md[[column]] <- as.character(identities[match(rownames(md), names(identities))])
        }
        table_export(md, name, "Seurat_metadata", paste0(id, "_metadata.tsv.gz"))
      }
      if ("assays" %in% sn) {
        for (an in names(slot(x, "assays"))) {
          ax <- slot(x, "assays")[[an]]
          layer_names <- if (isS4(ax) && "layers" %in% slotNames(ax)) names(slot(ax, "layers")) else if (isS4(ax)) intersect(c("counts", "data", "scale.data"), slotNames(ax)) else character()
          ad <- tryCatch(dim(ax), error = function(e) NULL)
          assays[[length(assays) + 1L]] <- data.frame(object = name, assay = an, class = paste(class(ax), collapse = ";"), dimensions = paste(ad, collapse = "x"), layers = paste(layer_names, collapse = ";"), active = "active.assay" %in% sn && identical(an, slot(x, "active.assay")), stringsAsFactors = FALSE)
        }
      }
      if ("commands" %in% sn) {
        commands <- slot(x, "commands")
        params <- lapply(commands, function(cmd) {
          if (!isS4(cmd)) return(list(class = class(cmd)))
          fields <- intersect(c("name", "time.stamp", "assay.used", "call.string", "params"), slotNames(cmd))
          setNames(lapply(fields, function(field) slot(cmd, field)), fields)
        })
        if (as.numeric(object.size(params)) <= max_table_mb * 1024^2) {
          file <- paste0(id, "_commands.R.txt")
          dput(params, file.path(outdir, file)); record(name, "Seurat_command_record", file, "exported", "Saved parameters; does not prove every historical operation was recorded")
        } else record(name, "Seurat_command_record", "", "skipped_size")
      }
    }
    # Saved ggplot input can identify the exact GO / dotplot rows displayed.
    if (inherits(x, "ggplot")) {
      plot_data <- tryCatch(x$data, error = function(e) NULL)
      table_export(plot_data, name, "ggplot_data", paste0(id, "_plot_data.tsv.gz"))
    }
    # Small gene/target vectors are useful; exclude matrices and functions.
    if (is.atomic(x) && is.null(dim(x)) && length(x) <= 100000L && as.numeric(object.size(x)) <= 1024^2) {
      file <- paste0(id, "_vector.R.txt")
      dput(x, file.path(outdir, file)); record(name, "atomic_vector", file, "exported")
    }
  }
  save_index <- function(rows, filename) {
    if (length(rows)) write.table(do.call(rbind, rows), file.path(outdir, filename), sep = "\t", quote = TRUE, row.names = FALSE)
  }
  save_index(inventory, "objects.tsv"); save_index(exports, "exports.tsv")
  save_index(assays, "assays.tsv"); save_index(bindings, "skipped_bindings.tsv")
  writeLines(capture.output(sessionInfo()), file.path(outdir, "inspection_sessionInfo.txt"))
  writeLines(c("Exported saved data without recomputing DE, assigning labels, joining assay layers or changing the source objects.", "sessionInfo describes this export session, not necessarily the session that produced the results.", "Metadata includes saved active identities; it does not establish earlier identities used by FindMarkers.", "No expression matrices or complete Seurat objects were exported."), file.path(outdir, "README.txt"))
  files <- list.files(outdir, full.names = TRUE)
  sums <- tools::md5sum(files)
  write.table(data.frame(file = basename(names(sums)), md5 = unname(sums)), file.path(outdir, "MD5SUMS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", normalizePath(outdir), "; inspect exports.tsv for skipped tables or errors")
  invisible(outdir)
}
