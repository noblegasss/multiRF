#!/usr/bin/env Rscript

# Static codebase review generator for multiRF.
# Usage:
#   Rscript tools/codebase_review.R

args <- commandArgs(trailingOnly = TRUE)
pkg_root <- if (length(args) >= 1L) normalizePath(args[1], mustWork = FALSE) else normalizePath(".", mustWork = FALSE)
r_dir <- file.path(pkg_root, "R")
out_dir <- file.path(pkg_root, "review")

if (!dir.exists(r_dir)) {
  stop("Cannot find R directory: ", r_dir, call. = FALSE)
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

r_files <- sort(list.files(r_dir, pattern = "[.]R$", full.names = TRUE))
if (length(r_files) == 0L) {
  stop("No R files found under: ", r_dir, call. = FALSE)
}

trim <- function(x) gsub("^\\s+|\\s+$", "", x)

file_role <- function(path) {
  nm <- basename(path)
  if (nm == "mrf3_workflow.R") return("Main end-to-end workflow orchestrator.")
  if (nm == "mrf3_workflow_simple.R") return("User-facing convenience wrappers for full pipeline runs.")
  if (nm == "mrf3_workflow_methods.R") return("S3 methods for summaries/printing of workflow results.")
  if (nm == "mrf3_init.R") return("Forest initialization and data preparation for multi-omics modeling.")
  if (nm == "mrf3_reconstr.R") return("Reconstruction and fusion of model-specific weight/similarity matrices.")
  if (nm == "mrf3_cl.R") return("Main clustering entrypoint for shared/specific labels.")
  if (nm == "mrf3_cl_prox.R") return("Proximity-based clustering variants (standard and enhanced).")
  if (nm == "mrf3_cl-util.R") return("Low-level tree/leaf/proximity helpers used by clustering.")
  if (nm == "mrf3-util.R") return("Core tree-network, importance, and weighting utilities.")
  if (nm == "mrf3_vs.R") return("Variable-selection routines and threshold utilities.")
  if (nm == "pairwise_imd.R") return("Pairwise IMD computation and related plotting/print helpers.")
  if (nm == "shared_fraction.R") return("Shared fraction scoring utilities.")
  if (nm == "shared_specific_weights.R") return("Shared/specific weight decomposition helpers.")
  if (nm == "signal_clustering.R") return("Similarity construction and clustering helpers for shared/specific signals.")
  if (nm == "findConnection.R") return("Connection discovery among omics blocks and quality scoring.")
  if (nm == "filter_omics.R") return("Input data filtering/normalization and connection enumeration helpers.")
  if (nm == "fit_rfsrc.R") return("Wrappers around randomForestSRC fitting for single and multi-block settings.")
  if (nm == "tune_t.R") return("Top-v tuning (model/fused) and tuning diagnostics.")
  if (nm == "tune_k_clusters.R") return("Cluster-count tuning utilities (PAM/spectral/gap/silhouette).")
  if (nm == "plot_tSNE.R") return("Visualization helpers (t-SNE/network/heatmap/circos/UMAP).")
  if (nm == "mrf3_tsne.R") return("t-SNE embedding wrappers and optimization helpers.")
  if (nm == "cluster_metrics.R") return("Clustering metric implementations and pairwise metric matrices.")
  if (nm == "workflow_summary_tools.R") return("Reporting/tidy summaries/stability tools for workflow outputs.")
  if (nm == "weights_preproc.R") return("Weight matrix preprocessing and validation utilities.")
  if (nm == "connect_utils.R") return("Connection list normalization/validation helpers.")
  if (nm == "KM_plot.R") return("Kaplan-Meier plotting helper.")
  if (nm == "sim_intersim_shared_specific.R") return("Simulation data generation helpers.")
  if (nm == "em-util.R") return("EM helper routines.")
  if (nm == "RcppExports.R") return("Auto-generated Rcpp bridge wrappers.")
  if (nm == "data.R") return("Package data object documentation stubs.")
  if (nm == "globals.R") return("Global variable declarations for non-standard evaluation.")
  if (nm == "multiRF-package.R") return("Package-level roxygen metadata.")
  "General utilities and helpers."
}

extract_comment_summary <- function(lines, start_line) {
  if (start_line <= 1L) return("")
  i <- start_line - 1L
  block <- character(0)
  while (i >= 1L) {
    ln <- lines[i]
    if (grepl("^\\s*$", ln)) break
    if (!grepl("^\\s*#", ln)) break
    block <- c(ln, block)
    i <- i - 1L
  }
  if (length(block) == 0L) return("")

  clean <- gsub("^\\s*#+'?\\s?", "", block)
  clean <- trim(clean)
  clean <- clean[nzchar(clean)]
  clean <- clean[!grepl("^@", clean)]
  clean <- clean[!grepl("^-+$", clean)]
  if (length(clean) == 0L) return("")
  clean[1]
}

is_fn_assign <- function(expr) {
  if (!is.call(expr)) return(FALSE)
  op <- as.character(expr[[1]])
  if (!op %in% c("<-", "=")) return(FALSE)
  if (length(expr) < 3L) return(FALSE)
  lhs <- expr[[2]]
  rhs <- expr[[3]]
  is.symbol(lhs) && is.call(rhs) && identical(as.character(rhs[[1]]), "function")
}

call_names_from_body <- function(body_expr) {
  nm <- all.names(body_expr, functions = TRUE, unique = TRUE)
  nm <- nm[nchar(nm) > 0L]
  nm <- setdiff(nm, c("{", "(", "[", "[[", "$", "@", "if", "for", "while", "repeat", "function", "return"))
  unique(nm)
}

functions <- list()
parse_errors <- character(0)

for (f in r_files) {
  exprs <- tryCatch(parse(f, keep.source = TRUE), error = function(e) e)
  if (inherits(exprs, "error")) {
    parse_errors <- c(parse_errors, sprintf("%s: %s", basename(f), conditionMessage(exprs)))
    next
  }
  lines <- readLines(f, warn = FALSE)
  for (expr in as.list(exprs)) {
    if (!is_fn_assign(expr)) next
    fn_name <- as.character(expr[[2]])
    fn_obj <- expr[[3]]
    fn_args <- names(as.pairlist(fn_obj[[2]]))
    if (length(fn_args) == 1L && is.na(fn_args)) fn_args <- character(0)
    fn_body <- fn_obj[[3]]
    sr <- attr(expr, "srcref")
    line <- if (!is.null(sr)) as.integer(sr[[1]]) else NA_integer_
    summary <- extract_comment_summary(lines, if (is.na(line)) 1L else line)
    calls <- call_names_from_body(fn_body)
    functions[[length(functions) + 1L]] <- data.frame(
      file = basename(f),
      file_path = f,
      function_name = fn_name,
      line = line,
      n_args = length(fn_args),
      args = paste(fn_args, collapse = ", "),
      summary = summary,
      calls_raw = paste(calls, collapse = ", "),
      stringsAsFactors = FALSE
    )
  }
}

if (length(parse_errors) > 0L) {
  cat("Parse errors detected:\n", paste0(" - ", parse_errors, collapse = "\n"), "\n", sep = "")
}

if (length(functions) == 0L) {
  stop("No function definitions found.", call. = FALSE)
}

fn_df <- do.call(rbind, functions)
fn_df <- fn_df[order(fn_df$file, fn_df$line, fn_df$function_name), ]
rownames(fn_df) <- NULL

all_fn <- unique(fn_df$function_name)
dup_fn <- sort(unique(all_fn[duplicated(all_fn)]))

edges <- do.call(rbind, lapply(seq_len(nrow(fn_df)), function(i) {
  caller <- fn_df$function_name[i]
  file <- fn_df$file[i]
  raw <- fn_df$calls_raw[i]
  if (!nzchar(raw)) return(NULL)
  calls <- unique(trim(strsplit(raw, ",", fixed = TRUE)[[1]]))
  calls <- calls[nzchar(calls)]
  internal <- intersect(calls, all_fn)
  if (length(internal) == 0L) return(NULL)
  data.frame(
    caller = caller,
    callee = internal,
    caller_file = file,
    stringsAsFactors = FALSE
  )
}))
if (is.null(edges)) {
  edges <- data.frame(caller = character(0), callee = character(0), caller_file = character(0), stringsAsFactors = FALSE)
} else {
  edges <- unique(edges)
  edges <- edges[order(edges$caller, edges$callee), ]
}

incoming <- table(edges$callee)
outgoing <- table(edges$caller)

fn_df$internal_calls <- vapply(fn_df$function_name, function(nm) {
  if (nm %in% names(outgoing)) as.integer(outgoing[[nm]]) else 0L
}, integer(1))
fn_df$internal_called_by <- vapply(fn_df$function_name, function(nm) {
  if (nm %in% names(incoming)) as.integer(incoming[[nm]]) else 0L
}, integer(1))

orphan_fns <- fn_df$function_name[fn_df$internal_calls == 0L & fn_df$internal_called_by == 0L]

write.csv(fn_df[, c("file", "function_name", "line", "n_args", "args", "summary", "internal_calls", "internal_called_by")],
          file.path(out_dir, "function_index.csv"), row.names = FALSE)
write.csv(edges, file.path(out_dir, "function_dependencies.csv"), row.names = FALSE)

md_path <- file.path(pkg_root, "CODEBASE_REVIEW.md")
con <- file(md_path, open = "wt")
on.exit(close(con), add = TRUE)

cat("# multiRF Codebase Review\n\n", file = con)
cat("Generated by `tools/codebase_review.R`.\n\n", file = con)
cat("## Scope\n\n", file = con)
cat("- Package root: `", pkg_root, "`\n", sep = "", file = con)
cat("- R files scanned: ", length(r_files), "\n", sep = "", file = con)
cat("- Functions found: ", nrow(fn_df), "\n", sep = "", file = con)
cat("- Internal dependency edges: ", nrow(edges), "\n\n", sep = "", file = con)

cat("## Health Checks\n\n", file = con)
if (length(parse_errors) == 0L) {
  cat("- Parse check: PASS\n", file = con)
} else {
  cat("- Parse check: FAIL\n", file = con)
  for (e in parse_errors) cat("  - ", e, "\n", sep = "", file = con)
}
if (length(dup_fn) == 0L) {
  cat("- Duplicate function names: none\n", file = con)
} else {
  cat("- Duplicate function names: ", paste(dup_fn, collapse = ", "), "\n", sep = "", file = con)
}
cat("- Orphan functions (no internal callers/callees): ", if (length(orphan_fns) == 0L) "none" else paste(orphan_fns, collapse = ", "), "\n\n", sep = "", file = con)

cat("## File Summaries\n\n", file = con)
for (f in basename(r_files)) {
  sub <- fn_df[fn_df$file == f, , drop = FALSE]
  cat("### ", f, "\n\n", sep = "", file = con)
  cat("- Role: ", file_role(f), "\n", sep = "", file = con)
  cat("- Functions: ", nrow(sub), "\n\n", sep = "", file = con)
}

cat("## Function Summaries\n\n", file = con)
for (f in unique(fn_df$file)) {
  sub <- fn_df[fn_df$file == f, , drop = FALSE]
  cat("### ", f, "\n\n", sep = "", file = con)
  for (i in seq_len(nrow(sub))) {
    row <- sub[i, ]
    summary <- if (nzchar(row$summary)) row$summary else "No nearby comment summary found."
    deps <- edges$callee[edges$caller == row$function_name]
    deps_txt <- if (length(deps) == 0L) "none" else paste(sort(unique(deps)), collapse = ", ")
    cat("- `", row$function_name, "` (line ", row$line, "): ", summary, "\n", sep = "", file = con)
    cat("  - Args: ", ifelse(nzchar(row$args), row$args, "none"), "\n", sep = "", file = con)
    cat("  - Depends on (internal): ", deps_txt, "\n", sep = "", file = con)
  }
  cat("\n", file = con)
}

cat("## Dependency Edges\n\n", file = con)
if (nrow(edges) == 0L) {
  cat("No internal function-to-function edges detected.\n", file = con)
} else {
  cat("| Caller | Callee | Caller File |\n", file = con)
  cat("|---|---|---|\n", file = con)
  for (i in seq_len(nrow(edges))) {
    cat("| `", edges$caller[i], "` | `", edges$callee[i], "` | ", edges$caller_file[i], " |\n", sep = "", file = con)
  }
}

cat("\n## Re-run Instructions\n\n", file = con)
cat("```bash\n", file = con)
cat("cd ", pkg_root, "\n", sep = "", file = con)
cat("Rscript tools/codebase_review.R\n", file = con)
cat("```\n\n", file = con)
cat("Generated artifacts:\n", file = con)
cat("- `CODEBASE_REVIEW.md`\n", file = con)
cat("- `review/function_index.csv`\n", file = con)
cat("- `review/function_dependencies.csv`\n", file = con)

cat("Wrote:\n")
cat(" - ", md_path, "\n", sep = "")
cat(" - ", file.path(out_dir, "function_index.csv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "function_dependencies.csv"), "\n", sep = "")

