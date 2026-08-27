script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
script_path <- if (length(script_argument) == 1L) {
  normalizePath(
    sub("^--file=", "", script_argument),
    mustWork = TRUE
  )
} else {
  normalizePath(
    file.path(
      "maze_preprocessing",
      "generate_structural_audit_report.R"
    ),
    mustWork = TRUE
  )
}
report_dir <- dirname(script_path)

read_audit_csv <- function(name) {
  read.csv(
    file.path(report_dir, name),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = ""
  )
}

html_escape <- function(value) {
  value <- as.character(value)
  value[is.na(value)] <- ""
  value <- gsub("&", "&amp;", value, fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  value <- gsub('"', "&quot;", value, fixed = TRUE)
  value
}

html_table <- function(data, caption) {
  header <- paste0(
    "<tr>",
    paste0("<th>", html_escape(names(data)), "</th>", collapse = ""),
    "</tr>"
  )
  rows <- apply(
    data,
    1L,
    function(row) {
      paste0(
        "<tr>",
        paste0("<td>", html_escape(row), "</td>", collapse = ""),
        "</tr>"
      )
    }
  )
  paste0(
    "<figure><figcaption>",
    html_escape(caption),
    "</figcaption><div class=\"table-wrap\"><table><thead>",
    header,
    "</thead><tbody>",
    paste0(rows, collapse = "\n"),
    "</tbody></table></div></figure>"
  )
}

source_manifest <- read_audit_csv("source_data_manifest.csv")
schema_manifest <- read_audit_csv("source_schema_manifest.csv")
session_audit <- read_audit_csv("analysis_session_audit.csv")
boundary_provenance <- read_audit_csv(
  "session_boundary_provenance.csv"
)
no_comparator <- read_audit_csv(
  "no_comparator_sensitivity.csv"
)
session_sensitivity <- read_audit_csv(
  "session_reconstruction_sensitivity.csv"
)
duplicate_provenance <- read_audit_csv(
  "duplicate_trajectory_provenance.csv"
)

duplicate_schema <- schema_manifest[
  schema_manifest$RawHeaderOccurrences > 1L,
  c(
    "Dataset",
    "ColumnPosition",
    "RawHeader",
    "NormalizedHeader",
    "RawHeaderOccurrences"
  ),
  drop = FALSE
]

stopifnot(
  nrow(source_manifest) == 19L,
  nrow(no_comparator) == 19L,
  sum(no_comparator$PrimaryNoComparatorRemoved) == 10047L,
  all(no_comparator$RetentionNoComparatorRemoved == 0L),
  all(session_audit$MetadataValidation == "NOT_AVAILABLE"),
  nrow(session_sensitivity) == 1L,
  nrow(duplicate_provenance) > 0L,
  !anyNA(duplicate_schema$RawHeaderOccurrences)
)

format_integer <- function(value) {
  format(value, big.mark = ",", scientific = FALSE, trim = TRUE)
}

summary_cards <- paste0(
  "<div class=\"cards\">",
  "<article><strong>19</strong><span>source CSVs checksummed</span></article>",
  "<article><strong>",
  format_integer(sum(no_comparator$PrimaryNoComparatorRemoved)),
  "</strong><span>structural removals in the primary policy</span></article>",
  "<article><strong>",
  format_integer(sum(no_comparator$FinalRTDifference)),
  "</strong><span>additional final RTs under retention</span></article>",
  "<article><strong>",
  round(max(no_comparator$MaximumCleanModelStandardizedShift), 3),
  "</strong><span>largest standardized diagnostic-model shift</span></article>",
  "</div>"
)

interpretation <- paste0(
  "<section class=\"finding warning\"><h2>What the sensitivity establishes</h2>",
  "<p>Retaining structurally non-comparable rows left the adaptive multiplier ",
  "unchanged in all datasets, but changed the Box–Cox lambda in ",
  sum(no_comparator$LambdaDifference != 0),
  " datasets and changed at least one non-intercept clean-RT diagnostic-model ",
  "coefficient sign in ",
  sum(no_comparator$CleanModelCoefficientSignChanges > 0),
  " datasets. The structural policy therefore cannot be described as ",
  "universally neutral. Final study-specific inferential models must carry ",
  "the same sensitivity comparison.</p></section>",
  "<section class=\"finding\"><h2>Reconstructed boundary result</h2>",
  "<p>No source file contains explicit run/session metadata, so every metadata ",
  "status is <code>NOT_AVAILABLE</code>. Orth A-Maze is the only dataset ",
  "requiring reconstruction after exact-copy removal. Leaving its source ",
  "label unsplit changed the final RT count by ",
  abs(
    session_sensitivity$UnsplitFinalRTs -
      session_sensitivity$ReconstructedFinalRTs
  ),
  ", changed no diagnostic-model coefficient signs, and had a maximum ",
  "standardized shift of ",
  round(
    session_sensitivity$MaximumCleanModelStandardizedShift,
    3
  ),
  ".</p></section>"
)

document <- paste0(
  "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
  "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
  "<title>Structural RT audit</title><style>",
  "body{font:15px/1.55 system-ui,sans-serif;margin:0;background:#f6f7f9;",
  "color:#18212b}main{max-width:1180px;margin:auto;padding:36px 22px 72px}",
  "h1{font-size:2.15rem;margin:.2rem 0}.lede{max-width:850px;color:#4b5563}",
  ".cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));",
  "gap:14px;margin:28px 0}.cards article,.finding,figure{background:white;",
  "border:1px solid #d9dee5;border-radius:10px;padding:18px}",
  ".cards strong{display:block;font-size:1.7rem}.cards span{color:#53606e}",
  ".finding{margin:18px 0}.warning{border-left:5px solid #a96300}",
  "h2{margin-top:0}code{background:#eef1f4;padding:2px 5px;border-radius:4px}",
  "figure{margin:24px 0}figcaption{font-weight:700;margin-bottom:12px}",
  ".table-wrap{overflow:auto;max-height:560px}table{border-collapse:collapse;",
  "font-size:12px;width:max-content;min-width:100%}th,td{border:1px solid ",
  "#d9dee5;padding:6px 8px;text-align:left;white-space:nowrap}th{background:",
  "#eef1f4;position:sticky;top:0}a{color:#1459a6}</style></head><body><main>",
  "<p>Reproducible analysis supplement</p><h1>Structural removals and ",
  "reconstructed-session audit</h1><p class=\"lede\">This report separates ",
  "structural no-comparator missingness from relative RT extremeness, validates ",
  "inferred boundaries where source metadata are unavailable, records ",
  "de-identified duplicate provenance, and publishes source identity and schema ",
  "checks.</p>",
  summary_cards,
  interpretation,
  "<p>Machine-readable files: <a href=\"source_data_manifest.csv\">source ",
  "checksums</a>, <a href=\"source_schema_manifest.csv\">schema manifest</a>, ",
  "<a href=\"no_comparator_sensitivity.csv\">no-comparator sensitivity</a>, ",
  "<a href=\"session_reconstruction_sensitivity.csv\">session sensitivity</a>, ",
  "and <a href=\"duplicate_trajectory_provenance.csv\">duplicate provenance",
  "</a>.</p>",
  html_table(no_comparator, "No-comparator retention sensitivity"),
  html_table(
    session_sensitivity,
    "Reconstructed versus unsplit session sensitivity"
  ),
  html_table(session_audit, "Session metadata and construction audit"),
  html_table(boundary_provenance, "Inferred boundary provenance"),
  html_table(duplicate_provenance, "De-identified duplicate trajectories"),
  html_table(duplicate_schema, "Explicitly normalized duplicate raw headers"),
  html_table(source_manifest, "Source-data SHA-256 and import manifest"),
  "</main></body></html>"
)

writeLines(
  document,
  file.path(report_dir, "structural_audit.html"),
  useBytes = TRUE
)

cat(
  "Wrote structural_audit.html with",
  nrow(no_comparator),
  "dataset sensitivity rows.\n"
)