# Source-file import, schema, and session-boundary audit helpers.
#
# These helpers intentionally preserve the raw CSV header before assigning
# deterministic analysis names. They do not use read.csv() name repair or
# make.unique(), so duplicate and blank source headers remain auditable.

is_experimental_condition <- function(condition) {
  condition_key <- tolower(trimws(as.character(condition)))
  normalized_key <- gsub("[^[:alnum:]]+", "", condition_key)

  nonexperimental_labels <- c(
    "",
    "consent",
    "debrief",
    "demographics",
    "exit",
    "intro",
    "intro1",
    "intro2",
    "intro3",
    "postpractice",
    "questionnaire",
    "brexit",
    "form",
    "questionalt",
    "subhtmlflash"
  )

  is_missing_or_administrative <-
    is.na(normalized_key) |
    normalized_key %in% nonexperimental_labels

  is_nonexperimental_item <- grepl(
    "^(filler|practice|burnin|intro)",
    normalized_key
  )

  !(is_missing_or_administrative | is_nonexperimental_item)
}

normalize_maze_headers <- function(raw_names) {
  raw_names <- as.character(raw_names)
  raw_names[is.na(raw_names)] <- ""
  base_names <- trimws(raw_names)
  base_names[base_names == ""] <- "unnamed"
  base_names <- gsub(
    "[^[:alnum:]_.]+",
    "_",
    base_names
  )
  base_names <- gsub("^_+|_+$", "", base_names)
  base_names[base_names == ""] <- "column"
  starts_with_number <- grepl("^[0-9]", base_names)
  base_names[starts_with_number] <- paste0(
    "column_",
    base_names[starts_with_number]
  )

  occurrence <- ave(
    seq_along(base_names),
    base_names,
    FUN = seq_along
  )
  totals <- table(base_names)
  duplicated_base <- as.integer(totals[base_names]) > 1L
  normalized <- base_names
  normalized[duplicated_base] <- paste0(
    base_names[duplicated_base],
    "__",
    sprintf("%02d", occurrence[duplicated_base])
  )

  if (anyDuplicated(normalized)) {
    stop("Explicit CSV header normalization did not produce unique names.")
  }

  normalized
}

read_maze_csv <- function(path) {
  dat <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  raw_header <- names(dat)
  normalized_header <- normalize_maze_headers(raw_header)
  names(dat) <- normalized_header

  required <- c("Subject", "Condition", "Item", "Word", "RT")
  missing_required <- setdiff(required, names(dat))
  if (length(missing_required) > 0L) {
    stop(
      basename(path),
      " is missing required columns after explicit header normalization: ",
      paste(missing_required, collapse = ", ")
    )
  }

  attr(dat, "maze_raw_header") <- raw_header
  attr(dat, "maze_normalized_header") <- normalized_header
  dat
}

sha256_file <- function(path) {
  sha256sum <- Sys.which("sha256sum")
  if (!nzchar(sha256sum)) {
    stop("The `sha256sum` utility is required to build the source manifest.")
  }

  output <- system2(
    sha256sum,
    shQuote(normalizePath(path, mustWork = TRUE)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("SHA-256 calculation failed for ", path, ": ", paste(output, collapse = " "))
  }

  sub("[[:space:]].*$", "", output[[1L]])
}

sha256_text <- function(text) {
  temporary_file <- tempfile("maze-header-", fileext = ".txt")
  on.exit(unlink(temporary_file), add = TRUE)
  writeLines(text, temporary_file, useBytes = TRUE)
  sha256_file(temporary_file)
}

find_session_metadata_columns <- function(dat) {
  raw_header <- attr(dat, "maze_raw_header")
  if (is.null(raw_header)) {
    raw_header <- names(dat)
  }
  normalized_key <- tolower(
    gsub("[^[:alnum:]]+", "", trimws(as.character(raw_header)))
  )
  accepted_keys <- c(
    "session",
    "sessionid",
    "run",
    "runid",
    "attempt",
    "attemptid",
    "recording",
    "recordingid",
    "submission",
    "submissionid",
    "batch",
    "batchid",
    "starttime",
    "timestamp",
    "datetime",
    "date"
  )

  names(dat)[normalized_key %in% accepted_keys]
}

make_analysis_subject_key <- function(
  dat,
  dataset_name
) {
  subject <- as.character(dat$Subject)
  experimental <- is_experimental_condition(dat$Condition)
  valid_subject <-
    experimental &
    !is.na(subject) &
    nzchar(trimws(subject))
  subject_counts <- table(subject[valid_subject])

  if (length(subject_counts) == 0L) {
    stop(dataset_name, " has no experimental rows with a participant label.")
  }

  metadata_columns <- find_session_metadata_columns(dat)
  analysis_subject <- subject
  provenance <- data.frame(
    Dataset = character(),
    BoundaryAuditGroup = character(),
    SessionToken = character(),
    SourceRowFirst = integer(),
    SourceRowLast = integer(),
    ExperimentalRows = integer(),
    SourceRowsContiguous = logical(),
    ExperimentalSequenceContiguous = logical(),
    DesignSchemaCompatible = logical(),
    stringsAsFactors = FALSE
  )

  if (length(metadata_columns) > 0L) {
    metadata_values <- do.call(
      paste,
      c(
        lapply(
          dat[metadata_columns],
          function(values) trimws(as.character(values))
        ),
        sep = "\034"
      )
    )
    metadata_present <-
      !is.na(metadata_values) &
      nzchar(metadata_values)
    if (any(valid_subject & !metadata_present)) {
      stop(
        dataset_name,
        " has incomplete source session metadata in experimental rows."
      )
    }

    source_session <- interaction(
      subject,
      metadata_values,
      drop = TRUE,
      lex.order = TRUE
    )
    analysis_subject[valid_subject] <- as.character(
      source_session[valid_subject]
    )
    source_sessions_by_subject <- split(
      as.character(source_session[valid_subject]),
      subject[valid_subject]
    )
    sessions_per_subject <- vapply(
      source_sessions_by_subject,
      function(values) length(unique(values)),
      integer(1L)
    )
    additional_sessions <- sum(pmax(0L, sessions_per_subject - 1L))
    analysis_counts <- table(analysis_subject[valid_subject])

    return(
      list(
        key = analysis_subject,
        unsplit_key = subject,
        source_participants = length(subject_counts),
        analysis_participants = length(analysis_counts),
        modal_rows = as.integer(
          names(which.max(table(as.integer(analysis_counts))))
        ),
        additional_sessions = additional_sessions,
        session_construction = "source_metadata",
        metadata_columns = paste(metadata_columns, collapse = ";"),
        metadata_validation = "PASS",
        boundary_validation = "PASS",
        boundary_provenance = provenance
      )
    )
  }

  count_frequencies <- table(as.integer(subject_counts))
  modal_candidates <- as.integer(
    names(count_frequencies)[
      count_frequencies == max(count_frequencies)
    ]
  )
  if (length(modal_candidates) != 1L) {
    stop(
      dataset_name,
      " does not have a unique modal experimental-row count ",
      "for identifying repeated complete sessions."
    )
  }

  modal_rows <- modal_candidates[[1L]]
  additional_sessions <- 0L
  boundary_group <- 0L
  boundary_schema_checks <- logical()
  experimental_row_position <- match(
    seq_len(nrow(dat)),
    which(valid_subject)
  )
  design_schema_signature <- function(rows) {
    paste(
      sort(
        paste(
          as.character(dat$Condition[rows]),
          as.character(dat$Item[rows]),
          as.character(dat$Word[rows]),
          sep = "\034"
        )
      ),
      collapse = "\035"
    )
  }
  complete_subjects <- names(subject_counts)[
    as.integer(subject_counts) == modal_rows
  ]
  complete_design_schemas <- vapply(
    complete_subjects,
    function(subject_id) {
      design_schema_signature(
        which(valid_subject & subject == subject_id)
      )
    },
    character(1L)
  )

  for (subject_id in names(subject_counts)) {
    subject_row_count <- as.integer(subject_counts[[subject_id]])
    is_complete_repeat <-
      subject_row_count > modal_rows &&
      subject_row_count %% modal_rows == 0L

    if (!is_complete_repeat) {
      next
    }

    subject_rows <- which(valid_subject & subject == subject_id)
    session_number <-
      ((seq_along(subject_rows) - 1L) %/% modal_rows) + 1L
    analysis_subject[subject_rows] <- paste0(
      subject_id,
      "::session",
      session_number
    )
    additional_sessions <-
      additional_sessions +
      max(session_number) - 1L
    boundary_group <- boundary_group + 1L

    session_rows <- split(subject_rows, session_number)
    design_schemas <- vapply(
      session_rows,
      design_schema_signature,
      character(1L)
    )
    design_schema_compatible <- if (
      length(complete_design_schemas) > 0L
    ) {
      design_schemas %in% complete_design_schemas
    } else {
      rep(
        length(unique(design_schemas)) == 1L,
        length(design_schemas)
      )
    }
    boundary_schema_checks <- c(
      boundary_schema_checks,
      design_schema_compatible
    )

    provenance <- rbind(
      provenance,
      do.call(
        rbind,
        lapply(
          seq_along(session_rows),
          function(index) {
            rows <- session_rows[[index]]
            data.frame(
              Dataset = dataset_name,
              BoundaryAuditGroup = sprintf(
                "B%03d",
                boundary_group
              ),
              SessionToken = sprintf(
                "S%03d",
                index
              ),
              SourceRowFirst = min(rows),
              SourceRowLast = max(rows),
              ExperimentalRows = length(rows),
              SourceRowsContiguous =
                identical(rows, seq.int(min(rows), max(rows))),
              ExperimentalSequenceContiguous = {
                positions <- experimental_row_position[rows]
                identical(
                  positions,
                  seq.int(min(positions), max(positions))
                )
              },
              DesignSchemaCompatible =
                design_schema_compatible[[index]],
              stringsAsFactors = FALSE
            )
          }
        )
      )
    )
  }

  analysis_counts <- table(analysis_subject[valid_subject])
  split_counts <- analysis_counts[
    grepl("::session[0-9]+$", names(analysis_counts))
  ]
  if (
    length(split_counts) > 0L &&
    any(split_counts != modal_rows)
  ) {
    stop(
      dataset_name,
      " produced an incomplete session while splitting a repeated ",
      "participant label."
    )
  }

  boundary_validation <- if (nrow(provenance) == 0L) {
    "NOT_APPLICABLE"
  } else if (
    all(provenance$ExperimentalSequenceContiguous) &&
    all(boundary_schema_checks)
  ) {
    "PASS"
  } else {
    "REVIEW"
  }

  list(
    key = analysis_subject,
    unsplit_key = subject,
    source_participants = length(subject_counts),
    analysis_participants = length(analysis_counts),
    modal_rows = modal_rows,
    additional_sessions = additional_sessions,
    session_construction = if (additional_sessions > 0L) {
      "row_count_reconstruction"
    } else {
      "source_label_only"
    },
    metadata_columns = "none",
    metadata_validation = "NOT_AVAILABLE",
    boundary_validation = boundary_validation,
    boundary_provenance = provenance
  )
}

build_source_manifests <- function(
  csv_files,
  datasets
) {
  if (
    length(csv_files) != length(datasets) ||
    !identical(
      tools::file_path_sans_ext(basename(csv_files)),
      names(datasets)
    )
  ) {
    stop("CSV paths and imported datasets are not aligned.")
  }

  checksum_rows <- vector("list", length(csv_files))
  schema_rows <- vector("list", length(csv_files))

  for (index in seq_along(csv_files)) {
    path <- csv_files[[index]]
    dat <- datasets[[index]]
    dataset_name <- names(datasets)[[index]]
    raw_header <- attr(dat, "maze_raw_header")
    normalized_header <- attr(dat, "maze_normalized_header")
    raw_header_occurrences <- vapply(
      raw_header,
      function(value) sum(raw_header == value),
      integer(1L)
    )
    unique_raw_header <- unique(raw_header)
    unique_raw_header_occurrences <- vapply(
      unique_raw_header,
      function(value) sum(raw_header == value),
      integer(1L)
    )
    duplicate_counts <- unique_raw_header_occurrences[
      unique_raw_header_occurrences > 1L
    ]
    experimental <- is_experimental_condition(dat$Condition)
    duplicate_experimental <-
      experimental & duplicated(dat)
    metadata_columns <- find_session_metadata_columns(dat)

    checksum_rows[[index]] <- data.frame(
      ManifestVersion = 1L,
      Dataset = dataset_name,
      RelativePath = file.path("data", basename(path)),
      Bytes = file.info(path)$size,
      SHA256 = sha256_file(path),
      HeaderSHA256 = sha256_text(
        paste(raw_header, collapse = "\034")
      ),
      Columns = ncol(dat),
      SourceRows = nrow(dat),
      ExperimentalRows = sum(experimental),
      NonexperimentalRows = sum(!experimental),
      ExactDuplicateExperimentalCopies =
        sum(duplicate_experimental),
      DuplicateRawHeaderNames = length(duplicate_counts),
      DuplicateRawHeaderExtraOccurrences = if (
        length(duplicate_counts) > 0L
      ) {
        sum(as.integer(duplicate_counts) - 1L)
      } else {
        0L
      },
      ExplicitSessionMetadataColumns = if (
        length(metadata_columns) > 0L
      ) {
        paste(metadata_columns, collapse = ";")
      } else {
        "none"
      },
      ImportDelimiter = "comma",
      ImportQuote = "double quote",
      ImportEncoding = "platform default",
      HeaderNormalization = "explicit-v1",
      stringsAsFactors = FALSE
    )

    schema_rows[[index]] <- data.frame(
      ManifestVersion = 1L,
      Dataset = dataset_name,
      ColumnPosition = seq_along(raw_header),
      RawHeader = ifelse(raw_header == "", "<blank>", raw_header),
      NormalizedHeader = normalized_header,
      RawHeaderOccurrences = raw_header_occurrences,
      Role = ifelse(
        normalized_header %in%
          c("Subject", "Condition", "Item", "Word", "RT"),
        "required-analysis",
        ifelse(
          normalized_header %in% metadata_columns,
          "session-metadata",
          "source-specific"
        )
      ),
      stringsAsFactors = FALSE
    )
  }

  list(
    checksum = do.call(rbind, checksum_rows),
    schema = do.call(rbind, schema_rows)
  )
}