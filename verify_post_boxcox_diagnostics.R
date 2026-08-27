# Independent verification for rt_process() post-Box-Cox diagnostics.
#
# Run from the project root with:
#   Rscript maze_preprocessing/verify_post_boxcox_diagnostics.R

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
script_path <- sub("^--file=", "", script_argument[[1L]])
script_directory <- dirname(normalizePath(script_path))
source(file.path(script_directory, "rt_process.R"))

reference_xapd <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  euler <- -digamma(1)
  z <- (x - mean(x)) /
    sqrt(var(x) * (n - 1) / n)
  z_nonzero <- z[z != 0]
  b2 <- sum(z^2 * sign(z)) / n
  k2 <- sum(z_nonzero^2 * log(abs(z_nonzero))) / n
  var_b2 <-
    (1 / n) *
    (3 - 8 / pi) *
    (1 - 1.9 / n)
  expected_net_k2 <-
    (
      (2 - log(2) - euler) / 2
    )^(1 / 3) *
    (1 - 1.026 / n)
  var_net_k2 <-
    (1 / n) *
    (1 / 72) *
    (
      (2 - log(2) - euler) / 2
    )^(-4 / 3) *
    (3 * pi^2 - 28) *
    (1 - 2.25 / n^0.8)
  net_k2 <- k2 - b2^2
  cube_root <- if (net_k2 == 0) {
    0
  } else {
    sign(net_k2) * abs(net_k2)^(1 / 3)
  }
  z_b2 <- b2 / sqrt(var_b2)
  z_net_k2 <-
    (cube_root - expected_net_k2) /
    sqrt(var_net_k2)
  statistic <- z_b2^2 + z_net_k2^2

  c(
    Z.B2 = z_b2,
    Z.net.K2 = z_net_k2,
    stat = statistic,
    p.value = pchisq(
      statistic,
      df = 2,
      lower.tail = FALSE
    )
  )
}

is_experimental_condition <- function(condition) {
  normalized <- gsub(
    "[^[:alnum:]]+",
    "",
    tolower(trimws(as.character(condition)))
  )
  administrative <- c(
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
  !(
    is.na(normalized) |
    normalized %in% administrative |
    grepl("^(filler|practice|burnin|intro)", normalized)
  )
}

make_analysis_subject <- function(data) {
  subject <- as.character(data$Subject)
  experimental <- is_experimental_condition(data$Condition)
  valid <- experimental &
    !is.na(subject) &
    nzchar(trimws(subject))
  counts <- table(subject[valid])
  frequencies <- table(as.integer(counts))
  modal_candidates <- as.integer(
    names(frequencies)[frequencies == max(frequencies)]
  )
  stopifnot(length(modal_candidates) == 1L)
  modal_rows <- modal_candidates[[1L]]
  result <- subject

  for (subject_id in names(counts)) {
    row_count <- as.integer(counts[[subject_id]])
    if (
      row_count <= modal_rows ||
      row_count %% modal_rows != 0L
    ) {
      next
    }
    rows <- which(valid & subject == subject_id)
    session <- ((seq_along(rows) - 1L) %/% modal_rows) + 1L
    result[rows] <- paste0(
      subject_id,
      "::session",
      session
    )
  }

  result
}

comparison_keys <- function(data, dataset_name) {
  condition <- data$Condition
  item <- data$Item
  region <- data$Word

  if (dataset_name %in% c(
    "fujita2026ex1",
    "fujita2026ex2",
    "fujita2026ex3"
  )) {
    condition <- data$Item
    item <- data$Condition
  }
  if (dataset_name %in% c(
    "hoeks2023exp1",
    "hoeks2023exp2",
    "hoeks2023exp3"
  )) {
    condition <- sub(
      "_[0-9]+$",
      "",
      as.character(data$Condition)
    )
  }
  if (dataset_name == "hoeks2023exp4") {
    region <- data$WordNo
  }

  list(condition = condition, item = item, region = region)
}

make_formula <- function(model_data) {
  predictors <- character()
  if (nlevels(model_data$Condition) > 1L) {
    predictors <- c(predictors, "Condition")
  }
  if (nlevels(model_data$Region) > 1L) {
    predictors <- c(predictors, "Region")
  }
  if (length(predictors) == 0L) {
    return(RT_BoxCox ~ 1)
  }
  reformulate(predictors, response = "RT_BoxCox")
}

direct_cluster_shift <- function(
  model,
  model_data,
  cluster,
  cluster_id
) {
  reduced <- model_data[cluster != cluster_id, , drop = FALSE]
  reduced_model <- lm(formula(model), data = reduced)
  full_coefficients <- coef(model)
  reduced_coefficients <- coef(reduced_model)
  coefficient_table <- coef(summary(model))
  coefficient_names <- setdiff(
    intersect(
      names(full_coefficients),
      names(reduced_coefficients)
    ),
    "(Intercept)"
  )
  standard_errors <-
    coefficient_table[coefficient_names, "Std. Error"]
  max(
    abs(
      reduced_coefficients[coefficient_names] -
      full_coefficients[coefficient_names]
    ) / standard_errors,
    na.rm = TRUE
  )
}

direct_cluster_oracle <- function(
  model,
  cluster,
  standardized_shift_threshold = 1
) {
  x <- model.matrix(model)
  y <- model.response(model.frame(model))
  full_coefficients <- coef(model)
  coefficient_table <- coef(summary(model))
  coefficient_names <- colnames(x)
  standard_errors <- rep(
    NA_real_,
    length(coefficient_names)
  )
  names(standard_errors) <- coefficient_names
  shared <- intersect(
    coefficient_names,
    rownames(coefficient_table)
  )
  standard_errors[shared] <-
    coefficient_table[shared, "Std. Error"]
  usable <-
    coefficient_names != "(Intercept)" &
    is.finite(standard_errors) &
    standard_errors > 0
  cluster <- as.character(cluster)
  cluster_levels <- unique(cluster)
  rank_failures <- 0L
  nonfinite_failures <- 0L
  assessed <- 0L
  maximum_shift <- -Inf
  maximum_cluster <- NA_character_
  maximum_coefficient <- NA_character_

  for (cluster_value in cluster_levels) {
    keep <- cluster != cluster_value
    if (sum(keep) <= ncol(x)) {
      rank_failures <- rank_failures + 1L
      next
    }
    reduced <- lm.fit(
      x = x[keep, , drop = FALSE],
      y = y[keep]
    )
    if (reduced$rank < ncol(x)) {
      rank_failures <- rank_failures + 1L
      next
    }
    if (!all(is.finite(reduced$coefficients))) {
      nonfinite_failures <- nonfinite_failures + 1L
      next
    }
    assessed <- assessed + 1L
    shifts <- abs(
      reduced$coefficients - full_coefficients
    ) / standard_errors
    shifts[!usable] <- NA_real_
    cluster_maximum <- max(shifts, na.rm = TRUE)
    if (cluster_maximum > maximum_shift) {
      maximum_shift <- cluster_maximum
      maximum_cluster <- cluster_value
      maximum_coefficient <- coefficient_names[
        which.max(replace(shifts, is.na(shifts), -Inf))
      ]
    }
  }

  if (!is.finite(maximum_shift)) {
    maximum_shift <- NA_real_
  }
  status <- if (
    rank_failures > 0L ||
    nonfinite_failures > 0L ||
    (
      is.finite(maximum_shift) &&
      maximum_shift >= standardized_shift_threshold
    )
  ) {
    "REVIEW"
  } else {
    "PASS"
  }

  list(
    status = status,
    assessed_clusters = assessed,
    rank_or_df_failures = rank_failures,
    nonfinite_failures = nonfinite_failures,
    maximum_standardized_shift = maximum_shift,
    most_influential_cluster = maximum_cluster,
    most_affected_coefficient = maximum_coefficient
  )
}

compare_cluster_result_to_oracle <- function(result, oracle) {
  stopifnot(
    result$status == oracle$status,
    result$assessed_clusters == oracle$assessed_clusters,
    result$rank_or_df_failures ==
      oracle$rank_or_df_failures,
    result$nonfinite_failures ==
      oracle$nonfinite_failures,
    identical(
      result$most_influential_cluster,
      oracle$most_influential_cluster
    ),
    identical(
      result$most_affected_coefficient,
      oracle$most_affected_coefficient
    ),
    isTRUE(all.equal(
      result$maximum_standardized_shift,
      oracle$maximum_standardized_shift,
      tolerance = 1e-7
    ))
  )
}

condition_verification_from_report <- function(
  report,
  data,
  keys,
  dataset_name,
  remove_no_comparator
) {
  rt <- suppressWarnings(
    as.numeric(trimws(as.character(data$RT)))
  )
  complete <- function(values) {
    !is.na(values) &
      nzchar(trimws(as.character(values)))
  }
  analysis_eligible <-
    is_experimental_condition(data$Condition) &
    is.finite(rt) &
    rt > 0 &
    complete(data$.AnalysisSubject) &
    complete(keys$condition) &
    complete(keys$item) &
    complete(keys$region)
  analysis_condition <- as.character(keys$condition)
  condition_levels <- sort(unique(
    analysis_condition[analysis_eligible]
  ))
  removals <- report$removed_rows
  removal_condition <- analysis_condition[removals$Row]
  do.call(
    rbind,
    lapply(
      condition_levels,
      function(condition_value) {
        in_condition <-
          analysis_eligible &
          analysis_condition == condition_value
        removed_in_condition <-
          removal_condition == condition_value
        reasons <- removals$Reason[removed_in_condition]
        cutoff_rows <- !is.na(
          removals$CutoffMs[removed_in_condition]
        )
        total_removed <- sum(removed_in_condition)
        data.frame(
          Dataset = dataset_name,
          RemoveNoComparator = remove_no_comparator,
          AnalysisCondition = condition_value,
          EligibleRTs = sum(in_condition),
          NoOtherItemComparator = sum(
            reasons == "no_other_item_comparison"
          ),
          ParticipantAndItemRatio = sum(
            reasons == "participant_and_item_ratio"
          ),
          ParticipantOnlyRatio = sum(
            reasons ==
              "participant_ratio_no_other_participant"
          ),
          CutoffRemoved = sum(cutoff_rows),
          TotalRemoved = total_removed,
          RetainedRTs = sum(in_condition) - total_removed,
          stringsAsFactors = FALSE
        )
      }
    )
  )
}

policy_verification_from_report <- function(
  report,
  dataset_name,
  remove_no_comparator,
  data,
  keys
) {
  reason_count <- function(reason) {
    counts <- report$iterative_removal_reason_counts
    if (!reason %in% names(counts)) {
      return(0L)
    }
    as.integer(counts[[reason]])
  }
  final_rts <- as.integer(
    report$post_boxcox_diagnostics$
      numerical_integrity$observations
  )
  list(
    row = data.frame(
      Dataset = dataset_name,
      RemoveNoComparator = remove_no_comparator,
      Multiplier = report$multiplier,
      EligibleRTs = final_rts + report$total_removed,
      NoOtherItemComparator = reason_count(
        "no_other_item_comparison"
      ),
      ParticipantAndItemRatio = reason_count(
        "participant_and_item_ratio"
      ),
      ParticipantOnlyRatio = reason_count(
        "participant_ratio_no_other_participant"
      ),
      CutoffRemoved = report$cutoff_removed,
      TotalRemoved = report$total_removed,
      FinalRTs = final_rts,
      Lambda = report$lambda,
      OverallStatus =
        report$post_boxcox_diagnostics$overall_status,
      ParticipantDeletion =
        report$post_boxcox_diagnostics$
          cluster_deletion$participant$status,
      ItemDeletion =
        report$post_boxcox_diagnostics$
          cluster_deletion$item$status,
      stringsAsFactors = FALSE
    ),
    clean_model_coefficients =
      report$clean_model_coefficients,
    condition_rows = condition_verification_from_report(
      report,
      data,
      keys,
      dataset_name,
      remove_no_comparator
    )
  )
}

verify_real_dataset <- function(dataset_name) {
  path <- file.path(
    script_directory,
    "data",
    paste0(dataset_name, ".csv")
  )
  data <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  experimental <- is_experimental_condition(data$Condition)
  data <- data[!(experimental & duplicated(data)), , drop = FALSE]
  data$.AnalysisSubject <- make_analysis_subject(data)
  keys <- comparison_keys(data, dataset_name)
  plot_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_path)
  output <- NULL
  invisible(
    capture.output(
      output <- rt_process(
        data = data,
        subject = data$.AnalysisSubject,
        condition = keys$condition,
        item = keys$item,
        region = keys$region,
        rt = data$RT,
        eligible =
          is_experimental_condition(data$Condition),
        remove_no_comparator = TRUE,
        rt_cutoff_ms = 9000
      )
    )
  )
  grDevices::dev.off()
  unlink(plot_path)

  report <- attr(output, "rt_process_report")
  stopifnot(
    !is.null(report$post_boxcox_diagnostics),
    !is.null(attr(output, "iterative_remove_report")$
      post_boxcox_diagnostics),
    report$total_removed == nrow(data) - nrow(output),
    nrow(report$removed_rows) == report$total_removed,
    nrow(report$iterative_removal_log) ==
      report$iterative_removed_total,
    identical(
      sort(report$iterative_removal_log$Row),
      report$iterative_removed_rows
    ),
    all(
      report$structural_removal_log$Reason ==
        "no_other_item_comparison"
    ),
    report$initial_multiplier == 1.5,
    report$multiplier_step == 0.5,
    report$max_reference_removal_percent == 5,
    report$multiplier >= report$initial_multiplier,
    report$percent_reference_removed <
      report$max_reference_removal_percent,
    all(
      abs(
        (
          report$removal_history$Multiplier -
            report$initial_multiplier
        ) / report$multiplier_step -
          round(
            (
              report$removal_history$Multiplier -
                report$initial_multiplier
            ) / report$multiplier_step
          )
      ) < 1e-12
    ),
    !report$lambda_profile_at_boundary,
    nrow(report$box_cox_profile) > 0L
  )

  output_keys <- comparison_keys(output, dataset_name)
  output_rt <- suppressWarnings(
    as.numeric(trimws(as.character(output$RT)))
  )
  complete <- function(x) {
    !is.na(x) & nzchar(trimws(as.character(x)))
  }
  analysis_rows <-
    is_experimental_condition(output$Condition) &
    is.finite(output_rt) &
    output_rt > 0 &
    complete(output$.AnalysisSubject) &
    complete(output_keys$condition) &
    complete(output_keys$item) &
    complete(output_keys$region)
  rt <- output_rt[analysis_rows]
  transformed <- output$RT_BoxCox[analysis_rows]
  lambda <- report$lambda

  stopifnot(
    length(transformed) ==
      report$post_boxcox_diagnostics$
        numerical_integrity$observations,
    all(is.finite(transformed)),
    is.finite(var(transformed)),
    var(transformed) > 0
  )

  inverse <- if (abs(lambda) < sqrt(.Machine$double.eps)) {
    exp(transformed)
  } else {
    (lambda * transformed + 1)^(1 / lambda)
  }
  direct_inverse_error <- max(
    abs(inverse - rt) / pmax(1, abs(rt))
  )
  stopifnot(
    isTRUE(all.equal(
      direct_inverse_error,
      report$post_boxcox_diagnostics$
        numerical_integrity$inverse_relative_error,
      tolerance = 1e-12
    ))
  )

  model_data <- data.frame(
    RT_BoxCox = transformed,
    Condition = factor(output_keys$condition[analysis_rows]),
    Region = factor(output_keys$region[analysis_rows])
  )
  model <- lm(make_formula(model_data), data = model_data)
  stopifnot(
    model$rank ==
      report$post_boxcox_diagnostics$model_estimability$rank,
    ncol(model.matrix(model)) ==
      report$post_boxcox_diagnostics$model_estimability$columns
  )

  direct_maximum_count <- sum(
    transformed == max(transformed)
  )
  stopifnot(
    direct_maximum_count ==
      report$post_boxcox_diagnostics$
        upper_boundary$maximum_count
  )

  clusters <- list(
    participant = output$.AnalysisSubject[analysis_rows],
    item = as.character(output_keys$item[analysis_rows])
  )
  for (cluster_name in names(clusters)) {
    result <- report$post_boxcox_diagnostics$
      cluster_deletion[[cluster_name]]
    if (
      !is.na(result$most_influential_cluster) &&
      result$rank_or_df_failures == 0L &&
      result$nonfinite_failures == 0L
    ) {
      direct_shift <- direct_cluster_shift(
        model = model,
        model_data = model_data,
        cluster = clusters[[cluster_name]],
        cluster_id = result$most_influential_cluster
      )
      stopifnot(
        isTRUE(all.equal(
          direct_shift,
          result$maximum_standardized_shift,
          tolerance = 1e-6
        ))
      )
    }
  }

  verification <- data.frame(
    Dataset = dataset_name,
    Rows = nrow(output),
    EligibleTransformed = length(transformed),
    Overall =
      report$post_boxcox_diagnostics$overall_status,
    StructuralRemoved =
      report$structural_removed_total,
    StructuralReferenceRemoved =
      report$structural_reference_removed,
    SelectedNoComparatorRemoved = sum(
      report$iterative_removal_log$Reason ==
        "no_other_item_comparison"
    ),
    IntegrityFailures = sum(
      report$post_boxcox_diagnostics$
        summary$Status == "FAIL"
    ),
    ReviewFindings = sum(
      report$post_boxcox_diagnostics$
        summary$Status == "REVIEW"
    ),
    ParticipantDeletionStatus =
      report$post_boxcox_diagnostics$
        cluster_deletion$participant$status,
    ParticipantDeletionShift =
      report$post_boxcox_diagnostics$
        cluster_deletion$participant$maximum_standardized_shift,
    ParticipantRankOrDfFailures =
      report$post_boxcox_diagnostics$
        cluster_deletion$participant$rank_or_df_failures,
    ParticipantNonfiniteFailures =
      report$post_boxcox_diagnostics$
        cluster_deletion$participant$nonfinite_failures,
    ItemDeletionStatus =
      report$post_boxcox_diagnostics$
        cluster_deletion$item$status,
    ItemDeletionShift =
      report$post_boxcox_diagnostics$
        cluster_deletion$item$maximum_standardized_shift,
    ItemDeletionCoefficient =
      report$post_boxcox_diagnostics$
        cluster_deletion$item$most_affected_coefficient,
    ItemDeletionFullEstimate =
      report$post_boxcox_diagnostics$
        cluster_deletion$item$most_affected_full_estimate,
    ItemDeletionEstimate =
      report$post_boxcox_diagnostics$
        cluster_deletion$item$most_affected_deletion_estimate,
    ItemDeletionSignChanged =
      report$post_boxcox_diagnostics$
        cluster_deletion$item$most_affected_sign_changed,
    ItemRankOrDfFailures =
      report$post_boxcox_diagnostics$
        cluster_deletion$item$rank_or_df_failures,
    ItemNonfiniteFailures =
      report$post_boxcox_diagnostics$
        cluster_deletion$item$nonfinite_failures,
    stringsAsFactors = FALSE
  )
  attr(
    verification,
    "policy_verification"
  ) <- policy_verification_from_report(
    report,
    dataset_name,
    TRUE,
    data,
    keys
  )
  verification
}

real_datasets <- c(
  "fujita2024memoryEx1",
  "fujita2024memoryEx2",
  "fujita2025lmaze",
  "fujita2026ex1",
  "fujita2026ex2",
  "fujita2026ex3",
  "fujitayoshida2024ex1",
  "fujitayoshida2024ex2",
  "fujitayoshida2024ex3",
  "fujitayoshida2024ex4",
  "hao2025exp1b",
  "hao2025exp2b",
  "hoeks2023exp1",
  "hoeks2023exp2",
  "hoeks2023exp3",
  "hoeks2023exp4",
  "orth2025",
  "orth2026lmaze",
  "orth2026amaze"
)
real_verifications <- lapply(
  real_datasets,
  verify_real_dataset
)
primary_policy_verifications <- lapply(
  real_verifications,
  attr,
  which = "policy_verification"
)
real_results <- do.call(rbind, real_verifications)

verify_retention_policy <- function(dataset_name) {
  path <- file.path(
    script_directory,
    "data",
    paste0(dataset_name, ".csv")
  )
  data <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  experimental <- is_experimental_condition(data$Condition)
  data <- data[!(experimental & duplicated(data)), , drop = FALSE]
  data$.AnalysisSubject <- make_analysis_subject(data)
  keys <- comparison_keys(data, dataset_name)
  plot_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_path)
  on.exit(
    {
      grDevices::dev.off()
      unlink(plot_path)
    },
    add = TRUE
  )
  output <- NULL
  invisible(capture.output(
    output <- rt_process(
      data = data,
      subject = data$.AnalysisSubject,
      condition = keys$condition,
      item = keys$item,
      region = keys$region,
      rt = data$RT,
      eligible = is_experimental_condition(data$Condition),
      remove_no_comparator = FALSE,
      rt_cutoff_ms = 9000
    )
  ))
  policy_verification_from_report(
    attr(output, "rt_process_report"),
    dataset_name,
    FALSE,
    data,
    keys
  )
}

retention_policy_verifications <- lapply(
  real_datasets,
  verify_retention_policy
)

# Freeze the publication-level sensitivity interpretation to this execution.
fujita_memory_1 <- real_results[
  real_results$Dataset == "fujita2024memoryEx1",
  ,
  drop = FALSE
]
stopifnot(
  sum(real_results$ParticipantDeletionStatus == "REVIEW") == 7L,
  sum(real_results$ItemDeletionStatus == "REVIEW") == 17L,
  all(real_results$ParticipantRankOrDfFailures == 0L),
  all(real_results$ParticipantNonfiniteFailures == 0L),
  all(real_results$ItemRankOrDfFailures == 0L),
  all(real_results$ItemNonfiniteFailures == 0L),
  nrow(fujita_memory_1) == 1L,
  round(fujita_memory_1$ItemDeletionShift, 3) == 6.852,
  identical(fujita_memory_1$ItemDeletionCoefficient, "Region7"),
  round(fujita_memory_1$ItemDeletionFullEstimate, 6) ==
    0.072689,
  round(fujita_memory_1$ItemDeletionEstimate, 6) ==
    0.064420,
  identical(fujita_memory_1$ItemDeletionSignChanged, FALSE)
)

canonical_manuscript_text <- paste(
  readLines(
    file.path(script_directory, "new_paper.tex"),
    warn = FALSE
  ),
  collapse = "\n"
)
canonical_findings_text <- paste(
  readLines(
    file.path(
      script_directory,
      "publication_s4_findings.tex"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)
one_column_position <- regexpr(
  "\\onecolumn",
  canonical_manuscript_text,
  fixed = TRUE
)[[1L]]
appendix_position <- regexpr(
  "\\appendix",
  canonical_manuscript_text,
  fixed = TRUE
)[[1L]]
stopifnot(
  one_column_position > 0L,
  appendix_position > one_column_position
)
fujita_shift_text <- sprintf(
  "%.3f",
  fujita_memory_1$ItemDeletionShift
)
stopifnot(
  grepl(
    paste0(
      "Item-deletion sensitivity & REVIEW & maximum shift ",
      fujita_shift_text,
      " SE"
    ),
    canonical_findings_text,
    fixed = TRUE
  ),
  grepl(
    "\\input{publication_s4_findings.tex}",
    canonical_manuscript_text,
    fixed = TRUE
  ),
  grepl(
    "screening model",
    canonical_manuscript_text,
    fixed = TRUE
  ),
  grepl(
    "Final-model sensitivity",
    canonical_manuscript_text,
    fixed = TRUE
  ),
  grepl(
    "is therefore explicitly",
    canonical_manuscript_text,
    fixed = TRUE
  ),
  grepl(
    "\\texttt{NOT\\_ASSESSED}",
    canonical_manuscript_text,
    fixed = TRUE
  ),
  !grepl("6.589", canonical_findings_text, fixed = TRUE)
)

policy_results <- read.csv(
  file.path(
    script_directory,
    "no_comparator_policy_results.csv"
  ),
  stringsAsFactors = FALSE,
  colClasses = "character",
  check.names = FALSE
)
condition_results <- read.csv(
  file.path(
    script_directory,
    "condition_exclusion_counts.csv"
  ),
  stringsAsFactors = FALSE,
  colClasses = "character",
  check.names = FALSE
)
policy_numeric_columns <- c(
  "EligibleRTs",
  "NoOtherItemComparator",
  "ParticipantAndItemRatio",
  "ParticipantOnlyRatio",
  "CutoffRemoved",
  "TotalRemoved",
  "FinalRTs"
)
condition_numeric_columns <- c(
  "EligibleRTs",
  "NoOtherItemComparator",
  "ParticipantAndItemRatio",
  "ParticipantOnlyRatio",
  "CutoffRemoved",
  "TotalRemoved",
  "RetainedRTs"
)
computed_policy_verifications <- c(
  primary_policy_verifications,
  retention_policy_verifications
)
computed_policy_rows <- do.call(
  rbind,
  lapply(
    computed_policy_verifications,
    function(verification) verification$row
  )
)
computed_policy_key <- paste(
  computed_policy_rows$Dataset,
  computed_policy_rows$RemoveNoComparator,
  sep = "::"
)
frozen_policy_key <- paste(
  policy_results$Dataset,
  policy_results$RemoveNoComparator,
  sep = "::"
)
computed_policy_rows <- computed_policy_rows[
  match(frozen_policy_key, computed_policy_key),
  ,
  drop = FALSE
]
stopifnot(
  !any(is.na(computed_policy_rows$Dataset)),
  identical(
    computed_policy_rows$Dataset,
    policy_results$Dataset
  ),
  identical(
    as.character(computed_policy_rows$RemoveNoComparator),
    policy_results$RemoveNoComparator
  )
)
computed_condition_rows <- do.call(
  rbind,
  lapply(
    computed_policy_verifications,
    function(verification) verification$condition_rows
  )
)
computed_condition_key <- paste(
  computed_condition_rows$Dataset,
  computed_condition_rows$RemoveNoComparator,
  computed_condition_rows$AnalysisCondition,
  sep = "::"
)
frozen_condition_key <- paste(
  condition_results$Dataset,
  condition_results$RemoveNoComparator,
  condition_results$AnalysisCondition,
  sep = "::"
)
computed_condition_rows <- computed_condition_rows[
  match(frozen_condition_key, computed_condition_key),
  ,
  drop = FALSE
]
stopifnot(
  nrow(computed_condition_rows) == nrow(condition_results),
  !any(is.na(computed_condition_rows$Dataset)),
  identical(
    computed_condition_rows$Dataset,
    condition_results$Dataset
  ),
  identical(
    as.character(
      computed_condition_rows$RemoveNoComparator
    ),
    condition_results$RemoveNoComparator
  ),
  identical(
    computed_condition_rows$AnalysisCondition,
    condition_results$AnalysisCondition
  )
)
for (column_name in condition_numeric_columns) {
  stopifnot(identical(
    as.integer(computed_condition_rows[[column_name]]),
    as.integer(condition_results[[column_name]])
  ))
}
computed_integer_columns <- c(
  "EligibleRTs",
  "NoOtherItemComparator",
  "ParticipantAndItemRatio",
  "ParticipantOnlyRatio",
  "CutoffRemoved",
  "TotalRemoved",
  "FinalRTs"
)
for (column_name in computed_integer_columns) {
  stopifnot(identical(
    as.integer(computed_policy_rows[[column_name]]),
    as.integer(policy_results[[column_name]])
  ))
}
stopifnot(
  isTRUE(all.equal(
    computed_policy_rows$Multiplier,
    as.numeric(policy_results$Multiplier),
    tolerance = 1e-12,
    check.attributes = FALSE
  )),
  isTRUE(all.equal(
    computed_policy_rows$Lambda,
    as.numeric(policy_results$Lambda),
    tolerance = 1e-12,
    check.attributes = FALSE
  )),
  identical(
    computed_policy_rows$OverallStatus,
    policy_results$OverallStatus
  ),
  identical(
    computed_policy_rows$ParticipantDeletion,
    policy_results$ParticipantDeletion
  ),
  identical(
    computed_policy_rows$ItemDeletion,
    policy_results$ItemDeletion
  )
)

compare_policy_coefficients <- function(
  primary_verification,
  retention_verification
) {
  comparison <- merge(
    primary_verification$clean_model_coefficients,
    retention_verification$clean_model_coefficients,
    by = "Term",
    suffixes = c("Primary", "Retention")
  )
  comparison <- comparison[
    comparison$Term != "(Intercept)",
    ,
    drop = FALSE
  ]
  c(
    PolicyShiftSE = max(
      abs(
        comparison$EstimateRetention -
          comparison$EstimatePrimary
      ) / comparison$StandardErrorPrimary,
      na.rm = TRUE
    ),
    PolicySignChanges = sum(
      sign(comparison$EstimatePrimary) !=
        sign(comparison$EstimateRetention)
    )
  )
}
computed_policy_sensitivity <- do.call(
  rbind,
  Map(
    compare_policy_coefficients,
    primary_policy_verifications,
    retention_policy_verifications
  )
)
computed_retention_rows <-
  policy_results$RemoveNoComparator == "FALSE"
computed_policy_sensitivity <- computed_policy_sensitivity[
  match(
    policy_results$Dataset[computed_retention_rows],
    real_datasets
  ),
  ,
  drop = FALSE
]
stopifnot(
  isTRUE(all.equal(
    round(computed_policy_sensitivity[, "PolicyShiftSE"], 3L),
    as.numeric(
      policy_results$PolicyShiftSE[computed_retention_rows]
    ),
    tolerance = 1e-12,
    check.attributes = FALSE
  )),
  identical(
    as.integer(
      computed_policy_sensitivity[, "PolicySignChanges"]
    ),
    as.integer(
      policy_results$PolicySignChanges[
        computed_retention_rows
      ]
    )
  )
)
stopifnot(
  nrow(policy_results) == 2L * length(real_datasets),
  nrow(condition_results) > nrow(policy_results),
  setequal(policy_results$Dataset, real_datasets),
  setequal(condition_results$Dataset, real_datasets),
  all(policy_results$Multiplier == "1.5"),
  all(policy_results$OverallStatus %in% c("PASS", "REVIEW", "FAIL")),
  all(policy_results$ParticipantDeletion %in%
    c("PASS", "REVIEW", "FAIL", "NOT_ASSESSABLE")),
  all(policy_results$ItemDeletion %in%
    c("PASS", "REVIEW", "FAIL", "NOT_ASSESSABLE")),
  all(grepl(
    "^-?[0-9]+\\.[0-9]{2}$",
    policy_results$LambdaTwoDecimals
  )),
  !any(policy_results$LambdaTwoDecimals == "-0.00"),
  all(
    table(
      policy_results$Dataset,
      policy_results$RemoveNoComparator
    ) == 1L
  ),
  sum(
    as.integer(
      policy_results$FinalRTs[
        policy_results$RemoveNoComparator == "FALSE"
      ]
    )
  ) -
    sum(
      as.integer(
        policy_results$FinalRTs[
          policy_results$RemoveNoComparator == "TRUE"
        ]
      )
    ) == 10059L
)

policy_pairs <- split(policy_results, policy_results$Dataset)
lambda_change_count <- sum(vapply(
  policy_pairs,
  function(rows) {
    length(unique(rows$LambdaTwoDecimals)) > 1L
  },
  logical(1L)
))
retention_policy_rows <-
  policy_results$RemoveNoComparator == "FALSE"
stopifnot(
  lambda_change_count == 2L,
  sum(
    as.integer(
      policy_results$PolicySignChanges[retention_policy_rows]
    ) > 0L
  ) == 4L,
  sum(
    as.integer(
      policy_results$PolicySignChanges[retention_policy_rows]
    )
  ) == 5L
)

for (row_number in seq_len(nrow(policy_results))) {
  policy_row <- policy_results[row_number, , drop = FALSE]
  condition_rows <- condition_results[
    condition_results$Dataset == policy_row$Dataset &
      condition_results$RemoveNoComparator ==
        policy_row$RemoveNoComparator,
    ,
    drop = FALSE
  ]
  stopifnot(nrow(condition_rows) > 0L)
  condition_sums <- vapply(
    condition_numeric_columns,
    function(column_name) {
      sum(as.integer(condition_rows[[column_name]]))
    },
    integer(1L)
  )
  names(condition_sums)[
    names(condition_sums) == "RetainedRTs"
  ] <- "FinalRTs"
  stopifnot(
    identical(
      unname(condition_sums[policy_numeric_columns]),
      as.integer(policy_row[1L, policy_numeric_columns])
    )
  )
}

required_supplement_inputs <- c(
  "publication_s2_source_audit.tex",
  "publication_s2_session_audit.tex",
  "publication_s2_session_sensitivity.tex",
  "publication_s3_adaptive_results.tex",
  "publication_s3_policy_results.tex",
  "publication_s3_condition_counts.tex",
  "publication_s4_status.tex",
  "publication_s4_findings.tex"
)
stopifnot(
  all(file.exists(
    file.path(script_directory, required_supplement_inputs)
  )),
  all(vapply(
    required_supplement_inputs,
    function(input_file) {
      grepl(
        paste0("\\input{", input_file, "}"),
        canonical_manuscript_text,
        fixed = TRUE
      )
    },
    logical(1L)
  ))
)

policy_table_text <- paste(
  readLines(
    file.path(
      script_directory,
      "publication_s3_policy_results.tex"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)
format_policy_integer <- function(value) {
  format(
    as.integer(value),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}
condition_table_text <- paste(
  readLines(
    file.path(
      script_directory,
      "publication_s3_condition_counts.tex"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)
escape_condition_latex <- function(value) {
  slash <- intToUtf8(92L)
  paste0(
    ifelse(
      strsplit(value, "", fixed = TRUE)[[1L]] %in%
        c("&", "%", "#", "_", "$"),
      paste0(
        slash,
        strsplit(value, "", fixed = TRUE)[[1L]]
      ),
      strsplit(value, "", fixed = TRUE)[[1L]]
    ),
    collapse = ""
  )
}
for (row_number in seq_len(nrow(condition_results))) {
  condition_row <- condition_results[
    row_number,
    ,
    drop = FALSE
  ]
  policy_label <- if (
    condition_row$RemoveNoComparator == "TRUE"
  ) {
    "Remove"
  } else {
    "Retain"
  }
  expected_condition_cells <- paste(
    c(
      policy_label,
      escape_condition_latex(
        condition_row$AnalysisCondition
      ),
      vapply(
        condition_numeric_columns,
        function(column_name) {
          format_policy_integer(
            condition_row[[column_name]]
          )
        },
        character(1L)
      )
    ),
    collapse = " & "
  )
  stopifnot(grepl(
    expected_condition_cells,
    condition_table_text,
    fixed = TRUE
  ))
}
for (row_number in seq_len(nrow(policy_results))) {
  policy_row <- policy_results[row_number, , drop = FALSE]
  policy_label <- if (
    policy_row$RemoveNoComparator == "TRUE"
  ) {
    "Remove"
  } else {
    "Retain"
  }
  shift_text <- if (
    is.na(suppressWarnings(
      as.numeric(policy_row$PolicyShiftSE)
    ))
  ) {
    "\\textemdash{}"
  } else {
    sprintf("%.3f", as.numeric(policy_row$PolicyShiftSE))
  }
  sign_change_text <- if (
    is.na(suppressWarnings(
      as.integer(policy_row$PolicySignChanges)
    ))
  ) {
    "\\textemdash{}"
  } else {
    as.character(as.integer(policy_row$PolicySignChanges))
  }
  expected_policy_cells <- paste(
    c(
      policy_label,
      policy_row$Multiplier,
      format_policy_integer(policy_row$EligibleRTs),
      format_policy_integer(
        policy_row$NoOtherItemComparator
      ),
      format_policy_integer(
        policy_row$ParticipantAndItemRatio
      ),
      format_policy_integer(
        policy_row$ParticipantOnlyRatio
      ),
      format_policy_integer(policy_row$CutoffRemoved),
      format_policy_integer(policy_row$FinalRTs),
      policy_row$LambdaTwoDecimals,
      policy_row$OverallStatus,
      policy_row$ParticipantDeletion,
      policy_row$ItemDeletion,
      shift_text,
      sign_change_text
    ),
    collapse = " & "
  )
  stopifnot(grepl(
    expected_policy_cells,
    policy_table_text,
    fixed = TRUE
  ))
}

run_rt_process_silently <- function(...) {
  plot_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_path)
  on.exit(
    {
      grDevices::dev.off()
      unlink(plot_path)
    },
    add = TRUE
  )
  result <- NULL
  invisible(
    capture.output(
      result <- rt_process(...)
    )
  )
  result
}

# Cross-check XAPD against a direct implementation of the published formulas.
xapd_inputs <- list(
  symmetric = qnorm(ppoints(101L)),
  right_skewed = exp(seq(-1, 1.5, length.out = 120L)),
  mixed = c(
    seq(-3, -0.1, length.out = 70L),
    seq(0.2, 4, length.out = 90L)
  )
)
for (xapd_input in xapd_inputs) {
  expected <- reference_xapd(xapd_input)
  observed <- Xapd.test(xapd_input)
  stopifnot(
    isTRUE(all.equal(
      unlist(
        observed[c(
          "Z.B2",
          "Z.net.K2",
          "stat",
          "p.value"
        )]
      ),
      expected,
      tolerance = 1e-12,
      check.attributes = FALSE
    ))
  )
}

# The public default retains rows lacking an other-item comparator; explicit
# structural removal removes them and records their original-row provenance.
structural_contract <- expand.grid(
  Subject = paste0("S", seq_len(8L)),
  Item = seq_len(6L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
structural_contract$Condition <- "A"
structural_contract$Word <- "1"
structural_subject_factor <- exp(seq(
  log(350),
  log(2800),
  length.out = 8L
))
structural_item_factor <- seq(0.9, 1.1, length.out = 6L)
structural_contract$RT <- round(
  structural_subject_factor[
    match(
      structural_contract$Subject,
      paste0("S", seq_len(8L))
    )
  ] *
    structural_item_factor[structural_contract$Item]
)
structural_contract <- rbind(
  structural_contract,
  data.frame(
    Subject = "isolated",
    Item = 99L,
    Condition = "B",
    Word = "1",
    RT = 700
  )
)
structural_default_output <- run_rt_process_silently(
  structural_contract
)
structural_default_report <- attr(
  structural_default_output,
  "rt_process_report"
)
structural_remove_output <- run_rt_process_silently(
  structural_contract,
  remove_no_comparator = TRUE
)
structural_remove_report <- attr(
  structural_remove_output,
  "rt_process_report"
)
stopifnot(
  nrow(structural_default_output) == nrow(structural_contract),
  any(structural_default_output$Subject == "isolated"),
  identical(
    structural_default_report$no_comparator_policy,
    "retain_rows_without_other_item_comparator"
  ),
  structural_default_report$structural_removed_total == 0L,
  nrow(structural_remove_output) == nrow(structural_contract) - 1L,
  !any(structural_remove_output$Subject == "isolated"),
  identical(
    structural_remove_report$no_comparator_policy,
    "remove_rows_without_other_item_comparator"
  ),
  identical(structural_remove_report$iterative_removed_rows, 49L),
  identical(structural_remove_report$removed_rows$Row, 49L),
  identical(
    structural_remove_report$removed_rows$Stage,
    "iterative_relative_rt"
  ),
  identical(
    structural_remove_report$removed_rows$Reason,
    "no_other_item_comparison"
  )
)

# NULL applies no absolute cutoff. An explicit cutoff removes only RTs strictly
# above its boundary, leaves excluded analysis rows untouched, and records
# dynamic cutoff provenance.
cutoff_contract <- expand.grid(
  Subject = paste0("S", seq_len(8L)),
  Item = seq_len(6L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
cutoff_contract$Condition <- "A"
cutoff_contract$Word <- "1"
cutoff_contract$RT <- round(
  structural_subject_factor[
    match(
      cutoff_contract$Subject,
      paste0("S", seq_len(8L))
    )
  ] *
    structural_item_factor[cutoff_contract$Item]
)
cutoff_contract$RT[cutoff_contract$Subject == "S8"] <- 9000
cutoff_contract$RT[[8L]] <- 10000
cutoff_contract <- rbind(
  cutoff_contract,
  data.frame(
    Subject = "",
    Item = 99L,
    Condition = "A",
    Word = "1",
    RT = 10000
  )
)
no_cutoff_contract_output <- run_rt_process_silently(
  cutoff_contract
)
no_cutoff_contract_report <- attr(
  no_cutoff_contract_output,
  "rt_process_report"
)
cutoff_contract_output <- run_rt_process_silently(
  cutoff_contract,
  rt_cutoff_ms = 9000
)
cutoff_contract_report <- attr(
  cutoff_contract_output,
  "rt_process_report"
)
stopifnot(
  nrow(no_cutoff_contract_output) == nrow(cutoff_contract),
  is.null(no_cutoff_contract_report$rt_cutoff_ms),
  identical(
    no_cutoff_contract_report$rt_cutoff_origin,
    "not_applied"
  ),
  no_cutoff_contract_report$cutoff_removed == 0L,
  identical(
    no_cutoff_contract_report$cutoff_removed_rows,
    integer()
  ),
  !any(grepl(
    "_ms_cutoff$",
    no_cutoff_contract_report$removed_rows$Stage
  )),
  all(is.na(no_cutoff_contract_report$removed_rows$CutoffMs)),
  nrow(cutoff_contract_output) ==
    nrow(cutoff_contract) - 1L,
  any(cutoff_contract_output$Item == 99L),
  any(cutoff_contract_output$RT == 9000),
  identical(cutoff_contract_report$rt_cutoff_ms, 9000),
  identical(cutoff_contract_report$rt_cutoff_origin, "user_supplied"),
  identical(
    cutoff_contract_report$cutoff_removed_rows,
    8L
  ),
  cutoff_contract_report$cutoff_removed == 1L,
  cutoff_contract_report$excluded_analysis_rows == 1L,
  identical(
    cutoff_contract_report$removed_rows$Row,
    8L
  ),
  identical(
    cutoff_contract_report$removed_rows$Stage,
    "9000_ms_cutoff"
  ),
  identical(
    cutoff_contract_report$removed_rows$Reason,
    "above_9000_ms"
  ),
  identical(
    names(cutoff_contract_report$software_versions),
    c("R", "MASS")
  ),
  !cutoff_contract_report$lambda_profile_at_boundary,
  nrow(cutoff_contract_report$box_cox_profile) > 0L
)

invalid_cutoffs <- list(
  0,
  -1,
  Inf,
  NaN,
  NA_real_,
  numeric(),
  c(9000, 10000),
  "9000"
)
for (invalid_cutoff in invalid_cutoffs) {
  cutoff_error <- try(
    run_rt_process_silently(
      cutoff_contract,
      rt_cutoff_ms = invalid_cutoff
    ),
    silent = TRUE
  )
  stopifnot(
    inherits(cutoff_error, "try-error"),
    grepl(
      "`rt_cutoff_ms` must be NULL or one finite positive number.",
      conditionMessage(attr(cutoff_error, "condition")),
      fixed = TRUE
    )
  )
}

# Iterative removals retain original row-level provenance.
iterative_contract <- expand.grid(
  Subject = paste0("S", seq_len(4L)),
  Item = seq_len(3L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
iterative_contract$Condition <- "A"
iterative_contract$Word <- "1"
iterative_contract$RT <- 500 + seq_len(nrow(iterative_contract))
iterative_contract$RT[[1L]] <- 8000
iterative_contract_output <- run_rt_process_silently(
  iterative_contract
)
iterative_contract_report <- attr(
  iterative_contract_output,
  "rt_process_report"
)
stopifnot(
  identical(
    iterative_contract_report$iterative_removed_rows,
    1L
  ),
  identical(
    iterative_contract_report$removed_rows$Row,
    1L
  ),
  identical(
    iterative_contract_report$removed_rows$Stage,
    "iterative_relative_rt"
  ),
  identical(
    iterative_contract_report$removed_rows$Reason,
    "participant_and_item_ratio"
  ),
  iterative_contract_report$total_removed ==
    nrow(iterative_contract) -
      nrow(iterative_contract_output)
)

# The public minimum must agree with XAPD's four-observation requirement.
three_clean_rows <- data.frame(
  Subject = c("S1", "S1", "S1", "S2"),
  Condition = "A",
  Item = c(1L, 2L, 3L, 1L),
  Word = "1",
  RT = c(300, 400, 500, 5000)
)
three_clean_error <- try(
  run_rt_process_silently(
    three_clean_rows,
    remove_no_comparator = TRUE
  ),
  silent = TRUE
)
stopifnot(
  inherits(three_clean_error, "try-error"),
  grepl(
    "At least four positive, varying RTs",
    conditionMessage(attr(three_clean_error, "condition")),
    fixed = TRUE
  )
)

# A profile maximum still on the outer search boundary is not an estimate.
boundary_contract <- expand.grid(
  Subject = paste0("S", seq_len(5L)),
  Item = seq_len(6L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
boundary_contract$Condition <- "A"
boundary_contract$Word <- "1"
boundary_contract$RT <-
  1000 +
  exp(
    seq(
      0,
      3,
      length.out = 6L
    )
  )[boundary_contract$Item]
boundary_error <- try(
  run_rt_process_silently(boundary_contract),
  silent = TRUE
)
stopifnot(
  inherits(boundary_error, "try-error"),
  grepl(
    "profile maximum remained at the boundary",
    conditionMessage(attr(boundary_error, "condition")),
    fixed = TRUE
  )
)

# Synthetic passing and review paths through the public function.
synthetic <- expand.grid(
  Subject = paste0("S", seq_len(20L)),
  Item = seq_len(50L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
synthetic$Condition <- ifelse(
  synthetic$Item <= 25L,
  "A",
  "B"
)
synthetic$Word <- "1"
set.seed(9182)
synthetic$RT <- pmin(
  850,
  pmax(
    250,
    round(
      exp(
        rnorm(
          nrow(synthetic),
          log(500),
          0.25
        )
      )
    )
  )
)
synthetic$RT[
  synthetic$Subject %in% paste0("S", seq_len(10L)) &
    synthetic$Item == 50L
] <- 900
synthetic$RT[synthetic$Subject == "S2"] <-
  synthetic$RT[synthetic$Subject == "S1"]

synthetic_plot <- tempfile(fileext = ".pdf")
grDevices::pdf(synthetic_plot)
synthetic_output <- NULL
synthetic_default_printout <- capture.output(
  synthetic_output <- rt_process(synthetic)
)
synthetic_verbose_output <- NULL
synthetic_verbose_printout <- capture.output(
  synthetic_verbose_output <- rt_process(
    synthetic,
    verbose = TRUE
  )
)
grDevices::dev.off()
unlink(synthetic_plot)
synthetic_report <- attr(
  synthetic_output,
  "rt_process_report"
)$post_boxcox_diagnostics
synthetic_duplicate_provenance <-
  synthetic_report$duplicate_sessions$duplicate_group_provenance
stopifnot(
  synthetic_report$numerical_integrity$finite_values,
  synthetic_report$numerical_integrity$inverse_recovery,
  synthetic_report$duplicate_sessions$duplicate_groups >= 1L,
  nrow(synthetic_duplicate_provenance) ==
    synthetic_report$duplicate_sessions$duplicated_sessions,
  identical(
    names(synthetic_duplicate_provenance),
    c(
      "GroupToken",
      "SessionToken",
      "Observations",
      "FirstSourceRow",
      "LastSourceRow",
      "SourceRowsContiguous"
    )
  ),
  all(grepl(
    "^duplicate_group_[0-9]{3}$",
    synthetic_duplicate_provenance$GroupToken
  )),
  all(grepl(
    "^duplicate_session_[0-9]{3}$",
    synthetic_duplicate_provenance$SessionToken
  )),
  all(synthetic_duplicate_provenance$Observations >= 20L),
  all(synthetic_duplicate_provenance$FirstSourceRow >= 1L),
  all(
    synthetic_duplicate_provenance$LastSourceRow <=
      nrow(synthetic)
  ),
  !any(
    c("Subject", "TrajectorySignature") %in%
      names(synthetic_duplicate_provenance)
  ),
  synthetic_report$summary$Status[
    synthetic_report$summary$Diagnostic ==
      "Upper-boundary concentration"
  ] == "REVIEW"
)
stopifnot(
  any(grepl(
    "Exact duplicated session trajectories",
    synthetic_default_printout,
    fixed = TRUE
  )),
  !any(grepl(
    "Finite transformed values",
    synthetic_default_printout,
    fixed = TRUE
  )),
  any(grepl(
    "Finite transformed values",
    synthetic_verbose_printout,
    fixed = TRUE
  )),
  identical(
    attr(synthetic_verbose_output, "rt_process_report")$
      post_boxcox_diagnostics$summary,
    synthetic_report$summary
  ),
  inherits(
    try(rt_process(synthetic, verbose = NA), silent = TRUE),
    "try-error"
  ),
  !any(grepl(
    "How the procedure works",
    synthetic_default_printout,
    fixed = TRUE
  )),
  !any(grepl(
    "Interpretation rules:",
    synthetic_default_printout,
    fixed = TRUE
  )),
  any(grepl(
    "Descriptive Xapd distribution-shape diagnostics",
    synthetic_default_printout,
    fixed = TRUE
  )),
  any(grepl(
    "planned final inferential model",
    synthetic_default_printout,
    fixed = TRUE
  )),
  any(grepl(
    "Absolute RT cutoff: not applied",
    synthetic_default_printout,
    fixed = TRUE
  )),
  !any(grepl(
    "detected departure from normality",
    synthetic_default_printout,
    fixed = TRUE
  ))
)
xapd_summary_lines <- grep(
  "^Xapd summary:",
  synthetic_default_printout,
  value = TRUE
)
stopifnot(
  length(xapd_summary_lines) == 1L,
  grepl("Initial RT", xapd_summary_lines, fixed = TRUE),
  grepl(
    "Clean RT (no absolute cutoff)",
    xapd_summary_lines,
    fixed = TRUE
  ),
  grepl(
    "End product (Box-Cox)",
    xapd_summary_lines,
    fixed = TRUE
  )
)

# Brute-force every synthetic participant and item deletion with direct QR fits.
synthetic_analysis_rows <- is.finite(synthetic_output$RT_BoxCox)
synthetic_model_data <- data.frame(
  RT_BoxCox =
    synthetic_output$RT_BoxCox[synthetic_analysis_rows],
  Condition = factor(
    synthetic_output$Condition[synthetic_analysis_rows]
  ),
  Region = factor(
    synthetic_output$Word[synthetic_analysis_rows]
  )
)
synthetic_model <- lm(
  make_formula(synthetic_model_data),
  data = synthetic_model_data
)
synthetic_clusters <- list(
  participant =
    synthetic_output$Subject[synthetic_analysis_rows],
  item = as.character(
    synthetic_output$Item[synthetic_analysis_rows]
  )
)
for (cluster_name in names(synthetic_clusters)) {
  oracle <- direct_cluster_oracle(
    model = synthetic_model,
    cluster = synthetic_clusters[[cluster_name]]
  )
  compare_cluster_result_to_oracle(
    synthetic_report$cluster_deletion[[cluster_name]],
    oracle
  )
}

# The one-standard-error threshold is inclusive.
synthetic_design_cell <- interaction(
  synthetic_model_data$Condition,
  synthetic_model_data$Region,
  drop = TRUE,
  lex.order = TRUE
)
participant_oracle <- direct_cluster_oracle(
  model = synthetic_model,
  cluster = synthetic_clusters$participant
)
threshold_review <- .post_boxcox_cluster_deletion(
  model = synthetic_model,
  cluster = synthetic_clusters$participant,
  cluster_name = "participant",
  design_cell = synthetic_design_cell,
  standardized_shift_threshold =
    participant_oracle$maximum_standardized_shift *
      (1 - 1e-10)
)
threshold_pass <- .post_boxcox_cluster_deletion(
  model = synthetic_model,
  cluster = synthetic_clusters$participant,
  cluster_name = "participant",
  design_cell = synthetic_design_cell,
  standardized_shift_threshold =
    participant_oracle$maximum_standardized_shift *
      (1 + 1e-10)
)
stopifnot(
  threshold_review$status == "REVIEW",
  threshold_pass$status == "PASS"
)

# Deleting the only cluster that supplies one condition must report rank loss.
rank_loss_data <- data.frame(
  y = c(1.0, 1.1, 0.9, 1.2, 2.0, 2.1),
  condition = factor(c("A", "A", "A", "A", "B", "B")),
  region = factor("1"),
  participant = c("S1", "S1", "S2", "S2", "S3", "S3")
)
rank_loss_model <- lm(
  y ~ condition,
  data = rank_loss_data
)
rank_loss_design_cell <- interaction(
  rank_loss_data$condition,
  rank_loss_data$region,
  drop = TRUE,
  lex.order = TRUE
)
rank_loss_result <- .post_boxcox_cluster_deletion(
  model = rank_loss_model,
  cluster = rank_loss_data$participant,
  cluster_name = "participant",
  design_cell = rank_loss_design_cell
)
rank_loss_oracle <- direct_cluster_oracle(
  model = rank_loss_model,
  cluster = rank_loss_data$participant
)
compare_cluster_result_to_oracle(
  rank_loss_result,
  rank_loss_oracle
)
stopifnot(
  rank_loss_result$status == "REVIEW",
  rank_loss_result$rank_or_df_failures == 1L
)

# Stable weighted-QR deletion must agree with direct QR in an ill-conditioned
# full-rank design.
set.seed(441)
ill_conditioned_data <- data.frame(
  x1 = seq(-1, 1, length.out = 240L)
)
ill_conditioned_data$x2 <-
  ill_conditioned_data$x1 +
  rnorm(nrow(ill_conditioned_data), sd = 1e-6)
ill_conditioned_data$y <-
  2 +
  0.5 * ill_conditioned_data$x1 -
  0.5 * ill_conditioned_data$x2 +
  rnorm(nrow(ill_conditioned_data), sd = 0.1)
ill_conditioned_data$cluster <- rep(
  paste0("C", seq_len(24L)),
  each = 10L
)
ill_conditioned_model <- lm(
  y ~ x1 + x2,
  data = ill_conditioned_data
)
ill_conditioned_result <- .post_boxcox_cluster_deletion(
  model = ill_conditioned_model,
  cluster = ill_conditioned_data$cluster,
  cluster_name = "cluster",
  design_cell = seq_len(nrow(ill_conditioned_data))
)
ill_conditioned_oracle <- direct_cluster_oracle(
  model = ill_conditioned_model,
  cluster = ill_conditioned_data$cluster
)
compare_cluster_result_to_oracle(
  ill_conditioned_result,
  ill_conditioned_oracle
)

# Synthetic hard-failure paths through the independent private diagnostic.
valid_rt <- seq(300, 499)
valid_transformed <- log(valid_rt)
condition <- factor(rep(c("A", "B"), each = 100))
region <- factor(rep(c("1", "2"), times = 100))
valid_model <- lm(
  valid_transformed ~ condition + region
)
corrupted <- valid_transformed
corrupted[[1L]] <- Inf
corrupted_result <- .post_boxcox_diagnostics(
  rt = valid_rt,
  transformed = corrupted,
  lambda = 0,
  model = valid_model,
  subject = rep(paste0("S", seq_len(10L)), each = 20L),
  condition = condition,
  item = rep(seq_len(20L), times = 10L),
  region = region,
  source_row = seq_along(valid_rt)
)
stopifnot(
  corrupted_result$overall_status == "FAIL",
  corrupted_result$summary$Status[
    corrupted_result$summary$Diagnostic ==
      "Finite transformed values"
  ] == "FAIL"
)

cat("Independent real-dataset verification passed:\n")
print(real_results, row.names = FALSE)
cat(
  "\nSynthetic pass/review and hard-failure paths passed.\n"
)