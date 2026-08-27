# Adaptive iterative relative RT removal, Box-Cox transformation,
# and descriptive Xapd distribution-shape diagnostics
#
# Public function:
#   rt_process()
#
# The public function starts at multiplier 1.5 and increases it in 0.5 steps
# until fewer than 5% of eligible finite RTs from 200 through 4,000 ms are
# removed. By default, rows without an other-item comparator are retained;
# `remove_no_comparator = TRUE` explicitly removes such rows, including when
# comparator availability changes during iterative removal.
# By default, no absolute RT cutoff is applied. Supplying `rt_cutoff_ms` removes
# remaining analysis-eligible finite RTs strictly above that cutoff. The
# function then estimates a Box-Cox lambda, adds the transformed RT as
# `RT_BoxCox`, reports descriptive Xapd shape diagnostics, and draws
# Raw-versus-End-product Q-Q and residual-versus-fitted plots.
#
# The private fixed-multiplier helper removes one currently qualifying finite
# eligible RT at a time.
# Within each independent condition-region group, the largest qualifying RT is
# removed first and all comparison maxima are then recalculated. Iteration
# stops at a fixed point.
#
# A focal RT m qualifies when:
#   1. no other-item comparison can be calculated for the participant and
#      `remove_no_comparator` is TRUE; or
#   2. m is greater than `multiplier` times every finite RT from the
#      participant's other items at that condition and region, and either:
#      a. m is greater than `multiplier` times every finite RT from other
#         participants at the same condition, item, and region; or
#      b. no other participant has a comparison value.
#
# Rows that are not eligible, and missing or non-finite RTs, are never removed.
# The public workflow also retains rows with nonpositive RTs or incomplete
# participant/condition/item/region keys; those rows are excluded from
# multiplier selection, lambda estimation, and diagnostic summaries.
#
# Return value:
#   The retained rows with the original columns plus `RT_BoxCox`. The selected
#   multiplier, Box-Cox lambda, removal history, descriptive statistics, and
#   conservative post-Box-Cox integrity/sensitivity diagnostics are stored in
#   the `rt_process_report` attribute. Original row positions for every removal
#   are retained in that report. Exact duplicated trajectories are represented
#   only by sequential group/session tokens and source-row provenance, never by
#   raw participant/session identifiers or trajectory signatures.
#
# Minimal use:
#   source("rt_process.R")
#   cleaned_data <- rt_process(data = my_data)
#   detailed_data <- rt_process(data = my_data, verbose = TRUE)

# XAPD formulas follow Desgagne and Lafaye de Micheaux (2018),
# Journal of Applied Statistics, 45, 2307-2327,
# https://doi.org/10.1080/02664763.2017.1415311.
.validate_xapd_input <- function(
  x,
  function_name = "Xapd"
) {
  x <- if (is.numeric(x)) {
    as.numeric(x)
  } else {
    suppressWarnings(
      as.numeric(trimws(as.character(x)))
    )
  }
  x <- x[is.finite(x)]

  if (length(x) < 4L) {
    stop(
      function_name,
      " requires at least four finite numeric observations."
    )
  }
  if (!is.finite(var(x)) || var(x) <= 0) {
    stop(
      function_name,
      " requires a distribution with positive finite variance."
    )
  }

  x
}

B2 <- function(x) {
  x <- .validate_xapd_input(x, "B2")
  n <- length(x)
  z <- (x - mean(x)) /
    sqrt(var(x) * (n - 1) / n)
  mean(z^2 * sign(z))
}

K2 <- function(x) {
  x <- .validate_xapd_input(x, "K2")
  n <- length(x)
  z <- (x - mean(x)) /
    sqrt(var(x) * (n - 1) / n)
  z <- z[!(z == 0)]
  sum(z^2 * log(abs(z))) / n
}

Z.B2 <- function(x) {
  x <- .validate_xapd_input(x, "Z.B2")
  n <- length(x)
  var.B2 <-
    (1 / n) *
    (3.0 - 8.0 / pi) *
    (1.0 - 1.9 / n)
  B2(x) / sqrt(var.B2)
}

Z.net.K2 <- function(x) {
  x <- .validate_xapd_input(x, "Z.net.K2")
  n <- length(x)
  euler <- -digamma(1)
  esp.net.K2 <-
    (
      (2.0 - log(2.0) - euler) / 2.0
    )^(1.0 / 3.0) *
    (1 - 1.026 / n)
  var.net.K2 <-
    (1.0 / n) *
    (1.0 / 72.0) *
    (
      (2.0 - log(2.0) - euler) / 2.0
    )^(-4.0 / 3.0) *
    (3.0 * pi^2 - 28.0) *
    (1.0 - 2.25 / n^0.8)
  net.K2 <- K2(x) - B2(x)^2
  real_cube_root <- if (net.K2 == 0) {
    0
  } else {
    sign(net.K2) * abs(net.K2)^(1.0 / 3.0)
  }
  (
    real_cube_root - esp.net.K2
  ) / sqrt(var.net.K2)
}

Xapd.test <- function(x) {
  x <- .validate_xapd_input(x, "Xapd")

  Z.B2.x <- Z.B2(x)
  Z.net.K2.x <- Z.net.K2(x)
  Xapd.stat <- Z.B2.x^2 + Z.net.K2.x^2
  p.value <- pchisq(
    Xapd.stat,
    2,
    lower.tail = FALSE
  )

  if (
    !all(
      is.finite(
        c(
          Z.B2.x,
          Z.net.K2.x,
          Xapd.stat,
          p.value
        )
      )
    )
  ) {
    stop(
      "Xapd could not calculate finite diagnostics for this ",
      "distribution."
    )
  }

  list(
    test =
      "2nd-power skewness and kurtosis omnibus normality test",
    Z.B2 = Z.B2.x,
    Z.net.K2 = Z.net.K2.x,
    stat = Xapd.stat,
    p.value = p.value
  )
}

.remove_iterative_rt_outliers_at_multiplier <- function(
  data,
  subject = data$Subject,
  condition = data$Condition,
  item = data$Item,
  region = data$Word,
  rt = data$RT,
  eligible = rep(TRUE, nrow(data)),
  multiplier,
  remove_no_comparator = TRUE
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  n <- nrow(data)
  inputs <- list(
    subject = subject,
    condition = condition,
    item = item,
    region = region,
    rt = rt,
    eligible = eligible
  )
  invalid_lengths <- vapply(
    inputs,
    length,
    integer(1)
  ) != n

  if (any(invalid_lengths)) {
    stop(
      "Every comparison vector must have one value per data row."
    )
  }

  if (
    length(multiplier) != 1L ||
    !is.numeric(multiplier) ||
    !is.finite(multiplier) ||
    multiplier <= 1
  ) {
    stop("`multiplier` must be one finite number greater than 1.")
  }
  if (
    !is.logical(remove_no_comparator) ||
      length(remove_no_comparator) != 1L ||
      is.na(remove_no_comparator)
  ) {
    stop("`remove_no_comparator` must be either TRUE or FALSE.")
  }

  subject <- as.character(subject)
  condition <- as.character(condition)
  item <- as.character(item)
  region <- as.character(region)
  rt_numeric <- suppressWarnings(
    as.numeric(trimws(as.character(rt)))
  )
  rt_numeric[!is.finite(rt_numeric)] <- NA_real_
  eligible <- as.logical(eligible)
  eligible[is.na(eligible)] <- FALSE

  remove <- rep(FALSE, n)
  removal_log <- data.frame(
    Row = integer(),
    Reason = character(),
    stringsAsFactors = FALSE
  )

  max_excluding_group <- function(values, groups) {
    result <- rep(NA_real_, length(values))

    for (indices in split(seq_along(values), groups)) {
      if (length(indices) <= 1L) {
        next
      }

      group_values <- values[indices]
      largest <- max(group_values)
      at_largest <- group_values == largest

      if (sum(at_largest) > 1L) {
        result[indices] <- largest
      } else {
        result[indices[!at_largest]] <- largest
        result[indices[at_largest]] <-
          max(group_values[!at_largest])
      }
    }

    result
  }

  evaluate_active_rows <- function(
    active,
    cell_id,
    group_rt,
    n_items
  ) {
    cell_max <- tapply(
      group_rt[active],
      cell_id[active],
      max
    )
    active_cells <- as.integer(names(cell_max))
    cell_max <- as.numeric(cell_max)
    cell_subject <-
      ((active_cells - 1L) %/% n_items) + 1L
    cell_item <-
      ((active_cells - 1L) %% n_items) + 1L

    within_cell_max <- max_excluding_group(
      cell_max,
      cell_subject
    )
    between_cell_max <- max_excluding_group(
      cell_max,
      cell_item
    )
    within_cell_count <-
      tabulate(cell_subject)[cell_subject] - 1L
    between_cell_count <-
      tabulate(cell_item)[cell_item] - 1L

    cell_position <- match(
      cell_id[active],
      active_cells
    )
    within_max <- within_cell_max[cell_position]
    within_count <- within_cell_count[cell_position]
    between_max <- between_cell_max[cell_position]
    between_count <-
      between_cell_count[cell_position]

    first_missing <-
      is.na(within_max) | within_count == 0L
    first_met <-
      !first_missing &
      within_max < (group_rt[active] / multiplier)
    second_missing <-
      is.na(between_max) | between_count == 0L
    second_met <-
      !second_missing &
      between_max < (group_rt[active] / multiplier)
    qualifies <-
      (remove_no_comparator & first_missing) |
      (
        first_met &
        (second_met | second_missing)
      )

    data.frame(
      LocalRow = active[qualifies],
      Reason = ifelse(
        first_missing[qualifies],
        "no_other_item_comparison",
        ifelse(
          second_missing[qualifies],
          "participant_ratio_no_other_participant",
          "participant_and_item_ratio"
        )
      ),
      stringsAsFactors = FALSE
    )
  }

  iterate_group <- function(rows) {
    group_subject <- subject[rows]
    group_item <- item[rows]
    group_rt <- rt_numeric[rows]
    subject_id <- match(
      group_subject,
      unique(group_subject)
    )
    item_id <- match(group_item, unique(group_item))
    n_items <- max(item_id)
    cell_id <-
      (subject_id - 1L) * n_items + item_id
    active <- seq_along(rows)
    removed <- data.frame(
      Row = integer(),
      Reason = character(),
      stringsAsFactors = FALSE
    )

    repeat {
      if (length(active) == 0L) {
        break
      }

      qualifying_rows <- evaluate_active_rows(
        active = active,
        cell_id = cell_id,
        group_rt = group_rt,
        n_items = n_items
      )

      if (nrow(qualifying_rows) == 0L) {
        break
      }

      chosen_candidate <- order(
        -group_rt[qualifying_rows$LocalRow],
        rows[qualifying_rows$LocalRow]
      )[1L]
      chosen_local_row <- qualifying_rows$LocalRow[
        chosen_candidate
      ]
      chosen_reason <- qualifying_rows$Reason[
        chosen_candidate
      ]
      chosen_global_row <- rows[chosen_local_row]
      removed <- rbind(
        removed,
        data.frame(
          Row = chosen_global_row,
          Reason = chosen_reason,
          stringsAsFactors = FALSE
        )
      )
      active <- active[active != chosen_local_row]
    }

    removed
  }

  finite_eligible <- eligible & !is.na(rt_numeric)
  complete_keys <-
    !is.na(subject) &
    !is.na(condition) &
    !is.na(item) &
    !is.na(region)
  unkeyed_rows <- which(
    finite_eligible & !complete_keys
  )
  valid_rows <- which(finite_eligible & complete_keys)

  if (
    remove_no_comparator &&
      length(unkeyed_rows) > 0L
  ) {
    remove[unkeyed_rows] <- TRUE
    removal_log <- rbind(
      removal_log,
      data.frame(
        Row = unkeyed_rows,
        Reason = "incomplete_comparison_key",
        stringsAsFactors = FALSE
      )
    )
  }

  if (length(valid_rows) > 0L) {
    condition_id <- match(
      condition[valid_rows],
      unique(condition[valid_rows])
    )
    region_id <- match(
      region[valid_rows],
      unique(region[valid_rows])
    )
    n_regions <- max(region_id)
    group_id <-
      (condition_id - 1L) * n_regions + region_id

    for (rows in split(valid_rows, group_id)) {
      if (
        length(unique(condition[rows])) != 1L ||
        length(unique(region[rows])) != 1L
      ) {
        stop(
          "A comparison group crossed conditions or regions."
        )
      }

      group_removals <- iterate_group(rows)
      if (nrow(group_removals) > 0L) {
        remove[group_removals$Row] <- TRUE
        removal_log <- rbind(
          removal_log,
          group_removals
        )
      }
    }
  }

  if (!remove_no_comparator) {
    stopifnot(
      !any(
        removal_log$Reason %in%
          c(
            "no_other_item_comparison",
            "incomplete_comparison_key"
          )
      )
    )
  }
  stopifnot(
    identical(
      sort(unique(removal_log$Row)),
      which(remove)
    )
  )
  result <- data[!remove, , drop = FALSE]
  attr(result, "iterative_removal_log") <- removal_log
  result
}


.post_boxcox_cluster_deletion <- function(
  model,
  cluster,
  cluster_name,
  design_cell,
  standardized_shift_threshold = 1
) {
  model_matrix <- model.matrix(model)
  response <- model.response(model.frame(model))
  coefficient_names <- colnames(model_matrix)
  coefficient_rows <- coefficient_names != "(Intercept)"
  cluster <- as.character(cluster)
  design_cell <- as.character(design_cell)
  cluster_levels <- unique(cluster)

  if (length(design_cell) != nrow(model_matrix)) {
    stop(
      "`design_cell` must have one value per transformed-model row."
    )
  }

  empty_result <- function(status, interpretation) {
    list(
      status = status,
      cluster_type = cluster_name,
      clusters = length(cluster_levels),
      assessed_clusters = 0L,
      rank_or_df_failures = 0L,
      nonfinite_failures = 0L,
      maximum_standardized_shift = NA_real_,
      threshold = standardized_shift_threshold,
      most_influential_cluster = NA_character_,
      most_affected_coefficient = NA_character_,
      most_affected_full_estimate = NA_real_,
      most_affected_deletion_estimate = NA_real_,
      most_affected_sign_changed = NA,
      interpretation = interpretation
    )
  }

  if (length(cluster_levels) < 2L) {
    return(
      empty_result(
        "NOT_ASSESSABLE",
        paste0(
          "Fewer than two ",
          cluster_name,
          " clusters were available."
        )
      )
    )
  }
  if (!any(coefficient_rows)) {
    return(
      empty_result(
        "NOT_ASSESSABLE",
        "The transformed diagnostic model contains only an intercept."
      )
    )
  }
  if (model$rank < ncol(model_matrix)) {
    return(
      empty_result(
        "NOT_ASSESSABLE",
        paste0(
          "The full transformed diagnostic model is rank-deficient, so ",
          cluster_name,
          " deletion coefficients are not uniquely estimable."
        )
      )
    )
  }

  coefficient_table <- coef(summary(model))
  full_standard_errors <- rep(
    NA_real_,
    length(coefficient_names)
  )
  names(full_standard_errors) <- coefficient_names
  shared_names <- intersect(
    rownames(coefficient_table),
    coefficient_names
  )
  full_standard_errors[shared_names] <-
    coefficient_table[shared_names, "Std. Error"]
  usable_coefficients <-
    coefficient_rows &
    is.finite(full_standard_errors) &
    full_standard_errors > 0

  if (!any(usable_coefficients)) {
    return(
      empty_result(
        "NOT_ASSESSABLE",
        "No non-intercept coefficient had a positive finite standard error."
      )
    )
  }

  cell_levels <- unique(design_cell)
  cell_id <- match(design_cell, cell_levels)
  cell_first_row <- match(cell_levels, design_cell)
  cell_model_matrix <- model_matrix[
    cell_first_row,
    ,
    drop = FALSE
  ]
  cell_counts <- tabulate(
    cell_id,
    nbins = length(cell_levels)
  )
  cell_response_sums <- as.numeric(
    rowsum(
      response,
      cell_id,
      reorder = FALSE
    )
  )
  represented_matrix <- cell_model_matrix[cell_id, , drop = FALSE]
  if (!isTRUE(all.equal(
    represented_matrix,
    model_matrix,
    check.attributes = FALSE
  ))) {
    stop(
      "The supplied design cells do not uniquely identify transformed-",
      "model rows."
    )
  }

  full_coefficients <- coef(model)
  cluster_rows <- split(
    seq_len(nrow(model_matrix)),
    cluster
  )
  assessed_clusters <- 0L
  rank_or_df_failures <- 0L
  nonfinite_failures <- 0L
  maximum_standardized_shift <- -Inf
  most_influential_cluster <- NA_character_
  most_affected_coefficient <- NA_character_
  most_affected_full_estimate <- NA_real_
  most_affected_deletion_estimate <- NA_real_

  for (cluster_value in names(cluster_rows)) {
    rows <- cluster_rows[[cluster_value]]
    retained_n <- nrow(model_matrix) - length(rows)
    removed_cell_counts <- tabulate(
      cell_id[rows],
      nbins = length(cell_levels)
    )
    removed_cell_sums <- numeric(length(cell_levels))
    cluster_sums <- rowsum(
      response[rows],
      cell_id[rows],
      reorder = FALSE
    )
    removed_cell_sums[
      as.integer(rownames(cluster_sums))
    ] <- as.numeric(cluster_sums)
    reduced_cell_counts <- cell_counts - removed_cell_counts
    reduced_cell_sums <-
      cell_response_sums - removed_cell_sums
    represented_cells <- reduced_cell_counts > 0L

    if (retained_n <= ncol(model_matrix)) {
      rank_or_df_failures <- rank_or_df_failures + 1L
      next
    }

    reduced_fit <- tryCatch(
      lm.wfit(
        x = cell_model_matrix[
          represented_cells,
          ,
          drop = FALSE
        ],
        y = reduced_cell_sums[represented_cells] /
          reduced_cell_counts[represented_cells],
        w = reduced_cell_counts[represented_cells]
      ),
      error = function(error) {
        NULL
      }
    )

    if (
      is.null(reduced_fit) ||
      reduced_fit$rank < ncol(model_matrix)
    ) {
      rank_or_df_failures <- rank_or_df_failures + 1L
      next
    }

    reduced_coefficients <- reduced_fit$coefficients
    if (!all(is.finite(reduced_coefficients))) {
      nonfinite_failures <- nonfinite_failures + 1L
      next
    }

    assessed_clusters <- assessed_clusters + 1L
    standardized_shift <- abs(
      reduced_coefficients - full_coefficients
    ) / full_standard_errors
    standardized_shift[!usable_coefficients] <- NA_real_
    cluster_maximum <- max(
      standardized_shift,
      na.rm = TRUE
    )

    if (
      is.finite(cluster_maximum) &&
      cluster_maximum > maximum_standardized_shift
    ) {
      affected_index <- which.max(
        replace(
          standardized_shift,
          is.na(standardized_shift),
          -Inf
        )
      )
      maximum_standardized_shift <- cluster_maximum
      most_influential_cluster <- cluster_value
      most_affected_coefficient <-
        coefficient_names[[affected_index]]
      most_affected_full_estimate <-
        full_coefficients[[affected_index]]
      most_affected_deletion_estimate <-
        reduced_coefficients[[affected_index]]
    }
  }

  if (!is.finite(maximum_standardized_shift)) {
    maximum_standardized_shift <- NA_real_
  }
  requires_review <-
    rank_or_df_failures > 0L ||
    nonfinite_failures > 0L ||
    (
      is.finite(maximum_standardized_shift) &&
      maximum_standardized_shift >=
        standardized_shift_threshold
    )
  status <- if (requires_review) "REVIEW" else "PASS"
  most_affected_sign_changed <- if (
    is.finite(most_affected_full_estimate) &&
      is.finite(most_affected_deletion_estimate)
  ) {
    sign(most_affected_full_estimate) !=
      sign(most_affected_deletion_estimate)
  } else {
    NA
  }
  interpretation <- if (requires_review) {
    paste0(
      "At least one ",
      cluster_name,
      " deletion caused non-estimability, nonfinite coefficients, or a ",
      "coefficient movement of at least one full-model standard error in ",
      "the additive screening model. Rerun the planned final inferential ",
      "model with and without the influential cluster and report whether ",
      "the focal estimate, uncertainty, and substantive conclusion change."
    )
  } else {
    paste0(
      "No single ",
      cluster_name,
      " deletion caused non-estimability, nonfinite coefficients, or a ",
      "coefficient movement of one full-model standard error in the ",
      "additive screening model. This does not replace leave-one-cluster ",
      "checks for the planned final inferential model."
    )
  }

  list(
    status = status,
    cluster_type = cluster_name,
    clusters = length(cluster_levels),
    assessed_clusters = assessed_clusters,
    rank_or_df_failures = rank_or_df_failures,
    nonfinite_failures = nonfinite_failures,
    maximum_standardized_shift =
      maximum_standardized_shift,
    threshold = standardized_shift_threshold,
    most_influential_cluster = most_influential_cluster,
    most_affected_coefficient = most_affected_coefficient,
    most_affected_full_estimate =
      most_affected_full_estimate,
    most_affected_deletion_estimate =
      most_affected_deletion_estimate,
    most_affected_sign_changed =
      most_affected_sign_changed,
    interpretation = interpretation
  )
}


.post_boxcox_diagnostics <- function(
  rt,
  transformed,
  lambda,
  model,
  subject,
  condition,
  item,
  region,
  source_row,
  minimum_session_length = 20L,
  boundary_minimum_count = 10L,
  boundary_minimum_proportion = 0.01,
  standardized_shift_threshold = 1
) {
  n <- length(rt)
  inputs <- list(
    transformed = transformed,
    subject = subject,
    condition = condition,
    item = item,
    region = region,
    source_row = source_row
  )
  invalid_lengths <- vapply(
    inputs,
    length,
    integer(1L)
  ) != n
  if (any(invalid_lengths)) {
    stop(
      "Post-Box-Cox diagnostic vectors must have equal lengths."
    )
  }
  if (
    !is.numeric(source_row) ||
      any(!is.finite(source_row)) ||
      any(source_row < 1) ||
      any(source_row > .Machine$integer.max) ||
      any(source_row != floor(source_row)) ||
      anyDuplicated(source_row)
  ) {
    stop(
      "`source_row` must contain unique positive whole-number original ",
      "row positions."
    )
  }
  source_row <- as.integer(source_row)

  finite_values <- all(is.finite(transformed))
  transformed_variance <- if (finite_values && n > 1L) {
    var(transformed)
  } else {
    NA_real_
  }
  positive_variation <-
    is.finite(transformed_variance) &&
    transformed_variance > 0

  ordered_rows <- order(rt)
  ordered_rt <- rt[ordered_rows]
  ordered_transformed <- transformed[ordered_rows]
  distinct_steps <- diff(ordered_rt) > 0
  transformed_steps <- diff(ordered_transformed)
  step_tolerance <-
    64 * .Machine$double.eps *
    pmax(
      1,
      abs(ordered_transformed[-length(ordered_transformed)]),
      abs(ordered_transformed[-1L])
    )
  monotonicity_violations <- sum(
    distinct_steps &
      (
        !is.finite(transformed_steps) |
        transformed_steps < -step_tolerance
      )
  )

  inverse_values <- rep(NA_real_, n)
  inverse_domain_valid <- finite_values
  if (inverse_domain_valid) {
    if (abs(lambda) < sqrt(.Machine$double.eps)) {
      inverse_values <- exp(transformed)
    } else {
      inverse_base <- lambda * transformed + 1
      inverse_domain_valid <- all(
        is.finite(inverse_base) &
        inverse_base > 0
      )
      if (inverse_domain_valid) {
        inverse_values <- inverse_base^(1 / lambda)
      }
    }
  }
  inverse_relative_error <- if (
    inverse_domain_valid &&
    all(is.finite(inverse_values))
  ) {
    max(
      abs(inverse_values - rt) /
        pmax(1, abs(rt))
    )
  } else {
    Inf
  }
  inverse_tolerance <- sqrt(.Machine$double.eps)
  inverse_recovery <-
    is.finite(inverse_relative_error) &&
    inverse_relative_error <= inverse_tolerance

  model_matrix <- model.matrix(model)
  model_rank <- model$rank
  model_columns <- ncol(model_matrix)
  model_estimable <-
    model_rank == model_columns &&
    df.residual(model) >= 1L &&
    all(is.finite(coef(model)))

  session_rows <- split(seq_len(n), as.character(subject))
  assessable_sessions <- session_rows[
    lengths(session_rows) >= minimum_session_length
  ]
  session_signatures <- if (length(assessable_sessions) > 0L) {
    vapply(
      assessable_sessions,
      function(rows) {
        row_signatures <- paste(
          condition[rows],
          item[rows],
          region[rows],
          formatC(
            transformed[rows],
            digits = 17,
            format = "g"
          ),
          sep = "\034"
        )
        paste(row_signatures, collapse = "\035")
      },
      character(1L)
    )
  } else {
    character()
  }
  duplicated_signature <-
    duplicated(session_signatures) |
    duplicated(session_signatures, fromLast = TRUE)
  duplicated_session_count <- sum(duplicated_signature)
  duplicate_group_count <- if (
    duplicated_session_count > 0L
  ) {
    length(unique(session_signatures[duplicated_signature]))
  } else {
    0L
  }
  duplicate_group_provenance <- data.frame(
    GroupToken = character(),
    SessionToken = character(),
    Observations = integer(),
    FirstSourceRow = integer(),
    LastSourceRow = integer(),
    SourceRowsContiguous = logical(),
    stringsAsFactors = FALSE
  )
  if (duplicated_session_count > 0L) {
    duplicate_indices <- which(duplicated_signature)
    duplicate_indices <- duplicate_indices[
      order(
        vapply(
          assessable_sessions[duplicate_indices],
          function(rows) source_row[rows[[1L]]],
          integer(1L)
        )
      )
    ]
    duplicate_signatures <- session_signatures[duplicate_indices]
    signature_order <- unique(duplicate_signatures)
    group_number <- match(duplicate_signatures, signature_order)
    duplicate_rows <- unname(
      assessable_sessions[duplicate_indices]
    )
    provenance_source_rows <- lapply(
      duplicate_rows,
      function(rows) source_row[rows]
    )
    duplicate_group_provenance <- data.frame(
      GroupToken = sprintf("duplicate_group_%03d", group_number),
      SessionToken = sprintf(
        "duplicate_session_%03d",
        seq_along(duplicate_indices)
      ),
      Observations = lengths(duplicate_rows),
      FirstSourceRow = vapply(
        provenance_source_rows,
        function(rows) rows[[1L]],
        integer(1L)
      ),
      LastSourceRow = vapply(
        provenance_source_rows,
        function(rows) rows[[length(rows)]],
        integer(1L)
      ),
      SourceRowsContiguous = vapply(
        provenance_source_rows,
        function(rows) {
          length(rows) <= 1L || all(diff(rows) == 1L)
        },
        logical(1L)
      ),
      stringsAsFactors = FALSE
    )
    stopifnot(
      nrow(duplicate_group_provenance) ==
        duplicated_session_count,
      length(unique(duplicate_group_provenance$GroupToken)) ==
        duplicate_group_count
    )
  }
  duplicate_status <- if (length(session_signatures) < 2L) {
    "NOT_ASSESSABLE"
  } else if (duplicate_group_count > 0L) {
    "REVIEW"
  } else {
    "PASS"
  }

  maximum_transformed <- max(transformed)
  maximum_count <- sum(transformed == maximum_transformed)
  maximum_proportion <- maximum_count / n
  boundary_assessable <- n >= 100L
  boundary_threshold_count <- max(
    boundary_minimum_count,
    ceiling(boundary_minimum_proportion * n)
  )
  boundary_status <- if (!boundary_assessable) {
    "NOT_ASSESSABLE"
  } else if (maximum_count >= boundary_threshold_count) {
    "REVIEW"
  } else {
    "PASS"
  }
  design_cell <- interaction(
    condition,
    region,
    drop = TRUE,
    lex.order = TRUE
  )

  participant_deletion <- .post_boxcox_cluster_deletion(
    model = model,
    cluster = subject,
    cluster_name = "participant",
    design_cell = design_cell,
    standardized_shift_threshold =
      standardized_shift_threshold
  )
  item_deletion <- .post_boxcox_cluster_deletion(
    model = model,
    cluster = item,
    cluster_name = "item",
    design_cell = design_cell,
    standardized_shift_threshold =
      standardized_shift_threshold
  )

  format_cluster_result <- function(result) {
    if (!is.finite(result$maximum_standardized_shift)) {
      return("not assessable")
    }

    estimate_change <- if (
      is.finite(result$most_affected_full_estimate) &&
        is.finite(result$most_affected_deletion_estimate)
    ) {
      paste0(
        "; ",
        result$most_affected_coefficient,
        " changed from ",
        format(result$most_affected_full_estimate, digits = 5),
        " to ",
        format(result$most_affected_deletion_estimate, digits = 5),
        if (isTRUE(result$most_affected_sign_changed)) {
          " (sign changed)"
        } else {
          ""
        }
      )
    } else {
      ""
    }

    paste0(
      "maximum shift ",
      format(result$maximum_standardized_shift, digits = 4),
      " SE",
      estimate_change
    )
  }

  summary <- data.frame(
    Diagnostic = c(
      "Finite transformed values",
      "Positive finite variance",
      "Box-Cox monotonicity",
      "Box-Cox inverse recovery",
      "Transformed-model estimability",
      "Exact duplicated session trajectories",
      "Upper-boundary concentration",
      "Participant-deletion sensitivity",
      "Item-deletion sensitivity"
    ),
    Type = c(
      rep("Integrity", 5L),
      "Integrity",
      "Sensitivity",
      "Sensitivity",
      "Sensitivity"
    ),
    Status = c(
      if (finite_values) "PASS" else "FAIL",
      if (positive_variation) "PASS" else "FAIL",
      if (monotonicity_violations == 0L) "PASS" else "FAIL",
      if (inverse_recovery) "PASS" else "FAIL",
      if (model_estimable) "PASS" else "FAIL",
      duplicate_status,
      boundary_status,
      participant_deletion$status,
      item_deletion$status
    ),
    Result = c(
      paste0(sum(is.finite(transformed)), " of ", n, " finite"),
      format(transformed_variance, digits = 6),
      paste0(monotonicity_violations, " order reversals"),
      paste0(
        "maximum relative error ",
        format(inverse_relative_error, digits = 4)
      ),
      paste0(
        "rank ",
        model_rank,
        " of ",
        model_columns,
        "; residual df ",
        df.residual(model)
      ),
      paste0(
        duplicate_group_count,
        " duplicate groups among ",
        length(session_signatures),
        " assessable sessions"
      ),
      paste0(
        maximum_count,
        " maximum ties (",
        sprintf("%.3f%%", 100 * maximum_proportion),
        ")"
      ),
      format_cluster_result(participant_deletion),
      format_cluster_result(item_deletion)
    ),
    Criterion = c(
      "All final eligible values must be finite.",
      "Final eligible values must have positive finite variance.",
      "Unequal RTs must not reverse order after transformation.",
      paste0(
        "Maximum relative inverse error <= ",
        format(inverse_tolerance, digits = 4),
        "."
      ),
      "Full rank, finite coefficients, and at least one residual df.",
      paste0(
        "Review only exact full trajectories with at least ",
        minimum_session_length,
        " observations."
      ),
      paste0(
        "Review when the maximum contains at least ",
        boundary_minimum_count,
        " values and ",
        100 * boundary_minimum_proportion,
        "% of observations."
      ),
      paste0(
        "Review rank/nonfinite failures or shifts >= ",
        standardized_shift_threshold,
        " SE."
      ),
      paste0(
        "Review rank/nonfinite failures or shifts >= ",
        standardized_shift_threshold,
        " SE."
      )
    ),
    stringsAsFactors = FALSE
  )
  overall_status <- if (any(summary$Status == "FAIL")) {
    "FAIL"
  } else if (any(summary$Status == "REVIEW")) {
    "REVIEW"
  } else {
    "PASS"
  }

  list(
    overall_status = overall_status,
    summary = summary,
    numerical_integrity = list(
      observations = n,
      finite_values = finite_values,
      variance = transformed_variance,
      positive_variation = positive_variation,
      monotonicity_violations = monotonicity_violations,
      inverse_relative_error = inverse_relative_error,
      inverse_tolerance = inverse_tolerance,
      inverse_recovery = inverse_recovery
    ),
    model_estimability = list(
      rank = model_rank,
      columns = model_columns,
      residual_df = df.residual(model),
      finite_coefficients = all(is.finite(coef(model))),
      estimable = model_estimable
    ),
    duplicate_sessions = list(
      minimum_session_length = minimum_session_length,
      assessable_sessions = length(session_signatures),
      duplicate_groups = duplicate_group_count,
      duplicated_sessions = duplicated_session_count,
      duplicate_group_provenance = duplicate_group_provenance
    ),
    upper_boundary = list(
      transformed_maximum = maximum_transformed,
      maximum_count = maximum_count,
      maximum_proportion = maximum_proportion,
      minimum_count = boundary_minimum_count,
      minimum_proportion = boundary_minimum_proportion,
      effective_count_threshold = boundary_threshold_count
    ),
    cluster_deletion = list(
      participant = participant_deletion,
      item = item_deletion
    )
  )
}


rt_process <- function(
  data,
  subject = if (
    ".AnalysisSubject" %in% names(data)
  ) {
    data$.AnalysisSubject
  } else {
    data$Subject
  },
  condition = data$Condition,
  item = data$Item,
  region = data$Word,
  rt = data$RT,
  eligible = rep(TRUE, nrow(data)),
  verbose = FALSE,
  remove_no_comparator = FALSE,
  rt_cutoff_ms = NULL
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  if (
    !is.logical(verbose) ||
      length(verbose) != 1L ||
      is.na(verbose)
  ) {
    stop("`verbose` must be either TRUE or FALSE.")
  }
  if (
    !is.logical(remove_no_comparator) ||
      length(remove_no_comparator) != 1L ||
      is.na(remove_no_comparator)
  ) {
    stop("`remove_no_comparator` must be either TRUE or FALSE.")
  }
  if (
    !is.null(rt_cutoff_ms) &&
      (
        !is.numeric(rt_cutoff_ms) ||
          length(rt_cutoff_ms) != 1L ||
          !is.finite(rt_cutoff_ms) ||
          rt_cutoff_ms <= 0
      )
  ) {
    stop(
      "`rt_cutoff_ms` must be NULL or one finite positive number."
    )
  }

  if ("RT_BoxCox" %in% names(data)) {
    stop(
      "`data` already contains an `RT_BoxCox` column. ",
      "Remove or rename it before running the function."
    )
  }

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop(
      "The MASS package is required for the Box-Cox step. ",
      "Install MASS and try again."
    )
  }

  n <- nrow(data)
  inputs <- list(
    subject = subject,
    condition = condition,
    item = item,
    region = region,
    rt = rt,
    eligible = eligible
  )
  invalid_lengths <- vapply(
    inputs,
    length,
    integer(1L)
  ) != n

  if (any(invalid_lengths)) {
    stop(
      "Every comparison vector must have one value per data row."
    )
  }

  subject <- as.character(subject)
  condition <- as.character(condition)
  item <- as.character(item)
  region <- as.character(region)
  rt_numeric <- suppressWarnings(
    as.numeric(trimws(as.character(rt)))
  )
  rt_numeric[!is.finite(rt_numeric)] <- NA_real_
  eligible <- as.logical(eligible)
  eligible[is.na(eligible)] <- FALSE

  key_is_present <- function(values) {
    !is.na(values) & nzchar(trimws(values))
  }
  analysis_eligible <-
    eligible &
    !is.na(rt_numeric) &
    rt_numeric > 0 &
    key_is_present(subject) &
    key_is_present(condition) &
    key_is_present(item) &
    key_is_present(region)
  reference_rows <-
    analysis_eligible &
    rt_numeric >= 200 &
    rt_numeric <= 4000
  reference_count <- sum(reference_rows)
  initial_multiplier <- 1.5
  multiplier_step <- 0.5
  max_reference_removal_percent <- 5
  raw_eligible_rows <- which(eligible)
  raw_subject_is_present <- key_is_present(subject)
  raw_subjects <- subject[
    eligible & raw_subject_is_present
  ]
  raw_rows_by_subject <- table(raw_subjects)
  raw_rows_are_balanced <-
    length(raw_rows_by_subject) > 0L &&
    length(unique(as.integer(raw_rows_by_subject))) == 1L

  if (reference_count == 0L) {
    stop(
      "The multiplier cannot be selected because there are no ",
      "eligible finite RTs from 200 through 4,000 ms."
    )
  }

  row_id_name <- ".rt_process_row_id"
  while (row_id_name %in% names(data)) {
    row_id_name <- paste0(row_id_name, "_")
  }

  working_data <- data
  working_data[[row_id_name]] <- seq_len(n)
  multiplier_history <- data.frame(
    Multiplier = numeric(),
    Removed = integer(),
    ReferenceRTs = integer(),
    PercentRemoved = numeric()
  )
  selected_data <- NULL
  selected_multiplier <- NA_real_

  analysis_rt <- rt_numeric[analysis_eligible]
  upper_multiplier <-
    max(
      initial_multiplier,
      ceiling(max(analysis_rt) / min(analysis_rt)) + 1
    )
  if (!is.finite(upper_multiplier)) {
    stop(
      "The RT range is too extreme to determine a stable multiplier."
    )
  }

  count_removed_reference_rows <- function(candidate_data) {
    retained_ids <- candidate_data[[row_id_name]]

    sum(
      reference_rows &
        !(seq_len(n) %in% retained_ids)
    )
  }

  upper_data <-
    .remove_iterative_rt_outliers_at_multiplier(
      data = working_data,
      subject = subject,
      condition = condition,
      item = item,
      region = region,
      rt = rt_numeric,
      eligible = analysis_eligible,
      multiplier = upper_multiplier,
      remove_no_comparator = remove_no_comparator
    )
  structural_removal_log <- attr(
    upper_data,
    "iterative_removal_log"
  )
  structural_removal_log$Multiplier <- rep(
    upper_multiplier,
    nrow(structural_removal_log)
  )
  if (
    nrow(structural_removal_log) > 0L &&
    any(
      structural_removal_log$Reason !=
        "no_other_item_comparison"
    )
  ) {
    stop(
      "The high-multiplier structural-removal audit produced an ",
      "unexpected ratio-based removal."
    )
  }
  minimum_reference_removed <-
    count_removed_reference_rows(upper_data)
  minimum_percent_removed <-
    100 * minimum_reference_removed / reference_count

  if (
    minimum_percent_removed >=
      max_reference_removal_percent
  ) {
    stop(
      "The under-5% rule cannot be reached for this dataset. ",
      sprintf(
        "%.2f%%",
        minimum_percent_removed
      ),
      " of eligible 200-4,000 ms RTs are removed even when the ",
      "multiplier is high enough to disable ratio-based removal. ",
      "This usually means too many participant-item groups lack ",
      "the comparison rows required by the iterative rule."
    )
  }

  candidate_multiplier <- initial_multiplier
  while (candidate_multiplier <= upper_multiplier) {
    candidate_data <-
      .remove_iterative_rt_outliers_at_multiplier(
        data = working_data,
        subject = subject,
        condition = condition,
        item = item,
        region = region,
        rt = rt_numeric,
        eligible = analysis_eligible,
        multiplier = candidate_multiplier,
        remove_no_comparator = remove_no_comparator
      )

    removed_reference_count <-
      count_removed_reference_rows(candidate_data)
    percent_removed <-
      100 * removed_reference_count / reference_count

    multiplier_history <- rbind(
      multiplier_history,
      data.frame(
        Multiplier = candidate_multiplier,
        Removed = removed_reference_count,
        ReferenceRTs = reference_count,
        PercentRemoved = percent_removed
      )
    )

    if (
      percent_removed <
        max_reference_removal_percent
    ) {
      selected_data <- candidate_data
      selected_multiplier <- candidate_multiplier
      break
    }

    candidate_multiplier <-
      candidate_multiplier + multiplier_step
  }

  if (is.null(selected_data)) {
    stop(
      "The multiplier search ended unexpectedly before reaching ",
      "the under-5% rule."
    )
  }

  iterative_retained_ids <-
    selected_data[[row_id_name]]
  selected_removal_log <- attr(
    selected_data,
    "iterative_removal_log"
  )
  iterative_removed_ids <- which(
    !(seq_len(n) %in% iterative_retained_ids)
  )
  if (!identical(
    sort(selected_removal_log$Row),
    iterative_removed_ids
  )) {
    stop(
      "The iterative-removal audit log did not match the removed rows."
    )
  }
  selected_removal_log$Multiplier <- rep(
    selected_multiplier,
    nrow(selected_removal_log)
  )
  iterative_removal_reason_counts <- table(
    selected_removal_log$Reason
  )
  iterative_removed_total <-
    length(iterative_removed_ids)
  cutoff_remove <- rep(FALSE, length(iterative_retained_ids))
  if (!is.null(rt_cutoff_ms)) {
    cutoff_remove <-
      analysis_eligible[iterative_retained_ids] &
      !is.na(rt_numeric[iterative_retained_ids]) &
      rt_numeric[iterative_retained_ids] > rt_cutoff_ms
  }
  cutoff_removed_ids <-
    iterative_retained_ids[cutoff_remove]
  cutoff_removed_count <-
    length(cutoff_removed_ids)
  selected_data <- selected_data[
    !cutoff_remove,
    ,
    drop = FALSE
  ]
  attr(selected_data, "iterative_removal_log") <- NULL
  retained_ids <- selected_data[[row_id_name]]
  selected_data[[row_id_name]] <- NULL
  cleaned_data <- selected_data

  cleaned_subject <- subject[retained_ids]
  cleaned_condition <- condition[retained_ids]
  cleaned_item <- item[retained_ids]
  cleaned_region <- region[retained_ids]
  cleaned_rt <- rt_numeric[retained_ids]
  cleaned_analysis_eligible <-
    analysis_eligible[retained_ids]

  make_model_data <- function(
    response,
    condition_values,
    region_values,
    rows,
    response_name
  ) {
    model_data <- data.frame(
      Response = response[rows],
      Condition = factor(condition_values[rows]),
      Region = factor(region_values[rows])
    )
    names(model_data)[1L] <- response_name
    model_data
  }

  make_model_formula <- function(model_data, response_name) {
    predictors <- character()

    if (nlevels(model_data$Condition) > 1L) {
      predictors <- c(predictors, "Condition")
    }
    if (nlevels(model_data$Region) > 1L) {
      predictors <- c(predictors, "Region")
    }

    if (length(predictors) == 0L) {
      return(
        as.formula(
          paste(response_name, "~ 1")
        )
      )
    }

    reformulate(
      predictors,
      response = response_name
    )
  }

  raw_model_rows <- analysis_eligible
  clean_model_rows <- cleaned_analysis_eligible

  if (
    sum(raw_model_rows) < 4L ||
    sum(clean_model_rows) < 4L ||
    length(unique(cleaned_rt[clean_model_rows])) < 2L
  ) {
    stop(
      "At least four positive, varying RTs with complete condition ",
      "and region values are required for Box-Cox and Xapd diagnostics."
    )
  }

  raw_model_data <- make_model_data(
    response = rt_numeric,
    condition_values = condition,
    region_values = region,
    rows = raw_model_rows,
    response_name = "RT"
  )
  clean_model_data <- make_model_data(
    response = cleaned_rt,
    condition_values = cleaned_condition,
    region_values = cleaned_region,
    rows = clean_model_rows,
    response_name = "RT"
  )

  raw_model <- lm(
    make_model_formula(raw_model_data, "RT"),
    data = raw_model_data,
    y = TRUE
  )
  clean_model <- lm(
    make_model_formula(clean_model_data, "RT"),
    data = clean_model_data,
    y = TRUE
  )
  clean_model_summary <- coef(summary(clean_model))
  clean_model_coefficients <- data.frame(
    Term = rownames(clean_model_summary),
    Estimate = clean_model_summary[, "Estimate"],
    StandardError = clean_model_summary[, "Std. Error"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (df.residual(clean_model) < 1L) {
    stop(
      "The cleaned data do not have enough residual degrees of ",
      "freedom to estimate a stable Box-Cox lambda."
    )
  }

  lambda_initial_limit <- 3
  lambda_maximum_limit <- 12
  lambda_step <- 0.05
  lambda_limit <- lambda_initial_limit
  repeat {
    box_cox_profile <- suppressWarnings(
      MASS::boxcox(
        clean_model,
        lambda = seq(
          -lambda_limit,
          lambda_limit,
          by = lambda_step
        ),
        plotit = FALSE
      )
    )
    valid_profile <-
      is.finite(box_cox_profile$x) &
      is.finite(box_cox_profile$y)

    if (!any(valid_profile)) {
      stop(
        "Box-Cox could not calculate a finite lambda profile for ",
        "the cleaned data."
      )
    }

    valid_x <- box_cox_profile$x[valid_profile]
    valid_y <- box_cox_profile$y[valid_profile]
    lambda_index <- which.max(valid_y)
    lambda <- valid_x[lambda_index]
    profile_is_at_boundary <-
      lambda_index %in% c(1L, length(valid_x))

    if (!profile_is_at_boundary) {
      break
    }

    if (lambda_limit >= lambda_maximum_limit) {
      stop(
        "The Box-Cox profile maximum remained at the boundary of ",
        "the searched lambda range [-",
        lambda_maximum_limit,
        ", ",
        lambda_maximum_limit,
        "]. A finite interior lambda estimate could not be established."
      )
    }

    lambda_limit <- min(
      lambda_limit * 2,
      lambda_maximum_limit
    )
  }
  box_cox_profile <- data.frame(
    Lambda = valid_x,
    LogLikelihood = valid_y
  )

  box_cox_transform <- function(values, lambda_value) {
    transformed <- rep(NA_real_, length(values))
    valid <- !is.na(values) & values > 0

    if (abs(lambda_value) < sqrt(.Machine$double.eps)) {
      transformed[valid] <- log(values[valid])
    } else {
      transformed[valid] <-
        (
          values[valid]^lambda_value - 1
        ) / lambda_value
    }

    transformed
  }

  cleaned_data$RT_BoxCox <-
    box_cox_transform(cleaned_rt, lambda)
  transformed_model_data <- make_model_data(
    response = cleaned_data$RT_BoxCox,
    condition_values = cleaned_condition,
    region_values = cleaned_region,
    rows = clean_model_rows,
    response_name = "RT_BoxCox"
  )
  transformed_model <- lm(
    make_model_formula(
      transformed_model_data,
      "RT_BoxCox"
    ),
    data = transformed_model_data,
    y = TRUE
  )

  initial_subjects <- raw_subjects
  final_subjects <- cleaned_subject[clean_model_rows]

  summarize_values <- function(values) {
    values <- values[is.finite(values)]

    c(
      N = length(values),
      Mean = mean(values),
      SD = sd(values),
      Median = median(values),
      IQR = IQR(values),
      Minimum = min(values),
      Maximum = max(values)
    )
  }

  initial_values <- rt_numeric[raw_model_rows]
  clean_values <- cleaned_rt[clean_model_rows]
  transformed_values <-
    cleaned_data$RT_BoxCox[clean_model_rows]
  post_boxcox_diagnostics <- .post_boxcox_diagnostics(
    rt = clean_values,
    transformed = transformed_values,
    lambda = lambda,
    model = transformed_model,
    subject = cleaned_subject[clean_model_rows],
    condition = cleaned_condition[clean_model_rows],
    item = cleaned_item[clean_model_rows],
    region = cleaned_region[clean_model_rows],
    source_row = retained_ids[clean_model_rows]
  )

  interpret_z_b2 <- function(value) {
    if (value < -4) {
      return("Large negative value: left-skewed.")
    }
    if (value > 4) {
      return("Large positive value: right-skewed.")
    }

    paste0(
      "Between -4 and 4: not strongly classified by this ",
      "skewness rule; inspect the Q-Q plot."
    )
  }

  interpret_z_net_k2 <- function(value) {
    if (value < -1) {
      return("Large negative value: shorter-than-normal tail.")
    }
    if (value > 1) {
      return("Large positive value: longer-than-normal tail.")
    }

    paste0(
      "Between -1 and 1: not informative; visually inspect ",
      "the Q-Q plot."
    )
  }

  interpret_xapd <- function(statistic, p_value) {
    reference_text <- if (p_value > 0.01) {
      paste0(
        "Under the independent-observation reference calculation, ",
        "p > 0.01."
      )
    } else {
      paste0(
        "Under the independent-observation reference calculation, ",
        "p <= 0.01."
      )
    }

    paste0(
      reference_text,
      " Repeated participant-by-item RTs are not independent, so this ",
      "p-value is only a descriptive shape reference; it is not evidence ",
      "for or against the residual-normality assumption of the final ",
      "inferential model. The Xapd statistic is ",
      format(statistic, digits = 5),
      "; smaller values indicate closer agreement with a normal-reference ",
      "marginal shape."
    )
  }

  cutoff_display <- if (is.null(rt_cutoff_ms)) {
    NULL
  } else {
    format(
      rt_cutoff_ms,
      trim = TRUE,
      scientific = FALSE,
      digits = 15,
      big.mark = ","
    )
  }
  clean_stage_label <- if (is.null(rt_cutoff_ms)) {
    "Clean RT (no absolute cutoff)"
  } else {
    paste0("Clean RT (<= ", cutoff_display, " ms)")
  }
  clean_statistics_label <- if (is.null(rt_cutoff_ms)) {
    "Clean (ms; no absolute cutoff)"
  } else {
    paste0("Clean (ms; <= ", cutoff_display, ")")
  }
  xapd_results <- lapply(
    list(
      "Initial RT" = initial_values,
      clean_values,
      "End product (Box-Cox)" = transformed_values
    ),
    Xapd.test
  )
  names(xapd_results)[[2L]] <- clean_stage_label
  xapd_diagnostics <- data.frame(
    Stage = names(xapd_results),
    z_b2 = vapply(
      xapd_results,
      function(result) result$Z.B2,
      numeric(1L)
    ),
    z_net_k2 = vapply(
      xapd_results,
      function(result) result$Z.net.K2,
      numeric(1L)
    ),
    xapd_statistic = vapply(
      xapd_results,
      function(result) result$stat,
      numeric(1L)
    ),
    p_value = vapply(
      xapd_results,
      function(result) result$p.value,
      numeric(1L)
    ),
    stringsAsFactors = FALSE
  )
  xapd_diagnostics$z_b2_interpretation <- vapply(
    xapd_diagnostics$z_b2,
    interpret_z_b2,
    character(1L)
  )
  xapd_diagnostics$z_net_k2_interpretation <- vapply(
    xapd_diagnostics$z_net_k2,
    interpret_z_net_k2,
    character(1L)
  )
  xapd_diagnostics$distribution_shape_interpretation <- mapply(
    interpret_xapd,
    statistic = xapd_diagnostics$xapd_statistic,
    p_value = xapd_diagnostics$p_value,
    USE.NAMES = FALSE
  )

  descriptive_statistics <- data.frame(
    Statistic = names(summarize_values(initial_values)),
    `Initial (ms)` =
      as.numeric(summarize_values(initial_values)),
    Clean =
      as.numeric(summarize_values(clean_values)),
    Transformed =
      as.numeric(summarize_values(transformed_values)),
    check.names = FALSE
  )
  names(descriptive_statistics)[[3L]] <- clean_statistics_label

  selected_percent <-
    tail(multiplier_history$PercentRemoved, 1L)
  selected_removed <-
    tail(multiplier_history$Removed, 1L)

  if (abs(lambda) < sqrt(.Machine$double.eps)) {
    transformation_text <- "natural logarithm: log(RT)"
  } else {
    transformation_text <- paste0(
      "(RT^",
      format(lambda, digits = 3),
      " - 1) / ",
      format(lambda, digits = 3)
    )
  }

  cat("\nRT cleaning overview\n")
  cat("--------------------\n")
  cat(
    "No-other-item-comparator policy: ",
    if (remove_no_comparator) {
      "remove rows lacking a comparator"
    } else {
      "retain rows lacking a comparator"
    },
    ".\n",
    sep = ""
  )
  cat(
    "Raw experimental rows supplied: ",
    length(raw_eligible_rows),
    ". No experimental rows were excluded before iterative ",
    "processing.\n",
    sep = ""
  )
  if (raw_rows_are_balanced) {
    rows_per_participant <-
      unique(as.integer(raw_rows_by_subject))
    cat(
      "Raw-row check: ",
      length(raw_rows_by_subject),
      " participants x ",
      rows_per_participant,
      " word positions per participant = ",
      length(raw_eligible_rows),
      " rows.\n",
      sep = ""
    )
  } else if (length(raw_rows_by_subject) > 0L) {
    cat(
      "Raw-row check: participants contributed ",
      min(raw_rows_by_subject),
      " to ",
      max(raw_rows_by_subject),
      " word rows each; their rows sum to ",
      sum(raw_rows_by_subject),
      ".\n",
      sep = ""
    )
  }
  cat(
    "Multiplier used: ",
    selected_multiplier,
    ". ",
    sep = ""
  )
  if (selected_multiplier == initial_multiplier) {
    cat(
      "Multiplier ",
      initial_multiplier,
      " was retained because it removed ",
      sprintf("%.2f%%", selected_percent),
      " of RTs from 200 through 4,000 ms, which is below ",
      max_reference_removal_percent,
      "%.\n",
      sep = ""
    )
  } else {
    cat(
      "The function started at ",
      initial_multiplier,
      " and increased the multiplier in steps of ",
      multiplier_step,
      " because lower settings removed ",
      max_reference_removal_percent,
      "% or more. Multiplier ",
      selected_multiplier,
      " removed ",
      sprintf("%.2f%%", selected_percent),
      ", which is below ",
      max_reference_removal_percent,
      "%.\n",
      sep = ""
    )
  }
  cat(
    "Multiplier checks: ",
    paste0(
      "m",
      multiplier_history$Multiplier,
      " = ",
      sprintf(
        "%.2f%%",
        multiplier_history$PercentRemoved
      ),
      collapse = "; "
    ),
    ".\n",
    sep = ""
  )
  cat(
    "Reference-range RTs removed: ",
    selected_removed,
    " of ",
    reference_count,
    ".\n",
    sep = ""
  )
  cat(
    "Total rows removed by iterative comparisons: ",
    iterative_removed_total,
    ".\n",
    sep = ""
  )
  if (is.null(rt_cutoff_ms)) {
    cat(
      "Absolute RT cutoff: not applied. No rows were removed by an ",
      "absolute cutoff before Box-Cox estimation.\n",
      sep = ""
    )
  } else {
    cat(
      "Additional eligible RTs removed strictly above the user-supplied ",
      cutoff_display,
      "-ms cutoff: ",
      cutoff_removed_count,
      ". The cutoff was applied after iterative removal and before ",
      "Box-Cox estimation.\n",
      sep = ""
    )
  }
  cat(
    "Box-Cox lambda: ",
    format(lambda, digits = 3),
    ". Transformation applied: ",
    transformation_text,
    ".\n",
    sep = ""
  )
  cat(
    "Participants: ",
    length(unique(initial_subjects)),
    " initially; ",
    length(unique(final_subjects)),
    " in the end product.\n\n",
    sep = ""
  )
  excluded_analysis_rows <- sum(eligible & !analysis_eligible)
  if (excluded_analysis_rows > 0L) {
    cat(
      "Analysis note: ",
      excluded_analysis_rows,
      " otherwise eligible rows lacked a positive finite RT or ",
      "complete participant, condition, item, or region keys. ",
      "They were retained in the returned data but were not used ",
      "to select the multiplier or estimate lambda.\n\n",
      sep = ""
    )
  }
  cat("Basic descriptive statistics\n")
  printable_statistics <- descriptive_statistics
  numeric_columns <- vapply(
    printable_statistics,
    is.numeric,
    logical(1L)
  )
  printable_statistics[numeric_columns] <- lapply(
    printable_statistics[numeric_columns],
    function(values) {
      round(values, 2)
    }
  )
  print(
    printable_statistics,
    row.names = FALSE
  )
  cat(
    "\nThe transformed values are on a new scale, so their means ",
    "and ranges should not be read as milliseconds.\n"
  )
  cat("\nPost-Box-Cox integrity and sensitivity diagnostics\n")
  cat("--------------------------------------------------\n")
  cat(
    "Overall status: ",
    post_boxcox_diagnostics$overall_status,
    ". PASS indicates that no severe issue was detected; REVIEW ",
    "identifies a conservative sensitivity finding that requires ",
    "explanation but does not invalidate or remove data; FAIL ",
    "indicates a mathematical or estimability violation.\n",
    sep = ""
  )
  printable_post_boxcox <- post_boxcox_diagnostics$summary
  if (!verbose) {
    printable_post_boxcox <- printable_post_boxcox[
      printable_post_boxcox$Status %in% c("REVIEW", "FAIL"),
      ,
      drop = FALSE
    ]
  }
  if (nrow(printable_post_boxcox) > 0L) {
    print(
      printable_post_boxcox[
        ,
        c("Diagnostic", "Type", "Status", "Result")
      ],
      row.names = FALSE
    )
  } else {
    cat("No REVIEW or FAIL diagnostic findings were detected.\n")
  }
  if (!verbose) {
    cat(
      "Only REVIEW and FAIL findings are shown. Use ",
      "`verbose = TRUE` to print every diagnostic result.\n"
    )
  }
  cat(
    "\nThese checks use only final analysis-eligible RT_BoxCox ",
    "values. They do not remove rows and are not omnibus tests of ",
    "normality, skewness, heteroskedasticity, or serial dependence.\n"
  )
  cat("\nDescriptive Xapd distribution-shape diagnostics\n")
  cat("-----------------------------------------------\n")
  cat(
    "These statistics describe marginal RT shape. Their reference ",
    "p-values assume independent observations and are not calibrated ",
    "for crossed repeated RTs. Evaluate normality and constant variance ",
    "from residuals of the planned final inferential model.\n"
  )
  printable_xapd <- xapd_diagnostics[
    ,
    c(
      "Stage",
      "z_b2",
      "z_net_k2",
      "xapd_statistic",
      "p_value"
    )
  ]
  printable_xapd$z_b2 <- round(
    printable_xapd$z_b2,
    3
  )
  printable_xapd$z_net_k2 <- round(
    printable_xapd$z_net_k2,
    3
  )
  printable_xapd$xapd_statistic <- round(
    printable_xapd$xapd_statistic,
    3
  )
  printable_xapd$p_value <- format.pval(
    printable_xapd$p_value,
    digits = 3,
    eps = .Machine$double.eps
  )
  print(
    printable_xapd,
    row.names = FALSE
  )
  xapd_stage_summaries <- vapply(
    seq_len(nrow(xapd_diagnostics)),
    function(row_index) {
      paste0(
        xapd_diagnostics$Stage[[row_index]],
        " (z_b2 = ",
        format(
          xapd_diagnostics$z_b2[[row_index]],
          digits = 5
        ),
        ", z_net_k2 = ",
        format(
          xapd_diagnostics$z_net_k2[[row_index]],
          digits = 5
        ),
        ", Xapd = ",
        format(
          xapd_diagnostics$xapd_statistic[[row_index]],
          digits = 5
        ),
        ", p = ",
        format.pval(
          xapd_diagnostics$p_value[[row_index]],
          digits = 3,
          eps = .Machine$double.eps
        ),
        ")"
      )
    },
    character(1L)
  )
  xapd_change <-
    xapd_diagnostics$xapd_statistic[[3L]] -
    xapd_diagnostics$xapd_statistic[[1L]]
  xapd_overall_summary <- if (xapd_change < 0) {
    paste0(
      "Overall, the final Xapd statistic is smaller than the ",
      "initial statistic, indicating closer agreement with a ",
      "normal-reference marginal shape by this descriptive measure."
    )
  } else if (xapd_change > 0) {
    paste0(
      "Overall, the final Xapd statistic is larger than the ",
      "initial statistic, indicating less agreement with a ",
      "normal-reference marginal shape by this descriptive measure."
    )
  } else {
    paste0(
      "Overall, the initial and final Xapd statistics are equal ",
      "by this descriptive shape measure."
    )
  }
  cat(
    "\nXapd summary: ",
    paste(xapd_stage_summaries, collapse = "; "),
    ". ",
    xapd_overall_summary,
    "\n",
    sep = ""
  )

  qq_values <- function(values, maximum_points = 5000L) {
    values <- values[is.finite(values)]

    if (length(values) <= maximum_points) {
      return(values)
    }

    sorted_values <- sort(values)
    sample_rows <- unique(
      round(
        seq(
          1,
          length(sorted_values),
          length.out = maximum_points
        )
      )
    )
    sorted_values[sample_rows]
  }

  plot_model_diagnostics <- function(
    model,
    title,
    maximum_points = 5000L
  ) {
    fitted_values <- fitted(model)
    residual_values <- residuals(model)
    row_count <- length(fitted_values)
    plot_rows <- if (row_count <= maximum_points) {
      seq_len(row_count)
    } else {
      unique(
        round(
          seq(
            1,
            row_count,
            length.out = maximum_points
          )
        )
      )
    }

    plot(
      fitted_values[plot_rows],
      residual_values[plot_rows],
      pch = 16,
      cex = 0.35,
      col = "steelblue",
      xlab = "Fitted value",
      ylab = "Residual",
      main = title
    )
    abline(h = 0, col = "red", lwd = 1.2)
  }

  old_parameters <- par(no.readonly = TRUE)
  on.exit(par(old_parameters), add = TRUE)
  par(
    mfrow = c(2, 2),
    mar = c(4, 4, 3, 1)
  )

  raw_qq <- qq_values(
    rt_numeric[raw_model_rows]
  )
  transformed_qq <- qq_values(
    cleaned_data$RT_BoxCox[
      clean_model_rows
    ]
  )
  qqnorm(
    raw_qq,
    pch = 16,
    cex = 0.45,
    col = "steelblue",
    main = "Raw RT: normal Q-Q",
    xlab = "Theoretical quantiles",
    ylab = "Observed RT"
  )
  qqline(raw_qq, col = "red", lwd = 1.2)
  qqnorm(
    transformed_qq,
    pch = 16,
    cex = 0.45,
    col = "steelblue",
    main = "End product: normal Q-Q",
    xlab = "Theoretical quantiles",
    ylab = "Transformed RT"
  )
  qqline(
    transformed_qq,
    col = "red",
    lwd = 1.2
  )
  plot_model_diagnostics(
    raw_model,
    "Raw RT: residuals vs fitted"
  )
  plot_model_diagnostics(
    transformed_model,
    "End product: residuals vs fitted"
  )
  cat(
    "\nResidual-versus-fitted panels use the function's additive ",
    "condition-and-region screening models. They do not establish the ",
    "assumptions of a study's final inferential model; inspect that ",
    "model's conditional residuals and influence diagnostics directly.\n"
  )

  cutoff_token <- if (is.null(rt_cutoff_ms)) {
    NULL
  } else {
    format(
      rt_cutoff_ms,
      trim = TRUE,
      scientific = FALSE,
      digits = 15
    )
  }
  cutoff_stage <- if (is.null(cutoff_token)) {
    character()
  } else {
    paste0(cutoff_token, "_ms_cutoff")
  }
  cutoff_reason <- if (is.null(cutoff_token)) {
    character()
  } else {
    paste0("above_", cutoff_token, "_ms")
  }
  report <- list(
    multiplier = selected_multiplier,
    no_comparator_policy = if (remove_no_comparator) {
      "remove_rows_without_other_item_comparator"
    } else {
      "retain_rows_without_other_item_comparator"
    },
    initial_multiplier = initial_multiplier,
    multiplier_step = multiplier_step,
    max_reference_removal_percent =
      max_reference_removal_percent,
    removal_history = multiplier_history,
    reference_range = c(200, 4000),
    reference_removed = selected_removed,
    reference_total = reference_count,
    percent_reference_removed = selected_percent,
    iterative_removed_total = iterative_removed_total,
    iterative_removed_rows = iterative_removed_ids,
    iterative_removal_log = selected_removal_log,
    iterative_removal_reason_counts =
      iterative_removal_reason_counts,
    structural_audit_multiplier = upper_multiplier,
    structural_removed_total =
      nrow(structural_removal_log),
    structural_removed_rows =
      structural_removal_log$Row,
    structural_reference_removed =
      minimum_reference_removed,
    structural_percent_reference_removed =
      minimum_percent_removed,
    structural_removal_log =
      structural_removal_log,
    rt_cutoff_ms = rt_cutoff_ms,
    rt_cutoff_origin = if (is.null(rt_cutoff_ms)) {
      "not_applied"
    } else {
      "user_supplied"
    },
    cutoff_removed = cutoff_removed_count,
    cutoff_removed_rows = cutoff_removed_ids,
    total_removed =
      iterative_removed_total + cutoff_removed_count,
    removed_rows = data.frame(
      Row = c(
        selected_removal_log$Row,
        cutoff_removed_ids
      ),
      Stage = c(
        rep(
          "iterative_relative_rt",
          nrow(selected_removal_log)
        ),
        rep(
          cutoff_stage,
          length(cutoff_removed_ids)
        )
      ),
      Reason = c(
        selected_removal_log$Reason,
        rep(
          cutoff_reason,
          length(cutoff_removed_ids)
        )
      ),
      Multiplier = c(
        selected_removal_log$Multiplier,
        rep(
          NA_integer_,
          length(cutoff_removed_ids)
        )
      ),
      CutoffMs = c(
        rep(
          NA_real_,
          nrow(selected_removal_log)
        ),
        rep(
          rt_cutoff_ms,
          length(cutoff_removed_ids)
        )
      ),
      stringsAsFactors = FALSE
    ),
    lambda = lambda,
    transformation = transformation_text,
    box_cox_profile = box_cox_profile,
    lambda_search_range = c(-lambda_limit, lambda_limit),
    lambda_search_step = lambda_step,
    lambda_profile_at_boundary = profile_is_at_boundary,
    software_versions = c(
      R = paste(
        R.version$major,
        R.version$minor,
        sep = "."
      ),
      MASS = as.character(
        utils::packageVersion("MASS")
      )
    ),
    initial_participants =
      length(unique(initial_subjects)),
    final_participants =
      length(unique(final_subjects)),
    raw_eligible_rows = length(raw_eligible_rows),
    raw_rows_by_participant = raw_rows_by_subject,
    raw_rows_are_balanced = raw_rows_are_balanced,
    excluded_analysis_rows = excluded_analysis_rows,
    clean_model_coefficients = clean_model_coefficients,
    descriptive_statistics = descriptive_statistics,
    xapd_diagnostics = xapd_diagnostics,
    post_boxcox_diagnostics = post_boxcox_diagnostics
  )
  attr(cleaned_data, "rt_process_report") <- report
  attr(cleaned_data, "iterative_remove_report") <- report

  invisible(cleaned_data)
}


remove_iterative_rt_outliers <- rt_process
iterative_remove <- rt_process