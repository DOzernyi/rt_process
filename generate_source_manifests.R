script_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", script_arguments, value = TRUE)
script_path <- if (length(file_argument) == 1L) {
  normalizePath(sub("^--file=", "", file_argument))
} else {
  normalizePath("maze_preprocessing/generate_source_manifests.R")
}
script_directory <- dirname(script_path)

source(file.path(script_directory, "source_audit_helpers.R"))

csv_files <- sort(
  list.files(
    path = file.path(script_directory, "data"),
    pattern = "\\.csv$",
    full.names = TRUE
  )
)
datasets <- lapply(csv_files, read_maze_csv)
names(datasets) <- tools::file_path_sans_ext(basename(csv_files))
manifests <- build_source_manifests(csv_files, datasets)

utils::write.csv(
  manifests$checksum,
  file.path(script_directory, "source_data_manifest.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  manifests$schema,
  file.path(script_directory, "source_schema_manifest.csv"),
  row.names = FALSE,
  na = ""
)

cat(
  "Wrote ",
  nrow(manifests$checksum),
  " source checksum rows and ",
  nrow(manifests$schema),
  " schema rows.\n",
  sep = ""
)