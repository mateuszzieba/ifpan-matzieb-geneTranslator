# preprocessing/R/src/hgnc_download.R

download_hgnc_complete_set <- function(
  save_path = "data/hgnc_complete_set.tsv",
  overwrite = FALSE
) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required.")
  }

  if (!dir.exists(dirname(save_path))) {
    dir.create(dirname(save_path), recursive = TRUE)
  }

  if (file.exists(save_path) && !overwrite) {
    message("✅ Plik już istnieje: ", save_path)
    return(invisible(save_path))
  }

  url <- "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"
  message("⬇️  Pobieranie HGNC complete set...")

  resp <- httr::GET(url)
  if (httr::status_code(resp) != 200) {
    stop("❌ Nie udało się pobrać pliku HGNC (HTTP ", httr::status_code(resp), ")")
  }

  writeBin(httr::content(resp, as = "raw"), save_path)
  message("📁 Zapisano jako: ", save_path)

  return(invisible(save_path))
}