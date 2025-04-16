# preprocessing/R/src/ncbi_gene_info.R

download_ncbi_gene_info <- function(
  organism_tax_id = 9606,
  save_path = "data/ncbi_gene_info.tsv",
  overwrite = FALSE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Pakiet 'data.table' jest wymagany.")
  }

  if (!dir.exists(dirname(save_path))) {
    dir.create(dirname(save_path), recursive = TRUE)
  }

  if (file.exists(save_path) && !overwrite) {
    message("✅ Plik już istnieje: ", save_path)
    return(invisible(save_path))
  }

  url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
  message("⬇️  Pobieranie pliku gene_info.gz z NCBI...")

  temp_file <- tempfile(fileext = ".gz")
  utils::download.file(url, destfile = temp_file, mode = "wb")

  message("📖 Wczytywanie i filtrowanie danych...")
  gene_info <- data.table::fread(temp_file, sep = "\t", header = TRUE)

  # Filtrowanie tylko dla Homo sapiens (tax_id = 9606)
  gene_info_human <- gene_info[`#tax_id` == organism_tax_id]

  data.table::fwrite(gene_info_human, file = save_path, sep = "\t")
  message("📁 Zapisano do: ", save_path)

  return(invisible(save_path))
}
