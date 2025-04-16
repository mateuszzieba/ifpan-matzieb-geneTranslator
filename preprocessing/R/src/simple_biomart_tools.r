library(tidyverse)
library(biomaRt)
library(magrittr)

download_biomart_genes <- function(
    version = "v110",
    species = "hsapiens",
    biotype = NULL,
    output_dir = "ensembl_gene_info",
    save_rds = FALSE,
    save_tsv = FALSE) {
  # Check if biomaRt package is available
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("biomaRt package is required.")
  }

  # Extract numeric version (e.g., "v110" -> 110)
  version_num <- as.numeric(sub("v", "", version))
  if (is.na(version_num)) stop("❌ Invalid version format. Use e.g. 'v110'.")

  # Create output directory if needed
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  message(sprintf("🔄 Connecting to Ensembl %s (version %d)...", version, version_num))

  # Try to download data
  result <- tryCatch(
    {
      ensembl <- biomaRt::useEnsembl(
        biomart = "genes",
        dataset = paste0(species, "_gene_ensembl"),
        version = version_num
      )

      attributes <- c(
        "ensembl_gene_id",
        "external_gene_name",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
        "gene_biotype"
      )

      if (is.null(biotype)) {
        biomaRt::getBM(attributes = attributes, mart = ensembl)
      } else {
        biomaRt::getBM(attributes = attributes, filters = "biotype", values = biotype, mart = ensembl)
      }
    },
    error = function(e) {
      stop(sprintf("❌ Failed to retrieve data from Ensembl %s: %s", version, e$message))
    }
  )

  if (nrow(result) == 0) {
    warning("⚠️ Retrieved 0 rows. Possibly no data matched the filters.")
  } else {
    message(sprintf("✅ Retrieved %d genes.", nrow(result)))
  }

  # Save as .rds if enabled
  if (save_rds) {
    rds_path <- file.path(output_dir, sprintf("gene_info_%s.rds", version))
    saveRDS(result, rds_path)
    message(sprintf("💾 Saved RDS to: %s", rds_path))
  }

  # Save as .tsv if enabled
  if (save_tsv) {
    tsv_path <- file.path(output_dir, sprintf("gene_info_%s.tsv", version))
    readr::write_tsv(result, tsv_path)
    message(sprintf("📄 Also saved TSV to: %s", tsv_path))
  }

  return(result)
}

get_gene_aliases <- function(gene_symbol, species = "hsapiens", version = NULL, retries = 3, wait_seconds = 2) {
  # This function gets gene aliases (synonyms) for one gene
  # It connects to Ensembl using biomaRt and looks for 'external_synonym'
  # It retries if there's a connection problem

  # Load the biomaRt package
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("biomaRt package is required.")
  }

  # Create the dataset name, like "hsapiens_gene_ensembl"
  dataset_name <- paste0(species, "_gene_ensembl")

  # Try to connect to Ensembl
  mart <- tryCatch(
    {
      if (!is.null(version)) {
        version_num <- as.numeric(sub("v", "", version))
        if (is.na(version_num)) stop("Invalid version format. Use for example 'v110'.")
        biomaRt::useEnsembl(biomart = "genes", dataset = dataset_name, version = version_num)
      } else {
        biomaRt::useEnsembl(biomart = "genes", dataset = dataset_name)
      }
    },
    error = function(e) {
      warning(paste("Failed to connect to Ensembl:", e$message))
      return(NULL)
    }
  )

  if (is.null(mart)) {
    return(NA)
  }

  # Try to download the data, with retries if it fails
  attempt <- 1
  while (attempt <= retries) {
    res <- tryCatch(
      {
        biomaRt::getBM(
          attributes = c("external_synonym", "hgnc_symbol"),
          filters = "hgnc_symbol",
          values = gene_symbol,
          mart = mart
        )
      },
      error = function(e) {
        warning(paste("Try", attempt, "failed for", gene_symbol, ":", e$message))
        return(NULL)
      }
    )

    if (!is.null(res)) break # worked!
    Sys.sleep(wait_seconds)
    attempt <- attempt + 1
  }

  # If results are OK, return the list of aliases
  if (!is.null(res) && nrow(res) > 0) {
    aliases <- unique(na.omit(res$external_synonym))
  } else {
    aliases <- NA
  }

  return(aliases)
}

get_aliases_df <- function(
    gene_symbols,
    species = "hsapiens",
    version = NULL,
    output_dir = NULL,
    save_rds = FALSE,
    save_tsv = FALSE,
    retries = 3,
    wait_seconds = 2) {
  # This function goes through a list of gene symbols and gets aliases for each one
  # It uses get_gene_aliases() inside a loop
  # You can also save the results to a file if you want

  if (!requireNamespace("progress", quietly = TRUE)) install.packages("progress")
  library(progress)
  library(dplyr)

  # Create a progress bar so we know it's working
  pb <- progress_bar$new(
    total = length(gene_symbols),
    format = "  Getting [:bar] :current/:total (:percent) :gene",
    clear = FALSE,
    width = 60
  )

  # Go through each gene and collect aliases
  results <- lapply(gene_symbols, function(gene) {
    pb$tick(tokens = list(gene = gene))
    aliases <- get_gene_aliases(
      gene_symbol = gene,
      species = species,
      version = version,
      retries = retries,
      wait_seconds = wait_seconds
    )
    alias_str <- if (all(is.na(aliases))) NA else paste(aliases, collapse = "|")
    data.frame(gene_symbol = gene, aliases = alias_str)
  })

  # Combine all rows into one table
  result_df <- bind_rows(results)

  # Save to files if needed
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    version_label <- ifelse(is.null(version), "latest", version)

    if (save_rds) {
      saveRDS(result_df, file.path(output_dir, paste0("gene_aliases_", version_label, ".rds")))
    }

    if (save_tsv) {
      if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
      readr::write_tsv(result_df, file.path(output_dir, paste0("gene_aliases_", version_label, ".tsv")))
    }
  }

  return(result_df)
}


translate_genes_to_human <- function(gene_symbols, from_species = "mouse", version = NULL) {
  # Supported species: "mouse", "rat"
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("❌ Package 'biomaRt' is required.")
  }

  # Map species to Ensembl dataset and attribute
  species_map <- list(
    mouse = list(dataset = "mmusculus_gene_ensembl", attr_symbol = "mgi_symbol"),
    rat   = list(dataset = "rnorvegicus_gene_ensembl", attr_symbol = "rgd_symbol")
  )

  if (!from_species %in% names(species_map)) {
    stop("❌ Only 'mouse' and 'rat' are supported.")
  }

  dataset_info <- species_map[[from_species]]

  # Connect to Ensembl (with or without version)
  mart <- tryCatch(
    {
      if (!is.null(version)) {
        version_num <- as.numeric(sub("v", "", version))
        biomaRt::useEnsembl("genes", dataset = dataset_info$dataset, version = version_num)
      } else {
        biomaRt::useEnsembl("genes", dataset = dataset_info$dataset)
      }
    },
    error = function(e) {
      stop(paste("❌ Failed to connect to Ensembl:", e$message))
    }
  )

  # Query orthologs
  result <- biomaRt::getBM(
    attributes = c(dataset_info$attr_symbol, "hsapiens_homolog_associated_gene_name"),
    filters = dataset_info$attr_symbol,
    values = gene_symbols,
    mart = mart
  )

  # Rename columns for clarity
  colnames(result) <- c("original_symbol", "human_symbol")

  # Remove duplicates
  result <- unique(result)

  return(result)
}

create_species_to_human_mapper <- function(
  from_species = "mouse",
  version = NULL,
  save_tsv = FALSE,
  output_path = NULL
) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) stop("Package 'biomaRt' is required.")

  # Map species to dataset
  dataset_map <- list(
    mouse = "mmusculus_gene_ensembl",
    rat = "rnorvegicus_gene_ensembl",
    zebrafish = "drerio_gene_ensembl",
    pig = "sscrofa_gene_ensembl"
  )

  if (!from_species %in% names(dataset_map)) {
    stop("❌ Supported species: ", paste(names(dataset_map), collapse = ", "))
  }

  # Connect to Ensembl
  dataset_name <- dataset_map[[from_species]]
  mart <- tryCatch({
    if (!is.null(version)) {
      version_num <- as.numeric(sub("v", "", version))
      biomaRt::useEnsembl("genes", dataset = dataset_name, version = version_num)
    } else {
      biomaRt::useEnsembl("genes", dataset = dataset_name)
    }
  }, error = function(e) {
    stop("❌ Failed to connect to Ensembl: ", e$message)
  })

  # Use ensembl_gene_id instead of species-specific symbol
  message("📥 Downloading ortholog mapping from ", from_species, " to human...")

  result <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_associated_gene_name"),
    mart = mart
  )

  colnames(result) <- c("ensembl_gene_id", "original_symbol", "human_symbol")
  result <- unique(result)
  result <- result[result$original_symbol != "" & result$human_symbol != "", ]

  # Save if needed
  if (save_tsv) {
    if (is.null(output_path)) {
      output_path <- paste0(from_species, "_to_human_mapper.tsv")
    }
    write.table(result, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("✅ Saved mapper to: ", output_path)
  }

  return(result)
}
create_species_to_human_mapper <- function(
  from_species = "mouse",
  version = NULL,
  save_tsv = FALSE,
  output_path = NULL,
  clean = FALSE
) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) stop("Package 'biomaRt' is required.")

  dataset_map <- list(
    mouse = "mmusculus_gene_ensembl",
    rat = "rnorvegicus_gene_ensembl",
    zebrafish = "drerio_gene_ensembl",
    pig = "sscrofa_gene_ensembl"
  )

  if (!from_species %in% names(dataset_map)) {
    stop("❌ Supported species: ", paste(names(dataset_map), collapse = ", "))
  }

  dataset_name <- dataset_map[[from_species]]
  version_label <- if (!is.null(version)) version else "latest"

  species_prefix <- paste0(from_species, "_gene") # np. mouse_gene
  ens_prefix <- paste0(from_species, "_ensembl_gene_id")
  biotype_col <- paste0(from_species, "_gene_biotype")
  symbol_col <- paste0(from_species, "_gene_symbol")

  # Connect to Ensembl (source species)
  mart_species <- tryCatch({
    if (!is.null(version)) {
      version_num <- as.numeric(sub("v", "", version))
      biomaRt::useEnsembl("genes", dataset = dataset_name, version = version_num)
    } else {
      biomaRt::useEnsembl("genes", dataset = dataset_name)
    }
  }, error = function(e) {
    stop("❌ Failed to connect to Ensembl: ", e$message)
  })

  message("📥 Downloading gene info for ", from_species, "...")
  gene_info <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
    mart = mart_species
  )

  message("📥 Downloading homolog info...")
  homolog_info <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_ensembl_gene"),
    mart = mart_species
  )

  species_data <- merge(gene_info, homolog_info, by = "ensembl_gene_id", all.x = TRUE)

  # Connect to latest human Ensembl
  mart_human <- tryCatch({
    biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
  }, error = function(e) {
    stop("❌ Failed to connect to Ensembl for human data: ", e$message)
  })

  message("📥 Downloading human biotype info...")
  human_biotypes <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "gene_biotype"),
    mart = mart_human
  )
  colnames(human_biotypes) <- c("hsapiens_homolog_ensembl_gene", "human_gene_biotype")

  merged_data <- merge(species_data, human_biotypes, by = "hsapiens_homolog_ensembl_gene", all.x = TRUE)

  # Rename columns
  colnames(merged_data) <- c(
    "human_ensembl_gene_id",
    ens_prefix,
    symbol_col,
    biotype_col,
    "human_gene_symbol",
    "human_gene_biotype"
  )

  # 🧹 Optional cleaning
  if (clean) {
    merged_data <- merged_data[merged_data[[symbol_col]] != "" & merged_data[["human_gene_symbol"]] != "", ]
  }

  merged_data <- merged_data %>%
    dplyr::select(
      tidyselect::any_of(c(
        ens_prefix,          # np. mouse_ensembl_gene_id
        symbol_col,          # np. mouse_gene_symbol
        biotype_col,         # np. mouse_gene_biotype
        "human_ensembl_gene_id",
        "human_gene_symbol",
        "human_gene_biotype"
      ))
    )  

  # Save if needed
  if (save_tsv) {
    if (is.null(output_path)) {
      output_path <- paste0(from_species, "_to_human_mapper_", version_label, if (clean) "_clean", ".tsv")
    }
    write.table(merged_data, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("✅ Saved mapper to: ", output_path)
  }

  return(merged_data)
}


# mapper_mouse <- create_species_to_human_mapper(
#   from_species = "mouse",
#   version = "v110",
#   save_tsv = TRUE,
#   clean = TRUE
# )

