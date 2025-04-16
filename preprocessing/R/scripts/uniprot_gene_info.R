require(dplyr)
require(tidyr)
require(purrr)
require(stringr)
require(magrittr)
require(tibble)

# preprocessing/R/src/uniprot_gene_info.R

get_uniprot_gene_info <- function(
  reviewed_only = TRUE,
  limit = 500,
  sleep = 3
) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Packages 'httr' and 'jsonlite' are required.")
  }

  base_url <- "https://rest.uniprot.org/uniprotkb/search"

  query <- "organism_id:9606"
  if (reviewed_only) {
    query <- paste(query, "AND reviewed:true")
  }

  fields <- paste(
    c("accession", "id", "gene_names", "organism_name", 
      "protein_name", "xref_ensembl", "xref_geneid"), collapse = ","
  )

  url <- paste0(
    base_url,
    "?query=", URLencode(query),
    "&format=json",
    "&fields=", URLencode(fields),
    "&size=", limit
  )

  message("⏳ Querying UniProt: ", url)

  resp <- httr::GET(url)
  if (httr::status_code(resp) != 200) {
    stop("Failed to fetch data from UniProt API")
  }

  content <- httr::content(resp, as = "text", encoding = "UTF-8")
  json <- jsonlite::fromJSON(content)

  if (!"results" %in% names(json)) {
    stop("No results found.")
  }

  # Przekształcenie danych
  df <- json$results
  out <- data.frame(
    uniprot_id = sapply(df$primaryAccession, as.character),
    uniprot_symbol = sapply(df$uniProtkbId, as.character),
    gene_names = sapply(df$genes, function(g) paste(unique(g$geneName$value), collapse = ",")),
    protein_name = sapply(df$proteinDescription$recommendedName$fullName$value, as.character),
    ensembl_ids = sapply(df$uniProtKBCrossReferences, function(x) {
      ens <- x[x$database == "Ensembl", "id"]
      if (length(ens) == 0) return(NA)
      return(paste(unique(ens), collapse = ","))
    }),
    ncbi_geneid = sapply(df$uniProtKBCrossReferences, function(x) {
      ncbi <- x[x$database == "GeneID", "id"]
      if (length(ncbi) == 0) return(NA)
      return(paste(unique(ncbi), collapse = ","))
    }),
    stringsAsFactors = FALSE
  )

  Sys.sleep(sleep)  # szacunek dla API :)
  return(json)
}

# preprocessing/R/src/uniprot_detailed_info.R
get_uniprot_detailed_info <- function(
  uniprot_ids,
  sleep = 1
) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Packages 'httr' and 'jsonlite' are required.")
  }

  base_url <- "https://rest.uniprot.org/uniprotkb/"
  result_list <- list()

  for (id in uniprot_ids) {
    url <- paste0(base_url, id, ".json")
    message("🔍 Fetching: ", id)

    resp <- httr::GET(url)
    if (httr::status_code(resp) != 200) {
      message("⚠️  Failed for ID: ", id)
      next
    }

    entry <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))

    # GENES
    gene_names <- NA
    if (!is.null(entry$genes)) {
      primary <- tryCatch(entry$genes$geneName$value[1], error = function(e) NA)
      aliases <- tryCatch(
        unlist(lapply(entry$genes$synonyms, function(s) s$value)),
        error = function(e) character(0)
      )
      gene_names <- paste(unique(na.omit(c(primary, aliases))), collapse = ", ")
    }

    # FUNCTION COMMENT
    function_comment <- NA
    if ("comments" %in% names(entry)) {
      fun_comments <- entry$comments[entry$comments$commentType == "FUNCTION", ]
      if (nrow(fun_comments) > 0 && !is.null(fun_comments$texts)) {
        function_comment <- paste(
          unlist(lapply(fun_comments$texts, function(x) x$value)),
          collapse = " "
        )
      }
    }

    # SUBCELLULAR LOCATION
    sub_location <- NA
    if ("comments" %in% names(entry)) {
      sub_comm <- entry$comments[entry$comments$commentType == "SUBCELLULAR LOCATION", ]
      if (nrow(sub_comm) > 0 && !is.null(sub_comm$subcellularLocations)) {
        locs <- unlist(lapply(sub_comm$subcellularLocations, function(x) {
          sapply(x, function(z) z$location$value)
        }))
        sub_location <- paste(unique(locs), collapse = ", ")
      }
    }

    # DOMAINS
    domains <- NA
    if ("features" %in% names(entry)) {
      doms <- entry$features[entry$features$type == "Domain", , drop = FALSE]
      if (nrow(doms) > 0 && "description" %in% names(doms)) {
        domains <- paste(unique(na.omit(doms$description)), collapse = ", ")
      }
    }

    result_list[[id]] <- data.frame(
      uniprot_id = id,
      gene_names = gene_names,
      function_comment = function_comment,
      subcellular_location = sub_location,
      domains = domains,
      stringsAsFactors = FALSE
    )

    Sys.sleep(sleep)
  }

  return(do.call(rbind, result_list))
}



# pobranie pełnych danych, później będę się zastanawiać nad tym, czy nie można tego zrobić w jednym zapytaniu

df <- get_uniprot_gene_info(limit = 50)

df$results %>% 
  head %>%
  colnames()



df$results %>% 
  head %>%
  select(
    "organism"
  ) %>% 
  as.tibble() %>% 
  unnest() %>% 
  unnest() %>% 
  unnest() %>% 
  unnest()


df$results %>% 
  head %>% 



df$results



ids <- c("P04637", "Q96GD4", "Q96RI1")
details_df <- get_uniprot_detailed_info(ids)
print(details_df)

details_df$domains
