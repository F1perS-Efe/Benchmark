#!/usr/bin/env Rscript
# SigBenchR ortam kurulumu.
# Depo kokunde:  Rscript setup.R
#
# 1. R >= 4.2 kontrolu
# 2. renv proje kutuphanesini acma
# 3. renv.lock geri yukleme
# 4. DESCRIPTION Suggests (Shiny, rapor, test)
# 5. SigBenchR kurulumunu dogrulama

find_sigbenchr_root <- function() {
  starts <- character()

  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (is.character(ofile) && length(ofile) == 1L && nzchar(ofile)) {
      starts <- c(starts, dirname(ofile))
    }
  }

  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg)) {
    starts <- c(starts, dirname(sub("^--file=", "", file_arg[[1L]])))
  }

  starts <- c(starts, getwd())

  is_root <- function(dir) {
    desc_path <- file.path(dir, "DESCRIPTION")
    if (!file.exists(desc_path)) {
      return(FALSE)
    }
    desc <- tryCatch(read.dcf(desc_path), error = function(e) NULL)
    if (is.null(desc) || !("Package" %in% colnames(desc))) {
      return(FALSE)
    }
    identical(as.character(desc[1, "Package"]), "SigBenchR")
  }

  for (start in unique(starts)) {
    cur <- tryCatch(
      normalizePath(start, winslash = "/", mustWork = TRUE),
      error = function(e) NA_character_
    )
    if (is.na(cur)) {
      next
    }
    for (i in seq_len(8L)) {
      if (is_root(cur)) {
        return(cur)
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) {
        break
      }
      cur <- parent
    }
  }

  stop(
    "SigBenchR kok dizini bulunamadi. Komutu depo kokunde calistirin:\n  Rscript setup.R",
    call. = FALSE
  )
}

parse_dep_field <- function(desc_path, field) {
  desc <- read.dcf(desc_path)
  if (!field %in% colnames(desc)) {
    return(character())
  }
  raw <- as.character(desc[1, field])
  if (is.na(raw) || !nzchar(raw)) {
    return(character())
  }
  parts <- strsplit(raw, ",", fixed = TRUE)[[1L]]
  parts <- trimws(gsub("\\([^)]*\\)", "", parts))
  parts[nzchar(parts)]
}

root <- find_sigbenchr_root()
setwd(root)
message("Proje: ", root)

r_ver <- getRversion()
if (r_ver < "4.2.0") {
  stop(
    sprintf("R >= 4.2.0 gerekli. Kurulu surum: %s", r_ver),
    call. = FALSE
  )
}
if (r_ver < "4.6.0") {
  message(sprintf(
    paste0(
      "Uyari: renv.lock R 4.6.1 ile uretildi (kurulu surum: %s). ",
      "Bazi paketler kaynak koddan derlenebilir. Windows'ta Rtools gerekebilir."
    ),
    r_ver
  ))
}

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  renv.consent = TRUE
)
Sys.setenv(RENV_PROJECT = root)

renv_install <- function(packages, project) {
  args <- list(packages = packages, project = project, prompt = FALSE)
  if ("lock" %in% names(formals(renv::install))) {
    args$lock <- FALSE
  }
  do.call(renv::install, args)
}

activate <- file.path(root, "renv", "activate.R")
if (!file.exists(activate)) {
  stop(
    "renv/activate.R bulunamadi. Depoyu renv klasoruyle birlikte indirin.",
    call. = FALSE
  )
}

message("renv proje kutuphanesi aciliyor...")
source(activate, local = FALSE)

if (!requireNamespace("renv", quietly = TRUE)) {
  stop(
    "renv yuklenemedi. CRAN erisimini ve R kurulumunu kontrol edin.",
    call. = FALSE
  )
}

message("renv.lock paketleri kuruluyor...")
renv::restore(project = root, prompt = FALSE)

suggests <- parse_dep_field(file.path(root, "DESCRIPTION"), "Suggests")
if (length(suggests)) {
  message(
    "Arayuz, rapor ve test paketleri kuruluyor: ",
    paste(suggests, collapse = ", ")
  )
  renv_install(suggests, root)
}

message("SigBenchR kuruluyor...")
renv_install(".", root)
renv::load(project = root, quiet = TRUE)

needed <- unique(c("SigBenchR", suggests))
ok <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) {
  stop(
    sprintf("Su paketler yuklenemedi: %s", paste(needed[!ok], collapse = ", ")),
    call. = FALSE
  )
}

ver <- as.character(utils::packageVersion("SigBenchR"))
message("")
message("Kurulum tamam. SigBenchR ", ver)
message("Bu klasorde bir R oturumu acin:")
message("  library(SigBenchR)")
message("  SigBenchR::run_app()")

if (requireNamespace("rmarkdown", quietly = TRUE) && !rmarkdown::pandoc_available()) {
  message(
    "HTML rapor icin Pandoc yok. RStudio ile gelen Pandoc veya ayri bir Pandoc kurulumu gerekir."
  )
}
