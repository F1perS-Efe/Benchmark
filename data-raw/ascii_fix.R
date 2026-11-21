## Bir kerelik: R kaynaklarindaki Turkce karakterleri ASCII'ye cevir
map <- c(
  "\u0131" = "i", "\u0130" = "I", "\u015f" = "s", "\u015e" = "S",
  "\u011f" = "g", "\u011e" = "G", "\u00fc" = "u", "\u00dc" = "U",
  "\u00f6" = "o", "\u00d6" = "O", "\u00e7" = "c", "\u00c7" = "C",
  "\u2013" = "-", "\u2014" = "-", "\u2019" = "'",
  "\u201c" = "\"", "\u201d" = "\""
)
files <- c(
  list.files("R", full.names = TRUE, pattern = "[.]R$"),
  list.files("inst/shiny", full.names = TRUE, pattern = "[.]R$", recursive = TRUE)
)
for (f in files) {
  txt <- readLines(f, warn = FALSE, encoding = "UTF-8")
  for (nm in names(map)) {
    txt <- gsub(nm, map[[nm]], txt, fixed = TRUE)
  }
  writeLines(txt, f)
  rem <- tools::showNonASCIIfile(f)
  cat(f, ": remaining=", length(rem), "\n")
}
