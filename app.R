# shinyapps.io / rsconnect entry point (not part of the built package).
# Install SigBenchR first, then deploy from the repository root.
# See README.md section "Dagitim (shinyapps.io)".
# Fail-closed: remote deploy never enables local file upload.

SigBenchR::run_app(allow_upload = FALSE)
