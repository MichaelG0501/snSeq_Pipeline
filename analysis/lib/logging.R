####################
# Analysis registry:
#   Status: active shared library; not a standalone workflow script
#   Script: analysis/lib/logging.R
#   Methodology: analysis/methodology/lib/shared_analysis_library_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Input: run metadata supplied by calling scripts
#   Output: run-summary log helper functions in the calling R session
####################

####################
# Shared run-summary logging helpers.
#
# Long-running scripts should write one log file per run in their logs/ tier.
####################

start_run_summary <- function(script, inputs = character(), outputs = character(),
                              parameters = list(), cached_reused = FALSE) {
  list(
    script = script,
    start_time = Sys.time(),
    end_time = NA,
    inputs = inputs,
    outputs = outputs,
    parameters = parameters,
    cached_reused = cached_reused,
    status = "started"
  )
}

finish_run_summary <- function(run_summary, status = "ok", error = NULL) {
  run_summary$end_time <- Sys.time()
  run_summary$status <- status
  run_summary$error <- error
  run_summary$duration_seconds <- as.numeric(
    difftime(run_summary$end_time, run_summary$start_time, units = "secs")
  )
  run_summary
}

write_run_summary <- function(run_summary, path, include_session = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- c(
    paste0("script: ", run_summary$script),
    paste0("status: ", run_summary$status),
    paste0("start_time: ", format(run_summary$start_time, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("end_time: ", format(run_summary$end_time, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("duration_seconds: ", round(run_summary$duration_seconds, 3)),
    paste0("cached_reused: ", isTRUE(run_summary$cached_reused)),
    "",
    "inputs:",
    paste0("  - ", run_summary$inputs),
    "",
    "outputs:",
    paste0("  - ", run_summary$outputs),
    "",
    "parameters:",
    capture.output(str(run_summary$parameters, give.attr = FALSE)),
    ""
  )
  if (!is.null(run_summary$error)) {
    lines <- c(lines, "error:", paste0("  ", run_summary$error), "")
  }
  if (include_session) {
    lines <- c(lines, "session_info:", paste0("  ", capture.output(sessionInfo())))
  }
  writeLines(lines, con = path)
  invisible(path)
}
