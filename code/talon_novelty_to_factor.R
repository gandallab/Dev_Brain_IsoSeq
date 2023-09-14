library(tidyverse)

talon_novelty_to_factor = function(df, split_ISMs = F, split_ISMs_include_both = F, make_other = T) {
  df = df %>%
    mutate(
      ISM_subtype = ISM_subtype %>% na_if("None")
    )
  levels_to_keep = c("Known", "ISM", "NIC", "NNC")
  if (split_ISMs) {
    df = df %>%
      mutate(
        transcript_novelty = if_else(
          is.na(ISM_subtype),
          transcript_novelty,
          str_c(transcript_novelty, ISM_subtype, sep = "_")
        )
      )
    if (split_ISMs_include_both | !make_other) {
      levels_to_keep = c("Known", "ISM_Prefix", "ISM_Suffix", "ISM_Both", "NIC", "NNC")
    } else {
      levels_to_keep = c("Known", "ISM_Prefix", "ISM_Suffix", "NIC", "NNC")
    }
  }
  if (make_other) {
    df = df %>%
      mutate(
        transcript_novelty = transcript_novelty %>%
          fct_other(keep = levels_to_keep) %>%
          fct_relevel(c(levels_to_keep, "Other"))
      )
  } else {
    df = df %>%
      mutate(
        transcript_novelty = transcript_novelty %>%
          fct_infreq() %>%
          fct_relevel(c(levels_to_keep))
      )
  }
  return(df)
}
