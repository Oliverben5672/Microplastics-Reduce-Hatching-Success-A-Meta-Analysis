library(dplyr)
library(kableExtra)

model_results <- data.frame(
  Class = c("Simple Model", "Phylogenetic Model", "Class Model", "Actinopterygii", "Copepoda", "Gastropoda"),
  Estimate = c(-1.24, -1.38, -1.23, -0.92, -0.46, -3.25),
  SE = c(0.31, 0.37, NA, 0.39, 0.77, 0.92),
  Significance = c("***", "***", "***", "*", "NS", "***"),
  n = c(29, 29, 25, 18, 4, 3)
)


model_results %>%
  kable(
    caption = "Mixed-Effects Model Results by Class",
    col.names = c("Class", "Estimate", "SE", "Significance", "n"),
    align = "lcccc"
  ) %>%
  kable_styling(full_width = FALSE, position = "center",
                bootstrap_options = c("striped", "hover", "condensed"))
