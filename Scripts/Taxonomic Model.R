library(dplyr)

# Count observations per Class
class_counts <- table(test$Class)

# Keep only classes with at least 3 observations
valid_classes <- names(class_counts[class_counts >= 3])

# Filter the dataset
test_Class <- test %>%
  filter(Class %in% valid_classes)

escalc_data_Class <- escalc(
  measure = "SMD",
  m1i = Treat_MHS,
  sd1i = Treat_SD,
  n1i = Treat_Reps,
  m2i = Control_MHS,
  sd2i = Control_SD,
  n2i = Control_Reps,
  data = test_Class
)

# Convert Class to factor
escalc_data_Class$Class <- as.factor(escalc_data_Class$Class)

# Fit random-effects model with Class as moderator
res_Class <- rma(yi, vi, mods = ~ Class - 1, data = escalc_data_Class, method = "REML")

library(dplyr)

res_by_Class <- escalc_data_Class %>%
  group_by(Class) %>%
  group_map(~ rma(yi, vi, data = .x, method = "REML"), .keep = TRUE)
escalc_data_Class <- escalc_data_Class %>%
  arrange(Class)

k <- length(escalc_data_Class$yi)
escalc_data_Class$Class <- as.character(escalc_data_Class$Class)

forest(escalc_data_Class$yi,
       vi = escalc_data_Class$vi,
       slab = escalc_data_Class$Paper_ID,
       ilab = escalc_data_Class$Class,
       ilab.xpos = -19,
       xlim = c(-25, 10),
       ylim = c(-6.25, k + 2),
       xlab = "Hedges' g",
       header = "Paper_ID",
       psize = 1.5,
       cex = 0.8)

text(-19, k + 1, "Class", font = 2)

# Find y-positions for each group
rows_per_class <- table(escalc_data_Class$Class)
end_rows <- cumsum(rows_per_class)

# Add subgroup summary diamonds at the right positions
addpoly(res_by_Class[[1]], row = -0.5, mlab = "Actinopterygii pooled effect")
addpoly(res_by_Class[[2]], row = -2, mlab = "Copepoda pooled effect")
addpoly(res_by_Class[[3]], row = -3.5, mlab = "Gastropoda pooled effect")
res_overall <- rma(yi, vi, data = escalc_data_Class, method = "REML")
addpoly(res_overall, row = -5, mlab = "Overall effect", cex = 0.9)


summary(res_Class)
  

