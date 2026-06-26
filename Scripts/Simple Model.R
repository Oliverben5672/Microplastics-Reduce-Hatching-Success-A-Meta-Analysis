#test
library(tidyverse)
library(metafor)

test <- MetaAnalysisTest[!is.na(MetaAnalysisTest$Control_MHS), ]

# Calculate Hedges' g (standardized mean difference) and variance
escalc_data <- escalc(
  measure = "SMD",                 # standardized mean difference
  m1i = Treat_MHS,                 # treatment mean
  sd1i = Treat_SD,                 # treatment SD
  n1i = Treat_Reps,                # treatment replicates
  m2i = Control_MHS,               # control mean
  sd2i = Control_SD,               # control SD
  n2i = Control_Reps,              # control replicates
  data = test
)

res <- rma(yi = escalc_data$yi,        # effect sizes
           vi = escalc_data$vi,        # corresponding variances
           method = "REML")            # restricted maximum likelihood estimator


par(mar = c(5, 6, 4, 4))  # reduce left margin slightly

# Reorder data by Class ---
escalc_data <- escalc_data[order(escalc_data$Class), ]

# Refit model using the reordered data 
res <- rma(yi, vi, data = escalc_data, method = "REML")

forest(res,
       xlim = c(-25, 10),               # bring the text area closer to forest
       xlab = "Hedges' g",
       mlab = "Overall effect",
       slab = escalc_data$Paper_ID,    # Paper_ID column
       ilab = escalc_data$Class,       # Class column
       ilab.xpos = -19,                 # move 'Class' closer to forest
       header = "Paper ID",
       psize = 1.5,
       cex = 0.8,
       col = "black") +
  text(-19, length(res$yi) + 2, "Class", font = 2)

summary(res)

#-- plot figure 1 for basic visualization of studies--------------------------------

library(ggplot2)
library(ggalt)  # for geom_dumbbell
library(dplyr)
library(ggthemes)

# Calculate change relative to control
dumbbell_data <- MetaAnalysisTest %>%
  select(Paper_ID, Control_Success, Hatch_Change, Class, Significance) %>%
  # Remove rows with missing values
  filter(!is.na(Hatch_Change), !is.na(Control_Success))

# Create the dumbbell plot
dumbbell_data <- dumbbell_data %>%
  arrange(Class) %>%
  mutate(Paper_ID = factor(Paper_ID, levels = Paper_ID))  # preserve order

# Plot
ggplot(dumbbell_data, aes(y = Paper_ID, colour = Class)) +
  geom_dumbbell(aes(x = 0, xend = Hatch_Change),
                size = 1.5,
                size_x = 2.5,
                size_xend = 2.5,
                dot_guide = FALSE,
                dot_guide_size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", alpha = 0.95) +
  geom_text(aes(x = Hatch_Change / 2,  # midpoint between start and end
                y = Paper_ID,
                label = Paper_ID),
            vjust = 0.45,   
            size = 3,
            color = "black") +
  xlab("Change in Hatching Success (%)") +
  ylab("Significance") +
  scale_y_discrete(labels = dumbbell_data$Significance) +
  scale_colour_manual(values = c(
    "#2A6FBD",  # Blue (3rd block)
    "#7E4FAE",  # Purple (3rd block)
    "#2A9A8F",  # Teal (3rd block)
    "#4C9A2A",  # Green (3rd block)
    "#D6B32A",   # Yellow (3rd block)
    "#E68A2E",  # Orange (3rd block)
    "#C73630"  # Red (3rd block)
  )) +
  coord_cartesian(xlim = c(-80, 10)) +
  scale_x_continuous(breaks = seq(-80, 10, by = 10)) +
  theme_linedraw() +
  ggtitle("Dumbbell Plot of Hatching % Change Relative to Control") +
  theme(plot.title = element_text(hjust = 0.5))




