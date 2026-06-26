#Simple Model ---------------------------------
# ── Cook's distance and studentised residuals ──────────────────────────────
cooks    <- cooks.distance(res)
rstudent <- rstudent(res)

# ── Plot to visualise ──────────────────────────────────────────────────────
par(mfrow = c(1, 2))

# Cook's distance - spikes above the threshold are influential
plot(cooks, type = "o", pch = 19, 
     xlab = "Observation", ylab = "Cook's Distance",
     main = "Cook's Distance")
abline(h = mean(cooks) + 3*sd(cooks), lty = 2, col = "red")  # 3SD threshold

# Studentised residuals - values beyond ±1.96 are potential outliers
plot(rstudent$z, type = "o", pch = 19,
     xlab = "Observation", ylab = "Studentised Residual",
     main = "Studentised Residuals")
abline(h =  1.96, lty = 2, col = "red")
abline(h = -1.96, lty = 2, col = "red")

# ── Identify which observations are flagged ────────────────────────────────
# Cook's distance outliers
cook_outliers <- which(cooks > mean(cooks) + 3*sd(cooks))
cat("Cook's distance outliers (row numbers):\n")
print(escalc_data[cook_outliers, c("Species", "Paper_ID", "yi", "vi")])

# Studentised residual outliers  
resid_outliers <- which(abs(rstudent$z) > 1.96)
cat("Studentised residual outliers (row numbers):\n")
print(escalc_data[resid_outliers, c("Species", "Paper_ID", "yi", "vi")])

# ── Combined: flagged by both ──────────────────────────────────────────────
both_outliers <- intersect(cook_outliers, resid_outliers)
cat("Flagged by both methods:\n")
print(escalc_data[both_outliers, c("Species", "Paper_ID", "yi", "vi")])


#Phylo model ---------------------------------
# ── Cook's distance and studentised residuals ──────────────────────────────
cooks    <- cooks.distance(res_phylo)
rstudent <- rstudent(res_phylo)

# ── Plot to visualise ──────────────────────────────────────────────────────
par(mfrow = c(1, 2))

# Cook's distance - spikes above the threshold are influential
plot(cooks, type = "o", pch = 19, 
     xlab = "Observation", ylab = "Cook's Distance",
     main = "Cook's Distance")
abline(h = mean(cooks) + 3*sd(cooks), lty = 2, col = "red")  # 3SD threshold

# Studentised residuals - values beyond ±1.96 are potential outliers
plot(rstudent$z, type = "o", pch = 19,
     xlab = "Observation", ylab = "Studentised Residual",
     main = "Studentised Residuals")
abline(h =  1.96, lty = 2, col = "red")
abline(h = -1.96, lty = 2, col = "red")

# ── Identify which observations are flagged ────────────────────────────────
# Cook's distance outliers
cook_outliers <- which(cooks > mean(cooks) + 3*sd(cooks))
cat("Cook's distance outliers (row numbers):\n")
print(escalc_data[cook_outliers, c("Species", "Paper_ID", "yi", "vi")])

# Studentised residual outliers  
resid_outliers <- which(abs(rstudent$z) > 1.96)
cat("Studentised residual outliers (row numbers):\n")
print(escalc_data[resid_outliers, c("Species", "Paper_ID", "yi", "vi")])

# ── Combined: flagged by both ──────────────────────────────────────────────
both_outliers <- intersect(cook_outliers, resid_outliers)
cat("Flagged by both methods:\n")
print(escalc_data[both_outliers, c("Species", "Paper_ID", "yi", "vi")])


#Class model ---------------------------------
# ── Cook's distance and studentised residuals ──────────────────────────────
cooks    <- cooks.distance(res_Class)
rstudent <- rstudent(res_Class)

# ── Plot to visualise ──────────────────────────────────────────────────────
par(mfrow = c(1, 2))

# Cook's distance - spikes above the threshold are influential
plot(cooks, type = "o", pch = 19, 
     xlab = "Observation", ylab = "Cook's Distance",
     main = "Cook's Distance")
abline(h = mean(cooks) + 3*sd(cooks), lty = 2, col = "red")  # 3SD threshold

# Studentised residuals - values beyond ±1.96 are potential outliers
plot(rstudent$z, type = "o", pch = 19,
     xlab = "Observation", ylab = "Studentised Residual",
     main = "Studentised Residuals")
abline(h =  1.96, lty = 2, col = "red")
abline(h = -1.96, lty = 2, col = "red")

# ── Identify which observations are flagged ────────────────────────────────
# Cook's distance outliers
cook_outliers <- which(cooks > mean(cooks) + 3*sd(cooks))
cat("Cook's distance outliers (row numbers):\n")
print(escalc_data[cook_outliers, c("Species", "Paper_ID", "yi", "vi")])

# Studentised residual outliers  
resid_outliers <- which(abs(rstudent$z) > 1.96)
cat("Studentised residual outliers (row numbers):\n")
print(escalc_data[resid_outliers, c("Species", "Paper_ID", "yi", "vi")])

# ── Combined: flagged by both ──────────────────────────────────────────────
both_outliers <- intersect(cook_outliers, resid_outliers)
cat("Flagged by both methods:\n")
print(escalc_data[both_outliers, c("Species", "Paper_ID", "yi", "vi")])
