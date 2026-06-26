# Basic funnel plot
funnel(res, 
       xlab = "Hedges' g", 
       ylab = "Standard Error",
       main = "Funnel Plot of Effect Sizes (Hedges' g)")

regtest(res, model = "rma")  # Egger’s regression test for funnel plot asymmetry


tf <- trimfill(res)
funnel(tf, 
       main = "Trim and Fill Funnel Plot",
       xlab = "Hedges' g")
summary(tf)

funnel(tf,
       xlab = "Hedges' g (Standardized Mean Difference)",
       ylab = "Standard Error",
       main = "Trim and Fill Funnel Plot (Publication Bias Assessment)",
       refline = tf$b,              # pooled effect line
       shade = c("gray95", "gray85", "gray75"),
       level = 95,                  # pseudo 95% confidence funnel
       back = "white")


#-class bias---------------------------------------------------------------
# Basic funnel plot
funnel(res_Class, 
       xlab = "Hedges' g (Standardized Mean Difference)",
       ylab = "Standard Error",
       main = "Trim and Fill Funnel Plot (Publication Bias Assessment)",
       refline = tf$b,              # pooled effect line
       shade = c("gray95", "gray85", "gray75"),
       level = 95,                  # pseudo 95% confidence funnel
       back = "white")

regtest(res_Class, model = "rma")  # Egger’s regression test for funnel plot asymmetry

