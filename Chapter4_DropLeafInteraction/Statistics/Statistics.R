library(lme4)
library(lmerTest)
library(DHARMa)
library(multcomp)

# Influence on drop impact location on bending and twisting angles
angle_data <- read.csv("AngleBending.csv")

# Perform a repeated measures mixed effects model
fittedModel <- lmer(log(Angle) ~ log(CorMLf)*Location + (1|Species), data = angle_data, REML=FALSE)
summary(fittedModel)
anova(fittedModel)

# Compare to model without interaction
fittedModel2 <- lmer(log(Angle) ~ log(CorMLf) + Location + (1|Species), data = angle_data, REML=FALSE)
anova(fittedModel, fittedModel2)

# Testing goodness of fit
simulationOutput <- simulateResiduals(fittedModel = fittedModel, plot = F)
plot(simulationOutput)

# Bending Angle n = 38
# Model: log(Angle) ~ log(CorMLf)*Location + (1|Species)
# Significant difference between impact locations
# using log transformed angles: F (2,76) = 5.5494, p-val = 0.005629

# Reporting anova of comparing two models 
# Model full: log(Angle) ~ log(CorMLf)*Location + (1|Species)
# Model null: log(Angle) ~ log(CorMLf) + Location + (1|Species)
# The interaction of Location and log(1/ML2f)(χ2(2)=6.3809, p=0.04115) affected the log(angle)
# For Lateral location and log(1/ML2f) log(angle) was lowered by -0.02569 degree +- 0.04386 degree
# For Tip location and log(1/ML2f) log(angle) was lowered by -0.10830 degree +- 0.04386 degree

angle_data <- read.csv("AngleTwisting.csv")

# Perform a repeated measures mixed effects model
fittedModel <- lmer(log(Angle) ~ log(CorIpf)*Location + (1|Species), data = angle_data, REML=FALSE)
summary(fittedModel)
anova(fittedModel)

# Compare to model without interaction
fittedModel2 <- lmer(log(Angle) ~ log(CorIpf) + Location + (1|Species), data = angle_data, REML=FALSE)
anova(fittedModel, fittedModel2)

# Testing goodness of fit
simulationOutput <- simulateResiduals(fittedModel = fittedModel, plot = F)
plot(simulationOutput)

# Twisting Angle n = 31
# Model: log(Angle) ~ log(CorIpf)*Location + (1|Species)
# Significant difference between impact locations
# using log transformed angles: F (2,62) = 8.1858, p-val = 0.0007002

# Reporting anova of comparing two models 
# Model full: log(Angle) ~ log(CorIpf)*Location + (1|Species)
# Model null: log(Angle) ~ log(CorIpf) + Location + (1|Species)
# The interaction of Location and log(1/Ipf)(χ2(2)=15.789, p=0.0003728) affected the log(angle)
# For Lateral location and log(1/Ipf) log(angle) was lowered by -0.21979 degree +- 0.11242 degree
# For Tip location and log(1/ML2f) log(angle) was lowered by -0.47626 degree +- 0.11242 degree


#-----------------------------------------------------------------------------------------------------------------
# Influence on drop impact location on damping coeffcients

# Switch between two datasets
damping_data <- read.csv("DampingBending.csv")
#damping_data <- read.csv("DampingTwisting.csv")

# Perform a repeated measures mixed effects model
fittedModel <- lmer(log(Damping) ~ Location + (1|Species), data = damping_data)
summary(fittedModel)
anova(fittedModel)

# Testing goodness of fit
simulationOutput <- simulateResiduals(fittedModel = fittedModel, plot = F)
plot(simulationOutput)

# Posthoc tukey test
summary(glht(fittedModel, linfct = mcp(Location = "Tukey")), test = adjusted("holm"))

# Notes on results:

# Damping coefficient in bending n=38
# Model: log(Damping) ~ Location + (1|Species)
# Significant difference between impact locations 
# using log transformed damping coeffs: F (2,74) = 10.215, p-val = 0.00012
# Log transform necessary for meeting assumptions of lmer to have normally distributed residuals, ect
# Posthoc tukey shows that 
# Lateral and Central do not differ  z = 1.380, p-val = 0.168
# Tip and Central differ z = 4.418, p-val < 0.001
# Tip and Lateral differ z = 3.038, p-val = 0.005

# Damping coefficient in twisting n=25
# Model: log(Damping) ~ Location + (1|Species)
# No difference between impact locations 
# using log transformed damping coeffs: F (2,48) = 0.9681, p-val = 0.3871
# Log transform necessary for meeting assumptions of lmer to have normally distributed residuals, ect
# Posthoc tukey shows that 
# Lateral and Central do not differ  z = -0.605, p-val = 0.867
# Tip and Central do not differ z = 0.783, p-val = 0.867
# Tip and Lateral do not differ z = 1.388, p-val = 0.496