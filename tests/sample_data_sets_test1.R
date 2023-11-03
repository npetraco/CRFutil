library(CRFutil)
library(MASS)
library(gRim)

# Brunner lizards
di <- generate_brunner_lizards(typ = 1)
di$lizards.contingency.table
di$lizards.samples
di$lizards.configuration.probabilities
di$lizards.configuration.frequencies

# gRbase lizard
di2 <- generate_gRbase_lizards()
di2$lizards.configuration.frequencies
lizardAGG

di2$lizards.samples
lizardRAW

di2$lizards.contingency.table
lizard

di2$lizards.configuration.probabilities

# Triangle model
tri <- generate_triangle_model(num.samples = 1000, Temp = 5, plot.sampleQ = F)
tri$tri.empirical.joint.distribution
tri$tri.exact.joint.distribution
plot_crf(tri$triangle.model)


# Schmidt small model
sch <- generate_schmidt_small_model(num.samples = 1000, Temp=-10, plot.sampleQ = T)
sch$schmidt.small.model
head(sch$schmidt.small.samples)
sch$schmidt.small.contingency.table
sch$schmidt.small.exact.joint.distribution

# Koller-Fiedman misconception model
msc <- generate_koller_misconception_model(num.samples = 1000, Temp=1, plot.sampleQ = T)
msc$koller.misconception.configuration.frequencies

# Star model
starm <- generate_star_model(num.samples = 1000, plot.sampleQ = T)
#starm$star.contingency.table
dump.crf(starm$star.model)
starm$star.model$par

# Tesseract model
tes <- generate_tesseract_model(num.samples = 1000, Temp=1, plot.sampleQ = T)

# Random graph
randm <- generate_random_model(num.samples = 1000, p.or.numedges = 12, type = "gnm", Temp=1, plot.sampleQ = T)
randm$random.contingency.table
