require(fence)

## test for lmer
data(iris)
full = Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + (1|Species)
test_af = fence.lmer(full, iris)
# test_naf = fence.lmer(full, iris, fence = "nonadaptive", cn = 12)

## test for NER
## same to lmer
# require(sae)
# data(cornsoybean)
# data(cornsoybeanmeans)
# cornsoybean = cornsoybean[-33,]
# cornsoybean$County = as.factor(cornsoybean$County)
# full = SoyBeansHec ~ CornPix + SoyBeansPix + CornPix ^ 2 + SoyBeansPix ^ 2 + CornPix : SoyBeansPix + (1 | County)
# test_af = fence.lmer(full, cornsoybean)

## test for FH

## test for glmer
## bin
# data = read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
# data$rank = as.factor(data$rank)
# full = admit~gre+gpa+(1|rank)
# test_af = fence.glmer(full, data, family = "b", 100)
# plot(test_af)
# test_af$sel_model

## poi
# data(grouseticks)
# data = grouseticks
# data$YEAR = as.factor(data$YEAR)
# full = TICKS~HEIGHT+YEAR+(1|BROOD)
# test_af = fence.glmer(full, data, family = "p", 50)
# plot(test_af)
# test_af$sel_model

## test for glmmadmb
# data(epil2)
# data = epil2
# rm(epil2)
# data$subject = as.factor(data$subject)
# full = y~Base+Age+Visit+(1|subject)
# 
# data = data.frame(ID = rep(1:50, each = 2), x = rnorm(100, 10))
# data$ID = as.factor(data$ID)
# data$y = rgamma(100, 5, scale = exp(1 + data$x + rep(rnorm(50), each = 2)) / 5)
# full = y~x+x^2+(1|ID)
# # sometimes throw error from glmmadmb
# # test_af = fence.glmmadmb(full, data, family = "g", 50)
# test_naf = fence.glmmadmb(full, data, family = "g", fence = "nonadaptive", cn = 10)

