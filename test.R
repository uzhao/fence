require(lme4)
data(iris)
samplem = lmer
#samplem = lm
samplef = Sepal.Length~Sepal.Width+Petal.Length+Petal.Width+(1|Species)+(1|Species2)
#samplef = Sepal.Length~Sepal.Width+Petal.Length+Petal.Width
sampled = iris
sampled$Species2 = sample(sampled$Species, length(sampled$Species))
full = samplem(samplef, sampled)
samplel1 = function(model) {
  sum(residuals(model)^2)
}
samplel2 = AIC

junk = adaptivefence(samplem, samplef,
                     list(Sepal.Length ~ Sepal.Width+(1|Species),
                          Sepal.Length ~ Petal.Width+(1|Species)),
                     sampled, samplel1, list(mode = "p", B = 10), sizef = size_lmerMod)
