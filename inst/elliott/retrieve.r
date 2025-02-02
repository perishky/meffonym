
load(url("https://github.com/MRCIEU/godmc/raw/refs/heads/master/resources/smoking/illig.RData"))

model = data.frame(
    pred.var=Illig_data$cpgs,
    weight=Illig_data$weights,
    dir=sign(Illig_data$all_effect),
    ref=Illig_data$reference_never_median_beta_all)

intercept = -sum(model$weight*model$dir*model$ref)

model = rbind(
    data.frame(pred.var="intercept", coef=intercept),
    data.frame(pred.var=model$pred.var, coef=model$dir*model$weight))

write.csv(model, file="coefs.csv",row.names=F)

