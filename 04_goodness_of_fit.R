#
# 04_goodness_of_fit.R
#

library(dplyr)

# check whether [2 sigmas; PCA + all others] model matches GLM estimates
df1 <- read.csv("data/individual_dat_2sigmas_pca_and_all_others_v2.csv")
result1 <- readRDS("data/results_2sigmas_pca_and_all_others_v2.RDS")

model1 <- facilitybirth ~ birthorder + age + married + pca + surface_var
fit1 <- glm(model1, family = binomial, data = df1)

coef(fit1)
result1$par


matrix(
  c(round(result1$par, digits = 4), round(sqrt(diag(solve(-1 * result1$hessian))), digits = 4)),
  ncol = 2
)



# check [2 sigmas + PCA] model
df2 <- read.csv("data/individual_dat_2sigmas_pca_v2.csv")
result2 <- readRDS("data/results_2sigmas_pca_v2.RDS")

model2 <- facilitybirth ~  pca + surface_var
fit2 <- glm(model2, family = binomial, data = df2)

coef(fit2)
result2$par


# SES model

fit3 <- glm(facilitybirth ~ pca, family = binomial, data = df2)
summary(fit3)


# ROC plots

library(ROCR)

df2$facilitybirth <- as.factor(df2$facilitybirth)

df_pred1 <- df1[, c("birthorder", "age", "married", "pca", "surface_var")]
df_pred2 <- df2[, c("pca", "surface_var")]
df_pred3 <- data.frame(pca = df2$pca)

pred1 <- predict(fit1, newdata = df_pred1, type = "response")
pred_result1 <- prediction(pred1, df1$facilitybirth)
performance1 <- performance(pred_result1, measure = "tpr", x.measure = "fpr")

pred2 <- predict(fit2, newdata = df_pred2, type = "response")
pred_result2 <- prediction(pred2, df2$facilitybirth)
performance2 <- performance(pred_result2, measure = "tpr", x.measure = "fpr")

pred3 <- predict(fit3, newdata = df_pred3, type = "response")
pred_result3 <- prediction(pred3, df2$facilitybirth)
performance3 <- performance(pred_result3, measure = "tpr", x.measure = "fpr")

plot(performance3, col = "green")
plot(performance2, col = "red", add = TRUE)
plot(performance1, col = "blue", add = TRUE)
# plot(performance1, lwd = 2)
abline(0,1, col = "gray")

auc1 <- performance(pred_result1, measure = "auc")
(auc1 <- auc1@y.values[[1]])

# legend(x = 0.45, y = 0.4,
#   legend = paste0("AUC = ", round(auc1, digits = 4)),
#   bty = "n", cex = 1.3
# )

auc2 <- performance(pred_result2, measure = "auc")
(auc2 <- auc2@y.values[[1]])

auc3 <- performance(pred_result3, measure = "auc")
(auc3 <- auc3@y.values[[1]])

legend(x = 0.3, y = 0.3,
  legend = c(
    paste0("Model 1; AUC = ", round(auc1, digits = 4)),
    paste0("Model 2; AUC = ", round(auc2, digits = 4)),
    paste0("Model 3; AUC = ", round(auc3, digits = 4)) ),
  lty = 1, col = c("blue", "red", "green"), bty = "n", cex = 1
)





# # cross validation
#
# loocv <- function (obj) {
#   data <- obj$data
#   m <- dim(data)[1]
#   form <- formula(obj)
#   fam <- obj$family$family
#   loo <- rep(NA, m)
#
#   for (i in 1:m) {
#     i.glm <- glm(form, data = data[-i, ], family = fam)
#     loo[i] <- predict(i.glm, newdata = data[i,], family = fam, type = "response")
#   }
#
#   loo
# }
#
# yhat.1.cv <- loocv(model1)
# yhat.2.cv <- loocv(model2)
#
#





