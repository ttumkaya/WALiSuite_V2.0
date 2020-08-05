library(caret)
library(boot)
library(MLmetrics)
library(ROCR)

df<-read.csv("/Users/Tayfun/Desktop/Codes/Odor_val_prediction_R/Transformed_mean_data_odor_prediction.csv")
df_minimized<-df[,c(3, 6:28)]
head(df_minimized)

table(df_minimized$Valence_class)
## Split train/test datasets
Train <- createDataPartition(df_minimized$Valence_class, p=0.8, list=FALSE)

training <- df_minimized[ Train, ]
testing <- df_minimized[ -Train, ]
dim(training)

## Handling class imbalance by Down- or Up-sampling 
   
  ### for training classes
'%ni%' <- Negate('%in%')

down_training <- downSample(x = training[, colnames(training) %ni% "Valence_class"],
                         y = training$Valence_class,
                         yname = 'Valence_class')

up_training <- upSample(x = training[, colnames(training) %ni% "Valence_class"],
                            y = training$Valence_class,
                            yname = 'Valence_class')

table(down_training$Valence_class)

  ### for full dataset

down_training <- downSample(x = df_minimized[, colnames(df_minimized) %ni% "Valence_class"],
                            y = df_minimized$Valence_class,
                            yname = 'Valence_class')

up_training <- upSample(x = df_minimized[, colnames(df_minimized) %ni% "Valence_class"],
                        y = df_minimized$Valence_class,
                        yname = 'Valence_class')

## See https://topepo.github.io/caret/model-training-and-tuning.html#metrics
f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[1])
  c(F1 = f1_val)
}

# Define train control for k fold cross validation
train_control <- trainControl(method="cv", number=10, summaryFunction = f1, savePredictions = TRUE)

# Fit Logistic Regression
model <- train(Valence_class~., data=training, trControl=train_control,
               preProcess = c("center", "scale"),
               metric="F1",
               tuneLength = 5,
               method="glm", family="binomial")

# Fit Decision Tree
model <- train(Valence_class~., data=training, trControl=train_control,
               #preProcess = c("center", "scale"),
               tuneLength = 5,
               method="ctree")

pred <- predict(model, newdata=testing)
confusionMatrix(mode="prec_recall", data=pred, testing$Valence_class)

# Summarise Results
varImp(model)
print(model)

## Plot ROC curve
prob <- predict(model1, newdata=testing, type="prob")
pred <- prediction(prob, testing$Valence_class)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

getModelInfo(model)
## install.packages('LiblineaR', repo = "https://mac.R-project.org")


### Playing around
featurePlot(x=df_minimized[-1])

### Linear Regression Model
training_data<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/Plots/LinearRegression/training_8LVs.csv")
target_data<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/Plots/LinearRegression/target.csv")
weights<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/Plots/LinearRegression/weights.csv")
names(weights)[1]<-'weights'

regressionDF<-training_data
regressionDF$valence<-target_data$Valence
regressionDF$weight<-weights$weights

ols_model<- lm(valence ~ ., data=regressionDF,
               x=TRUE, y = TRUE)

ols_model_w<- lm(valence ~ ., data=regressionDF,
               x=TRUE, y = TRUE, weights = weights$weights)

summary(ols_model)
summary(ols_model_w)

layout(matrix(c(1,2,3,4),2,2))
plot(ols_model_w)

library(lmvar)
cv.lm(ols_model, k = 10,seed=333)

cv.lm(ols_model_w, k = 10,seed=333)

set.seed(107)
data_ctrl <- trainControl(method = "LOOCV", number = 10)
model_caret <- train(valence ~ .,   # model to fit
                     data = regressionDF,                        
                     trControl = data_ctrl,              # folds
                     weights = weights$weights,
                     method = "lm",                      # specifying regression model
                     na.action = na.pass)

model_caret
capture.output(model_caret, file = "/Users/Tayfun/Desktop/Codes/OdorValPrediction/Plots/LinearRegression/MLRegression_Weighted_LOOCV.txt")

## Manually LOOCV test
#lwo <- rep(0,nrow(regressionDF))

#for (j in 1:nrow(regressionDF))
#{
#  fj <- lm(valence~X0+X1+X2+X3+X4+X5+X6+X7, family=quasibinomial, data=regressionDF[-j,])
#  lwo[j] <- predict(fj, regressionDF[j,], type="response")
#}

#caret::postResample(lwo, regressionDF$valence)

#rss <- sum((lwo - regressionDF$valence) ^ 2)  ## residual sum of squares
#tss <- sum((regressionDF$valence - mean(regressionDF$valence)) ^ 2)  ## total sum of squares
#rsq <- 1 - rss/tss
