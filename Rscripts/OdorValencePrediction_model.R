library(caret)
#library(boot)
#library(MLmetrics)
#library(ROCR)

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

### Linear Regression Model
training_data_4LVs<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/PLSDA_components/PLSDA_4LVs.csv")
training_data_6LVs<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/PLSDA_components/PLSDA_6LVs.csv")
training_data_8LVs<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/Plots/LinearRegression/training_8LVs.csv")
training_data_10LVs<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/PLSDA_components/PLSDA_10LVs.csv")
training_data_12LVs<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/PLSDA_components/PLSDA_12LVs.csv")
training_data_all<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/SSR_dataset.csv")

target_data<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/Plots/LinearRegression/target.csv")
weights<-read.csv("/Users/Tayfun/Desktop/Codes/OdorValPrediction/Plots/LinearRegression/weights.csv")
names(weights)[1]<-'weights'

regressionDF<-training_data_8LVs

regressionDF_scaled<-as.data.frame(scale(training_data_8LVs, center = TRUE, scale = TRUE))
regressionDF_scaled$valence<-target_data$Valence

head(regressionDF)

## Outlier detection
boxplot(regressionDF)$out

regressionDF_cleaned <- remove_outliers(regressionDF)

regressionDF$valence<-target_data$Valence

## Models
ols_model<- lm(valence ~ ., data=regressionDF_scaled,
               x=TRUE, y = TRUE)

ols_model_w<- lm(valence ~ ., data=regressionDF_scaled,
               x=TRUE, y = TRUE, weights = weights$weights)

summary(ols_model)
summary(ols_model_w)

layout(matrix(c(1,2,3,4),2,2))
plot(ols_model_w)

#### SVM MODELS

## tuning hyperparameters
svmGrid <-  expand.grid(C = 10^(-3:2),
                        degree = seq(1, 5, 1),
                        scale = seq(0.1, 1, 0.1))
                        #sigma = seq(0.1, 1, 0.1))

set.seed(29510)
algo <- "lm"

data_ctrl <- trainControl(method = "cv", p = 0.8, number = 10)

model_caret <- train(form = valence ~ .,   # model to fit
                     data = regressionDF, 
                     trControl = data_ctrl,              # folds
                     method = algo,                      # specifying regression model
                     #tuneGrid = svmGrid,
                     na.action = na.pass)

model_caret_w <- train(form = valence ~ .,   # model to fit
                     data = regressionDF,   
                     trControl = data_ctrl,              # folds
                     weights = weights$weights,
                     method = algo,                      # specifying regression model
                     #tuneGrid = svmGrid,
                     na.action = na.pass)

resamp <- resamples(list(model_caret, model_caret_w))
summary(resamp)

plot(model_caret_w)
model_caret_w$finalModel

capture.output(model_caret, file = "/Users/Tayfun/Desktop/Codes/OdorValPrediction/Plots/LinearRegression/MLRegression_Weighted_LOOCV.txt")

### TRY IN PYTHONNN WEIGHTINGG
varImp(model_caret)
varImp(model_caret_w)

### Learning curve
set.seed(29510)

lda_data <- learning_curve_dat(dat = regressionDF, 
                              outcome = "valence",
                              test_prop = 1/5,
                              ## `train` arguments:
                              method = "svmRadial", 
                              metric = "RMSE",
                              trControl = data_ctrl)

library(dplyr)
target <- c('Testing', 'Training')
lda_data <- filter(lda_data, Data %in% target)

ggplot(lda_data, aes(x = Training_Size, y = RMSE, color = Data)) + 
  geom_line() +
  #geom_smooth(method = loess, span = .2) + 
  ylim(0, 0.8) +
  xlab ('Training size') +
  theme_classic()



