####################################################################################################
# First author: Bernhard Ulm
# Second author: Dimitrij Tschodu
# Code status: Still implementing :) 
####################################################################################################




library("tidyr")
library("dplyr")
library("purrr")
library("broom")
library("ipred")

#### Global constants ##############################################################################

stime = 120 # survival time in months



#### Helper functions ##############################################################################
gl = compose(get, load)


#### Load METABRIC data ############################################################################
# You have get access in order to download METABRIC (MB) data from the official  website.
# All results in the paper can be replicated using the training and test sets:
# patient_ids* provide the patients ids we used.
####################################################################################################
X = gl(file= "Metabric.data")
patient_ids.training = gl(file = "patient_ids.training")
patient_ids.test =     gl(file = "patient_ids.test")

# Find patients who died due to disease. Among these patients, count patients with survival > 120 
# months as non-events.   
X.train = X %>%  
  filter(Pat %in% patient_ids.training) %>%
  mutate(Status = replace(!((Vital.Status != "Died of Disease") | (Overall.Survival..Months. > stime)), 0, 1),
         Stime = replace(Overall.Survival..Months., Overall.Survival..Months. > stime, stime)) %>%
  tibble()

X.train$Status == auswahl.pro.jahr$status.pro.jahr

X.test =  X %>% 
  filter(Pat %in% patient_ids.test) %>%
  mutate(Status = replace(!((Vital.Status != "Died of Disease") | (Overall.Survival..Months. > stime)), 0, 1),
         Stime = replace(Overall.Survival..Months., Overall.Survival..Months. > stime, stime)) %>%
  tibble()

randomSignatures = gl(file = "randomSignatures")

#### Get Cox regression on training data ###########################################################
# These are lists of gene selections for each signature
selections = gl(file = "selections")

### Helper functions 

#Select variables in training and test sets for each signature
sv = function(df, selection) {return(df %>% select(Stime, Status, selection))}

# Fits Cox model for each signature
signature_model = function(df){
  lhs = "Surv(Stime, Status)"
  rhs = paste(names(df[3:length(df)]), collapse =  "+")
  form = as.formula(paste(lhs, "~", rhs))
  coxph(form,  data = df, x=T, y=T)
}

auc = function(df, model_prediction){
  #Compute AUC over time
  auc_i = function(time_i){
    #compute AUC at time point time_i
    eventROC(Stime  = df$Stime, 
             status = df$Status,
             marker = model_prediction,
             predict.time = time_i,
             method='iso')$AUC
  }
  map(sort(df$Stime[df$Status==1]), auc_i) #time points up to the stime
}


dt_train = function(df, cox_prediction){
  # Train Decision Trees
  tmp = cbind(df, cox_prediction)
  ctree(Surv(Stime, Status) ~ cox_prediction, data = tmp)
}

dt_predict = function(df, dt_model, cox_prediction){
  # Predict using Decision Trees
  df$cox_prediction = cox_prediction
  predict(dt_model, type ="node", newdata = df)
}

c_index = function(model){return(model$concordance[[6]])}

r2 = function(model){
  # Nagelkerkes version of R2
  m = model
  l0<-        summary(m)$loglik[1] 
  lm =        summary(m)$loglik[2] 
  n =         summary(m)$n
  rsq     <-  (1 - exp(2*(l0-lm)/n)) / (1 - exp(2*(l0)/n)) 
  return(rsq)  
}

logrank = function(model){summary(model)$sctest[[1]]}
surv_probs = function(model, df){as.numeric(predictSurvProb(model, df, times = stime ))}
brier_score = function(df, surv_probs){sbrier(Surv(df$Stime, df$Status), surv_probs, btime=stime)[[1]]}
sss = function(bs){
  #Signature Skill Score
  bs.avg = mean(unlist(randomSignatures$BS)) 
  ((bs.avg - bs)/bs.avg)
}

# Main dataframe with everything inside: data for each signature, Cox models, everything you want.
df = 
  tibble(signature = c("Hybrid", "EndoPredict", "OncotypeDx", "NPI", "Random", "RandomNPI"),
                  training_data = list(    tibble(sv(X.train, selections$hybrid)),
                                  tibble(sv(X.train, selections$endo)),
                                  tibble(sv(X.train, selections$otdx)),
                                  tibble(sv(X.train, "Nottingham.prognostic.index")),
                                  tibble(sv(X.train, selections$random)),
                                  tibble(sv(X.train, selections$randomnpi)) 
                                  ),
                  test_data = list(tibble(sv(X.test, selections$hybrid)),
                                  tibble(sv(X.test, selections$endo)),
                                  tibble(sv(X.test, selections$otdx)),
                                  tibble(sv(X.test, "Nottingham.prognostic.index")),
                                  tibble(sv(X.test, selections$random)),
                                  tibble(sv(X.test, selections$randomnpi))) 
  ) %>%
  mutate(model = map(training_data, signature_model)) %>%
  mutate(cox_prediction = map2(model, test_data, predict)) %>%
  mutate(auc = map2(test_data, cox_prediction,  roc))

tmp = df %>%
  mutate(
    model = map(training_data, signature_model),
    cox_prediction_on_training = map2(model, training_data, predict),
    dtree_model = map2(training_data, cox_prediction_on_training, dt_train),
    dtree_preds   = pmap(list(test_data, dtree_model, cox_prediction), dt_predict),
    Cindex = map(model, c_index),
    R2 = map(model, r2),
    LogRank = map(model, logrank),
    survProbs = map2(model, test_data, surv_probs),
    brierScore = map2(test_data, survProbs, brier_score),
    SSS =map(brierScore, sss)
  )

iauc.comp(tmp$auc[[1]])
unlist(tmp$auc[[1]])





#### Plots ######################################################################################### 


R = data.frame(Cindex = unlist(randomSignatures$Cindex),
               R2 = unlist(randomSignatures$R2),
               LR = unlist(randomSignatures$LR),
               BS = unlist(randomSignatures$BS))
               

P = pal_lancet("lanonc", alpha = 1)(6)
Ptable = pal_lancet("lanonc", alpha = 0.5)(6)
FONT = "Helvitica"
FONTSIZE = 16
ps = 2.5
ls = 1.0

results = t.test(R$Cindex,conf.level = 0.95)
low =  round(results$conf.int[1], 4)
high = round(results$conf.int[2], 4)
mean = round(results$estimate[[1]], 4)

p.cindex = ggplot(R, aes(x=Cindex)) + geom_histogram(binwidth = 0.006, color="black", fill=P[6]) +
  theme_tufte() +
  theme_linedraw() +
  ylab("Counts") +
  xlab("C-index") +
  theme(text = element_text(size = FONTSIZE, family=FONT)) + 
  annotate("text", x = 0.66, y = 11, label = paste("Mean = ", mean, "\n", 
                                                   "95% CI: ", low, "-", high), 
           size=6 )+
  annotate("rect", xmin = low, xmax = high, ymin = 0, ymax = 18,alpha = .2)


results = t.test(R$BS,conf.level = 0.95)
low =  round(results$conf.int[1], 4)
high = round(results$conf.int[2], 4)
mean = round(results$estimate[[1]], 4)

p.bs = ggplot(R, aes(x=BS)) + geom_histogram(binwidth = 0.0025, color="black", fill=P[6]) +
  theme_tufte() +
  theme_linedraw() +
  ylab("Counts") +
  xlab("Brier Score") +
  theme(text = element_text(size = FONTSIZE, family=FONT)) + 
  annotate("text", x = 0.18, y = 11, label = paste("Mean = ", mean, "\n", 
                                                   "95% CI: ", low, "-", high), 
           size=6 )+
  annotate("rect", xmin = low, xmax = high, ymin = 0, ymax = 17,alpha = .2)


results = t.test(R$R2,conf.level = 0.95)
low =  round(results$conf.int[1], 4)
high = round(results$conf.int[2], 4)
mean = round(results$estimate[[1]], 4)

p.r2 = ggplot(R, aes(x=R2)) + geom_histogram(binwidth = 0.005, color="black", fill=P[6]) +
  theme_tufte() +
  theme_linedraw() +
  ylab("Counts") +
  xlab("R2") +
  theme(text = element_text(size = FONTSIZE, family=FONT)) + 
  annotate("text", x = 0.1, y = 11, label = paste("Mean = ", mean, "\n", 
                                                   "95% CI: ", low, "-", high), 
           size=6 )+
  annotate("rect", xmin = low, xmax = high, ymin = 0, ymax = 16,alpha = .2)


results = t.test(R$LR,conf.level = 0.95)
low =  round(results$conf.int[1], 1)
high = round(results$conf.int[2], 1)
mean = round(results$estimate[[1]], 1)

p.lr = ggplot(R, aes(x=LR)) + geom_histogram(binwidth = 4, color="black", fill=P[6]) +
  theme_tufte() +
  theme_linedraw() +
  ylab("Counts") +
  xlab("LogRank statistic") +
  theme(text = element_text(size = FONTSIZE, family=FONT)) + 
  annotate("text", x = 75, y = 11, label = paste("Mean = ", mean, "\n", 
                                                   "95% CI: ", low, "-", high), 
           size=6 )+
  annotate("rect", xmin = low, xmax = high, ymin = 0, ymax = 16,alpha = .2)


grid.arrange(p.cindex, p.r2, p.bs, p.lr, nrow = 2)
















