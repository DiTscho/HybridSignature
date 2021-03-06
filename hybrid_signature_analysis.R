####################################################################################################
# First author: Bernhard Ulm
# Second author: Dimitrij Tschodu
####################################################################################################
##### Load libraries ###############################################################################
library("tidyr")
library("dplyr")
library("purrr")
library("broom")
library("ipred")
library("pec")
#### Define global constants #######################################################################

stime = 120 # survival time in months

#### Load METABRIC data ############################################################################
# You have to get access in order to download METABRIC (MB) data from the official  website.
# Main results in the paper can be replicated using the training and test sets:
# patient_ids* provides the patients ids we used.

gl = compose(get, load) # syntactic sugar

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

#### Define functions ##############################################################################
#Select variables in training and test sets for each signature
sv = function(df, selection) {return(df %>% select(Stime, Status, selection))}
# Fit Cox model for each signature
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
iauc.compute = function(signature){
  # Computes the integrated AUC for each signature and returns corresponding p-values, iauc1, iauc2
  # as they are used in the "iauc.comp" function from the "survcomp" package.
  auc = tmp$auc
  test_data = tmp$test_data
  
  inds=1:length(tmp$signature)
  i = inds[tmp$signature==signature] # index of signature which the iauc is computed for
  d = size(test_data)[2]-1 #dimension without the first iauc
  auc.1 = unlist(auc[i])
  auc.other = auc[-i] 
  time = sort(test_data[[i]]$Stime[test_data[[i]]$Status==1]) 
  auc.integrate = partial(iauc.comp, auc1 = auc.1, time = time) #syntactic sugar    
  IAUCs = vector("list", d) 
  for (a in seq_along(auc.other)) {IAUCs[[a]] = auc.integrate(unlist(auc.other[[a]]))}
  IAUCs = data.frame(matrix(unlist(IAUCs), nrow=length(IAUCs), byrow=TRUE))
  colnames(IAUCs) = c("pvalue", "iauc1", "iauc2")
  rownames(IAUCs) = tmp$signature[-i]
  return(IAUCs)
}

### Do analysis ####################################################################################
# Load lists of gene selections for each signature
selections = gl(file = "selections")
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
                                  tibble(sv(X.test, selections$randomnpi)))) %>%
              mutate(model = map(training_data, signature_model),
                     cox_prediction = map2(model, test_data, predict),
                     cox_prediction_on_training = map2(model, training_data, predict),
                     auc = map2(test_data, cox_prediction,  auc),
                     dtree_model = map2(training_data, cox_prediction_on_training, dt_train),
                     dtree_preds   = pmap(list(test_data, dtree_model, cox_prediction), dt_predict),
                     Cindex = map(model, c_index),
                     R2 = map(model, r2),
                     LogRank = map(model, logrank),
                     survProbs = map2(model, test_data, surv_probs),
                     brierScore = map2(test_data, survProbs, brier_score),
                     SSS =map(brierScore, sss),
                     IAUCs = map(signature, iauc.compute)
                    )
