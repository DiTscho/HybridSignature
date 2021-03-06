#### Plots for 100 random signatures################################################################

load(file = "randomSignatures")

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
















