#rm(list = ls())

source("TMA4900_Data.R")

library(INLA)
library(ggplot2)


# Model 3: Run
################################################################################
inla.model3.replicate.pred = function(complete_data){
  
  # number of datasets
  N = 16
  
  # length of each dataset
  n = nrow(complete_data[[1]])
  m = n-1
  
  # define response matrix
  Y1 = Y2 = Y3 = Y4 = Y5 = wzb1_ = wzb2_ = wzb3_ = c()
  for(i in c(1,2,3,4,6,8,9,10,12,14,15,17)){
    Y1 = c(Y1, rep(10, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y2 = c(Y2, rep(NA, n), complete_data[[i]]$h, rep(NA, n), rep(NA, m), rep(NA, m))
    Y3 = c(Y3, rep(NA, n), rep(NA, n), complete_data[[i]]$g, rep(NA, m), rep(NA, m))
    Y4 = c(Y4, rep(NA, n), rep(NA, n), rep(NA, n), rep(0, m), rep(NA, m))
    Y5 = c(Y5, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(0, m))
    wzb1_ = c(wzb1_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$T[1:m])
    wzb2_ = c(wzb2_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$M[1:m])
    wzb3_ = c(wzb3_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$n[1:m]/4900)
  }
  
  for(i in c(5,7,11,16)){
    Y1 = c(Y1, rep(10, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y2 = c(Y2, rep(NA, n), complete_data[[i]]$h, rep(NA, n), rep(NA, m), rep(NA, m))
    Y3 = c(Y3, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y4 = c(Y4, rep(NA, n), rep(NA, n), rep(NA, n), rep(0, m), rep(NA, m))
    Y5 = c(Y5, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(0, m))
    wzb1_ = c(wzb1_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$T[1:m])
    wzb2_ = c(wzb2_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$M[1:m])
    wzb3_ = c(wzb3_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -(1/4900)*complete_data[[i]]$n[1:m])
  }
  
  Y = cbind(Y1, Y2, Y3, Y4, Y5)
  
  # replicate object
  rep = c()
  for(i in 1:N){
    rep = c(rep, rep(i, 3*n + 2*m))
  }
  
  # digression: 1 kg O2 gives 1.067 kg CO2. 1 kg CO2 requires 0.937 kg O2
  # scaling constant 0.937 for CO2
  
  cscale = 0.937
  
  # define indices and weights of latent fields
  
  # latent z(t)
  idz = rep(c(1:n, rep(NA, n), rep(NA, n), 2:n, rep(NA, m)), N)
  wz = rep(c(rep(1, n), rep(NA, n), rep(NA, n), rep(1, m), rep(NA, m)), N)
  
  # latent h(t)
  idh = rep(c(rep(NA, n), 1:n, rep(NA, n), 2:n, rep(NA, m)), N)
  wh = rep(c(rep(NA, n), rep(1, n), rep(NA, n), rep(-1, m), rep(NA, m)), N)
  
  # latent g(t)
  idg = rep(c(rep(NA, n), rep(NA, n), 1:n, 2:n, 2:n), N)
  wg = rep(c(rep(NA, n), rep(NA, n), rep(1, n), rep(cscale, m), rep(cscale, m)), N)
  
  # latent z(t-1)
  idz_ = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), 1:m, rep(NA, m)), N)
  wz_ = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(-1, m), rep(NA, m)), N)
  
  # latent b0*z(t-1)
  idzb0_ = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), 1:m), N)
  wzb0_ = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(-1, m)), N)
  
  # latent b1*T*z(t-1) 
  idzb1_ = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), 1:m), N)
  
  # latent b2*M*z(t-1) 
  idzb2_ = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), 1:m), N)
  
  # latent b3*Dens*z(t-1) 
  idzb3_ = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), 1:m), N)
  
  
  
  
  # define residuals
  ide0 = rep(c(1:n, rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m)), N)
  ideh = rep(c(rep(NA, n), 1:n, rep(NA, n), rep(NA, m), rep(NA, m)), N)
  ideg = rep(c(rep(NA, n), rep(NA, n), 1:n, rep(NA, m), rep(NA, m)), N)
  idem = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), 1:m, rep(NA, m)), N)
  idea = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), 1:m), N)
  
  nrow(Y)
  length(rep)
  # define inla data
  inla.data = data.frame("Y" = Y, "idh" = idh, "wh" = wh, "idg" = idg, "wg" = wg, 
                         "idz" = idz, "wz" = wz, "idz_" = idz_, "wz_" = wz_, 
                         "idzb0_" = idzb0_, "wzb0_", wzb0_,
                         "idzb1_" = idzb1_, "wzb1_", wzb1_,
                         "idzb2_" = idzb2_, "wzb2_", wzb2_,
                         "idzb3_" = idzb3_, "wzb3_", wzb3_,
                         "ide0" = ide0, "ideh" = ideh, "ideg" = ideg, 
                         "idem" = idem, "idea" = idea, "replicate" = rep)
  
  # priors for likelihood variances
  pc.prec.0 = list(prior = "pc.prec", param = c(5, 0.1), initial = 0.5, fixed = TRUE)
  pc.prec.h = list(prior = "pc.prec", param = c(0.01, 0.5))
  pc.prec.g = list(prior = "pc.prec", param = c(0.01, 0.5))
  pc.prec.m = list(prior = "pc.prec", param = c(1e-2, 1e-3))
  pc.prec.a = list(prior = "pc.prec", param = c(0.01, 0.5))
  
  # prior for random walk hyperparameter in z
  pc.prec.z = list(prior = "pc.prec", param = c(0.01, 0.001))
  
  # define formula object
  formula = Y ~ -1 + 
    f(idz, wz, model = "rw1", constr = F, hyper = list(prec = pc.prec.z), replicate = rep) +
    f(idh, wh, model = "rw1", constr = F, replicate = rep) +
    f(idg, wg, model = "rw1", constr = F, replicate = rep) +
    f(idz_, wz_, copy = "idz", hyper = list(beta = list(fixed = TRUE))) +
    f(idzb0_, wzb0_, copy = "idz", hyper = list(beta = list(fixed = FALSE, prior = "gaussian", param = c(0, 1e2)))) +
    #f(idzb1_, wzb1_, copy = "idz", hyper = list(beta = list(fixed = FALSE, prior = "gaussian", param = c(0, 1e4)))) +
    f(idzb2_, wzb2_, copy = "idz", hyper = list(beta = list(fixed = FALSE, prior = "gaussian", param = c(0, 1e2))))
  #f(idzb3_, wzb3_, copy = "idz", hyper = list(beta = list(fixed = FALSE, prior = "gaussian", param = c(0, 1e6)))) #+
  #f(ide0, model = "iid", hyper = list(prec = pc.prec.0), replicate = rep) +
  #f(ideh, model = "iid", hyper = list(prec = pc.prec.h), replicate = rep) +
  #f(ideg, model = "iid", hyper = list(prec = pc.prec.g), replicate = rep) +
  #f(idem, model = "iid", hyper = list(prec = pc.prec.m), replicate = rep) +
  #f(idea, model = "iid", hyper = list(prec = pc.prec.a), replicate = rep) 
  
  set.seed(1)
  r = inla(formula, data = inla.data, family = rep("gaussian", 5),
           control.family = list(list(hyper = list(prec = pc.prec.0)),
                                 list(hyper = list(prec = pc.prec.h)),
                                 list(hyper = list(prec = pc.prec.g)),
                                 list(hyper = list(prec = pc.prec.m)),
                                 list(hyper = list(prec = pc.prec.a))),
           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
  
  return(r)
  
}
r3.pred = inla.model3.replicate.pred(complete_data)
plot(r3.pred$summary.random$idg$mean, type = "l")
r3.pred$summary.hyperpar

summary(r3.pred)

-sum(r3.pred$cpo$cpo[!is.na(r3.pred$cpo$cpo)])

(MSE_m3t1 = mean((complete_data[[5]]$g - r3.pred$summary.random$idg$mean[1:288 + 12*288])^2))
(MSE_m3t2 = mean((complete_data[[7]]$g - r3.pred$summary.random$idg$mean[1:288 + 13*288])^2))
(MSE_m3t3 = mean((complete_data[[11]]$g - r3.pred$summary.random$idg$mean[1:288 + 14*288])^2))
(MSE_m3t4 = mean((complete_data[[16]]$g - r3.pred$summary.random$idg$mean[1:288 + 15*288])^2))
mean(c(MSE_m3t1, MSE_m3t2, MSE_m3t3, MSE_m3t4))

################################################################################

# Model 3: Posterior distribution of a
################################################################################
post.b0.ggdf = data.frame("x" = r3.pred$marginals.hyperpar$`Beta for idzb0_`[,1],
                         "y" = r3.pred$marginals.hyperpar$`Beta for idzb0_`[,2])

(post.b0.plot = ggplot(aes(x = x, y = y), data = post.b0.ggdf) + 
    geom_line(lwd = 1.5, col = "brown4") +
    geom_area(fill = "brown4", alpha = 0.5) +
    scale_y_continuous("Density") + 
    scale_x_continuous(bquote(beta[0]),
                       breaks = c(0.03, 0.04)) + theme_bw() +
    theme(axis.title = element_text(size = 35),
          axis.text = element_text(size = 20)))

post.b2.ggdf = data.frame("x" = r3.pred$marginals.hyperpar$`Beta for idzb2_`[,1],
                          "y" = r3.pred$marginals.hyperpar$`Beta for idzb2_`[,2])

(post.b2.plot = ggplot(aes(x = x, y = y), data = post.b2.ggdf) + 
    geom_line(lwd = 1.5, col = "brown4") +
    geom_area(fill = "brown4", alpha = 0.5) +
    scale_y_continuous("Density") + 
    scale_x_continuous(bquote(beta[2]),
                       breaks = c(0.0030, 0.0035)) + theme_bw() +
    theme(axis.title = element_text(size = 35),
          axis.text = element_text(size = 20)))

ggarrange(post.b0.plot, post.b2.plot + rremove("ylab"), ncol = 2, align = "h")

################################################################################

# Model 3: Plot posterior standard deviations
################################################################################

marg.variance.h.r3 <- inla.tmarginal(function(x) 1/x,
                                     r3.pred$marginals.hyperpar$`Precision for the Gaussian observations[2]`)
sigma2.h.mod3.df <- data.frame(marg.variance.h.r3)

(post.sigma.h.mod3 = ggplot(sigma2.h.mod3.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "brown4") +
    geom_area(alpha = .5, fill = "brown4") +
    scale_x_continuous(bquote(sigma[h]^2),
                       breaks = c(0.0187,0.0193)) + 
    scale_y_continuous("Density", limits = c(0, 4100)) + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))


marg.variance.g.r3 <- inla.tmarginal(function(x) 1/x,
                                     r3.pred$marginals.hyperpar$`Precision for the Gaussian observations[3]`)
sigma2.g.mod3.df <- data.frame(marg.variance.g.r3)

(post.sigma.g.mod3 = ggplot(sigma2.g.mod3.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "brown4") +
    geom_area(alpha = .5, fill = "brown4") +
    scale_x_continuous(bquote(sigma[g]^2), breaks = c(0.0175, 0.0185)) + 
    scale_y_continuous("Density", limits = c(0, 3000)) + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))

marg.variance.z.r3 <- inla.tmarginal(function(x) 1/x,
                                     r3.pred$marginals.hyperpar$`Precision for the Gaussian observations[4]`)
sigma2.z.mod3.df <- data.frame(marg.variance.z.r3)

(post.sigma.z.mod3 = ggplot(sigma2.z.mod3.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "brown4") +
    geom_area(alpha = .5, fill = "brown4") +
    scale_x_continuous(bquote(sigma[z]^2),
                       breaks = c(0.017, 0.019)) + 
    scale_y_continuous("Density", limits = c(0, 3000)) + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))

marg.variance.z_.r3 <- inla.tmarginal(function(x) 1/x,
                                      r3.pred$marginals.hyperpar$`Precision for the Gaussian observations[5]`)
sigma2.z_.mod3.df <- data.frame(marg.variance.z_.r3)

(post.sigma.z_.mod3 = ggplot(sigma2.z_.mod3.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "brown4") +
    geom_area(alpha = .5, fill = "brown4") +
    scale_x_continuous(bquote(sigma[z_]^2),
                       breaks = c(0.0175, 0.0185)) + 
    scale_y_continuous("Density", limits = c(0, 3000)) + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))



ggarrange(post.sigma.h.mod3, 
          post.sigma.g.mod3 + rremove("ylab"), 
          post.sigma.z.mod3 + rremove("ylab"), 
          post.sigma.z_.mod3 + rremove("ylab"),
          common.legend = T, ncol = 4, nrow = 1, align = "h")


################################################################################


# Model 3: Plotting predictions
################################################################################
r3.pred.sd = sqrt(r3.pred$summary.random$idg$sd^2 +
                    (1/r3.pred$summary.hyperpar$mean[2]))  



# plot predictions vs observations
r3.pred.df = data.frame(x = complete_data[[5]]$t,
                        obs1 = complete_data[[5]]$g,
                        obs2 = complete_data[[7]]$g,
                        obs3 = complete_data[[11]]$g,
                        obs4 = complete_data[[16]]$g,
                        g1 = r3.pred$summary.random$idg$mean[1:288 + 12*288],
                        g1.lq = r3.pred$summary.random$idg$mean[1:288 + 12*288] + qnorm(0.025)*r3.pred.sd[1:288 + 12*288],
                        g1.uq = r3.pred$summary.random$idg$mean[1:288 + 12*288] + qnorm(0.975)*r3.pred.sd[1:288 + 12*288],
                        g2 = r3.pred$summary.random$idg$mean[1:288 + 13*288],
                        g2.lq = r3.pred$summary.random$idg$mean[1:288 + 13*288] + qnorm(0.025)*r3.pred.sd[1:288 + 13*288],
                        g2.uq = r3.pred$summary.random$idg$mean[1:288 + 13*288] + qnorm(0.975)*r3.pred.sd[1:288 + 13*288],
                        g3 = r3.pred$summary.random$idg$mean[1:288 + 14*288],
                        g3.lq = r3.pred$summary.random$idg$mean[1:288 + 14*288] + qnorm(0.025)*r3.pred.sd[1:288 + 14*288],
                        g3.uq = r3.pred$summary.random$idg$mean[1:288 + 14*288] + qnorm(0.975)*r3.pred.sd[1:288 + 14*288],
                        g4 = r3.pred$summary.random$idg$mean[1:288 + 15*288],
                        g4.lq = r3.pred$summary.random$idg$mean[1:288 + 15*288] + qnorm(0.025)*r3.pred.sd[1:288 + 15*288],
                        g4.uq = r3.pred$summary.random$idg$mean[1:288 + 15*288] + qnorm(0.975)*r3.pred.sd[1:288 + 15*288])

# plot predicted vs observed 
(mod1.pred1 = ggplot(aes(x = x), data = r3.pred.df) +
    geom_line(aes(y = g1, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs1, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g1.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g1.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g1.uq, ymin = g1.lq), fill = "brown4", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "brown4", "obs" = "black"),
                       labels = c("pred" = "Model 3 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0, 1)) +
    xlab("t") + theme_bw() + ggtitle("Test 1") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

(mod1.pred2 = ggplot(aes(x = x), data = r3.pred.df) +
    geom_line(aes(y = g2, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs2, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g2.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g2.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g2.uq, ymin = g2.lq), fill = "brown4", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "brown4", "obs" = "black"),
                       labels = c("pred" = "Model 3 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw() + ggtitle("Test 2") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

(mod1.pred3 = ggplot(aes(x = x), data = r3.pred.df) +
    geom_line(aes(y = g3, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs3, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g3.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g3.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g3.uq, ymin = g3.lq), fill = "brown4", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "brown4", "obs" = "black"),
                       labels = c("pred" = "Model 3 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw() + ggtitle("Test 3") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 


(mod1.pred4 = ggplot(aes(x = x), data = r3.pred.df) +
    geom_line(aes(y = g4, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs4, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g4.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g4.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g4.uq, ymin = g4.lq), fill = "brown4", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "brown4", "obs" = "black"),
                       labels = c("pred" = "Model 3 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw() + ggtitle("Test 4") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

# Combine by ggarrange
mod1.pred.plot = ggarrange(mod1.pred1 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank()),
                           NULL,
                           mod1.pred2 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                           mod1.pred3 + rremove("xlab") + rremove("ylab"), 
                           NULL,
                           mod1.pred4 + rremove("xlab") + rremove("ylab") + theme(axis.text.y = element_blank()),
                           ncol = 3, nrow = 2,
                           widths = c(1,0.05,1,1,0.05,1), common.legend = TRUE)

annotate_figure(mod1.pred.plot,
                left = text_grob(bquote(g(t)), size = 35, rot = 90), 
                bottom = text_grob("t", size = 35))
################################################################################


# Model 3: Plot posterior z(t)
################################################################################

mod1.z.df = data.frame(x = complete_data[[1]]$t, 
                       y2m = r3.pred$summary.random$idz$mean[1:288 + 288],
                       y2uq = r3.pred$summary.random$idz$`0.975quant`[1:288 + 288],
                       y2lq = r3.pred$summary.random$idz$`0.025quant`[1:288 + 288],
                       y8m = r3.pred$summary.random$idz$mean[1:288 + 7*288],
                       y8uq = r3.pred$summary.random$idz$`0.975quant`[1:288 + 7*288],
                       y8lq = r3.pred$summary.random$idz$`0.025quant`[1:288 + 7*288])


min(r3.pred$summary.random$idz$`0.025quant`)
max(r3.pred$summary.random$idz$`0.975quant`)


# Define plot objects for all posterior idz




(mod1.z2 = ggplot(aes(x = x), data = mod1.z.df) + 
    geom_line(aes(y = y2m), lwd = 1.5, col = "brown4") +
    geom_line(aes(y = y2uq), lwd = 0.5, col = "brown4") + 
    geom_line(aes(y = y2lq), lwd = 0.5, col = "brown4") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y2uq, ymin = y2lq), fill = "brown4",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.5, 11)) +
    xlab("t") + theme_bw()+ ggtitle("Train 2") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)


(mod1.z8 = ggplot(aes(x = x), data = mod1.z.df) + 
    geom_line(aes(y = y8m), lwd = 1.5, col = "brown4") +
    geom_line(aes(y = y8uq), lwd = 0.5, col = "brown4") + 
    geom_line(aes(y = y8lq), lwd = 0.5, col = "brown4") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y8uq, ymin = y8lq), fill = "brown4",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.5, 11)) +
    xlab("t") + theme_bw()+ ggtitle("Train 8") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

annotate_figure(ggarrange(mod1.z2 + rremove("xlab") + rremove("ylab"), 
                          mod1.z8 + rremove("xlab") + rremove("ylab") + theme(axis.text.y = element_blank()), 
                          ncol = 2, align = "h"),
                bottom = text_grob("t", size = 35),
                left = text_grob(bquote(eta^z~(t)), size = 35, rot = 90))


mod1.h.df = data.frame(x = complete_data[[1]]$t, 
                       y2m = r3.pred$summary.random$idh$mean[1:288 + 288],
                       y2uq = r3.pred$summary.random$idh$`0.975quant`[1:288 + 288],
                       y2lq = r3.pred$summary.random$idh$`0.025quant`[1:288 + 288],
                       y8m = r3.pred$summary.random$idh$mean[1:288 + 7*288],
                       y8uq = r3.pred$summary.random$idh$`0.975quant`[1:288 + 7*288],
                       y8lq = r3.pred$summary.random$idh$`0.025quant`[1:288 + 7*288])


min(r3.pred$summary.random$idh$`0.025quant`)
max(r3.pred$summary.random$idh$`0.975quant`)


# Define plot objects for all posterior idh




(mod1.h2 = ggplot(aes(x = x), data = mod1.h.df) + 
    geom_line(aes(y = y2m), lwd = 1.5, col = "brown4") +
    geom_line(aes(y = y2uq), lwd = 0.5, col = "brown4") + 
    geom_line(aes(y = y2lq), lwd = 0.5, col = "brown4") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y2uq, ymin = y2lq), fill = "brown4",alpha = 0.2) +
    scale_y_continuous(bquote(eta^h~(t)), limits = c(0, 1)) +
    xlab("t") + theme_bw()+ ggtitle("Train 2") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)


(mod1.h8 = ggplot(aes(x = x), data = mod1.h.df) + 
    geom_line(aes(y = y8m), lwd = 1.5, col = "brown4") +
    geom_line(aes(y = y8uq), lwd = 0.5, col = "brown4") + 
    geom_line(aes(y = y8lq), lwd = 0.5, col = "brown4") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y8uq, ymin = y8lq), fill = "brown4",alpha = 0.2) +
    scale_y_continuous(bquote(eta^h~(t)), limits = c(0, 1)) +
    xlab("t") + theme_bw()+ ggtitle("Train 8") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

annotate_figure(ggarrange(mod1.h2 + rremove("xlab") + rremove("ylab"), 
                          mod1.h8 + rremove("xlab") + rremove("ylab") + theme(axis.text.y = element_blank()), 
                          ncol = 2, align = "h"),
                bottom = text_grob("t", size = 35),
                left = text_grob(bquote(eta^h~(t)), size = 35, rot = 90))




mod1.g.df = data.frame(x = complete_data[[1]]$t, 
                       y2m = r3.pred$summary.random$idg$mean[1:288 + 288],
                       y2uq = r3.pred$summary.random$idg$`0.975quant`[1:288 + 288],
                       y2lq = r3.pred$summary.random$idg$`0.025quant`[1:288 + 288],
                       y8m = r3.pred$summary.random$idg$mean[1:288 + 7*288],
                       y8uq = r3.pred$summary.random$idg$`0.975quant`[1:288 + 7*288],
                       y8lq = r3.pred$summary.random$idg$`0.025quant`[1:288 + 7*288])


min(r3.pred$summary.random$idg$`0.025quant`)
max(r3.pred$summary.random$idg$`0.975quant`)


# Define plot objects for all posterior idg




(mod1.g2 = ggplot(aes(x = x), data = mod1.g.df) + 
    geom_line(aes(y = y2m), lwd = 1.5, col = "brown4") +
    geom_line(aes(y = y2uq), lwd = 0.5, col = "brown4") + 
    geom_line(aes(y = y2lq), lwd = 0.5, col = "brown4") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y2uq, ymin = y2lq), fill = "brown4",alpha = 0.2) +
    scale_y_continuous(bquote(eta^g~(t)), limits = c(0.35, 0.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 2") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)


(mod1.g8 = ggplot(aes(x = x), data = mod1.g.df) + 
    geom_line(aes(y = y8m), lwd = 1.5, col = "brown4") +
    geom_line(aes(y = y8uq), lwd = 0.5, col = "brown4") + 
    geom_line(aes(y = y8lq), lwd = 0.5, col = "brown4") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y8uq, ymin = y8lq), fill = "brown4",alpha = 0.2) +
    scale_y_continuous(bquote(eta^g~(t)), limits = c(0.35, 0.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 8") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

annotate_figure(ggarrange(mod1.g2 + rremove("xlab") + rremove("ylab"), 
                          mod1.g8 + rremove("xlab") + rremove("ylab") + theme(axis.text.y = element_blank()), 
                          ncol = 2, align = "h"),
                bottom = text_grob("t", size = 35),
                left = text_grob(bquote(eta^g~(t)), size = 35, rot = 90))


