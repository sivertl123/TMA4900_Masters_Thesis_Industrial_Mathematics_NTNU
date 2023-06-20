#rm(list = ls())

source("TMA4900_Data.R")

library(INLA)
library(ggplot2)


# Model 2: a = b0 + b1*T - Prediction
################################################################################

# Model 2: Running model
################################################################################
inla.model2.replicate.pred = function(complete_data){
  
  # number of datasets
  N = 16
  
  # length of each dataset
  n = nrow(complete_data[[1]])
  m = n-1
  
  # define response matrix
  Y1 = Y2 = Y3 = Y4 = Y5 = wzb1_ = c()
  for(i in c(1,2,3,4,6,8,9,10,12,14,15,17)){
    Y1 = c(Y1, rep(10, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y2 = c(Y2, rep(NA, n), complete_data[[i]]$h, rep(NA, n), rep(NA, m), rep(NA, m))
    Y3 = c(Y3, rep(NA, n), rep(NA, n), complete_data[[i]]$g, rep(NA, m), rep(NA, m))
    Y4 = c(Y4, rep(NA, n), rep(NA, n), rep(NA, n), rep(0, m), rep(NA, m))
    Y5 = c(Y5, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(0, m))
    wzb1_ = c(wzb1_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$T[1:m])
  }
  
  for(i in c(5,7,11,16)){
    Y1 = c(Y1, rep(10, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y2 = c(Y2, rep(NA, n), complete_data[[i]]$h, rep(NA, n), rep(NA, m), rep(NA, m))
    Y3 = c(Y3, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y4 = c(Y4, rep(NA, n), rep(NA, n), rep(NA, n), rep(0, m), rep(NA, m))
    Y5 = c(Y5, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(0, m))
    wzb1_ = c(wzb1_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$T[1:m])
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
  
  
  # define residuals
  ide0 = rep(c(1:n, rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m)), N)
  ideh = rep(c(rep(NA, n), 1:n, rep(NA, n), rep(NA, m), rep(NA, m)), N)
  ideg = rep(c(rep(NA, n), rep(NA, n), 1:n, rep(NA, m), rep(NA, m)), N)
  idem = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), 1:m, rep(NA, m)), N)
  idea = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), 1:m), N)
  
  
  # define inla data
  inla.data = data.frame("Y" = Y, "idh" = idh, "wh" = wh, "idg" = idg, "wg" = wg, 
                         "idz" = idz, "wz" = wz, "idz_" = idz_, "wz_" = wz_, 
                         "idzb0_" = idzb0_, "wzb0_", wzb0_,
                         "idzb1_" = idzb1_, "wzb1_", wzb1_,
                         "ide0" = ide0, "ideh" = ideh, "ideg" = ideg, 
                         "idem" = idem, "idea" = idea, "rep" = rep)
  
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
    f(idzb1_, wzb1_, copy = "idz", hyper = list(beta = list(fixed = FALSE, prior = "gaussian", param = c(0, 1e4)))) #+
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
r2.pred = inla.model2.replicate.pred(complete_data)

# check results
plot(r2.pred$summary.random$idg$mean, type = "l")
r2.pred$summary.hyperpar


################################################################################

# Model 2: Plot posterior alpha
################################################################################
post.b0.ggdf.2 = data.frame("x" = r2.pred$marginals.hyperpar$`Beta for idzb0_`[,1],
                            "y" = r2.pred$marginals.hyperpar$`Beta for idzb0_`[,2])

(post.b0.plot.2 = ggplot(aes(x = x, y = y), data = post.b0.ggdf.2) + 
    geom_line(lwd = 1.5, col = "darkblue") +
    geom_area(alpha = .5, fill = "darkblue") +
    scale_y_continuous("Density") + 
    scale_x_continuous(bquote(beta[0]),
                       breaks = c(0.090, 0.094)) + 
    theme_bw() +
    theme(axis.title = element_text(size = 35),
          axis.text = element_text(size = 20)))

post.b1.ggdf.2 = data.frame("x" = r2.pred$marginals.hyperpar$`Beta for idzb1_`[,1],
                            "y" = r2.pred$marginals.hyperpar$`Beta for idzb1_`[,2])

(post.b1.plot.2 = ggplot(aes(x = x, y = y), data = post.b1.ggdf.2) + 
    geom_line(lwd = 1.5, col = "darkblue") +
    geom_area(alpha = .5, fill = "darkblue") +
    scale_y_continuous("Density") + 
    scale_x_continuous(bquote(beta[1]),
                       breaks = c(-0.0045,-0.0041)) +
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))

ggarrange(post.b0.plot.2, post.b1.plot.2 + rremove("ylab"),
          align = "h")



################################################################################

# Model 2: Plot posterior standard deviations
################################################################################

marg.variance.h.r2 <- inla.tmarginal(function(x) 1/x,
                                     r2.pred$marginals.hyperpar$`Precision for the Gaussian observations[2]`)
sigma2.h.mod2.df <- data.frame(marg.variance.h.r2)

(post.sigma.h.mod2 = ggplot(sigma2.h.mod2.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "darkblue") +
    geom_area(alpha = .5, fill = "darkblue") +
    scale_x_continuous(bquote(sigma[h]^2),
                       breaks = c(0.0174, 0.0184)) + 
    scale_y_continuous("Density", limits = c(0, 2500)) + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))


marg.variance.g.r2 <- inla.tmarginal(function(x) 1/x,
                                     r2.pred$marginals.hyperpar$`Precision for the Gaussian observations[3]`)
sigma2.g.mod2.df <- data.frame(marg.variance.g.r2)

(post.sigma.g.mod2 = ggplot(sigma2.g.mod2.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "darkblue") +
    geom_area(alpha = .5, fill = "darkblue") +
    scale_x_continuous(bquote(sigma[g]^2),
                       breaks = c(0.016, 0.017)) + 
    scale_y_continuous("Density", limits = c(0, 2500)) + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))

marg.variance.z.r2 <- inla.tmarginal(function(x) 1/x,
                                     r2.pred$marginals.hyperpar$`Precision for the Gaussian observations[4]`)
sigma2.z.mod2.df <- data.frame(marg.variance.z.r2)

(post.sigma.z.mod2 = ggplot(sigma2.z.mod2.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "darkblue") +
    geom_area(alpha = .5, fill = "darkblue") +
    scale_x_continuous(bquote(sigma[z]^2),
                       breaks = c(0.0172, 0.0182)) + 
    scale_y_continuous("Density", limits = c(0, 2500)) + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))

marg.variance.z_.r2 <- inla.tmarginal(function(x) 1/x,
                                      r2.pred$marginals.hyperpar$`Precision for the Gaussian observations[5]`)
sigma2.z_.mod2.df <- data.frame(marg.variance.z_.r2)

(post.sigma.z_.mod2 = ggplot(sigma2.z_.mod2.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "darkblue") +
    geom_area(alpha = .5, fill = "darkblue") +
    scale_x_continuous(bquote(sigma[z_]^2),
                       breaks = c(0.0175, 0.0185)) + 
    scale_y_continuous("Density", limits = c(0, 2500)) + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))




ggarrange(post.sigma.h.mod2, 
          post.sigma.g.mod2 + rremove("ylab"), 
          post.sigma.z.mod2 + rremove("ylab"), 
          post.sigma.z_.mod2 + rremove("ylab"),
          common.legend = T, ncol = 4, nrow = 1, align = "h")

################################################################################



# Model 2: Plot predictions
################################################################################
r2.pred.sd = sqrt(r2.pred$summary.random$idg$sd^2 +
                    (1/r2.pred$summary.hyperpar$mean[2]))  



# plot predictions vs observations
r2.pred.df = data.frame(x = complete_data[[5]]$t,
                        obs1 = complete_data[[5]]$g,
                        obs2 = complete_data[[7]]$g,
                        obs3 = complete_data[[11]]$g,
                        obs4 = complete_data[[16]]$g,
                        g1 = r2.pred$summary.random$idg$mean[1:288 + 12*288],
                        g1.lq = r2.pred$summary.random$idg$mean[1:288 + 12*288] + qnorm(0.025)*r2.pred.sd[1:288 + 12*288],
                        g1.uq = r2.pred$summary.random$idg$mean[1:288 + 12*288] + qnorm(0.975)*r2.pred.sd[1:288 + 12*288],
                        g2 = r2.pred$summary.random$idg$mean[1:288 + 13*288],
                        g2.lq = r2.pred$summary.random$idg$mean[1:288 + 13*288] + qnorm(0.025)*r2.pred.sd[1:288 + 13*288],
                        g2.uq = r2.pred$summary.random$idg$mean[1:288 + 13*288] + qnorm(0.975)*r2.pred.sd[1:288 + 13*288],
                        g3 = r2.pred$summary.random$idg$mean[1:288 + 14*288],
                        g3.lq = r2.pred$summary.random$idg$mean[1:288 + 14*288] + qnorm(0.025)*r2.pred.sd[1:288 + 14*288],
                        g3.uq = r2.pred$summary.random$idg$mean[1:288 + 14*288] + qnorm(0.975)*r2.pred.sd[1:288 + 14*288],
                        g4 = r2.pred$summary.random$idg$mean[1:288 + 15*288],
                        g4.lq = r2.pred$summary.random$idg$mean[1:288 + 15*288] + qnorm(0.025)*r2.pred.sd[1:288 + 15*288],
                        g4.uq = r2.pred$summary.random$idg$mean[1:288 + 15*288] + qnorm(0.975)*r2.pred.sd[1:288 + 15*288])


# plot predicted vs observed 
(mod2.pred1 = ggplot(aes(x = x), data = r2.pred.df) +
    geom_line(aes(y = g1, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs1, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g1.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g1.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g1.uq, ymin = g1.lq), fill = "darkblue", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "darkblue", "obs" = "black"),
                       labels = c("pred" = "Model 2 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw() + ggtitle("Test 1") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

(mod2.pred2 = ggplot(aes(x = x), data = r2.pred.df) +
    geom_line(aes(y = g2, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs2, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g2.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g2.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g2.uq, ymin = g2.lq), fill = "darkblue", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "darkblue", "obs" = "black"),
                       labels = c("pred" = "Model 2 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw() + ggtitle("Test 2") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

(mod2.pred3 = ggplot(aes(x = x), data = r2.pred.df) +
    geom_line(aes(y = g3, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs3, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g3.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g3.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g3.uq, ymin = g3.lq), fill = "darkblue", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "darkblue", "obs" = "black"),
                       labels = c("pred" = "Model 2 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw() + ggtitle("Test 3") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 


(mod2.pred4 = ggplot(aes(x = x), data = r2.pred.df) +
    geom_line(aes(y = g4, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs4, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g4.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g4.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g4.uq, ymin = g4.lq), fill = "darkblue", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "darkblue", "obs" = "black"),
                       labels = c("pred" = "Model 2 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw()+ ggtitle("Test 4") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

# Combine by ggarrange
mod2.pred.plot = ggarrange(mod2.pred1 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank()),
                           NULL,
                           mod2.pred2 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                           mod2.pred3 + rremove("xlab") + rremove("ylab"), 
                           NULL,
                           mod2.pred4 + rremove("xlab") + rremove("ylab") + theme(axis.text.y = element_blank()),
                           ncol = 3, nrow = 2,
                           widths = c(1,0.05,1,1,0.05,1), common.legend = TRUE)

annotate_figure(mod2.pred.plot,
                left = text_grob("g(t)", size = 35, rot = 90), 
                bottom = text_grob("t", size = 35))
################################################################################

# Model 2: Plot posterior z
################################################################################
mod2.z.df = data.frame(x = complete_data[[1]]$t, 
                       y1m = r2.pred$summary.random$idz$mean[1:288],
                       y1uq = r2.pred$summary.random$idz$`0.975quant`[1:288],
                       y1lq = r2.pred$summary.random$idz$`0.025quant`[1:288],
                       y2m = r2.pred$summary.random$idz$mean[1:288 + 288],
                       y2uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 288],
                       y2lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 288],
                       y3m = r2.pred$summary.random$idz$mean[1:288 + 2*288],
                       y3uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 2*288],
                       y3lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 2*288],
                       y4m = r2.pred$summary.random$idz$mean[1:288 + 3*288],
                       y4uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 3*288],
                       y4lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 3*288],
                       y5m = r2.pred$summary.random$idz$mean[1:288 + 4*288],
                       y5uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 4*288],
                       y5lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 4*288],
                       y6m = r2.pred$summary.random$idz$mean[1:288 + 5*288],
                       y6uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 5*288],
                       y6lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 5*288],
                       y7m = r2.pred$summary.random$idz$mean[1:288 + 6*288],
                       y7uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 6*288],
                       y7lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 6*288],
                       y8m = r2.pred$summary.random$idz$mean[1:288 + 7*288],
                       y8uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 7*288],
                       y8lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 7*288],
                       y9m = r2.pred$summary.random$idz$mean[1:288 + 8*288],
                       y9uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 8*288],
                       y9lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 8*288],
                       y10m = r2.pred$summary.random$idz$mean[1:288 + 9*288],
                       y10uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 9*288],
                       y10lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 9*288],
                       y11m = r2.pred$summary.random$idz$mean[1:288 + 10*288],
                       y11uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 10*288],
                       y11lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 10*288],
                       y12m = r2.pred$summary.random$idz$mean[1:288 + 11*288],
                       y12uq = r2.pred$summary.random$idz$`0.975quant`[1:288 + 11*288],
                       y12lq = r2.pred$summary.random$idz$`0.025quant`[1:288 + 11*288])


min(r2.pred$summary.random$idz$`0.025quant`) #9.35
max(r2.pred$summary.random$idz$`0.975quant`) #10.44
# Define plot objects for all posterior idz


(mod2.z1 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y1m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y1uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y1lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y1uq, ymin = y1lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 1") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)


(mod2.z2 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y2m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y2uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y2lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y2uq, ymin = y2lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 2") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

(mod2.z3 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y3m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y3uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y3lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y3uq, ymin = y3lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 3") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

(mod2.z4 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y4m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y4uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y4lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y4uq, ymin = y4lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 4") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

(mod2.z5 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y5m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y5uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y5lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y5uq, ymin = y5lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 5") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

(mod2.z6 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y6m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y6uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y6lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y6uq, ymin = y6lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 6") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

(mod2.z7 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y7m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y7uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y7lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y7uq, ymin = y7lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 7") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

(mod2.z8 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y8m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y8uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y8lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y8uq, ymin = y8lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 8") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

(mod2.z9 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y9m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y9uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y9lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y9uq, ymin = y9lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 9") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

(mod2.z10 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y10m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y10uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y10lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y10uq, ymin = y10lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 10") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

(mod2.z11 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y11m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y11uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y11lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y11uq, ymin = y11lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 11") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)

(mod2.z12 = ggplot(aes(x = x), data = mod2.z.df) + 
    geom_line(aes(y = y12m), lwd = 1.5, col = "darkblue") +
    geom_line(aes(y = y12uq), lwd = 0.5, col = "darkblue") + 
    geom_line(aes(y = y12lq), lwd = 0.5, col = "darkblue") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y12uq, ymin = y12lq), fill = "darkblue",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.35, 10.45)) +
    xlab("t") + theme_bw()+ ggtitle("Train 12") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)





# Plot all idz in a 4*3 grid
mod2.z.plot = ggarrange(mod2.z1 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank()),
                        NULL,
                        mod2.z2 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                        NULL,
                        mod2.z3 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                        NULL,
                        mod2.z4 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()),
                        mod2.z5 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank()), 
                        NULL,
                        mod2.z6 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()),
                        NULL,
                        mod2.z7 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()),
                        NULL, 
                        mod2.z8 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()),
                        mod2.z9 + rremove("ylab") + rremove("xlab"), 
                        NULL,
                        mod2.z10 + rremove("ylab") + rremove("xlab") + theme(axis.text.y = element_blank()), 
                        NULL,
                        mod2.z11 + rremove("ylab") + rremove("xlab") + theme(axis.text.y = element_blank()),
                        NULL,
                        mod2.z12 + rremove("ylab") + rremove("xlab") + theme(axis.text.y = element_blank()),
                        ncol = 7, nrow = 3,
                        widths = c(1,0.05,1,0.05,1,0.05,1,1,0.05,1,0.05,1,0.05,1,1,0.05,1,0.05,1,0.05,1))

annotate_figure(mod2.z.plot,
                left = text_grob(bquote(eta^z~(t)), size = 40, rot = 90), 
                bottom = text_grob("t", size = 40))



################################################################################

