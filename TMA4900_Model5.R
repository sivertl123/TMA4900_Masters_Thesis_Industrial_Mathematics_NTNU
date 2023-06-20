#rm(list = ls())

source("TMA4900_Data.R")
source("TMA4900_Model1.R")
library(INLA)
library(ggplot2)


# Model 4: a = b0 + b1*T + b2*M
################################################################################

# Model 4: Run Model
################################################################################

inla.model5.replicate.pred = function(complete_data){
  
  # number of datasets
  N = 16
  
  # length of each dataset
  n = nrow(complete_data[[1]])
  m = n-1
  
  # define response matrix
  Y1 = Y2 = Y3 = Y4 = Y5 = wzb1_ = wzb2_ = c()
  for(i in c(1,2,3,4,6,8,9,10,12,14,15,17)){
    Y1 = c(Y1, rep(10, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y2 = c(Y2, rep(NA, n), complete_data[[i]]$h, rep(NA, n), rep(NA, m), rep(NA, m))
    Y3 = c(Y3, rep(NA, n), rep(NA, n), complete_data[[i]]$g, rep(NA, m), rep(NA, m))
    Y4 = c(Y4, rep(NA, n), rep(NA, n), rep(NA, n), rep(0, m), rep(NA, m))
    Y5 = c(Y5, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(0, m))
    wzb1_ = c(wzb1_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$T[1:m])
    wzb2_ = c(wzb2_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$M[1:m])
  }
  
  for(i in c(5,7,11,16)){
    Y1 = c(Y1, rep(10, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y2 = c(Y2, rep(NA, n), complete_data[[i]]$h, rep(NA, n), rep(NA, m), rep(NA, m))
    Y3 = c(Y3, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y4 = c(Y4, rep(NA, n), rep(NA, n), rep(NA, n), rep(0, m), rep(NA, m))
    Y5 = c(Y5, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(0, m))
    wzb1_ = c(wzb1_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$T[1:m])
    wzb2_ = c(wzb2_, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), -1*complete_data[[i]]$M[1:m])
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
  
  # latent b1*T*z(t-1) 
  idzb2_ = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), 1:m), N)
  
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
    f(idzb1_, wzb1_, copy = "idz", hyper = list(beta = list(fixed = FALSE, prior = "gaussian", param = c(0, 1e4)))) +
    f(idzb2_, wzb2_, copy = "idz", hyper = list(beta = list(fixed = FALSE, prior = "gaussian", param = c(0, 1e2)))) #+
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
r5.pred = inla.model5.replicate.pred(complete_data)

# check results
plot(r5.pred$summary.random$idg$mean, type = "l")
plot(r5.pred$summary.random$idz$mean, type = "l")
r5.pred$summary.hyperpar

################################################################################

# Model 4: Plot posterior alpha
################################################################################
post.b0.ggdf.3 = data.frame("x" = r5.pred$marginals.hyperpar$`Beta for idzb0_`[,1],
                            "y" = r5.pred$marginals.hyperpar$`Beta for idzb0_`[,2])

(post.b0.plot.3 = ggplot(aes(x = x, y = y), data = post.b0.ggdf.3) + 
    geom_line(lwd = 1.5, col = "darkseagreen") +
    geom_area(fill = "darkseagreen", alpha = 0.5) +
    scale_y_continuous("Density") + 
    scale_x_continuous(bquote(beta[0]),
                       breaks = c(0.064, 0.072)) + 
    theme_bw() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 35)))

post.b1.ggdf.3 = data.frame("x" = r5.pred$marginals.hyperpar$`Beta for idzb1_`[,1],
                            "y" = r5.pred$marginals.hyperpar$`Beta for idzb1_`[,2])

(post.b1.plot.3 = ggplot(aes(x = x, y = y), data = post.b1.ggdf.3) + 
    geom_line(lwd = 1.5, col = "darkseagreen") +
    geom_area(fill = "darkseagreen", alpha = 0.5) +
    scale_y_continuous("Density", lim = c(0,23000)) + 
    scale_x_continuous(bquote(beta[1]),
                       breaks = c(-0.0022, -0.0018)) + 
    theme_bw() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 35)))

post.b2.ggdf.3 = data.frame("x" = r5.pred$marginals.hyperpar$`Beta for idzb2_`[,1],
                            "y" = r5.pred$marginals.hyperpar$`Beta for idzb2_`[,2])

(post.b2.plot.3 = ggplot(aes(x = x, y = y), data = post.b2.ggdf.3) + 
    geom_line(lwd = 1.5, col = "darkseagreen") +
    geom_area(fill = "darkseagreen", alpha = 0.5) +
    scale_y_continuous("Density", lim = c(0, 23000)) + 
    scale_x_continuous(bquote(beta[2]),
                       breaks = c(-0.0042, -0.0040)) + 
    theme_bw() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 35)))

ggarrange(post.b0.plot.3, 
          post.b1.plot.3 + rremove("ylab"),
          post.b2.plot.3 + rremove("ylab"), 
          ncol = 3, nrow = 1, align = "h")





################################################################################

# Model 4: Plot posterior standard deviations
################################################################################

marg.variance.h.r5 <- inla.tmarginal(function(x) 1/x,
                                     r5.pred$marginals.hyperpar$`Precision for the Gaussian observations[2]`)
sigma2.h.mod5.df <- data.frame(marg.variance.h.r5)

(post.sigma.h.mod5 = ggplot(sigma2.h.mod5.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "darkseagreen") +
    geom_area(alpha = .5, fill = "darkseagreen") +
    scale_x_continuous(bquote(sigma[h]^2),
                       breaks = c(0.015, 0.016)) + 
    scale_y_continuous("Density") + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))


marg.variance.g.r5 <- inla.tmarginal(function(x) 1/x,
                                     r5.pred$marginals.hyperpar$`Precision for the Gaussian observations[3]`)
sigma2.g.mod5.df <- data.frame(marg.variance.g.r5)

(post.sigma.g.mod5 = ggplot(sigma2.g.mod5.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "darkseagreen") +
    geom_area(alpha = .5, fill = "darkseagreen") +
    scale_x_continuous(bquote(sigma[g]^2),
                       breaks = c(0.0171, 0.0175)) + 
    scale_y_continuous("Density") + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))

marg.variance.z.r5 <- inla.tmarginal(function(x) 1/x,
                                     r5.pred$marginals.hyperpar$`Precision for the Gaussian observations[4]`)
sigma2.z.mod5.df <- data.frame(marg.variance.z.r5)

(post.sigma.z.mod5 = ggplot(sigma2.z.mod5.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "darkseagreen") +
    geom_area(alpha = .5, fill = "darkseagreen") +
    scale_x_continuous(bquote(sigma[z]^2),
                       breaks = c(0.0159, 0.0162)) + 
    scale_y_continuous("Density") + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))

marg.variance.z_.r5 <- inla.tmarginal(function(x) 1/x,
                                      r5.pred$marginals.hyperpar$`Precision for the Gaussian observations[5]`)
sigma2.z_.mod5.df <- data.frame(marg.variance.z_.r5)

(post.sigma.z_.mod5 = ggplot(sigma2.z_.mod5.df, aes(x, y)) +
    geom_line(lwd = 1.5, color = "darkseagreen") +
    geom_area(alpha = .5, fill = "darkseagreen") +
    scale_x_continuous(bquote(sigma[z_]^2),
                       breaks = c(0.0155, 0.0160)) + 
    scale_y_continuous("Density") + 
    theme_bw() + theme(axis.title = element_text(size = 35),
                       axis.text = element_text(size = 20)))






ggarrange(post.sigma.h.mod5, 
          post.sigma.g.mod5 + rremove("ylab"), 
          post.sigma.z.mod5 + rremove("ylab"), 
          post.sigma.z_.mod5 + rremove("ylab"),
          common.legend = T, ncol = 4, nrow = 1,
          align = "h")

################################################################################


# Model 4: Plot predictions
################################################################################

r5.pred.sd = sqrt(r5.pred$summary.random$idg$sd^2 +
                    (1/r5.pred$summary.hyperpar$mean[2]))  



# plot predictions vs observations
r5.pred.df = data.frame(x = complete_data[[5]]$t,
                        obs1 = complete_data[[5]]$g,
                        obs2 = complete_data[[7]]$g,
                        obs3 = complete_data[[11]]$g,
                        obs4 = complete_data[[16]]$g,
                        g1 = r5.pred$summary.random$idg$mean[1:288 + 12*288],
                        g1.lq = r5.pred$summary.random$idg$mean[1:288 + 12*288] + qnorm(0.025)*r5.pred.sd[1:288 + 12*288],
                        g1.uq = r5.pred$summary.random$idg$mean[1:288 + 12*288] + qnorm(0.975)*r5.pred.sd[1:288 + 12*288],
                        g2 = r5.pred$summary.random$idg$mean[1:288 + 13*288],
                        g2.lq = r5.pred$summary.random$idg$mean[1:288 + 13*288] + qnorm(0.025)*r5.pred.sd[1:288 + 13*288],
                        g2.uq = r5.pred$summary.random$idg$mean[1:288 + 13*288] + qnorm(0.975)*r5.pred.sd[1:288 + 13*288],
                        g3 = r5.pred$summary.random$idg$mean[1:288 + 14*288],
                        g3.lq = r5.pred$summary.random$idg$mean[1:288 + 14*288] + qnorm(0.025)*r5.pred.sd[1:288 + 14*288],
                        g3.uq = r5.pred$summary.random$idg$mean[1:288 + 14*288] + qnorm(0.975)*r5.pred.sd[1:288 + 14*288],
                        g4 = r5.pred$summary.random$idg$mean[1:288 + 15*288],
                        g4.lq = r5.pred$summary.random$idg$mean[1:288 + 15*288] + qnorm(0.025)*r5.pred.sd[1:288 + 15*288],
                        g4.uq = r5.pred$summary.random$idg$mean[1:288 + 15*288] + qnorm(0.975)*r5.pred.sd[1:288 + 15*288])


# plot predicted vs observed 
(mod3.pred1 = ggplot(aes(x = x), data = r5.pred.df) +
    geom_line(aes(y = g1, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs1, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g1.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g1.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g1.uq, ymin = g1.lq), fill = "darkseagreen", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "darkseagreen", "obs" = "black"),
                       labels = c("pred" = "Model 4 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(-0.1, 1)) +
    xlab("t") + theme_bw() + ggtitle("Test 1") +
    theme(axis.text = element_text(size = 25), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

(mod3.pred2 = ggplot(aes(x = x), data = r5.pred.df) +
    geom_line(aes(y = g2, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs2, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g2.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g2.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g2.uq, ymin = g2.lq), fill = "darkseagreen", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "darkseagreen", "obs" = "black"),
                       labels = c("pred" = "Model 4 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(-0.1, 1)) +
    xlab("t") + theme_bw() + ggtitle("Test 2") +
    theme(axis.text = element_text(size = 25), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

(mod3.pred3 = ggplot(aes(x = x), data = r5.pred.df) +
    geom_line(aes(y = g3, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs3, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g3.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g3.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g3.uq, ymin = g3.lq), fill = "darkseagreen", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "darkseagreen", "obs" = "black"),
                       labels = c("pred" = "Model 4 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(-0.1, 1)) +
    xlab("t") + theme_bw() + ggtitle("Test 3") +
    theme(axis.text = element_text(size = 25), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 


(mod3.pred4 = ggplot(aes(x = x), data = r5.pred.df) +
    geom_line(aes(y = g4, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs4, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g4.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g4.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g4.uq, ymin = g4.lq), fill = "darkseagreen", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "darkseagreen", "obs" = "black"),
                       labels = c("pred" = "Model 4 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(-0.1, 1)) +
    xlab("t") + theme_bw() + ggtitle("Test 4") +
    theme(axis.text = element_text(size = 25), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

# Combine by ggarrange
mod3.pred.plot = ggarrange(mod3.pred1 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank()),
                           NULL,
                           mod3.pred2 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                           mod3.pred3 + rremove("xlab") + rremove("ylab"), 
                           NULL,
                           mod3.pred4 + rremove("xlab") + rremove("ylab") + theme(axis.text.y = element_blank()),
                           ncol = 3, nrow = 2,
                           widths = c(1,0.05,1,1,0.05,1), common.legend = TRUE)

annotate_figure(mod3.pred.plot,
                left = text_grob(bquote(g(t)), size = 35, rot = 90), 
                bottom = text_grob("t", size = 35))

################################################################################

# Model 4: Plot posterior z(t)
################################################################################

mod3.z.df = data.frame(x = complete_data[[1]]$t, 
                       y2m = r5.pred$summary.random$idz$mean[1:288 + 288],
                       y2uq = r5.pred$summary.random$idz$`0.975quant`[1:288 + 288],
                       y2lq = r5.pred$summary.random$idz$`0.025quant`[1:288 + 288],
                       y8m = r5.pred$summary.random$idz$mean[1:288 + 7*288],
                       y8uq = r5.pred$summary.random$idz$`0.975quant`[1:288 + 7*288],
                       y8lq = r5.pred$summary.random$idz$`0.025quant`[1:288 + 7*288])


min(r5.pred$summary.random$idz$`0.025quant`) #9.47
max(r5.pred$summary.random$idz$`0.975quant`) #10.37
# Define plot objects for all posterior idz




(mod3.z2 = ggplot(aes(x = x), data = mod3.z.df) + 
    geom_line(aes(y = y2m), lwd = 1.5, col = "darkseagreen") +
    geom_line(aes(y = y2uq), lwd = 0.5, col = "darkseagreen") + 
    geom_line(aes(y = y2lq), lwd = 0.5, col = "darkseagreen") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y2uq, ymin = y2lq), fill = "darkseagreen",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.5, 11)) +
    xlab("t") + theme_bw()+ ggtitle("Train 2") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)


(mod3.z8 = ggplot(aes(x = x), data = mod3.z.df) + 
    geom_line(aes(y = y8m), lwd = 1.5, col = "darkseagreen") +
    geom_line(aes(y = y8uq), lwd = 0.5, col = "darkseagreen") + 
    geom_line(aes(y = y8lq), lwd = 0.5, col = "darkseagreen") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y8uq, ymin = y8lq), fill = "darkseagreen",alpha = 0.2) +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(9.5, 11)) +
    xlab("t") + theme_bw()+ ggtitle("Train 8") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)


# Plot all idz in a 4*3 grid
mod3.z.plot = ggarrange(mod1.z2 + rremove("xlab") + rremove("ylab"),
                        mod3.z2 + rremove("xlab") + rremove("ylab"),
                        mod1.z8 + rremove("xlab") + rremove("ylab"),
                        mod3.z8 + rremove("xlab") + rremove("ylab"),  
                        ncol = 2, nrow = 2, align = "hv")
                        
annotate_figure(mod3.z.plot,
                left = text_grob(bquote(eta^z~(t)), size = 35, rot = 90), 
                bottom = text_grob("t", size = 35))


################################################################################

mod3.h.df = data.frame(x = complete_data[[1]]$t, 
                       y2m = r5.pred$summary.random$idh$mean[1:288 + 288],
                       y2uq = r5.pred$summary.random$idh$`0.975quant`[1:288 + 288],
                       y2lq = r5.pred$summary.random$idh$`0.025quant`[1:288 + 288],
                       y8m = r5.pred$summary.random$idh$mean[1:288 + 7*288],
                       y8uq = r5.pred$summary.random$idh$`0.975quant`[1:288 + 7*288],
                       y8lq = r5.pred$summary.random$idh$`0.025quant`[1:288 + 7*288])


min(r5.pred$summary.random$idh$`0.025quant`) #9.47
max(r5.pred$summary.random$idh$`0.975quant`) #10.37
# Define plot objects for all posterior idh




(mod3.h2 = ggplot(aes(x = x), data = mod3.h.df) + 
    geom_line(aes(y = y2m), lwd = 1.5, col = "darkseagreen") +
    geom_line(aes(y = y2uq), lwd = 0.5, col = "darkseagreen") + 
    geom_line(aes(y = y2lq), lwd = 0.5, col = "darkseagreen") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y2uq, ymin = y2lq), fill = "darkseagreen",alpha = 0.2) +
    scale_y_continuous(bquote(eta^h~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw()+ ggtitle("Train 2") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)


(mod3.h8 = ggplot(aes(x = x), data = mod3.h.df) + 
    geom_line(aes(y = y8m), lwd = 1.5, col = "darkseagreen") +
    geom_line(aes(y = y8uq), lwd = 0.5, col = "darkseagreen") + 
    geom_line(aes(y = y8lq), lwd = 0.5, col = "darkseagreen") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y8uq, ymin = y8lq), fill = "darkseagreen",alpha = 0.2) +
    scale_y_continuous(bquote(eta^h~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw()+ ggtitle("Train 8") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)


# Plot all idh in a 4*3 grid
mod3.h.plot = ggarrange(mod1.h2 + rremove("xlab") + rremove("ylab"),
                        mod3.h2 + rremove("xlab") + rremove("ylab"),
                        mod1.h8 + rremove("xlab") + rremove("ylab"),
                        mod3.h8 + rremove("xlab") + rremove("ylab"),  
                        ncol = 2, nrow = 2, align = "hv")

annotate_figure(mod3.h.plot,
                left = text_grob(bquote(eta^h~(t)), size = 35, rot = 90), 
                bottom = text_grob("t", size = 35))

################################################################################

mod3.g.df = data.frame(x = complete_data[[1]]$t, 
                       y2m = r5.pred$summary.random$idg$mean[1:288 + 288],
                       y2uq = r5.pred$summary.random$idg$`0.975quant`[1:288 + 288],
                       y2lq = r5.pred$summary.random$idg$`0.025quant`[1:288 + 288],
                       y8m = r5.pred$summary.random$idg$mean[1:288 + 7*288],
                       y8uq = r5.pred$summary.random$idg$`0.975quant`[1:288 + 7*288],
                       y8lq = r5.pred$summary.random$idg$`0.025quant`[1:288 + 7*288])


min(r5.pred$summary.random$idg$`0.025quant`) #9.47
max(r5.pred$summary.random$idg$`0.975quant`) #10.37
# Define plot objects for all posterior idg




(mod3.g2 = ggplot(aes(x = x), data = mod3.g.df) + 
    geom_line(aes(y = y2m), lwd = 1.5, col = "darkseagreen") +
    geom_line(aes(y = y2uq), lwd = 0.5, col = "darkseagreen") + 
    geom_line(aes(y = y2lq), lwd = 0.5, col = "darkseagreen") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y2uq, ymin = y2lq), fill = "darkseagreen",alpha = 0.2) +
    scale_y_continuous(bquote(eta^g~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw()+ ggtitle("Train 2") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)


(mod3.g8 = ggplot(aes(x = x), data = mod3.g.df) + 
    geom_line(aes(y = y8m), lwd = 1.5, col = "darkseagreen") +
    geom_line(aes(y = y8uq), lwd = 0.5, col = "darkseagreen") + 
    geom_line(aes(y = y8lq), lwd = 0.5, col = "darkseagreen") +
    scale_x_datetime(date_labels = "%H:%M") +
    geom_ribbon(aes(ymax = y8uq, ymin = y8lq), fill = "darkseagreen",alpha = 0.2) +
    scale_y_continuous(bquote(eta^g~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw()+ ggtitle("Train 8") +
    theme(axis.text = element_text(size = 20),
          plot.title = element_text(size = 30, face = "bold"))
)


mod1.g2 = mod1.g2 + scale_y_continuous(bquote(eta^g~(t)), limits = c(0,1))
mod1.g8 = mod1.g8 + scale_y_continuous(bquote(eta^g~(t)), limits = c(0,1))

# Plot all idg in a 4*3 grid
mod3.g.plot = ggarrange(mod1.g2 + rremove("xlab") + rremove("ylab"),
                        mod3.g2 + rremove("xlab") + rremove("ylab"),
                        mod1.g8 + rremove("xlab") + rremove("ylab"),
                        mod3.g8 + rremove("xlab") + rremove("ylab"),  
                        ncol = 2, nrow = 2, align = "hv")

annotate_figure(mod3.g.plot,
                left = text_grob(bquote(eta^g~(t)), size = 35, rot = 90), 
                bottom = text_grob("t", size = 35))
