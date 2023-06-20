# Script for semi-informed predictions of Model 1



inla.model1.replicate.ipred = function(complete_data){
  # number of datasets
  N = 16
  
  # length of each dataset
  n = nrow(complete_data[[1]])
  m = n-1
  
  # define response matrix
  Y1 = Y2 = Y3 = Y4 = Y5 = c()
  for(i in c(1,2,3,4,6,8,9,10,12,14,15,17)){
    Y1 = c(Y1, rep(10, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y2 = c(Y2, rep(NA, n), complete_data[[i]]$h, rep(NA, n), rep(NA, m), rep(NA, m))
    Y3 = c(Y3, rep(NA, n), rep(NA, n), complete_data[[i]]$g, rep(NA, m), rep(NA, m))
    Y4 = c(Y4, rep(NA, n), rep(NA, n), rep(NA, n), rep(0, m), rep(NA, m))
    Y5 = c(Y5, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(0, m))
  }
  
  ng = n%/%2 
  
  for(i in c(5,7,11,16)){
    Y1 = c(Y1, rep(10, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m))
    Y2 = c(Y2, rep(NA, n), complete_data[[i]]$h, rep(NA, n), rep(NA, m), rep(NA, m))
    Y3 = c(Y3, rep(NA, n), rep(NA, n), complete_data[[i]]$g[1:ng], rep(NA, n-ng), rep(NA, m), rep(NA, m))
    Y4 = c(Y4, rep(NA, n), rep(NA, n), rep(NA, n), rep(0, m), rep(NA, m))
    Y5 = c(Y5, rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(0, m))
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
  idza_ = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), 1:m), N)
  wza_ = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), rep(-1, m)), N)
  
  
  
  
  # define residuals
  ide0 = rep(c(1:n, rep(NA, n), rep(NA, n), rep(NA, m), rep(NA, m)), N)
  ideh = rep(c(rep(NA, n), 1:n, rep(NA, n), rep(NA, m), rep(NA, m)), N)
  ideg = rep(c(rep(NA, n), rep(NA, n), 1:n, rep(NA, m), rep(NA, m)), N)
  idem = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), 1:m, rep(NA, m)), N)
  idea = rep(c(rep(NA, n), rep(NA, n), rep(NA, n), rep(NA, m), 1:m), N)
  
  
  # define inla data
  inla.data = data.frame("Y" = Y, "idh" = idh, "wh" = wh, "idg" = idg, "wg" = wg, 
                         "idz" = idz, "wz" = wz, "idz_" = idz_, "wz_" = wz_, 
                         "idza_" = idza_, "wza_", wza_,
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
    f(idz, wz, model = "rw1", constr = F, replicate = rep, hyper = list(prec = pc.prec.z)) +
    f(idh, wh, model = "rw1", constr = F, replicate = rep) +
    f(idg, wg, model = "rw1", constr = F, replicate = rep) +
    f(idz_, wz_, copy = "idz", hyper = list(beta = list(fixed = TRUE))) +
    f(idza_, wza_, copy = "idz", hyper = list(beta = list(fixed = FALSE, prior = "gaussian", param = c(0,1e2)))) #+
  #f(ide0, model = "iid", replicate = rep, hyper = list(prec = pc.prec.0)) +
  #f(ideh, model = "iid", replicate = rep, hyper = list(prec = pc.prec.h)) +
  #f(ideg, model = "iid", replicate = rep, hyper = list(prec = pc.prec.g)) +
  #f(idem, model = "iid", replicate = rep, hyper = list(prec = pc.prec.m)) +
  #f(idea, model = "iid", replicate = rep, hyper = list(prec = pc.prec.a)) 
  
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
r1.ipred = inla.model1.replicate.ipred(complete_data)

(ipredMSE1 = mean((complete_data[[5]]$g[145:288] - r1.ipred$summary.random$idg$mean[145:288 + 12*288])^2))
(ipredMSE2 = mean((complete_data[[7]]$g[145:288] - r1.ipred$summary.random$idg$mean[145:288 + 13*288])^2))
(ipredMSE3 = mean((complete_data[[11]]$g[145:288] - r1.ipred$summary.random$idg$mean[145:288 + 14*288])^2))
(ipredMSE4 = mean((complete_data[[16]]$g[145:288] - r1.ipred$summary.random$idg$mean[145:288 + 15*288])^2))
(mean.ipred.MSE = 0.25*(ipredMSE1 + ipredMSE2 + ipredMSE3 + ipredMSE4))


r1.ipred.sd = sqrt(r1.ipred$summary.random$idg$sd^2 +
                    (1/r1.ipred$summary.hyperpar$mean[2]))  



# plot predictions vs observations
r1.ipred.df = data.frame(x = complete_data[[5]]$t[145:288],
                        obs1 = complete_data[[5]]$g[145:288],
                        obs2 = complete_data[[7]]$g[145:288],
                        obs3 = complete_data[[11]]$g[145:288],
                        obs4 = complete_data[[16]]$g[145:288],
                        g1 = r1.ipred$summary.random$idg$mean[145:288 + 12*288],
                        g1.lq = r1.ipred$summary.random$idg$mean[145:288 + 12*288] + qnorm(0.025)*r1.ipred.sd[145:288 + 12*288],
                        g1.uq = r1.ipred$summary.random$idg$mean[145:288 + 12*288] + qnorm(0.975)*r1.ipred.sd[145:288 + 12*288],
                        g2 = r1.ipred$summary.random$idg$mean[145:288 + 13*288],
                        g2.lq = r1.ipred$summary.random$idg$mean[145:288 + 13*288] + qnorm(0.025)*r1.ipred.sd[145:288 + 13*288],
                        g2.uq = r1.ipred$summary.random$idg$mean[145:288 + 13*288] + qnorm(0.975)*r1.ipred.sd[145:288 + 13*288],
                        g3 = r1.ipred$summary.random$idg$mean[145:288 + 14*288],
                        g3.lq = r1.ipred$summary.random$idg$mean[145:288 + 14*288] + qnorm(0.025)*r1.ipred.sd[145:288 + 14*288],
                        g3.uq = r1.ipred$summary.random$idg$mean[145:288 + 14*288] + qnorm(0.975)*r1.ipred.sd[145:288 + 14*288],
                        g4 = r1.ipred$summary.random$idg$mean[145:288 + 15*288],
                        g4.lq = r1.ipred$summary.random$idg$mean[145:288 + 15*288] + qnorm(0.025)*r1.ipred.sd[145:288 + 15*288],
                        g4.uq = r1.ipred$summary.random$idg$mean[145:288 + 15*288] + qnorm(0.975)*r1.ipred.sd[145:288 + 15*288])

# plot predicted vs observed 
(mod1.ipred1 = ggplot(aes(x = x), data = r1.ipred.df) +
    geom_line(aes(y = g1, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs1, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g1.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g1.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g1.uq, ymin = g1.lq), fill = "salmon3", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "salmon3", "obs" = "black"),
                       labels = c("pred" = "Model 1 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw() + ggtitle("Test 1") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

(mod1.ipred2 = ggplot(aes(x = x), data = r1.ipred.df) +
    geom_line(aes(y = g2, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs2, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g2.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g2.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g2.uq, ymin = g2.lq), fill = "salmon3", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "salmon3", "obs" = "black"),
                       labels = c("pred" = "Model 1 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw() + ggtitle("Test 2") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

(mod1.ipred3 = ggplot(aes(x = x), data = r1.ipred.df) +
    geom_line(aes(y = g3, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs3, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g3.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g3.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g3.uq, ymin = g3.lq), fill = "salmon3", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "salmon3", "obs" = "black"),
                       labels = c("pred" = "Model 1 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw() + ggtitle("Test 3") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 


(mod1.ipred4 = ggplot(aes(x = x), data = r1.ipred.df) +
    geom_line(aes(y = g4, col = "pred"), lwd = 1.5) +
    geom_line(aes(y = obs4, col = "obs"), lwd = 1.5) +
    geom_line(aes(y = g4.uq, col = "pred"), lwd = 0.5) +
    geom_line(aes(y = g4.lq, col = "pred"), lwd = 0.5) +
    geom_ribbon(aes(ymax = g4.uq, ymin = g4.lq), fill = "salmon3", alpha = 0.2) +
    scale_color_manual(values = c("pred" = "salmon3", "obs" = "black"),
                       labels = c("pred" = "Model 1 predictions", "obs" = "Observed    "),
                       name = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(eta^z~(t)), limits = c(0,1)) +
    xlab("t") + theme_bw() + ggtitle("Test 4") +
    theme(axis.text = element_text(size = 20), 
          legend.text = element_text(size = 35),
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top")) 

# Combine by ggarrange
mod1.ipred.plot = ggarrange(mod1.ipred1 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank()),
                           NULL,
                           mod1.ipred2 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                           mod1.ipred3 + rremove("xlab") + rremove("ylab"), 
                           NULL,
                           mod1.ipred4 + rremove("xlab") + rremove("ylab") + theme(axis.text.y = element_blank()),
                           ncol = 3, nrow = 2,
                           widths = c(1,0.05,1,1,0.05,1), common.legend = TRUE)

annotate_figure(mod1.ipred.plot,
                left = text_grob(bquote(g(t)), size = 35, rot = 90), 
                bottom = text_grob("t", size = 35))

