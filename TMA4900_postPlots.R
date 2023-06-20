# Selected posterior plots from 4 considered models

source("TMA4900_Model1.R")
source("TMA4900_Model2.R")
source("TMA4900_Model3.R")
source("TMA4900_Model5.R")

# Plot alpha as a function of T

x = seq(8,14,by = 1)

plot.aT.df = data.frame("x" = x,
                        "y1" = rep(r1.pred$summary.hyperpar$mean[8], length(x)),
                        "y2" = r2.pred$summary.hyperpar$mean[8] + 
                          x*r2.pred$summary.hyperpar$mean[9],
                        "y4M0.5" = r5.pred$summary.hyperpar$mean[8] + 
                          x*r5.pred$summary.hyperpar$mean[9] +
                          0.5*r5.pred$summary.hyperpar$mean[10],
                        "y4M2.5" = r5.pred$summary.hyperpar$mean[8] + 
                          x*r5.pred$summary.hyperpar$mean[9] +
                          2.5*r5.pred$summary.hyperpar$mean[10],
                        "y4M4.5" = r5.pred$summary.hyperpar$mean[8] + 
                          x*r5.pred$summary.hyperpar$mean[9] +
                          4.5*r5.pred$summary.hyperpar$mean[10])

ggplot(plot.aT.df, aes(x = x)) +
  geom_line(aes(y = y1, col = "mod1"), lwd = 1.5) + 
  geom_line(aes(y = y2, col = "mod2"), lwd = 1.5) + 
  geom_line(aes(y = y4M0.5, col = "mod4.M0.5"), lwd = 1.5) +  
  geom_line(aes(y = y4M2.5, col = "mod4.M2.5"), lwd = 1.5) + 
  geom_line(aes(y = y4M4.5, col = "mod4.M4.5"), lwd = 1.5) +
  scale_colour_manual(values = c("mod1" = "salmon3", "mod2" = "darkblue", 
                                 "mod4.M0.5" = "darkseagreen",
                                 "mod4.M2.5" = "seagreen",
                                 "mod4.M4.5" = "aquamarine3"),
                      labels = c("mod1" = "Model 1", "mod2" = "Model 2", 
                                 "mod4.M0.5" = "Model 4, M = 0.5 kg",
                                 "mod4.M2.5" = "Model 4, M = 2.5 kg",
                                 "mod4.M4.5" = "Model 4, M = 4.5 kg"),
                      name = "") +
  scale_x_continuous("T [\u00B0C]") +
  scale_y_continuous(bquote(alpha(T))) +
  theme_bw() + theme(axis.text = element_text(size = 20), 
                     axis.title = element_text(size = 35), 
                     legend.text = element_text(size = 35))
  


# Plot alpha as a function of M

x = seq(0.1, 6, by = 0.1)

plot.aM.df = data.frame("x" = x,
                        "y1" = rep(r1.pred$summary.hyperpar$mean[8], length(x)),
                        "y3" = r3.pred$summary.hyperpar$mean[8] + 
                          x*r3.pred$summary.hyperpar$mean[9],
                        "y4T8" = r5.pred$summary.hyperpar$mean[8] + 
                          8*r5.pred$summary.hyperpar$mean[9] +
                          x*r5.pred$summary.hyperpar$mean[10],
                        "y4T10" = r5.pred$summary.hyperpar$mean[8] + 
                          10*r5.pred$summary.hyperpar$mean[9] +
                          x*r5.pred$summary.hyperpar$mean[10],
                        "y4T12" = r5.pred$summary.hyperpar$mean[8] + 
                          12*r5.pred$summary.hyperpar$mean[9] +
                          x*r5.pred$summary.hyperpar$mean[10])

ggplot(plot.aM.df, aes(x = x)) +
  geom_line(aes(y = y1, col = "mod1"), lwd = 1.5) + 
  geom_line(aes(y = y3, col = "mod3"), lwd = 1.5) + 
  geom_line(aes(y = y4T8, col = "mod4.T8"), lwd = 1.5) +  
  geom_line(aes(y = y4T10, col = "mod4.T10"), lwd = 1.5) + 
  geom_line(aes(y = y4T12, col = "mod4.T12"), lwd = 1.5) +
  scale_colour_manual(values = c("mod1" = "salmon3", "mod3" = "brown4", 
                                 "mod4.T8" = "darkseagreen",
                                 "mod4.T10" = "seagreen",
                                 "mod4.T12" = "aquamarine3"),
                      labels = c("mod1" = "Model 1", "mod3" = "Model 3", 
                                 "mod4.T8" = "Model 4, T = 8 \u00B0C",
                                 "mod4.T10" = "Model 4, T = 10 \u00B0C",
                                 "mod4.T12" = "Model 4, T = 12 \u00B0C"),
                      name = "") +
  scale_x_continuous("M [kg]") +
  scale_y_continuous(bquote(alpha(M))) +
  theme_bw() + theme(axis.text = element_text(size = 20), 
                     axis.title = element_text(size = 35), 
                     legend.text = element_text(size = 35))
