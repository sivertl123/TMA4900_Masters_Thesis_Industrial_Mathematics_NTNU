source("TMA4900_Data.R")

library(ggplot2)


## MONTHLY DATA - A21.2 as SAMPLE

# Only datapoints with daily M and n
A21.2 = all_data$A21.1
A21.2 = A21.2[which(!is.na(A21.2$M))]
A21.2$id = 1:length(A21.2$t)

head(A21.2)

# Plot M(t)
ggplot(data = A21.2, aes(x = t)) +
  geom_line(aes(x = id, y = M), lwd = 1) +
  scale_y_continuous("Average mass - M(t) - [g]") +
  xlab("Time (days after insert)") +
  #scale_x_datetime(date_labels = "%d") +
  theme_bw() + theme(text = element_text(size = 25)) 

# Plot n(t)
ggplot(data = A21.2, aes(x = t)) +
  geom_line(aes(x = id, y = n), lwd = 1) +
  scale_y_continuous("Number of salmon - n(t)") +
  xlab("Time (days after insert)") +
  #scale_x_datetime(date_labels = "%d") +
  theme_bw() + theme(text = element_text(size = 25)) 






## DAILY DATA - A11 as SAMPLE

df.A11 = complete_data$A11.221010

head(df.A11)



# Plot o2 added
ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = l), lwd = 1) + 
  scale_y_continuous(bquote(O[2]~"added, \u2113(t) [kg/min]")) +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M") +
  theme_bw() + theme(text = element_text(size = 25)) 

# Plot water flows
ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = q_f, colour = "f"), lwd = 1) + 
  geom_line(aes(y = q_r, colour = "r"), lwd = 1) + 
  scale_y_continuous("Water flows, q(t) [l/min]") +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M") +
  scale_colour_manual("", breaks = c("f", "r"),
                      values = c("f" = "black", "r" = "darkgrey"),
                      labels = c("f" = expression(q[f](t)), "r" = expression(q[r](t)))) + 
  theme_bw() +
  theme(text = element_text(size = 25), legend.position = "top") 

# Plot temperature 
ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = T), lwd = 1) +
  scale_y_continuous("Temperature, T(t) [\u00B0C]") +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M") + 
  theme_bw() + theme(text = element_text(size = 25)) 

# Plot O2 concentrations
ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = c_O2_center, colour = "c"), lwd = 1.5) + 
  geom_line(aes(y = c_O2_edge, colour = "e"), lwd = 1.5) + 
  scale_y_continuous(bquote(O[2]~"concentration, "~bold(c)[O[2]](t)~" [mg/l]")) +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M") +
  scale_colour_manual("", breaks = c("c", "e"),
                      values = c("c" = "black", "e" = "darkgrey"),
                      labels = c("c" = expression(c[paste(C, ",", O[2])](t)), 
                                 "e" = expression(c[paste(E, ",", O[2])](t)))) + 
  theme_bw() +
  theme(text = element_text(size = 25), legend.position = "top") 
  

# Plot carbon dioxide concentration

ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = c_CO2), lwd = 1) +
  scale_y_continuous(bquote(CO[2]~"concentration, "~c[CO[2]](t)~" [mg/l]")) +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M") + 
  theme_bw() + theme(text = element_text(size = 25)) 
  

# Plot feeding intensity

df.mplot = data.frame(t = seq(from=as.POSIXct("2020-01-01 00:00"),to=as.POSIXct("2020-01-01 23:55"),by="5 min"),
                      m = df.A11$m)

ggplot(data = df.A11, aes(x = t)) +
  geom_line(data = df.A11, aes(y = m), lwd = 1.5) +
  scale_y_continuous("Feeding intensity, m(t) [kg/5 min]") +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M") + 
  theme_bw() + theme(text = element_text(size = 25)) 
  
  
# plot carbon dioxide production

ggplot(data = df.A11, aes(x = t)) +
  geom_line(aes(y = g), lwd = 1) +
  scale_y_continuous(bquote(CO[2]~"produced - g(t) - [kg/min]")) +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M:%S") + 
  theme_bw() + theme(text = element_text(size = 25)) 
  

# plot oxygen consumption

ggplot(data = df.A11, aes(x = t)) +
  geom_line(aes(y = h), lwd = 1) +
  scale_y_continuous(bquote(O[2]~"consumed - h(t) - [kg/min]")) +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M:%S") + 
  theme_bw() + theme(text = element_text(size = 25)) 


# Cross-Autocorrelation h, g
ccf.hg = ccf(complete_data[[1]]$h, complete_data[[1]]$g,lag.max = 2000)
df.ccf.hg = data.frame("ccf1" = ccf.hg$acf, "lags" = ccf.hg$lag)
for(i in 2:17){
  ccf.hg = ccf(complete_data[[i]]$h, complete_data[[i]]$g,lag.max = 2000)
  df.ccf.hg[paste("ccf", i, sep = "")] = ccf.hg$acf
}

df.ccf.hg["mean"] = (1/17)*(df.ccf.hg$ccf1 + df.ccf.hg$ccf2 + df.ccf.hg$ccf3 + 
                            df.ccf.hg$ccf4 + df.ccf.hg$ccf5 + df.ccf.hg$ccf6 +
                            df.ccf.hg$ccf7 + df.ccf.hg$ccf8 + df.ccf.hg$ccf9 +
                            df.ccf.hg$ccf10 + df.ccf.hg$ccf11 + df.ccf.hg$ccf12 +
                            df.ccf.hg$ccf13 + df.ccf.hg$ccf14 + df.ccf.hg$ccf15 +
                            df.ccf.hg$ccf16 + df.ccf.hg$ccf17)

ggplot(data = df.ccf.hg, aes(x = lags)) +
  geom_line(aes(y = ccf1), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf2), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf3), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf4), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf5), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf6), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf7), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf8), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf9), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf10), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf11), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf12), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf13), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf14), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf15), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf16), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf17), lwd = 1, col = "grey") +
  geom_line(aes(y = mean), lwd = 2, col = "black") +
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = qnorm((1 + 0.95)/2)/sqrt(ccf.hg$n.used), 
             linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = - qnorm((1+0.95)/2)/sqrt(ccf.hg$n.used), 
             linetype = "dashed", lwd = 1) +
  scale_y_continuous("CCF(h, g)",
                     limits = c(-0.5, 1)) +
  xlab("Lag") + theme_bw() + theme(text = element_text(size = 35))



# Cross-Autocorrelation m, g

ccf.mg = ccf(df.A11$m, df.A11$g,lag.max = 2000)
df.ccf.mg = data.frame("ccf" = ccf.mg$acf, "lags" = ccf.mg$lag)

ggplot(data = df.ccf.mg, aes(x = lags)) +
  geom_line(aes(y = ccf), lwd = 1) +
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = qnorm((1 + 0.95)/2)/sqrt(ccf.mg$n.used), 
             linetype = "dashed") +
  geom_hline(yintercept = - qnorm((1+0.95)/2)/sqrt(ccf.mg$n.used), 
             linetype = "dashed") +
  scale_y_continuous("CACF(m, g)",
                     limits = c(-0.5, 0.5)) +
  xlab("Lag") + theme_bw() + theme(text = element_text(size = 25))


ccf.mg = ccf(complete_data[[1]]$m, complete_data[[1]]$g,lag.max = 2000)
df.ccf.mg = data.frame("ccf1" = ccf.mg$acf, "lags" = ccf.mg$lag)
for(i in 2:17){
  ccf.mg = ccf(complete_data[[i]]$m, complete_data[[i]]$g,lag.max = 2000)
  df.ccf.mg[paste("ccf", i, sep = "")] = ccf.mg$acf
}

df.ccf.mg["mean"] = (1/17)*(df.ccf.mg$ccf1 + df.ccf.mg$ccf2 + df.ccf.mg$ccf3 + 
                              df.ccf.mg$ccf4 + df.ccf.mg$ccf5 + df.ccf.mg$ccf6 +
                              df.ccf.mg$ccf7 + df.ccf.mg$ccf8 + df.ccf.mg$ccf9 +
                              df.ccf.mg$ccf10 + df.ccf.mg$ccf11 + df.ccf.mg$ccf12 +
                              df.ccf.mg$ccf13 + df.ccf.mg$ccf14 + df.ccf.mg$ccf15 +
                              df.ccf.mg$ccf16 + df.ccf.mg$ccf17)

ggplot(data = df.ccf.mg, aes(x = lags)) +
  geom_line(aes(y = ccf1), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf2), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf3), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf4), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf5), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf6), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf7), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf8), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf9), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf10), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf11), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf12), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf13), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf14), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf15), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf16), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf17), lwd = 1, col = "grey") +
  geom_line(aes(y = mean), lwd = 2, col = "black") +
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = qnorm((1 + 0.95)/2)/sqrt(ccf.mg$n.used), 
             linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = - qnorm((1+0.95)/2)/sqrt(ccf.mg$n.used), 
             linetype = "dashed", lwd = 1) +
  scale_y_continuous("CCF(m, g)",
                     limits = c(-0.7, 0.7)) +
  xlab("Lag") + theme_bw() + theme(text = element_text(size = 35))





# Cross-Autocorrelation m, h

ccf.mh = ccf(df.A11$feed, df.A11$h, lag.max = 1000)
df.ccf.mh = data.frame("ccf" = ccf.mh$acf, "lags" = ccf.mh$lag)

ggplot(data = df.ccf.mh, aes(x = lags)) +
  geom_line(aes(y = ccf), lwd = 1) +
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = qnorm((1 + 0.95)/2)/sqrt(ccf.mh$n.used), 
             linetype = "dashed") +
  geom_hline(yintercept = - qnorm((1+0.95)/2)/sqrt(ccf.mh$n.used), 
             linetype = "dashed") +
  scale_y_continuous("CACF(m, h)",
                     limits = c(-0.5, 0.5)) +
  xlab("Lag") + theme_bw() + theme(text = element_text(size = 25))

ccf.mh = ccf(complete_data[[1]]$m, complete_data[[1]]$h,lag.max = 2000)
df.ccf.mh = data.frame("ccf1" = ccf.mh$acf, "lags" = ccf.mh$lag)
for(i in 2:17){
  ccf.mh = ccf(complete_data[[i]]$m, complete_data[[i]]$h,lag.max = 2000)
  df.ccf.mh[paste("ccf", i, sep = "")] = ccf.mh$acf
}

df.ccf.mh["mean"] = (1/17)*(df.ccf.mh$ccf1 + df.ccf.mh$ccf2 + df.ccf.mh$ccf3 + 
                              df.ccf.mh$ccf4 + df.ccf.mh$ccf5 + df.ccf.mh$ccf6 +
                              df.ccf.mh$ccf7 + df.ccf.mh$ccf8 + df.ccf.mh$ccf9 +
                              df.ccf.mh$ccf10 + df.ccf.mh$ccf11 + df.ccf.mh$ccf12 +
                              df.ccf.mh$ccf13 + df.ccf.mh$ccf14 + df.ccf.mh$ccf15 +
                              df.ccf.mh$ccf16 + df.ccf.mh$ccf17)

ggplot(data = df.ccf.mh, aes(x = lags)) +
  geom_line(aes(y = ccf1), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf2), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf3), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf4), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf5), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf6), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf7), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf8), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf9), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf10), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf11), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf12), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf13), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf14), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf15), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf16), lwd = 1, col = "grey") +
  geom_line(aes(y = ccf17), lwd = 1, col = "grey") +
  geom_line(aes(y = mean), lwd = 2, col = "black") +
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = qnorm((1 + 0.95)/2)/sqrt(ccf.mh$n.used), 
             linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = - qnorm((1+0.95)/2)/sqrt(ccf.mh$n.used), 
             linetype = "dashed", lwd = 1) +
  scale_y_continuous("CCF(m, h)",
                     limits = c(-.7, .8)) +
  xlab("Lag") + theme_bw() + theme(text = element_text(size = 35))





A11.ccf = ccf(df.A11.std$feed, df.A11.std$g,lag.max = 20000)
A11.ccf.df = data.frame("ccf" = A11.ccf$acf, "lag" = A11.ccf$lag)
ggplot(data = A11.ccf.df, aes(x = lag)) + geom_line(aes(y = ccf))
A11.ccf.df$lag[which.max(A11.ccf.df$ccf)]
ggplot(data = df.A11.std, aes(x = t)) + geom_line(aes(y = feed), col = "red") + 
  geom_line(aes(y = g), col = "blue")



ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = h, colour = "h(t)"), lwd = 1.5) + 
  geom_line(aes(y = g_tot, colour = "g(t)"), lwd = 1.5) +
  geom_vline(xintercept = data$time[which.min(na.omit(data$h_tot))], col = "blue", lwd = 1.5) +
  geom_vline(xintercept = data$time[which.min(na.omit(data$g_tot))], col = "red", lwd = 1.5) +
  scale_colour_manual("", breaks = c("h(t)", "g(t)"),
                      values = c("h(t)" = "blue", "g(t)" = "red")) +
  scale_y_continuous(bquote(h[tot](t) ~ "and" ~ g[tot](t) ~ "[kg/min]"))  +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M:%S") +
  theme(text = element_text(size = 40))

ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = T), lwd = 1.5) +
  scale_y_continuous(bquote(T(t)~"[*C]"), limits = c(8, 13)) +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M:%S") +
  theme(text = element_text(size = 20))

ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = c_CO2), lwd = 1.5) + 
  scale_y_continuous(bquote(c[CO[2]](t)~"[mg/l]"), limits = c(6, 10)) +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M:%S") +
  theme(text = element_text(size = 20))

ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = g), lwd = 1.5) + 
  scale_y_continuous(bquote(g(t)~"[kg/5 min]"), limits = c(1.4,2.2)) +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M:%S") +
  theme(text = element_text(size = 20))

ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = feed), lwd = 1.5) + 
  geom_hline(aes(yintercept = mean(df.A11$feed)), lwd = 1.5) +
  scale_y_continuous(bquote(m(t)~"[kg/5 min]")) +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M") +
  theme(text = element_text(size = 20))

ggplot(data = df.A11, aes(x = t)) + 
  geom_line(aes(y = h), lwd = 1.5) + 
  scale_y_continuous(bquote(h(t)~"[kg/5 min]"), limits = c(1.4,2.2)) +
  xlab("Time") +
  scale_x_datetime(date_labels = "%H:%M:%S") +
  theme(text = element_text(size = 20))


# plot data used in validation

# Feeding data for test sets
test1 = complete_data[[5]]

(feed.test1 = ggplot(aes(x = t), data = test1) +
  geom_line(aes(y = m), lwd = 1.5, color = "black") +
  scale_x_datetime(date_labels = "%H:%M") +
  scale_y_continuous(bquote(m(t)), limits = c(0, 5)) +
  xlab("t") + theme_bw() + ggtitle("Test 1") +
  theme(axis.text = element_text(size = 30), 
        plot.title = element_text(size = 30, face = "bold"), 
        legend.position = "top"))

test2 = complete_data[[7]]

(feed.test2 = ggplot(aes(x = t), data = test2) +
    geom_line(aes(y = m), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(m(t)), limits = c(0, 5)) +
    xlab("t") + theme_bw() + ggtitle("Test 2") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

test3 = complete_data[[11]]

(feed.test3 = ggplot(aes(x = t), data = test3) +
    geom_line(aes(y = m), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(m(t)), limits = c(0, 5)) +
    xlab("t") + theme_bw() + ggtitle("Test 3") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

test4 = complete_data[[16]]

(feed.test4 = ggplot(aes(x = t), data = test4) +
    geom_line(aes(y = m), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(m(t)), limits = c(0, 5)) +
    xlab("t") + theme_bw() + ggtitle("Test 4") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

feed.test.plot = ggarrange(feed.test1 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank()),
                           NULL,
                           feed.test2 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                           feed.test3 + rremove("xlab") + rremove("ylab"), 
                           NULL,
                           feed.test4 + rremove("xlab") + rremove("ylab") + theme(axis.text.y = element_blank()),
                           ncol = 3, nrow = 2,
                           widths = c(1,0.05,1,1,0.05,1), common.legend = TRUE)

annotate_figure(feed.test.plot,
                left = text_grob(bquote(m(t)), size = 35, rot = 90), 
                bottom = text_grob("t", size = 35))

# oxygen data for test sets
test1 = complete_data[[5]]

(o2.test1 = ggplot(aes(x = t), data = test1) +
    geom_line(aes(y = h), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(h(t)), limits = c(0.26, 0.64)) +
    xlab("t") + theme_bw() + ggtitle("Test 1") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

test2 = complete_data[[7]]

(o2.test2 = ggplot(aes(x = t), data = test2) +
    geom_line(aes(y = h), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(h(t)), limits = c(0.26, 0.64)) +
    xlab("t") + theme_bw() + ggtitle("Test 2") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

test3 = complete_data[[11]]

(o2.test3 = ggplot(aes(x = t), data = test3) +
    geom_line(aes(y = h), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(h(t)), limits = c(0.26, 0.64)) +
    xlab("t") + theme_bw() + ggtitle("Test 3") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

test4 = complete_data[[16]]

(o2.test4 = ggplot(aes(x = t), data = test4) +
    geom_line(aes(y = h), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(h(t)), limits = c(0.26, 0.64)) +
    xlab("t") + theme_bw() + ggtitle("Test 4") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

o2.test.plot = ggarrange(o2.test1 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank()),
                           NULL,
                           o2.test2 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                           o2.test3 + rremove("xlab") + rremove("ylab"), 
                           NULL,
                           o2.test4 + rremove("xlab") + rremove("ylab") + theme(axis.text.y = element_blank()),
                           ncol = 3, nrow = 2,
                           widths = c(1,0.05,1,1,0.05,1), common.legend = TRUE)

annotate_figure(o2.test.plot,
                left = text_grob(bquote(h(t)), size = 35, rot = 90), 
                bottom = text_grob("t", size = 35))

# oxygen data for test sets
test1 = complete_data[[5]]

(temp.test1 = ggplot(aes(x = t), data = test1) +
    geom_line(aes(y = T), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(T(t)), limits = c(8.5, 14)) +
    xlab("t") + theme_bw() + ggtitle("Test 1") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

test2 = complete_data[[7]]

(temp.test2 = ggplot(aes(x = t), data = test2) +
    geom_line(aes(y = T), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(T(t)), limits = c(8.5, 14)) +
    xlab("t") + theme_bw() + ggtitle("Test 2") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

test3 = complete_data[[11]]

(temp.test3 = ggplot(aes(x = t), data = test3) +
    geom_line(aes(y = T), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(T(t)), limits = c(8.5,14)) +
    xlab("t") + theme_bw() + ggtitle("Test 3") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

test4 = complete_data[[16]]

(temp.test4 = ggplot(aes(x = t), data = test4) +
    geom_line(aes(y = T), lwd = 1.5, color = "black") +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_y_continuous(bquote(T(t)), limits = c(8.5, 14)) +
    xlab("t") + theme_bw() + ggtitle("Test 4") +
    theme(axis.text = element_text(size = 30), 
          plot.title = element_text(size = 30, face = "bold"), 
          legend.position = "top"))

temp.test.plot = ggarrange(temp.test1 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank()),
                         NULL,
                         temp.test2 + rremove("xlab") + rremove("ylab") + theme(axis.text.x = element_blank(), axis.text.y = element_blank()), 
                         temp.test3 + rremove("xlab") + rremove("ylab"), 
                         NULL,
                         temp.test4 + rremove("xlab") + rremove("ylab") + theme(axis.text.y = element_blank()),
                         ncol = 3, nrow = 2,
                         widths = c(1,0.05,1,1,0.05,1), common.legend = TRUE)

annotate_figure(temp.test.plot,
                left = text_grob(bquote(T(t)), size = 35, rot = 90), 
                bottom = text_grob("t", size = 35))

Tank = rep(NA, nrow(full_data))
for(i in 1:nrow(full_data)){
  if(full_data$G1[i] == 1){Tank[i] = "1"}
  if(full_data$G2[i] == 1){Tank[i] = "2"}
  if(full_data$G3[i] == 1){Tank[i] = "3"}
  if(full_data$G4[i] == 1){Tank[i] = "4"}
}

full_data$Tank = Tank

(dplot_temp = ggplot(data = full_data) +
  geom_density(aes(x = T, fill = Tank), color = "black", 
                 lwd = .5, alpha = 0.7) +
  scale_fill_manual(values = c("salmon4", "darkblue", "darkseagreen4", "deeppink4")) +
  ylab("Density") + xlab("T [\u00B0C]") +
  theme_bw() + theme(text = element_text(size = 35)))

(dplot_m = ggplot(data = full_data) +
  geom_density(aes(x = m, fill = Tank), color = "black", 
               lwd = .5, alpha = 0.7) +
  scale_fill_manual(values = c("salmon4", "darkblue", "darkseagreen4", "deeppink4")) +
  ylab("Density") + xlab("m [kg/min]") +
  theme_bw() + theme(text = element_text(size = 35)))
  
(dplot_h = ggplot(data = full_data) +
  geom_density(aes(x = h, fill = Tank), color = "black", 
               lwd = .5, alpha = 0.7) +
  scale_fill_manual(values = c("salmon4", "darkblue", "darkseagreen4", "deeppink4")) +
  scale_y_continuous("Density", limits = c(0,9)) +
  scale_x_continuous("h [kg/min]", limits = c(0,1)) +
  theme_bw() + theme(text = element_text(size = 33))
)

(dplot_g = ggplot(data = full_data) +
  geom_density(aes(x = g, fill = Tank), color = "black", 
               lwd = .5, alpha = 0.7) +
  scale_fill_manual(values = c("salmon4", "darkblue", "darkseagreen4", "deeppink4")) +
  scale_y_continuous("Density", limits = c(0, 9)) + 
  scale_x_continuous("g [kg/min]", limits = c(0,1)) +
  theme_bw() + theme(text = element_text(size = 33))
)

hplot_M = ggplot(data = full_data) +
  geom_histogram(aes(x = M, fill = Tank), color = "black", 
               lwd = .5, alpha = 0.7, position = "dodge", binwidth = 1) +
  scale_fill_manual(values = c("salmon4", "darkblue", "darkseagreen4", "deeppink4")) +
  scale_y_continuous("Count", limits = c(0, 1500)) + xlab("M [kg]") +
  theme_bw() + theme(text = element_text(size = 35))

hplot_d = ggplot(data = full_data) +
  geom_histogram(aes(x = n/4900, fill = Tank), color = "black", 
               lwd = .5, alpha = 0.7, binwidth = 10, position = "dodge") +
  scale_fill_manual(values = c("salmon4", "darkblue", "darkseagreen4", "deeppink4")) +
  scale_y_continuous("Count", limits = c(0, 1500)) + xlab(bquote("d ["~(m^3)^{-1}~"]")) +
  theme_bw() + theme(text = element_text(size = 35))


ggarrange(dplot_h, dplot_g + rremove("ylab"), ncol = 2, align = "h", common.legend = T)


d = rep(NA, 17)
M = rep(NA, 17)
Tank = c(rep("1", 5), rep("2", 4), rep("3", 3), rep("4", 4)) 

for(i in 1:17){
  d[i] = mean(complete_data[[i]]$n)/4900
  M[i] = mean(complete_data[[i]]$M)
}
d = d[-c(13)]
M = M[-c(13)]

hplot_df = data.frame("d" = d, "M" = M, "Tank" = Tank)

hplot_df

(hplot_M = ggplot(data = hplot_df) +
  geom_histogram(aes(x = M, fill = Tank), color = "black", 
                 lwd = .5, alpha = 0.7, position = "dodge", binwidth = 1) +
  scale_fill_manual(values = c("salmon4", "darkblue", "darkseagreen4", "deeppink4")) +
  scale_y_continuous("Count", limits = c(0,5)) + 
  scale_x_continuous("M [kg]", breaks = c(0,1,2,3,4,5)) +
  theme_bw() + theme(text = element_text(size = 35)))

(hplot_d = ggplot(data = hplot_df) +
  geom_histogram(aes(x = d, fill = Tank), color = "black", 
                 lwd = .5, alpha = 0.7, binwidth = 10, position = "dodge") +
  scale_fill_manual(values = c("salmon4", "darkblue", "darkseagreen4", "deeppink4")) +
  scale_y_continuous("Count", limits = c(0,5)) + xlab(bquote("d ["~1~"/"~m^3~"]")) +
  theme_bw() + theme(text = element_text(size = 35))
)

ggarrange(hplot_M, hplot_d + rremove("ylab"), common.legend = T, ncol = 2, align = "h")

sd(full_data$n/4900)


# plot training set observations

# plot h, test 2
(htest2 = ggplot(complete_data[[2]], aes(x = t, y = h)) +
  geom_line(lwd = 1.5) + 
  scale_y_continuous(bquote(h(t)~" [kg/5 min]"), limits = c(-0.1,1)) +
  scale_x_datetime(date_labels = "%H:%M") + ggtitle("Train 2") +
  theme_bw() + theme(axis.text = element_text(size = 20),
                     axis.title = element_text(size = 35),
                     plot.title = element_text(size = 30, face = "bold")))


# plot g, test 2
(gtest2 = ggplot(complete_data[[2]], aes(x = t, y = g)) +
  geom_line(lwd = 1.5) + 
  scale_y_continuous(bquote(g(t)~" [kg/5 min]"), limits = c(0.3,.7)) +
  scale_x_datetime(date_labels = "%H:%M") +
  theme_bw() + theme(axis.text = element_text(size = 20),
                     axis.title = element_text(size = 35)))

# plot h, test 8
(htest8 = ggplot(complete_data[[10]], aes(x = t, y = h)) +
  geom_line(lwd = 1.5) + 
  scale_y_continuous(bquote(h(t)~" [kg/5 min]"), limits = c(-0.1,1)) +
  scale_x_datetime(date_labels = "%H:%M") + ggtitle("Train 8") +
  theme_bw() + theme(axis.text = element_text(size = 20),
                     axis.title = element_text(size = 35),
                     plot.title = element_text(size = 30, face = "bold")))
# plot g, test 8
(gtest8 = ggplot(complete_data[[10]], aes(x = t, y = g)) +
  geom_line(lwd = 1.5) + 
  scale_y_continuous(bquote(g(t)~" [kg/5 min]"), limits = c(0.3,0.7)) +
  scale_x_datetime(date_labels = "%H:%M") +
  theme_bw() + theme(axis.text = element_text(size = 20),
                     axis.title = element_text(size = 35)))

ggarrange(htest2 + rremove("xlab"), htest8 + rremove("ylab") + rremove("xlab"),
          gtest2, gtest8 + rremove("ylab"), 
          nrow = 2, ncol = 2, align = "hv")


?inla.pc.dprec
library(INLA)

# plot PC(u = 0.01, a = 0.5)
x = seq(1, 10000, by = 1)
y = inla.pc.dprec(x, u = 0.01, alpha = 0.5)
ymid <- data.frame(inla.tmarginal(function(x) 1/x, data.frame(x, y)))
nx = c(ymid$x, 0.02)
ny = c(ymid$y, 0)
pc1 = data.frame("x" = nx, "y" = ny)

pcplot1 = ggplot(data = data.frame(pc1), aes(x = x, y = y)) +
  geom_line(lwd = 1.5) + 
  scale_y_continuous(bquote(pi(sigma^2))) +
  scale_x_continuous(bquote(sigma^2), breaks = c(0.00, 0.01, 0.02)) +
  ggtitle(bquote(P(sigma > 0.01)~" = 0.5")) + theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 35),
        plot.title = element_text(size = 35))

# plot PC(u = 0.01, a = 0.001)
x = seq(1, 1000000, by = 1)
y = inla.pc.dprec(x, u = 0.01, alpha = 0.001)
ymid <- data.frame(inla.tmarginal(function(x) 1/x, data.frame(x, y)))
nx = c(ymid$x, 0.02)
ny = c(ymid$y, 0)
pc2 = data.frame("x" = nx, "y" = ny)

pcplot2 = ggplot(data = data.frame(pc2), aes(x = x, y = y)) +
  geom_line(lwd = 1.5) + 
  scale_y_continuous(bquote(pi(sigma^2))) +
  scale_x_continuous(bquote(sigma^2), breaks = c(0.00, 0.01, 0.02)) +
  ggtitle(bquote(P(sigma > 0.01)~" = 0.001")) + theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 35),
        plot.title = element_text(size = 35))

ggarrange(pcplot1, pcplot2, align = "h")


# plot PC(u = 0.5, a = 0.5)
x = seq(1, 1000, by = 1)
y = inla.pc.dprec(x, u = 0.5, alpha = 0.5)
ymid <- data.frame(inla.tmarginal(function(x) 1/x, data.frame(x, y)))
nx = c(ymid$x)
ny = c(ymid$y)
pc1 = data.frame("x" = nx, "y" = ny)

(pcplot1 = ggplot(data = data.frame(pc1), aes(x = x, y = y)) +
  geom_line(lwd = 1.5) + 
  scale_y_continuous(bquote(pi(sigma^2))) +
  scale_x_continuous(bquote(sigma^2)) +
  ggtitle(bquote(P(sigma > 0.5)~" = 0.5")) + theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 35),
        plot.title = element_text(size = 35))
)
# plot PC(u = 0.001, a = 0.0001)
x = seq(1, 10000000, by = 1)
y = inla.pc.dprec(x, u = 0.001, alpha = 0.0001)
ymid <- data.frame(inla.tmarginal(function(x) 1/x, data.frame(x, y)))
nx = c(ymid$x, 0.02)
ny = c(ymid$y, 0)
pc2 = data.frame("x" = nx, "y" = ny)

(pcplot2 = ggplot(data = data.frame(pc2), aes(x = x, y = y)) +
  geom_line(lwd = 1.5) + 
  scale_y_continuous(bquote(pi(sigma^2))) +
  scale_x_continuous(bquote(sigma^2), breaks = c(0.00, 0.01, 0.02)) +
  ggtitle(bquote(P(sigma > 0.001)~" = 0.0001")) + theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 35),
        plot.title = element_text(size = 35))
)
ggarrange(pcplot1, pcplot2, align = "h")




