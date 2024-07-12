library(tidyverse)
library(ggpubr)


se.ratio <- function(phi){
  r <- sqrt((1-phi)/(4-phi))
  return(r)
}

var.ratio <- function(phi){
  r <- (1-phi)/(4-phi)
  return(r)
}

pd.samp.se <- function(N_FE, K, phi){
  phi_v <- rep(phi, length(N_FE))
  N_PD <- sqrt(4-phi_v)*(N_FE-K-1) / (sqrt(1-phi_v)) + K + 1
  return(N_PD)
}

pd.samp.var <- function(N_FE, K, phi){
  phi_v <- rep(phi, length(N_FE))
  N_PD <- (4-phi_v)*(N_FE-K-1) / ((1-phi_v)) + K + 1
  return(N_PD)
}

phi <- seq(0,1, 0.01)

d <- tibble(
  phi = phi,
  ratio = se.ratio(phi))

N_fe <- seq(1, 5001, 100)

d2 <- tibble(
  N_FE = rep(N_fe, 4),
  Phi = as.factor(c(rep(0.01, length(N_fe)),
          rep(0.1, length(N_fe)),
          rep(0.3, length(N_fe)),
          rep(0.5, length(N_fe)))),
  N_PD = c(pd.samp.se(N_fe, 0, 0.01),
           pd.samp.se(N_fe, 0, 0.1),
           pd.samp.se(N_fe, 0, 0.3),
           pd.samp.se(N_fe, 0, 0.5))
)

p1 <- ggplot(d, aes(x = phi, y = ratio)) +
  geom_line(color='#008BBC') +
  labs(x = expression(phi),
       y = 'FE/PD Estimator SE Ratio',
       title = expression('Estimator Standard Error Ratio by'~phi)) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,0.5)) +
  theme_minimal()

p2 <- ggplot(d2, aes(x = N_PD, y=N_FE, color=Phi)) +
  geom_line() +
  geom_abline(aes(intercept=0, slope=0.5), lty=2, alpha=0.4) +
  coord_cartesian(xlim = c(0,5000), ylim = c(0, 2500)) +
  labs(x = 'N PD Sibling Pairs',
       y = 'N FE Sibling Pairs',
       title = 'FE Sample Size Required to Match PD Precision',
       color = expression(phi)) +
  scale_color_manual(values=c('0.01'='#37FF8B',
                              '0.1'='#008BBC',
                              '0.3'='#A06B9A',
                              '0.5'='#F9ADA0')) +
  theme_minimal()

comb1 <- ggarrange(p1, p2, ncol=2, labels = c('[A]', '[B]'))




rho.se.rat <- function(beta, se, rho, M){
  big_var <- se^2 + (beta^2*(1+rho)^2)/((M-3))
  se_rat <- sqrt(big_var/(se^2))
  return(se_rat)
}

M <- seq(500, 5000, length.out=100)

d3 <- tibble(M = M)
d3$BMI <- rho.se.rat(0.415, 0.075, 0.5, M)
d3$Height <- rho.se.rat(0.610, 0.091, 0.61, M)
d3$`Educational Attainment` <- rho.se.rat(0.336, 0.053, 0.52, M)
d3$`Self-Rated Health` <- rho.se.rat(0.103, 0.079, 0.54, M)
d3$`Depressive Symptoms` <- rho.se.rat(0.126, 0.072, 0.52, M)

p3 <- d3 |>
  pivot_longer(-M, names_to='PGS', values_to='SE Ratio') |>
  ggplot(aes(x=M, y=`SE Ratio`, color=PGS)) +
  geom_line() +
  geom_vline(aes(xintercept=2107), lty=2) +
  scale_color_manual(values=c('BMI'='#37FF8B',
                              'Depressive Symptoms'='#4A5043',
                              'Educational Attainment'='#008BBC',
                              'Height'='#A06B9A',
                              'Self-Rated Health'='#F9ADA0')) +
  labs(title=expression('Inflation of Standard Errors for PGSs'),
       x=expression('N'~rho~'Sibling Pairs'),
       y='SE Inflation Factor') +
  ylim(c(1, 1.1)) +
  theme_minimal()

d4 <- tibble(M = rep(M,2))
d4$z_1.96 <- c(rho.se.rat(1.96, 1, 0.5, M),
               rho.se.rat(1.96, 1, 0.65, M))
d4$z_2.58 <- c(rho.se.rat(2.58, 1, 0.5, M),
               rho.se.rat(2.58, 1, 0.65, M))
d4$z_3.29 <- c(rho.se.rat(3.29, 1, 0.5, M),
               rho.se.rat(3.29, 1, 0.65, M))
d4$rho <- rep(c(0.5, 0.65), each=length(M))

p4 <- d4 |>
  pivot_longer(c(-M,-rho), names_to='z-Score', values_to='Inflated SE', names_prefix = 'z_') |>
  ggplot(aes(x=M, y=`Inflated SE`, color=`z-Score`, lty=as.character(rho))) +
  geom_line() +
  #geom_vline(aes(xintercept=2107), lty=2) +
  scale_color_manual(values=c('1.96'='#008BBC',
                              '2.58'='#A06B9A',
                              '3.29'='#F9ADA0')) +
  labs(title=expression('Inflation of Standard Errors by t-Statistic'),
       x=expression('N'~rho~'Sibling Pairs'),
       y='SE Inflation Factor',
       color=expression('t-Statistic for'~hat(beta)^PD),
       lty=expression(rho)) +
  ylim(c(1, 1.1)) +
  theme_minimal()


comb2 <- ggarrange(p4, p3, ncol=2, labels = c('[A]', '[B]'))


M_max <- 5000

d5 <- data.frame(M= rep(100:M_max, 2),
                 rho = rep(c(0.5, 0.65), each=length(100:M_max))) |>
  mutate(var = (1-rho^2)^2/(M-3)) |>
  mutate(bias = (1-rho)^2/(var + (1-rho)^2))


p5 <- ggplot(d5, aes(x = M, y = bias, color=as.factor(rho))) +
  geom_line() +
  labs(x=expression('N'~rho~'Sibling Pairs'),
       y = 'Bias',
       title = expression('Bias in Estimates of'~beta^PD~'by'~phi~'and Sample Size'),
       color = expression(rho)) +
  scale_color_manual(values=c('0.5'='#008BBC',
                              '0.65'='#A06B9A')) +
  theme_minimal()

ggsave('figures/phi-plots.png', comb1, height=6, width=10)
ggsave('figures/rho-plots.png', comb2, height=6, width=10)
ggsave('figures/bias-plot.png', p5, height=6, width=10)



