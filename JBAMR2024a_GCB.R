########################################################################
# PAPER 
# Journal: GCB Bioenergy
# Title:# Did the entry of the corn ethanol industry in Brazil 
        # affect the relationship between domestic and 
        # international corn prices?
# Authors: M. Justus, L. Bachioni, S. Arantes, M. Moreira, L. Rodrigues 
# Updated on July 1, 2024
########################################################################

################################
##### PART I: Modelagem VAR-VEC
###############################

rm(list=ls()) 

# Carregando os pacotes usados na parte 1

#install.packages("tidyverse")
#install.packages("readxl")
#install.packages("forecast")
#install.packages("vars")
#install.packages("urca")
#install.packages("lmtest")
#install.packages("ggplot2")
#install.packages("gridExtra")

library(tidyverse)
library(readxl)
library(forecast)
library(vars)
library(fUnitRoots)
library(urca)
library(lmtest)
library(ggplot2)
library(gridExtra)

# Definindo o diretório

#getwd()

setwd()

# Lendo os dados

data <- read_excel("JBAMR2024a_GCB.xlsx",
                    sheet = "data")

# Conhecendo os dados

class(data)
names(data)
head(data)
dim(data)

# Transformando em TS

PIMI <- ts((data$PIMI), start=c(2005,6), frequency=12)
PMIMT <- ts((data$PMIMT), start=c(2005,6), frequency=12)
D1 <- ts((data$D1), start=c(2005,6), frequency=12)
D2 <- ts((data$D2), start=c(2005,6), frequency=12)
D3 <- ts((data$D3), start=c(2005,6), frequency=12)
D4 <- ts((data$D4), start=c(2005,6), frequency=12)
D5 <- ts((data$D5), start=c(2005,6), frequency=12)
D6 <- ts((data$D6), start=c(2005,6), frequency=12)

###
### Análises de correlação
###

cor.test(log(PIMI), log(PMIMT), method = "pearson", 
         alternative = "two.sided", conf.level = 0.95)

###
### Olhando para a distribuição dos dados
###

# Densidades de probabilidade

par(mfrow=c(1,2))

plot(density(log(PIMI), kernel=c("gaussian")), 
     main="Distribuição do log(PIMI)",
     col=2)

plot(density(log(PMIMT), kernel=c("gaussian")), 
     main="Distribuição do log(PMIMT)",
     col=2)

### 
### Checando a estacionariedade das séries
###

## Aplicando análises gráficas e testes de RU

par(mfrow=c(1,2))

# PIMI

Acf(log(PIMI), main="log(PIMI)")
teste1_PIMI <- ur.df(log(PIMI), type="trend", selectlags="BIC")
summary(teste1_PIMI)

teste2_PIMI <- ur.df(log(PIMI), type="drift", selectlags="BIC")
summary(teste2_PIMI)

teste3_PIMI <- ur.df(log(PIMI), type="none", selectlags="BIC")
summary(teste3_PIMI)

ndiffs(log(PIMI),alpha=0.05,test=c("adf"),type=c("trend"),max.d=2)
ndiffs(log(PIMI),alpha=0.05,test=c("pp"),type=c("trend"),max.d=2)
ndiffs(log(PIMI),alpha=0.05,test=c("kpss"),type=c("trend"),max.d=2)

# PMIMT

Acf(log(PMIMT), main="log(PMIMT)")

teste1_PMIMT <- ur.df(log(PMIMT), type="trend", selectlags="BIC")
summary(teste1_PMIMT)

teste2_PMIMT <- ur.df(log(PMIMT), type="drift", selectlags="BIC")
summary(teste2_PMIMT)

teste3_PMIMT <- ur.df(log(PMIMT), type="none", selectlags="BIC")
summary(teste3_PMIMT)

ndiffs(log(PMIMT),alpha=0.05,test=c("adf"),type=c("trend"),max.d=2)
ndiffs(log(PMIMT),alpha=0.05,test=c("pp"),type=c("trend"),max.d=2)
ndiffs(log(PMIMT),alpha=0.05,test=c("kpss"),type=c("trend"),max.d=2)# PMIPR

###
### Modelagem VAR
###

## Construíndo os vetores 

YMT <- cbind(log(PMIMT), log(PIMI))

D <- cbind(D1,D2,D3,D4,D5,D6)

head(YMT)
head(D)
tail(D)

## Selecionado o VAR

VARselect(YMT, lag.max=4, type="cons", season=12, 
          exogen=D) 

VARselect(YMT, lag.max=4, type="both", season=12, 
          exogen=D) 

VARselect(YMT, lag.max=4, type="trend", season=12, 
          exogen=D) 

## Estimando o VAR(1) para Delta_1

summary(reg1 <- VAR(diff(YMT), p=1, type="cons", season=12, exogen=diff(D)))

## Diagnósticos nos resíduos

par(mfrow=c(1,2))

plot(density(reg1$varresult$log.PMIMT.$residuals, kernel=c("gaussian")), 
     main="Resíduos da equação em primeira diferença do log(PMIMT)",
     col=2)

plot(density(reg1$varresult$log.PIMI.$residuals, kernel=c("gaussian")), 
     main="Resíduos da equação em primeira diferença do log(PIMI)",
     col=2)

serial.test(reg1,lags.pt=6, type = "PT.adjusted")

arch.test(reg1, lags.single=6, lags.multi=6)

normality.test(reg1)

## Análise de causalidade de granger

(causality(reg1,cause="log.PIMI.", boot=TRUE, boot.runs=1000))
# PIMI -> PMIMT

(causality(reg1,cause = "log.PMIMT.", boot=TRUE, boot.runs=1000))
# PMIMT -> PIMI 

## Funções de impulso resposta 

# PIMI -> PMIMT

plot(irf(reg1, impulse="log.PIMI.",
          response="log.PMIMT.", n.ahead=5,
          ortho=T, cumulative=F, boot=TRUE,
          ci=0.95, runs=1000, seed=NULL))

irf1 <- irf(reg1, impulse="log.PIMI.",
         response="log.PMIMT.", n.ahead=5,
         ortho=T, cumulative=F, boot=TRUE,
         ci=0.95, runs=1000, seed=NULL)

print(irf1)

##
## Figura 3A
##

# Organizando os dados

impulse_response_irf1 <- data.frame(
  Time = 1:6,
  Response = c(0.0000000000, 0.0239023603, 0.0103720342, 0.0039708139, 0.0014931022, 0.0005598668),
  Lower = c(0.000000e+00, 1.138607e-02, 4.345042e-03, 1.224680e-03, 2.707097e-04, 4.268667e-05),
  Upper = c(0.000000000, 0.033450649, 0.016336900, 0.007732496, 0.003719934, 0.001784125)
)

# Plotando 

fig3a <- ggplot(impulse_response_irf1, aes(x = Time)) +
  geom_line(aes(y = Response), size = 1, show.legend = FALSE, color = "blue") + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "95% Bootstrap CI, 1000 runs"), alpha = 0.2) + 
  geom_hline(yintercept = 0, color = "red", linetype = "solid") + 
  labs(
    title = "Orthogonal Impulse Responses from International Price of Corn (log)",
    x = "Months", 
    y = "Domestic Price of Corn (log)",
    fill = "" 
  ) +
  scale_fill_manual(values = c("95% Bootstrap CI, 1000 runs" = "blue")) + 
  scale_x_continuous(breaks = 1:6) + 
  theme_minimal() +
  theme(
    legend.position = c(0.8, 0.8), 
    axis.title = element_text(size = 12), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12), 
    plot.title = element_text(size = 12), 
    panel.background = element_rect(fill = "gray95", color = NA), 
    panel.grid.major = element_line(size = 0.5, color = "white"), 
    panel.grid.minor = element_line(size = 0.5, color = "white") 
  )

# PMIMT -> PIMI 

plot(irf(reg1, impulse="log.PMIMT.",
         response="log.PIMI.", n.ahead=5,
         ortho=T, cumulative=F, boot=TRUE,
         ci=0.95, runs=1000, seed=NULL))

irf2 <- irf(reg1, impulse="log.PMIMT.",
         response="log.PIMI.", n.ahead=5,
         ortho=T, cumulative=F, boot=TRUE,
         ci=0.95, runs=1000, seed=NULL)

print(irf2)

##
## Figura 3B
##

impulse_response_irf2 <- data.frame(
  Time = 1:6,
  Response = c(0.0242492487, 0.0096793244, 0.0036625328, 0.0013746873, 0.0005153178, 0.0001931343),
  Lower = c(1.362622e-02, 1.547892e-03, -5.174481e-04, -3.631046e-04, -1.034340e-04, -2.884493e-05),
  Upper = c(0.0337003081, 0.0165107652, 0.0077647707, 0.0037255741, 0.0018534616, 0.0009085941)
)

fig3b <- ggplot(impulse_response_irf2, aes(x = Time)) +
  geom_line(aes(y = Response), size = 1, show.legend = FALSE, color = "blue") + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "95% Bootstrap CI, 1000 runs"), alpha = 0.2) + 
  geom_hline(yintercept = 0, color = "red", linetype = "solid") + 
  labs(
    title = "Orthogonal Impulse Responses from Domestic Price of Corn (log)",
    x = "Months", 
    y = "International Price of Corn (log)",
    fill = "" 
  ) +
  scale_fill_manual(values = c("95% Bootstrap CI, 1000 runs" = "blue")) + 
  scale_x_continuous(breaks = 1:6) + 
  theme_minimal() +
  theme(
    legend.position = c(0.8, 0.8), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 12), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12), 
    plot.title = element_text(size = 12), 
    panel.background = element_rect(fill = "gray95", color = NA), 
    panel.grid.major = element_line(size = 0.5, color = "white"), 
    panel.grid.minor = element_line(size = 0.5, color = "white") 
  )

### 
### Figura 3 do texto
###

# Configurando para salvar a imagem com 300 dpi

png("figure3_300dpi.png", width = 10 * 300, height = 10 * 300, res = 300)

grid.arrange(fig3a, fig3b, nrow = 2, ncol = 1)

dev.off()

###
### Análise de cointegração Johansen
###

## Testando a existência de cointegração

# Teste do traço: r=1

summary(joh1A <- ca.jo(YMT, type="trace", ecdet="none", K=2,
                 spec="transitory",
                 season=12, dumvar=D))

summary(joh2A <- ca.jo(YMT, type="trace", ecdet="cons", K=2,
                       spec="transitory",
                       season=12, dumvar=D))

summary(joh3A <- ca.jo(YMT, type="trace", ecdet="trend", K=2,
                       spec="transitory", 
                       season=12, dumvar=D))

# Testes do máximo autovalor

summary(joh4A <- ca.jo(YMT, type="eigen", ecdet="none", K=2,
                       spec="transitory",
                       season=12, dumvar=D))

summary(joh5A <- ca.jo(YMT, type="eigen", ecdet="cons", K=2,
                       spec="transitory",
                       season=12, dumvar=D))

summary(joh6A <- ca.jo(YMT, type="eigen",
                       ecdet="trend", K=2,
                       spec="transitory", season=12, dumvar=D))

## Estimando um VECM(1) 

vecm_MT <- ca.jo(YMT, type="trace", 
                 ecdet=c("trend"), K=2,
                 spec=c("transitory"),
                 season=12, dumvar=D)
summary(vecm_MT)

## Transformando o VECM em VAR em níveis

vec2_MT <- vec2var(vecm_MT, r=1)
list(vec2_MT)

plot(irf(vec2_MT, impulse="log.PIMI.",
         response="log.PMIMT.", n.ahead=11,
         ortho=T, cumulative=F, boot=TRUE,
         ci=0.95, runs=1000, seed=NULL))

plot(irf(vec2_MT, impulse="log.PMIMT.",
         response="log.PIMI.", n.ahead=11,
         ortho=T, cumulative=F, boot=TRUE,
         ci=0.95, runs=1000, seed=NULL))

# PIMI -> PMIMT

## Checando os resíduos

# Testes de autocorrelação univariados

res_PMIMT <- vec2_MT$resid[,1]

res_PIMI <- vec2_MT$resid[,2]

Box.test(res_PMIMT,lag=8, fitdf=2, type="Ljung-Box")
Box.test(res_PMIMT,lag=6, fitdf=2, type="Ljung-Box")
Box.test(res_PMIMT,lag=4, fitdf=2, type="Ljung-Box")

Box.test(res_PIMI,lag=8, fitdf=2, type="Ljung-Box")
Box.test(res_PIMI,lag=6, fitdf=2, type="Ljung-Box")
Box.test(res_PIMI,lag=4, fitdf=2, type="Ljung-Box")

# Testes de autocorrelação mulivariados

serial.test(vec2_MT, lags.pt=6)
serial.test(vec2_MT, lags.pt=4)

# Testes de heteroced. ARCH-LM 

arch.test(vec2_MT, lags.single=6, lags.multi=6,
          multivariate.only=F)

arch.test(vec2_MT, lags.single=4, lags.multi=4, multivariate.only=F)

## Testes de normalidade   

norm <- normality.test(vec2_MT, multivariate.only=F)

norm$jb.uni

norm$jb.mul

plot(density(res_PMIMT, kernel=c("gaussian")), 
     main="Resíduos da equação log(PMIMT) no VEC(1))",
     col=2)

plot(density(res_PIMI, kernel=c("gaussian")), 
     main="Resíduos da equação log(PIMi) no VEC(1))",
     col=2)

plot(density(vec2_MT$resid, kernel=c("gaussian")), 
     main="Resíduos multivariado do VEC(1))",
     col=2)

## Análise nas relações de coint. (beta) e
## Coef. de ajustamento (alpha)

summary(vecm_MT)

# Betas

vecm.r1_MT <- cajorls(vecm_MT, r=1)

vecm.r1_MT 

names(vecm.r1_MT)

summary(vecm.r1_MT$rlm)

vecm.r1_MT$beta

# Alphas

alpha <- coef(vecm.r1_MT$rlm)[1, ]

names(alpha) <- c("PMIMT", "PIMI")

alpha

resids <- resid(vecm.r1_MT$rlm)

N <- nrow(resids)

sigma <- crossprod(resids) / N

beta <- vecm.r1_MT$beta; beta

alpha.se <- sqrt(solve(crossprod(
  cbind(vecm_MT@ZK %*% beta, vecm_MT@Z1)))[1, 1] * diag(sigma))
alpha.se
names(alpha.se) <-  c("PMIMT", "PIMI")
(alpha.t <- alpha/alpha.se)

# Mais sobre os Betas

beta <- vecm.r1_MT$beta

beta

beta.se <- sqrt(diag(kronecker(
  solve(crossprod(vecm_MT@RK[, -1])),
  solve(t(alpha) %*% solve(sigma) %*% alpha))))

beta.t <- c(NA, beta[-1]/beta.se)

names(beta.t) <- rownames(vecm.r1_MT$beta)

beta.se

beta.t

## Testes de restrições

# H0:TREND=0 

HD0 <- matrix(c(1,0,0, 0,0,1), c(3,2))
HD0
summary(blrtest(joh3A, H=HD0, r=1))

# H0:PIMI=0 

HD0 <- matrix(c(1,0,0, 0,1,0), c(3,2))
HD0
summary(blrtest(joh3A, H=HD0, r=1))

# Alpha1=0 (na equação do PMIMT)

DA1 <- matrix(0, nrow = 2, ncol=1)
DA1[2,1] <-1
DA1
summary(alrtest(joh3A, A=DA1, r=1))

##############################
##### PART II: GRÁFICOS EXTRAS
##############################

##
## Figura 1 do texto
##

# Carregando pacotes 

#install.packages("ggplot2")
#install.packages("ggplot2")
#install.packages("lubridate")
#install.packages("scales")

library(ggplot2)
library(lubridate)
library(scales) 
library(readxl)

# Lendo os dados

#setwd( )

data_fig1 <- read_excel("JBAMR2024a_GCB.xlsx",
                           sheet = "figure1")

# Convertendo ano e mês em data

data_fig1$date <- make_date(data_fig1$year, data_fig1$month, 1)

# Configurando para salvar a imagem com 300 dpi

png("figure1_300dpi.png", width = 10 * 300, height = 5 * 300, res = 300)

# Definindo os tamanhos de fonte

tamanho_titulo <- 12
tamanho_eixo <- 10
tamanho_legenda <- 12
tamanho_linha <- 1  

ggplot(data_fig1, aes(x = date, y = price, color = Market)) +
  geom_line(linewidth = 0.5) +
  geom_smooth(aes(group = Market, linetype = "Smoothed by loess"), size=0.5, se = FALSE, method = "loess", color = "black") + 
  scale_color_manual(values = c("International" = "blue", "Mato Grosso" = "red")) +
  scale_linetype_manual(values = c("Smoothed by loess" = "dashed"), name = "") + 
  labs(title = '',
       x = 'Time', y = 'US$/mt', color = 'Market',
       caption = "") + 
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = tamanho_titulo, hjust = 0.5), 
        axis.text = element_text(size = tamanho_eixo), 
        axis.title = element_text(size = tamanho_eixo), 
        legend.text = element_text(size = tamanho_legenda), 
        plot.caption = element_text(size = tamanho_eixo, hjust = 0),
        panel.background = element_rect(fill = "gray95", color = NA), 
        panel.grid.major = element_line(size = 1, color = "white"), 
        panel.grid.minor = element_line(size = 1, color = "white")) + 
  scale_x_date(limits = as.Date(c("2005-01-01", "2023-12-31")), 
               date_breaks = "1 year", 
               date_labels = "%Y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_vline(xintercept = as.numeric(as.Date("2007-11-01")), linetype = "dashed", size = 0.6, color = "80grey") +
  annotate("text", x = as.Date("2007-11-01") + 15, y = max(data_fig1$price), label = "D1", size = 4, vjust = 1, hjust = 0, color = "black") +
  geom_vline(xintercept = as.numeric(as.Date("2008-09-01")), linetype = "dashed", size = 0.6, color = "80grey") +
  annotate("text", x = as.Date("2008-09-01") + 15, y = max(data_fig1$price), label = "D2", size = 4, vjust = 1, hjust = 0, color = "black") +
  geom_vline(xintercept = as.numeric(as.Date("2011-02-01")), linetype = "dashed", size = 0.6, color = "80grey") +
  annotate("text", x = as.Date("2011-02-01") + 15, y = max(data_fig1$price), label = "D3", size = 4, vjust = 1, hjust = 0, color = "black") +
  geom_vline(xintercept = as.numeric(as.Date("2016-04-01")), linetype = "dashed", size = 0.6, color = "80grey") +
  annotate("text", x = as.Date("2016-04-01") + 15, y = max(data_fig1$price), label = "D4", size = 4, vjust = 1, hjust = 0, color = "black") +
  geom_vline(xintercept = as.numeric(as.Date("2020-07-01")), linetype = "dashed", size = 0.6, color = "80grey") +
  annotate("text", x = as.Date("2020-07-01") + 15, y = max(data_fig1$price), label = "D6", size = 4, vjust = 1, hjust = 0, color = "black") +
  geom_vline(xintercept = as.numeric(as.Date("2017-05-01")), linetype = "dashed", size = 0.6, color = "80grey") +
  annotate("text", x = as.Date("2017-05-01") + 15, y = max(data_fig1$price), label = "D5", size = 4, vjust = 1, hjust = 0, color = "black")

dev.off()

##
## Figura 2
##

# Carregando pacotes 

#install.packages("tsibble")
#install.packages("dplyr")
#install.packages("gridExtra")
#install.packages("tidyr")
#install.packages("readxl")

library(tsibble) 
library(dplyr) 
library(gridExtra)
library(tidyr)
library(readxl) 

# Lendo os dados

# setwd()

data_fig2 <- read_excel("JBAMR2024a_GCB.xlsx",
                    sheet = "figure2")

prices <- data_fig2  %>%
  mutate(Month = as.Date(Tempo))

names(prices)

prices <- prices %>% 
  mutate(Month = yearmonth(Month)) 

prices <- prices %>% 
  as_tsibble(index  = Month)

print(prices)

prices$log_International <- log(prices$PIMI)
prices$log_MT <- log(prices$PMIMT)

prices |>
  ggplot(aes(x = log_International, y = log_MT)) +
  labs(title = "Domestic Corn Prices (MT, Brazil) vs International Corn rices and the Fitted Regression Line",
       subtitle = "Logarithm of US dollars per ton, 2005:M6 to 2023:M8",
       y = "Domestic Price (log)", x = "Internacional Price (log)") +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = T)

## Separando os períodos

# Definindo o ponto de corte

data_corte <- as.Date("2017-05-01")

# Dividindo os dados em dois subconjuntos

dados_antes <- subset(prices, Tempo < data_corte)

dados_depois <- subset(prices, Tempo >= data_corte)

# OLS para cada conjunto de dados

modelo_total  <- lm(log_MT ~ log_International, data = prices)

modelo_antes  <- lm(log_MT ~ log_International, data = dados_antes)

modelo_depois <- lm(log_MT ~ log_International, data = dados_depois)

# Extraindo R^2

(r2_total <- summary(modelo_total)$r.squared)
sqrt(0.6665243)

(r2_antes <- summary(modelo_antes)$r.squared)
sqrt(0.61970)

(r2_depois <- summary(modelo_depois)$r.squared)
sqrt(0.7642642)

# Plotando e adicionando R^2

g1 <- ggplot(prices, aes(x = log_International, y = log_MT)) + labs(title = "Full sample (2005:M6 to 2023:M8)",
       y = "Domestic Price (log)", x = "Internacional Price (log)") + geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, 
       color = "blue") + annotate("text", x = Inf, y = Inf, label = sprintf("R-squared = %.2f", r2_total), hjust = 6, vjust = 2, size = 5)

g2 <- ggplot(dados_antes, aes(x = log_International, y = log_MT)) + labs(title = "First period (2005:M6 to 2017:M4)",
       y = "Domestic Price (log)", x = "Internacional Price (log)") + geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, 
       color = "blue") + annotate("text", x = Inf, y = Inf,label = sprintf("R-squared = %.2f", r2_antes), hjust = 2.5, vjust = 2, size = 5)

g3 <- ggplot(dados_depois, aes(x = log_International, y = log_MT)) + labs(title = "Second period (2017:M5 to 2023:M8)",
      y = "Domestic Price (log)", x = "Internacional Price (log)") + geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, 
      color = "blue") + annotate("text", x = Inf, y = Inf, label = sprintf("R-squared = %.2f", r2_depois), 
      hjust = 2.5, vjust = 2, size = 5)

# Configurando para salvar a imagem com 300 dpi

png("figure2_300dpi.png", width = 10 * 300, height = 10 * 300, res = 300)

# Organizando os gráficos 

grid.arrange(grobs = list(g1, g2, g3),
  layout_matrix = rbind(c(1, 1),
                        c(2, 3)))

dev.off()

##
## Figura 4 do texto
##

#install.packages("tidyr")

library(tidyr)

# Lendo os dados 

data_fig4 <- read_excel("JBAMR2024a_GCB.xlsx",
                   sheet = "figure4")

# Calculando a proporção

data_fig4$proporcao_milho_soja <- (data_fig4$corn2_areaMT / data_fig4$soy_areaMT) * 100

# Fator de escala adequado para a proporção

fator_escala <- max(data_fig4$soy_areaMT, data_fig4$corn2_areaMT) / 100

# Identificando o valor da proporção para o último ano

ultimo_ano <- tail(data_fig4, 1)

valor_proporcao <- ultimo_ano$proporcao_milho_soja

ano <- as.character(ultimo_ano$crop)

# Transformando os dados para um formato longo

data_fig4 <- data_fig4 %>%
  pivot_longer(cols = c(soy_areaMT, corn2_areaMT), 
               names_to = "Cultura", 
               values_to = "Area")

# Ajustando os nomes das culturas para o gráfico

data_fig4$Cultura <- factor(data_fig4$Cultura, labels = c("Second-Crop Corn Area", "Soybean Area"))

# Definindo uma paleta de cores mais sóbria

cores_sobrias <- c("Soybean Area" = "#fd8d3c", "Second-Crop Corn Area" = "#6baed6", "Corn to Soybean Area Ratio" = "#74c476")

# Criando o gráfico

png("figure4_300dpi.png", width = 10 * 300, height = 5 * 300, res = 300)

ggplot(data_fig4, aes(x = crop)) +
  geom_bar(aes(y = Area, fill = Cultura), 
           stat = "identity", position = position_dodge()) +
  geom_line(data = data_fig4, aes(y = proporcao_milho_soja * fator_escala, group = 1, color = "Corn to Soybean Area Ratio"), size = 1) +
  geom_point(data = data_fig4, aes(y = proporcao_milho_soja * fator_escala, color = "Corn to Soybean Area Ratio"), size = 2) +
  annotate("text", x = ano, y = valor_proporcao * fator_escala, label = sprintf("%.0f%%", valor_proporcao), vjust = -0.5, color = "black", size=3) +
  scale_fill_manual(values = cores_sobrias) +
  scale_color_manual(values = cores_sobrias) +
  scale_y_continuous(
    name = "Area (1,000 ha)",
    sec.axis = sec_axis(~ . / fator_escala, name = "Corn to Soybean Area Ratio (%)")
  ) +
  labs(title = "State of Mato Grosso, Brazil",
       x = "Crop Year",
       y = "",
       fill = "",
       color = "",
       caption = "") +
  theme_minimal() +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(size = 12, hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y.right = element_text(angle = 90),
        panel.background = element_rect(fill = "gray90", color = NA), 
        panel.grid.major = element_line(size = 1, color = "white"), 
        panel.grid.minor = element_line(size = 1, color = "white"))  

dev.off()

### 
### FIM
###





