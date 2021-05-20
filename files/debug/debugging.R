library(knitr)
library(ggplot2)
library(reshape2)
library(plyr)
library(Rsolnp)
library(nloptr)

# Path to exchange folder
path_to_cryptos = "../files/csvs/libra/crypto"
files <- list.files(path_to_cryptos)
names <- gsub("\\ .*","", files)
csv_list <- list()

for(i in 1:length(files)){
  csv_list[[names[i]]] <- read.csv(file = paste0(path_to_cryptos, "/", files[i]))[c(1,2)]
}


df_fx_crypto <- merge(csv_list[[1]], csv_list[[2]], by = "Date")
colnames(df_fx_crypto)[c(2,3)] <- names[c(1,2)]

for(i in 3:length(csv_list)){
  
  df_fx_crypto <- merge(df_fx_crypto, csv_list[[i]], by = "Date")
  colnames(df_fx_crypto)[i + 1] <- names[i]
  
}

# Drop Doge coin too many spikes
df_fx_crypto$DOGE_USD <- NULL

df_fx_crypto$Date <- as.Date(df_fx_crypto$Date, "%b %d, %Y")
df_fx_crypto <- df_fx_crypto[order(df_fx_crypto$Date),]
df_fx_flipped <- df_fx[order(df_fx_crypto$Date, decreasing = F),]

# Currency names
curr_names <- c("USD", "BTC", "ETH", "LTC", "XRP")

# Get the x-rates matrix
srates <- c(1, as.numeric(df_fx_crypto[1, 2:5]))
xrates <- (srates) %*% t(1/srates)
rownames(xrates) <- curr_names
colnames(xrates) <- curr_names

nval <- xrates / apply(xrates, 2, function(x) (prod(x[x!=0]))^(1/sum(x!=0)))

df_nfx <- cbind(1, df_fx_crypto[,2:5])
df_nfx <- df_nfx/apply(df_nfx, 1, function(x) (prod(x[x!=0]))^(1/sum(x!=0)))

df_rnfx_crypto <- df_nfx / c(df_nfx[1,])

colnames(df_rnfx_crypto) <- curr_names
df_rnfx_crypto <- cbind(df_fx_crypto$Date, df_rnfx_crypto)
colnames(df_rnfx_crypto)[1] <- "Date"


optimWeights3 <- function(RNVals, target = 1, init_vals){
  
  
  target_curr <- RNVals[,target]
  df_rest <- RNVals[,-target]
  
  CovTarg <- cov(df_rest, target_curr)
  CovRest <- cov(df_rest)
  
  obj <- function(w, CovTarg, CovRest){
    
    x <- (w %*% CovTarg / sqrt((w) %*% CovRest %*% (w)))
    return(x)
  }
  
  equal <- function(w, CovTarg, CovRest) {
    sum(w) - 1 
  }
  
  # Lower and upper bounds
  lb <- c(0,0,0,0)
  ub <- c(1,1,1,1)
  #initial values
  n <- nrow(CovRest)
  z <- CovRest
  y <- CovTarg
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0,
                "check_derivatives" = TRUE,
                "print_options_doc" = TRUE)
  
  res <- nloptr (x0 = init_vals,
                 eval_f = obj,
                 lb = lb,
                 ub = ub,
                 eval_g_eq = equal,
                 opts = opts,
                 CovTarg = y,
                 CovRest = z)
  
  return(res)
}


Weights_crypto <- optimWeights3(df_rnfx_crypto[,2:6], target=1, init_vals = c(0.25,0.25,0.25,0.25))

####################################################################
####################################################################
##  Experimenting with FX
####################################################################
####################################################################

# Clear WS
rm(list = ls())

# Path to exchange folder
path_to_exrates = "../files/csvs/libra/fx"
files <- list.files(path_to_exrates)
names <- gsub("\\ .*","", files)
csv_list <- list()

for(i in 1:length(files)){
  csv_list[[names[i]]] <- read.csv(file = paste0(path_to_exrates, "/", files[i]))[c(1,2)]
}

df_fx <- merge(csv_list[[1]], csv_list[[2]], by = "Date")
colnames(df_fx)[c(2,3)] <- names[c(1,2)]

for(i in 3:length(csv_list)){
  
  df_fx <- merge(df_fx, csv_list[[i]], by = "Date")
  colnames(df_fx)[i + 1] <- names[i]
  
}


df_fx$Date <- as.Date(df_fx$Date, "%b %d, %Y")
df_fx <- df_fx[order(df_fx$Date),]

# Cut down the years 

df_fx <- df_fx[df_fx$Date <= as.Date("2003-01-01"),]

# Drop some currs INR for now
df_fx["USD_CNY"] <- NULL
df_fx["USD_INR"] <- NULL

df_fx_flipped <- df_fx[order(df_fx$Date, decreasing = T),]

# Currency names
curr_names <- c("USD", "AUD", "CAD", "CHF", "EUR", "GBP", "JPY")

# Get the x-rates matrix
srates <- c(1, as.numeric(df_fx[1,2:7]))
xrates <- (srates) %*% t(1/srates)
rownames(xrates) <- curr_names
colnames(xrates) <- curr_names

nval <- xrates / apply(xrates, 2, function(x) (prod(x[x!=0]))^(1/sum(x!=0)))

## Print the cross rates at t= 0
print(xrates)

# Get the RNVAls
df_nfx <- cbind(1, df_fx[,2:7])
df_nfx <- df_nfx/apply(df_nfx, 1, function(x) (prod(x[x!=0]))^(1/sum(x!=0)))

df_rnfx <- df_nfx / c(df_nfx[1,])

colnames(df_rnfx) <- curr_names
df_rnfx <- cbind(df_fx$Date, df_rnfx)
colnames(df_rnfx)[1] <- "Date"

kable(head(df_rnfx))
kable(cor(df_rnfx[,2:7]))


### Minimisation problem 

optimWeights3 <- function(RNVals, target = 1, init_vals, lb, ub, linear_cons = T){
  
  
  target_curr <- RNVals[,target]
  df_rest <- RNVals[,-target]
  
  CovTarg <- cov(df_rest, target_curr)
  CovRest <- cov(df_rest)
  
  # obj <- function(w, target_curr, df_rest){
  #   
  #   x <- -cor(target_curr, as.vector(w %*% t(df_rest)))
  #   return(x)
  #   
  # }
  
  obj <- function(w, CovTarg, CovRest){

    xx <- as.numeric(-1 * ((w %*% (CovTarg) / sqrt(((w) %*% CovRest %*% (w))))))
    return(xx)
    
  }
  
  equal_c <- function(w, CovTarg, CovRest) {
    
   eq_c <- list("constraints" = sum(w) - 1,
                "jacobian" = rep(1,length(w)))
   return(eq_c)
  }
  
  obj_grad <- function(w, CovTarg, CovRest){
    
    u <- as.numeric(w %*% (CovTarg))
    u_ <- CovTarg
    v <- as.numeric(t(w) %*% CovRest %*% (w))
    v_ <- 2 * t(w) %*% CovRest
    x <- ((u_*v - t(v_*u)) / v^2) * -1
    return(x)
    
  }
  
  equal_g <- function(w, CovTarg, CovRest) {
    return(c(0,0,0,0))
  }
  
  # Lower and upper bounds
  lb <- lb
  ub <- ub
  #initial values
  n <- nrow(CovRest)
  z <- CovRest
  y <- CovTarg
  
  if(linear_cons == TRUE){
    
    local_opts <- list("algorithm" = "NLOPT_LD_TNEWTON", "xtol_rel" = 1.0e-5)
    opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                  "xtol_abs"= 1.0e-5,
                  "maxeval"= 200000,
                  "local_opts" = local_opts,
                  "tol_constraints_eq" = 1e-3,
                  "print_level" = 3,
                  "check_derivatives" = TRUE,
                  "check_derivatives_print" = "all",
                  "print_options_doc" = TRUE)
    
    res <- nloptr (x0 = init_vals,
                   eval_f = obj,
                   eval_g_eq = equal_c, 
                   eval_grad_f = obj_grad,
                   lb = lb,
                   ub = ub,
                   opts = opts,
                   CovTarg = y,
                   CovRest = z)
    
  } else {
    
    res <- optimx(par = init_vals, fn = obj, gr = obj_grad,
                  method = "CG",
                  CovTarg = y,
                  CovRest = z,
                  control = list(trace = 9, parscale=c(rep(0.1,6))))
    
  }
  
  return(res)
}


Weights_fx <- optimWeights3(df_rnfx[,c(2:8)],
                            target=1,
                            init_vals = c(rep(1/6,6)),
                            lb = c(0,0,0,0,0,0),
                            ub = c(1,1,1,1,1,1),
                            linear_cons = T)

Weights_fx <- Weights_fx$solution

# Quantities
Q_fx <- Weights_fx/xrates[1,2:7]
USDFX <- as.vector(Weights_fx %*% t(df_rnfx[,3:8]))
df_rnfx_plus_usdfx <- cbind(df_rnfx, USDFX)
kable(round(cor(df_rnfx_plus_usdfx[2:9]), 3))
# 
# 
#Weights_HSK <- c(0.1519, 0.2064, 0.1390, 0.1442, 0.2040, 0.1545)
#USDFX <- as.vector(Weights_HSK %*% t(df_rnfx[,3:8]))

df_toplot <- melt(df_rnfx_plus_usdfx, id.vars = "Date")

ggplot(data = df_toplot, aes(y = value, x = Date, colour = variable)) +
  geom_line() +
  ylab("Reduced Normalised Value") +
  xlab("Date") +
  theme_classic() +
  scale_colour_brewer(type = "div", palette = "Dark2") +
  labs(colour = "Currency") +
  theme(legend.position="top",
        legend.text = element_text(size = 8))

USDFX_USD <- as.data.frame(t(Q_fx %*% t(cbind(1, 1/df_fx[,3:7]))))
colnames(USDFX_USD) <- "USDFX_USD"
df_fx_sac <- cbind(df_fx, USDFX_USD)

max(df_fx_sac$USDFX_USD)

df_fx_sac_plot <- melt(df_fx_sac[,c(1,2,3,6,8)], id.vars = "Date")

ggplot(data = df_fx_sac_plot, aes(y = value, x = Date, colour = variable)) + 
  geom_line() +
  ylab("Exchange Rate") +
  xlab("Date") +
  theme_classic() +
  scale_colour_brewer(type = "div", palette = "Dark2", labels = c("USD_FX")) +
  labs(colour = "Currency") +
  theme(legend.position="top",
        legend.text = element_text(size = 8))



cov_obj <- function(w, CovTarg, CovRest){
  
  x <- ((w %*% (CovTarg) / sqrt((w) %*% CovRest %*% (w))))
  return(x)
  
}

# Get the RNVAls
df_nfx_2 <- cbind(1, df_fx[2:3])
df_nfx_2 <- df_nfx_2/apply(df_nfx_2, 1, function(x) (prod(x[x!=0]))^(1/sum(x!=0)))

df_rnfx_2 <- df_nfx_2 / c(df_nfx_2[1,])
curr_names_2 <- c("USD", "AUD", "CAD")
colnames(df_rnfx_2) <- curr_names_2
  
a <- cov(df_rnfx_2[,c(2,3)], df_rnfx_2[,1])
W <- cov(df_rnfx_2[,c(2,3)])

w1 <- as.vector(seq(from = 0.001,to = 1, length.out = 100))
w2 <- as.vector(seq(from = 0.001,to = 1, length.out = 100))

z <- matrix(nrow = length(w1), ncol = length(w2))

for(i in 1:length(w1)){
  for(j in 1:length(w2)){
    ww <- c(w1[i],w2[j])
    z[i,j] <- cov_obj(ww, a, W)
  }
}

persp(z = z, x = w1, y = w2, theta=45, phi=10, r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, col="cyan")

Weights_fx <- optimWeights3(df_rnfx_2, 
                            target=1, 
                            init_vals = c(rep(1/2, 2)),
                            lb = c(0,0),
                            ub = c(1,1))


