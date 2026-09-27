### Adaptive MCMC to fit "shrinking and regression" model simultaneously to resource supply and periodic 
### starvation experiments using Biomphalaria glabrata and Schistosoma mansoni

library("adaptMCMC") # needed for MCMC
library("deSolve") # needed to simulate the models
library("tidyverse") # needed for plotting
#library("mvtnorm")
#library("LaplacesDemon")
#library("coda")
#library("GGally")

# The model definition is written in C, so we need to use Rtools and load the model as a dll
rtools <- "C:/rtools45/usr/bin"
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

# # compile my model from C definition
setwd("C:/Users/dcivite/OneDrive - Emory/Emory/Projects/DEBDung")
#setwd("C:/RData/Models_for_Dung_Fitting")

try(dyn.unload("IndividualModel_Shrink_P_dilute.dll")) # unload dll
system("R CMD SHLIB IndividualModel_Shrink_P_dilute.c")
dyn.load("IndividualModel_Shrink_P_dilute.dll") # Load dll

try(dyn.unload("SizeCompModel_Shrink.dll")) # unload dll
system("R CMD SHLIB SizeCompModel_Shrink.c")
dyn.load("SizeCompModel_Shrink.dll") # Load dll

try(dyn.unload("DEBDung_MoA_ingest.dll")) # unload dll
system("R CMD SHLIB DEBDung_MoA_ingest.c")
dyn.load("DEBDung_MoA_ingest.dll") # Load dll


#### Fixed information ####

# 1 - Data
setwd("C:/Users/dcivite/OneDrive - Emory/RData")

data = read.csv("ResourceSupplyExp.csv")
data2 = read.csv("PeriodicStarveExp.csv")
data3 = read.csv("Size_comp_LT.csv")
data4 = read.csv("Dung_LT_vetted.csv")

data = list(t = data$Date, L = data$Length, Negg = data$C_Eggs, Nworms = data$C_Worms, Alive=data$Last_Alive,
            L2 = data2$Length, Negg2 = data2$C_Eggs, Nworms2 = data2$C_Worms, Alive2=data2$Last_Alive,
            L3 = data3$Length_F, L3C = data3$Length_C, Negg3 = data3$C_Eggs, Nworms3 = data3$C_Cercs, Alive3=data3$F_last_alive,
            L4 = data4$Length, Negg4 = data4$C_Eggs, Nworms4 = data4$C_Worms, Alive4 = data4$Last_Alive)

# Dung experiment summary
dung_summary = data4 %>% group_by(Week, FoodType, Infection_Status, Dung_160ml) %>%
  summarise(mean_L = mean(Length, na.rm=T), SE_L = sd(Length, na.rm=T)/sqrt(n()),
            mean_E = mean(C_Eggs, na.rm=T), SE_E = sd(C_Eggs, na.rm=T)/sqrt(n()),
            mean_W = mean(C_Worms, na.rm=T), SE_W = sd(C_Worms, na.rm=T)/sqrt(n()),
            Survival = mean(Alive),
            Day = mean(Week)*7)

dung_summary %>% filter(Dung_160ml == 32 & FoodType =="High" & Infection_Status == 1) %>% select(Week, Survival)

# 2 - Initial conditions
setinits.Food<-function(F0 = 16.5, L0=0.7, e0=0.9, D0 = 0, RH0 = 0, P0 = 0, RP0 = 0, DAM0=0, HAZ0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0)
  return(inits)
}

setinits.Starve<-function(F0 = 10.74, L0=0.7, e0=0.9, D0 = 0, RH0 = 0, P0 = 0, RP0 = 0, DAM0=0, HAZ0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0)
  return(inits)
}

setinits.Size<-function(F0 = 10.76, L0=3.98, e0=0.9, D0 = 0, RH0 = 0, P0 = 2.85e-5, RP0 = 00, DAM0=0, HAZ0=0,
                        LC0=0, eC0=0.9, DC0=0, RHC0=0, DAMC0=0, HAZC0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0,
           LC = LC0, eC = eC0, DC = DC0, RHC=RHC0, DAMC = DAMC0, HAZC = HAZC0)
  return(inits)
}

setinits.Dung<-function(F0 = 10.76, L0=11.23, e0=0.9, D0 = 0, RH0 = 0, P0 = 2.85e-5, RP0 = 00, DAM0=0, DAM20=0, HAZ0=0,
                        Dung0=100){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, DAM2= DAM20, HAZ=HAZ0,
           Dung = Dung0)
  return(inits)
}


samps = readRDS("Dung_fitting2.RDA")
pars = samps$samples[which.max(samps$log.p),]


DEB_parameter_trans = function(x){
  log_par_names = c("ph", "alpha", "iPM", "EM", "DR", "Fh", "muD", "kR", "delta0", "hdelta", "theta", "mR", "hb", "rho", "kR2", "mR2", "kk2", "d02", "kkM", "d0M")
  x[log_par_names] = exp(x[log_par_names])
  x
}

DEB_parameter_backtrans = function(x){
  log_par_names = c("ph", "alpha", "iPM", "EM", "DR", "Fh", "muD", "kR", "delta0", "hdelta", "theta", "mR", "hb", "rho", "kR2", "mR2", "kk2", "d02", "kkM", "d0M")
  x[log_par_names] = log(x[log_par_names])
  x
}

# 3 - Functions for feeding events (these match the timing of feeding events in the actual experiments)
# For the dung experiment
Dung.events = function(initial.food, initial.dung){
  
  time.F = c(sort(c((1:16)*7,((1:16)*7 - 2),((1:16)*7 - 4)))) # gets the third, fifth, and seventh day of every week
  time.D = (1:16)*7 # gets the seventh day of every week
  
  var.F = rep("F", times = length(time.F))
  var.D = rep("Dung", times = length(time.D))
  
  value.F = rep(as.numeric(initial.food), times = length(time.F))
  value.D = rep(as.numeric(initial.dung), times = length(time.D))
  
  method.F = rep("replace", times = length(time.F))
  method.D = rep("replace", times = length(time.D))
  
  out = data.frame(var = c(var.F, var.D), time = c(time.F, time.D), value = c(value.F, value.D), method = c(method.F, method.D))
  out[order(out$time),]
}

# For the size competition experiment
size.events = function(initial.food){
  time.F = c(sort(c((5:19)*7,((5:19)*7 - 3))), 56) # gets the fourth and seventh day of every week
  
  var.F = c(rep("F", times = length(time.F)-1),  "HAZ")
  
  value.F = c(rep(as.numeric(initial.food), times = length(time.F)-1), 0)
  
  method = rep("replace", times = length(time.F))
  out = data.frame(var = var.F, time = time.F, value = value.F, method = method)
  out[order(out$time),]
}

# For the starvation experiment
periodic.starvation.events = function(initial.food, starve.period, infection.date, weeks.duration){
  infection.value = 2.85e-5
  all.fed.dates = infection.date + c(1,4) # Snails were infected on a Thursday, All fed on Friday and Monday, then treatments
  potential.feedings = infection.date + 5 + sort(c((1:weeks.duration)*7 - 6,((1:weeks.duration)*7 - 4),((1:weeks.duration)*7 - 2)))
  #print(0:(weeks.duration/starve.period - 1)*6 + 1)
  potential.feedings.vals = rep(initial.food, times = length(potential.feedings))
  # work out starvation
  if( starve.period == 2){
    potential.feedings.vals[(potential.feedings - (infection.date + 5))  %% 14 <= 7] = 0
  }
  if (starve.period == 3){
    potential.feedings.vals[(potential.feedings - (infection.date + 5)) %% 21 <= 14] = 0
  }
  if (starve.period == 4){
    potential.feedings.vals[(potential.feedings - (infection.date + 5)) %% 28 <= 21] = 0
  }
  # Here the final event date is to later condition survival of infecteds on surviving to being diagnosed
  event.dates = c(infection.date, all.fed.dates, potential.feedings, infection.date, infection.date+28)
  event.values = c(infection.value, initial.food, initial.food, potential.feedings.vals, 0, 0)
  methods = rep("replace", times=length(event.dates))
  vars = c("P", rep("F", times = length(event.dates)-3), "HAZ", "HAZ")
  data.frame(var=vars, time= event.dates, value = event.values, method= methods)
}

# For the resource supply experiment
feeding.events <- function(dates, var="F", value, method="replace", Infected=0){
  # Assemble all data
  events = length(dates) # number of events
  vars = rep(var, times = events)
  values = rep(value, times = events)
  methods = rep(method, times = events)
  
  #build data.frame
  result = data.frame(var = vars, time = dates, value = values, method = methods)
  if(Infected > 0){
    result = rbind(result, data.frame(var="P", time=Infected, value=2.85e-5, method="replace"))
  }
  result = rbind(result, data.frame(var="HAZ", time=63, value=0, method="replace"))
}

dur.S = 133 # This experiment starts on t = 28, ends on t = 133
dur.R = 245 # Starts on t = 0
dur.P = 140 # Starts on t = 0
dur.D = 112 # Starts on t = 0

in.S = setinits.Size()
in.R = setinits.Food()
in.P = setinits.Starve()
#in.D = setinits.Dung(D0 = as.numeric(params["DR"]))


Feed.R = feeding.events(dates = sort(c((1:35)*7,((1:35)*7 - 3))), var="F", in.R[1], method="replace", Infected=28)
Feed.S = size.events(10.76)
Feed.D = Dung.events(initial.food = 10.76, initial.dung = 100)

DEB_dung_vis = function(pars){
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=0, Dung0=0)
  dung0 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                            initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                            params.t[1:34],  rtol=1e-6, atol=1e-6,
                            events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung0
  
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=0, Dung0=40)
  dung40 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung40
  
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=0, Dung0=80)
  dung80 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung80
  
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=0, Dung0=120)
  dung120 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung120
  
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=0, Dung0=160)
  dung160 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung160
  
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=0, Dung0=200)
  dung200 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                           initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                           params.t[1:34],  rtol=1e-6, atol=1e-6,
                           events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung200
  
  control_df = dung_summary %>% filter(Dung_160ml == 0, Infection_Status==0, FoodType=="High")
  dung_40 = dung_summary %>% filter(Dung_160ml == 6.4, Infection_Status==0, FoodType=="High")
  dung_80 = dung_summary %>% filter(Dung_160ml == 12.8, Infection_Status==0, FoodType=="High")
  dung_120 = dung_summary %>% filter(Dung_160ml == 19.2, Infection_Status==0, FoodType=="High")
  dung_160 = dung_summary %>% filter(Dung_160ml == 25.6, Infection_Status==0, FoodType=="High")
  dung_200 = dung_summary %>% filter(Dung_160ml == 32, Infection_Status==0, FoodType=="High")
  
  
  par(mfrow=c(2, 2))
  # plot(F ~ time, data=dung0, type="l", lty=2, col="turquoise1", ylab="Food abundance")
  # lines(F ~ time, data=dung40, type="l", col="chocolate1")
  # lines(F ~ time, data=dung80, type="l", col="chocolate2")
  # lines(F ~ time, data=dung120, type="l", col="chocolate3")
  # lines(F ~ time, data=dung160, type="l", col="chocolate4")
  # lines(F ~ time, data=dung200, type="l", col="black")
  # 
  # plot(Dung ~ time, data=dung0, type="l", col="turquoise1", ylab="Dung abundance", ylim=c(0,1.1*inits["Dung"]))
  # lines(Dung ~ time, data=dung40, type="l", col="chocolate1")
  # lines(Dung ~ time, data=dung80, type="l", col="chocolate2")
  # lines(Dung ~ time, data=dung120, type="l", col="chocolate3")
  # lines(Dung ~ time, data=dung160, type="l", col="chocolate4")
  # lines(Dung ~ time, data=dung200, type="l", col="black")
  
  plot(LG ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylim=c(10, 20), ylab="Snail shell length")
  lines(LG ~ time, data=dung40, type="l", col="chocolate1")
  lines(LG ~ time, data=dung80, type="l", col="chocolate2")
  lines(LG ~ time, data=dung120, type="l", col="chocolate3")
  lines(LG ~ time, data=dung160, type="l", col="chocolate4")
  lines(LG ~ time, data=dung200, type="l", col="black")
  lines(L ~ time, data=dung0, type="l", lwd=3, lty=2, col="turquoise1")
  lines(L ~ time, data=dung40, type="l", lty=2, col="chocolate1")
  lines(L ~ time, data=dung80, type="l", lty=2, col="chocolate2")
  lines(L ~ time, data=dung120, type="l", lty=2, col="chocolate3")
  lines(L ~ time, data=dung160, type="l", lty=2, col="chocolate4")
  lines(L ~ time, data=dung200, type="l", lty=2, col="black")
  points(mean_L ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_L ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_L ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_L ~ Day, dung_120, pch=21, bg="chocolate3")
  points(mean_L ~ Day, dung_160, pch=21, bg="chocolate4")
  points(mean_L ~ Day, dung_200, pch=21, bg="black")
  
  # plot(e ~ time, data=dung0, type="l", col="turquoise1", ylim=c(0, 1.1*max(dung200$e, dung0$e)), ylab="Snail reserve density")
  # lines(e ~ time, data=dung40, type="l", col="chocolate1")
  # lines(e ~ time, data=dung80, type="l", col="chocolate2")
  # lines(e ~ time, data=dung120, type="l", col="chocolate3")
  # lines(e ~ time, data=dung160, type="l", col="chocolate4")
  # lines(e ~ time, data=dung200, type="l", col="black")
  
  plot(RH/0.015 ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylab="Cumulative host reproduction", ylim=c(0, 1.1*max(dung0$RH/0.015, control_df$mean_E)))
  lines(RH/0.015 ~ time, data=dung40, type="l", col="chocolate1")
  lines(RH/0.015 ~ time, data=dung80, type="l", col="chocolate2")
  lines(RH/0.015 ~ time, data=dung120, type="l", col="chocolate3")
  lines(RH/0.015 ~ time, data=dung160, type="l", col="chocolate4")
  lines(RH/0.015 ~ time, data=dung200, type="l", col="black")
  points(mean_E ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_E ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_E ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_E ~ Day, dung_120, pch=21, bg="chocolate3")
  points(mean_E ~ Day, dung_160, pch=21, bg="chocolate4")
  points(mean_E ~ Day, dung_200, pch=21, bg="black")
  
  plot(RP/4e-5 ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylab="Cumulative cercariae release")
  lines(RP/4e-5 ~ time, data=dung40, type="l", col="chocolate1")
  lines(RP/4e-5 ~ time, data=dung80, type="l", col="chocolate2")
  lines(RP/4e-5 ~ time, data=dung120, type="l", col="chocolate3")
  lines(RP/4e-5 ~ time, data=dung160, type="l", col="chocolate4")
  lines(RP/4e-5 ~ time, data=dung200, type="l", col="black")
  points(mean_W ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_W ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_W ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_W ~ Day, dung_120, pch=21, bg="chocolate3")
  points(mean_W ~ Day, dung_160, pch=21, bg="chocolate4")
  points(mean_W ~ Day, dung_200, pch=21, bg="black")
  
  # plot(DAM ~ time, data=dung0, type="l", col="turquoise1", ylab="Damage type 1", ylim=c(0, 1.1*max(dung$DAM, dung0$DAM)))
  # lines(DAM ~ time, data=dung40, type="l", col="chocolate1")
  # lines(DAM ~ time, data=dung80, type="l", col="chocolate2")
  # lines(DAM ~ time, data=dung120, type="l", col="chocolate3")
  # lines(DAM ~ time, data=dung160, type="l", col="chocolate4")
  # lines(DAM ~ time, data=dung200, type="l", col="black")
  # abline(a=params.t["delta0"], b=0, lty=2, col="red")
  # 
  # plot(DAM2 ~ time, data=dung0, type="l", col="turquoise1", ylab="Damage type 2", ylim=c(0, 1.1*max(dung200$DAM2, dung0$DAM2,params.t["d0M"])))
  # lines(DAM2 ~ time, data=dung40, type="l", col="chocolate1")
  # lines(DAM2 ~ time, data=dung80, type="l", col="chocolate2")
  # lines(DAM2 ~ time, data=dung120, type="l", col="chocolate3")
  # lines(DAM2 ~ time, data=dung160, type="l", col="chocolate4")
  # lines(DAM2 ~ time, data=dung200, type="l", col="black")
  # abline(a=params.t["d0M"], b=0, lty=2, col="purple")
  # abline(a=params.t["d02"], b=0, lty=2, col="red")
  
  plot(Survival ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylim=c(0,1.1), ylab="Probability of host survival")
  lines(Survival ~ time, data=dung40, type="l", col="chocolate1")
  lines(Survival ~ time, data=dung80, type="l", col="chocolate2")
  lines(Survival ~ time, data=dung120, type="l", col="chocolate3")
  lines(Survival ~ time, data=dung160, type="l", col="chocolate4")
  lines(Survival ~ time, data=dung200, type="l", col="black")
  points(Survival ~ Day, data = control_df, pch=21, bg="turquoise1")
  points(Survival ~ Day, dung_40, pch=21, bg="chocolate1", col="chocolate1", typ="b")
  points(Survival ~ Day, dung_80, pch=21, bg="chocolate2", col="chocolate2", typ="b")
  points(Survival ~ Day, dung_120, pch=21, bg="chocolate3", col="chocolate3", typ="b")
  points(Survival ~ Day, dung_160, pch=21, bg="chocolate4", col="chocolate4", typ="b")
  points(Survival ~ Day, dung_200, pch=21, bg="black", col="black", typ="b")
  
  ### Infected ###
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=0)
  dung0 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                            initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                            params.t[1:34],  rtol=1e-6, atol=1e-6,
                            events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung0
  
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=40)
  dung40 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                             initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                             params.t[1:34],  rtol=1e-6, atol=1e-6,
                             events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung40
  
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=80)
  dung80 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                             initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                             params.t[1:34],  rtol=1e-6, atol=1e-6,
                             events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung80
  
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=120)
  dung120 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung120
  
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=160)
  dung160 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung160
  
  inits = setinits.Dung(D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=200)
  dung200 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 10.76, initial.dung = as.numeric(inits["Dung"])))))
  dung200
  
  control_df = dung_summary %>% filter(Dung_160ml == 0, Infection_Status==1, FoodType=="High")
  dung_40 = dung_summary %>% filter(Dung_160ml == 6.4, Infection_Status==1, FoodType=="High")
  dung_80 = dung_summary %>% filter(Dung_160ml == 12.8, Infection_Status==1, FoodType=="High")
  dung_120 = dung_summary %>% filter(Dung_160ml == 19.2, Infection_Status==1, FoodType=="High")
  dung_160 = dung_summary %>% filter(Dung_160ml == 25.6, Infection_Status==1, FoodType=="High")
  dung_200 = dung_summary %>% filter(Dung_160ml == 32, Infection_Status==1, FoodType=="High")
  
  
  par(mfrow=c(2, 2))
  # plot(F ~ time, data=dung0, type="l", lty=2, col="turquoise1", ylab="Food abundance")
  # lines(F ~ time, data=dung40, type="l", col="chocolate1")
  # lines(F ~ time, data=dung80, type="l", col="chocolate2")
  # lines(F ~ time, data=dung120, type="l", col="chocolate3")
  # lines(F ~ time, data=dung160, type="l", col="chocolate4")
  # lines(F ~ time, data=dung200, type="l", col="black")
  # 
  # plot(Dung ~ time, data=dung0, type="l", col="turquoise1", ylab="Dung abundance", ylim=c(0,1.1*inits["Dung"]))
  # lines(Dung ~ time, data=dung40, type="l", col="chocolate1")
  # lines(Dung ~ time, data=dung80, type="l", col="chocolate2")
  # lines(Dung ~ time, data=dung120, type="l", col="chocolate3")
  # lines(Dung ~ time, data=dung160, type="l", col="chocolate4")
  # lines(Dung ~ time, data=dung200, type="l", col="black")
  
  plot(LG ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylim=c(10, 20), ylab="Snail shell length")
  lines(LG ~ time, data=dung40, type="l", col="chocolate1")
  lines(LG ~ time, data=dung80, type="l", col="chocolate2")
  lines(LG ~ time, data=dung120, type="l", col="chocolate3")
  lines(LG ~ time, data=dung160, type="l", col="chocolate4")
  lines(LG ~ time, data=dung200, type="l", col="black")
  lines(L ~ time, data=dung0, type="l", lty=2, col="turquoise1")
  lines(L ~ time, data=dung40, type="l", lty=2, col="chocolate1")
  lines(L ~ time, data=dung80, type="l", lty=2, col="chocolate2")
  lines(L ~ time, data=dung120, type="l", lty=2, col="chocolate3")
  lines(L ~ time, data=dung160, type="l", lty=2, col="chocolate4")
  lines(L ~ time, data=dung200, type="l", lty=2, col="black")
  points(mean_L ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_L ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_L ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_L ~ Day, dung_120, pch=21, bg="chocolate3")
  points(mean_L ~ Day, dung_160, pch=21, bg="chocolate4")
  points(mean_L ~ Day, dung_200, pch=21, bg="black")
  
  # plot(e ~ time, data=dung0, type="l", col="turquoise1", ylim=c(0, 1.1*max(dung200$e, dung0$e)), ylab="Snail reserve density")
  # lines(e ~ time, data=dung40, type="l", col="chocolate1")
  # lines(e ~ time, data=dung80, type="l", col="chocolate2")
  # lines(e ~ time, data=dung120, type="l", col="chocolate3")
  # lines(e ~ time, data=dung160, type="l", col="chocolate4")
  # lines(e ~ time, data=dung200, type="l", col="black")
  
  plot(RH/0.015 ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylab="Cumulative host reproduction", ylim=c(0, 1.1*max(dung0$RH/0.015, control_df$mean_E)))
  lines(RH/0.015 ~ time, data=dung40, type="l", col="chocolate1")
  lines(RH/0.015 ~ time, data=dung80, type="l", col="chocolate2")
  lines(RH/0.015 ~ time, data=dung120, type="l", col="chocolate3")
  lines(RH/0.015 ~ time, data=dung160, type="l", col="chocolate4")
  lines(RH/0.015 ~ time, data=dung200, type="l", col="black")
  points(mean_E ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_E ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_E ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_E ~ Day, dung_120, pch=21, bg="chocolate3")
  points(mean_E ~ Day, dung_160, pch=21, bg="chocolate4")
  points(mean_E ~ Day, dung_200, pch=21, bg="black")
  
  plot(RP/4e-5 ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylim=c(0,20000), ylab="Cumulative cercariae release")
  lines(RP/4e-5 ~ time, data=dung40, type="l", col="chocolate1")
  lines(RP/4e-5 ~ time, data=dung80, type="l", col="chocolate2")
  lines(RP/4e-5 ~ time, data=dung120, type="l", col="chocolate3")
  lines(RP/4e-5 ~ time, data=dung160, type="l", col="chocolate4")
  lines(RP/4e-5 ~ time, data=dung200, type="l", col="black")
  points(mean_W ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_W ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_W ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_W ~ Day, dung_120, pch=21, bg="chocolate3")
  points(mean_W ~ Day, dung_160, pch=21, bg="chocolate4")
  points(mean_W ~ Day, dung_200, pch=21, bg="black")
  
  # plot(DAM ~ time, data=dung0, type="l", col="turquoise1", ylab="Damage type 1", ylim=c(0, 1.1*max(dung$DAM, dung0$DAM)))
  # lines(DAM ~ time, data=dung40, type="l", col="chocolate1")
  # lines(DAM ~ time, data=dung80, type="l", col="chocolate2")
  # lines(DAM ~ time, data=dung120, type="l", col="chocolate3")
  # lines(DAM ~ time, data=dung160, type="l", col="chocolate4")
  # lines(DAM ~ time, data=dung200, type="l", col="black")
  # abline(a=params.t["delta0"], b=0, lty=2, col="red")
  # 
  # plot(DAM2 ~ time, data=dung0, type="l", col="turquoise1", ylab="Damage type 2", ylim=c(0, 1.1*max(dung200$DAM2, dung0$DAM2,params.t["d0M"])))
  # lines(DAM2 ~ time, data=dung40, type="l", col="chocolate1")
  # lines(DAM2 ~ time, data=dung80, type="l", col="chocolate2")
  # lines(DAM2 ~ time, data=dung120, type="l", col="chocolate3")
  # lines(DAM2 ~ time, data=dung160, type="l", col="chocolate4")
  # lines(DAM2 ~ time, data=dung200, type="l", col="black")
  # abline(a=params.t["d0M"], b=0, lty=2, col="purple")
  # abline(a=params.t["d02"], b=0, lty=2, col="red")
  
  plot(Survival ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylim=c(0,1.1), ylab="Probability of host survival")
  lines(Survival ~ time, data=dung40, type="l", col="chocolate1")
  lines(Survival ~ time, data=dung80, type="l", col="chocolate2")
  lines(Survival ~ time, data=dung120, type="l", col="chocolate3")
  lines(Survival ~ time, data=dung160, type="l", col="chocolate4")
  lines(Survival ~ time, data=dung200, type="l", col="black")
  points(Survival ~ Day, data = control_df, pch=21, bg="turquoise1")
  points(Survival ~ Day, dung_40, pch=21, bg="chocolate1", col="chocolate1", typ="b")
  points(Survival ~ Day, dung_80, pch=21, bg="chocolate2", col="chocolate2", typ="b")
  points(Survival ~ Day, dung_120, pch=21, bg="chocolate3", col="chocolate3", typ="b")
  points(Survival ~ Day, dung_160, pch=21, bg="chocolate4", col="chocolate4", typ="b")
  points(Survival ~ Day, dung_200, pch=21, bg="black", col="black", typ="b")
  
  ########## Low Food ##########
  
  inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=0, Dung0=0)
  dung0 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                            initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                            params.t[1:34],  rtol=1e-6, atol=1e-6,
                            events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  dung0
  
  inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=0, Dung0=40)
  dung40 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                             initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                             params.t[1:34],  rtol=1e-6, atol=1e-6,
                             events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  dung40
  
  inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=0, Dung0=80)
  dung80 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                             initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                             params.t[1:34],  rtol=1e-6, atol=1e-6,
                             events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  dung80
  
  inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=0, Dung0=120)
  dung120 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  dung120
  
  inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=0, Dung0=160)
  dung160 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  dung160
  
  inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=0, Dung0=200)
  dung200 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  dung200
  
  control_df = dung_summary %>% filter(Dung_160ml == 0, Infection_Status==0, FoodType=="Low")
  dung_40 = dung_summary %>% filter(Dung_160ml == 6.4, Infection_Status==0, FoodType=="Low")
  dung_80 = dung_summary %>% filter(Dung_160ml == 12.8, Infection_Status==0, FoodType=="Low")
  dung_120 = dung_summary %>% filter(Dung_160ml == 19.2, Infection_Status==0, FoodType=="Low")
  dung_160 = dung_summary %>% filter(Dung_160ml == 25.6, Infection_Status==0, FoodType=="Low")
  dung_200 = dung_summary %>% filter(Dung_160ml == 32, Infection_Status==0, FoodType=="Low")
  
  
  par(mfrow=c(2, 2))
  # plot(F ~ time, data=dung0, type="l", lty=2, col="turquoise1", ylab="Food abundance")
  # lines(F ~ time, data=dung40, type="l", col="chocolate1")
  # lines(F ~ time, data=dung80, type="l", col="chocolate2")
  # lines(F ~ time, data=dung120, type="l", col="chocolate3")
  # lines(F ~ time, data=dung160, type="l", col="chocolate4")
  # lines(F ~ time, data=dung200, type="l", col="black")
  # 
  # plot(Dung ~ time, data=dung0, type="l", col="turquoise1", ylab="Dung abundance", ylim=c(0,1.1*inits["Dung"]))
  # lines(Dung ~ time, data=dung40, type="l", col="chocolate1")
  # lines(Dung ~ time, data=dung80, type="l", col="chocolate2")
  # lines(Dung ~ time, data=dung120, type="l", col="chocolate3")
  # lines(Dung ~ time, data=dung160, type="l", col="chocolate4")
  # lines(Dung ~ time, data=dung200, type="l", col="black")
  
  plot(LG ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylim=c(10, 20), ylab="Snail shell length")
  lines(LG ~ time, data=dung40, type="l", col="chocolate1")
  lines(LG ~ time, data=dung80, type="l", col="chocolate2")
  lines(LG ~ time, data=dung120, type="l", col="chocolate3")
  lines(LG ~ time, data=dung160, type="l", col="chocolate4")
  lines(LG ~ time, data=dung200, type="l", col="black")
  lines(L ~ time, data=dung0, type="l", lwd=3, lty=2, col="turquoise1")
  lines(L ~ time, data=dung40, type="l", lty=2, col="chocolate1")
  lines(L ~ time, data=dung80, type="l", lty=2, col="chocolate2")
  lines(L ~ time, data=dung120, type="l", lty=2, col="chocolate3")
  lines(L ~ time, data=dung160, type="l", lty=2, col="chocolate4")
  lines(L ~ time, data=dung200, type="l", lty=2, col="black")
  points(mean_L ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_L ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_L ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_L ~ Day, dung_120, pch=21, bg="chocolate3")
  points(mean_L ~ Day, dung_160, pch=21, bg="chocolate4")
  points(mean_L ~ Day, dung_200, pch=21, bg="black")
  
  # plot(e ~ time, data=dung0, type="l", col="turquoise1", ylim=c(0, 1.1*max(dung200$e, dung0$e)), ylab="Snail reserve density")
  # lines(e ~ time, data=dung40, type="l", col="chocolate1")
  # lines(e ~ time, data=dung80, type="l", col="chocolate2")
  # lines(e ~ time, data=dung120, type="l", col="chocolate3")
  # lines(e ~ time, data=dung160, type="l", col="chocolate4")
  # lines(e ~ time, data=dung200, type="l", col="black")
  
  plot(RH/0.015 ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylab="Cumulative host reproduction", ylim=c(0, 1.1*max(dung0$RH/0.015, control_df$mean_E)))
  lines(RH/0.015 ~ time, data=dung40, type="l", col="chocolate1")
  lines(RH/0.015 ~ time, data=dung80, type="l", col="chocolate2")
  lines(RH/0.015 ~ time, data=dung120, type="l", col="chocolate3")
  lines(RH/0.015 ~ time, data=dung160, type="l", col="chocolate4")
  lines(RH/0.015 ~ time, data=dung200, type="l", col="black")
  points(mean_E ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_E ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_E ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_E ~ Day, dung_120, pch=21, bg="chocolate3")
  points(mean_E ~ Day, dung_160, pch=21, bg="chocolate4")
  points(mean_E ~ Day, dung_200, pch=21, bg="black")
  
  plot(RP/4e-5 ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylab="Cumulative cercariae release")
  lines(RP/4e-5 ~ time, data=dung40, type="l", col="chocolate1")
  lines(RP/4e-5 ~ time, data=dung80, type="l", col="chocolate2")
  lines(RP/4e-5 ~ time, data=dung120, type="l", col="chocolate3")
  lines(RP/4e-5 ~ time, data=dung160, type="l", col="chocolate4")
  lines(RP/4e-5 ~ time, data=dung200, type="l", col="black")
  points(mean_W ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_W ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_W ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_W ~ Day, dung_120, pch=21, bg="chocolate3")
  points(mean_W ~ Day, dung_160, pch=21, bg="chocolate4")
  points(mean_W ~ Day, dung_200, pch=21, bg="black")
  
  # plot(DAM ~ time, data=dung0, type="l", col="turquoise1", ylab="Damage type 1", ylim=c(0, 1.1*max(dung$DAM, dung0$DAM)))
  # lines(DAM ~ time, data=dung40, type="l", col="chocolate1")
  # lines(DAM ~ time, data=dung80, type="l", col="chocolate2")
  # lines(DAM ~ time, data=dung120, type="l", col="chocolate3")
  # lines(DAM ~ time, data=dung160, type="l", col="chocolate4")
  # lines(DAM ~ time, data=dung200, type="l", col="black")
  # abline(a=params.t["delta0"], b=0, lty=2, col="red")
  # 
  # plot(DAM2 ~ time, data=dung0, type="l", col="turquoise1", ylab="Damage type 2", ylim=c(0, 1.1*max(dung200$DAM2, dung0$DAM2,params.t["d0M"])))
  # lines(DAM2 ~ time, data=dung40, type="l", col="chocolate1")
  # lines(DAM2 ~ time, data=dung80, type="l", col="chocolate2")
  # lines(DAM2 ~ time, data=dung120, type="l", col="chocolate3")
  # lines(DAM2 ~ time, data=dung160, type="l", col="chocolate4")
  # lines(DAM2 ~ time, data=dung200, type="l", col="black")
  # abline(a=params.t["d0M"], b=0, lty=2, col="purple")
  # abline(a=params.t["d02"], b=0, lty=2, col="red")
  
  plot(Survival ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylim=c(0,1.1), ylab="Probability of host survival")
  lines(Survival ~ time, data=dung40, type="l", col="chocolate1")
  lines(Survival ~ time, data=dung80, type="l", col="chocolate2")
  lines(Survival ~ time, data=dung120, type="l", col="chocolate3")
  lines(Survival ~ time, data=dung160, type="l", col="chocolate4")
  lines(Survival ~ time, data=dung200, type="l", col="black")
  points(Survival ~ Day, data = control_df, pch=21, bg="turquoise1")
  points(Survival ~ Day, dung_40, pch=21, bg="chocolate1", col="chocolate1", typ="b")
  points(Survival ~ Day, dung_80, pch=21, bg="chocolate2", col="chocolate2", typ="b")
  points(Survival ~ Day, dung_120, pch=21, bg="chocolate3", col="chocolate3", typ="b")
  points(Survival ~ Day, dung_160, pch=21, bg="chocolate4", col="chocolate4", typ="b")
  points(Survival ~ Day, dung_200, pch=21, bg="black", col="black", typ="b")
  
  ### Infected ###
  inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=0)
  dung0 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                            initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                            params.t[1:34],  rtol=1e-6, atol=1e-6,
                            events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  dung0
  
  inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=40)
  dung40 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                             initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                             params.t[1:34],  rtol=1e-6, atol=1e-6,
                             events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  dung40
  
  inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=80)
  dung80 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                             initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                             params.t[1:34],  rtol=1e-6, atol=1e-6,
                             events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  dung80
  
  inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=120)
  dung120 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
                              initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
                              params.t[1:34],  rtol=1e-6, atol=1e-6,
                              events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  dung120
  
  # inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=160)
  # dung160 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
  #                             initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
  #                             params.t[1:34],  rtol=1e-6, atol=1e-6,
  #                             events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  # dung160
  # 
  # inits = setinits.Dung(F0 = 0.538, D0 = as.numeric(params.t["DR"]), P0=2.85e-5, Dung0=200)
  # dung200 <- data.frame(lsoda(inits, 0:dur.D, func = "derivs", dllname = "DEBDung_MoA_ingest",
  #                             initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=5e5,
  #                             params.t[1:34],  rtol=1e-6, atol=1e-6,
  #                             events = list(data = Dung.events(initial.food = 0.538, initial.dung = as.numeric(inits["Dung"])))))
  # dung200
  
  control_df = dung_summary %>% filter(Dung_160ml == 0, Infection_Status==1, FoodType=="Low")
  dung_40 = dung_summary %>% filter(Dung_160ml == 6.4, Infection_Status==1, FoodType=="Low")
  dung_80 = dung_summary %>% filter(Dung_160ml == 12.8, Infection_Status==1, FoodType=="Low")
  dung_120 = dung_summary %>% filter(Dung_160ml == 19.2, Infection_Status==1, FoodType=="Low")
  #dung_160 = dung_summary %>% filter(Dung_160ml == 25.6, Infection_Status==1, FoodType=="Low")
  #dung_200 = dung_summary %>% filter(Dung_160ml == 32, Infection_Status==1, FoodType=="Low")
  
  
  par(mfrow=c(2, 2))
  # plot(F ~ time, data=dung0, type="l", lty=2, col="turquoise1", ylab="Food abundance")
  # lines(F ~ time, data=dung40, type="l", col="chocolate1")
  # lines(F ~ time, data=dung80, type="l", col="chocolate2")
  # lines(F ~ time, data=dung120, type="l", col="chocolate3")
  # lines(F ~ time, data=dung160, type="l", col="chocolate4")
  # lines(F ~ time, data=dung200, type="l", col="black")
  # 
  # plot(Dung ~ time, data=dung0, type="l", col="turquoise1", ylab="Dung abundance", ylim=c(0,1.1*inits["Dung"]))
  # lines(Dung ~ time, data=dung40, type="l", col="chocolate1")
  # lines(Dung ~ time, data=dung80, type="l", col="chocolate2")
  # lines(Dung ~ time, data=dung120, type="l", col="chocolate3")
  # lines(Dung ~ time, data=dung160, type="l", col="chocolate4")
  # lines(Dung ~ time, data=dung200, type="l", col="black")
  
  plot(LG ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylim=c(10, 20), ylab="Snail shell length")
  lines(LG ~ time, data=dung40, type="l", col="chocolate1")
  lines(LG ~ time, data=dung80, type="l", col="chocolate2")
  lines(LG ~ time, data=dung120, type="l", col="chocolate3")
  #lines(LG ~ time, data=dung160, type="l", col="chocolate4")
  #lines(LG ~ time, data=dung200, type="l", col="black")
  lines(L ~ time, data=dung0, type="l",  lwd=3, lty=2, col="turquoise1")
  lines(L ~ time, data=dung40, type="l", lty=2, col="chocolate1")
  lines(L ~ time, data=dung80, type="l", lty=2, col="chocolate2")
  lines(L ~ time, data=dung120, type="l", lty=2, col="chocolate3")
  #lines(L ~ time, data=dung160, type="l", lty=2, col="chocolate4")
  #lines(L ~ time, data=dung200, type="l", lty=2, col="black")
  points(mean_L ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_L ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_L ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_L ~ Day, dung_120, pch=21, bg="chocolate3")
  #points(mean_L ~ Day, dung_160, pch=21, bg="chocolate4")
  #points(mean_L ~ Day, dung_200, pch=21, bg="black")
  
  # plot(e ~ time, data=dung0, type="l", col="turquoise1", ylim=c(0, 1.1*max(dung200$e, dung0$e)), ylab="Snail reserve density")
  # lines(e ~ time, data=dung40, type="l", col="chocolate1")
  # lines(e ~ time, data=dung80, type="l", col="chocolate2")
  # lines(e ~ time, data=dung120, type="l", col="chocolate3")
  # lines(e ~ time, data=dung160, type="l", col="chocolate4")
  # lines(e ~ time, data=dung200, type="l", col="black")
  
  plot(RH/0.015 ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylab="Cumulative host reproduction")
  lines(RH/0.015 ~ time, data=dung40, type="l", col="chocolate1")
  lines(RH/0.015 ~ time, data=dung80, type="l", col="chocolate2")
  lines(RH/0.015 ~ time, data=dung120, type="l", col="chocolate3")
  #lines(RH/0.015 ~ time, data=dung160, type="l", col="chocolate4")
  #lines(RH/0.015 ~ time, data=dung200, type="l", col="black")
  points(mean_E ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_E ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_E ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_E ~ Day, dung_120, pch=21, bg="chocolate3")
  #points(mean_E ~ Day, dung_160, pch=21, bg="chocolate4")
  #points(mean_E ~ Day, dung_200, pch=21, bg="black")
  
  plot(RP/4e-5 ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylim=c(0,20000), ylab="Cumulative cercariae release")
  lines(RP/4e-5 ~ time, data=dung40, type="l", col="chocolate1")
  lines(RP/4e-5 ~ time, data=dung80, type="l", col="chocolate2")
  lines(RP/4e-5 ~ time, data=dung120, type="l", col="chocolate3")
  #lines(RP/4e-5 ~ time, data=dung160, type="l", col="chocolate4")
  #lines(RP/4e-5 ~ time, data=dung200, type="l", col="black")
  points(mean_W ~ Day, control_df, pch=21, bg="turquoise1")
  points(mean_W ~ Day, dung_40, pch=21, bg="chocolate1")
  points(mean_W ~ Day, dung_80, pch=21, bg="chocolate2")
  points(mean_W ~ Day, dung_120, pch=21, bg="chocolate3")
  #points(mean_W ~ Day, dung_160, pch=21, bg="chocolate4")
  #points(mean_W ~ Day, dung_200, pch=21, bg="black")
  
  # plot(DAM ~ time, data=dung0, type="l", col="turquoise1", ylab="Damage type 1", ylim=c(0, 1.1*max(dung$DAM, dung0$DAM)))
  # lines(DAM ~ time, data=dung40, type="l", col="chocolate1")
  # lines(DAM ~ time, data=dung80, type="l", col="chocolate2")
  # lines(DAM ~ time, data=dung120, type="l", col="chocolate3")
  # lines(DAM ~ time, data=dung160, type="l", col="chocolate4")
  # lines(DAM ~ time, data=dung200, type="l", col="black")
  # abline(a=params.t["delta0"], b=0, lty=2, col="red")
  # 
  # plot(DAM2 ~ time, data=dung0, type="l", col="turquoise1", ylab="Damage type 2", ylim=c(0, 1.1*max(dung200$DAM2, dung0$DAM2,params.t["d0M"])))
  # lines(DAM2 ~ time, data=dung40, type="l", col="chocolate1")
  # lines(DAM2 ~ time, data=dung80, type="l", col="chocolate2")
  # lines(DAM2 ~ time, data=dung120, type="l", col="chocolate3")
  # lines(DAM2 ~ time, data=dung160, type="l", col="chocolate4")
  # lines(DAM2 ~ time, data=dung200, type="l", col="black")
  # abline(a=params.t["d0M"], b=0, lty=2, col="purple")
  # abline(a=params.t["d02"], b=0, lty=2, col="red")
  
  plot(Survival ~ time, data=dung0, type="l", lwd=3, col="turquoise1", ylim=c(0,1.1), ylab="Probability of host survival")
  lines(Survival ~ time, data=dung40, type="l", col="chocolate1")
  lines(Survival ~ time, data=dung80, type="l", col="chocolate2")
  lines(Survival ~ time, data=dung120, type="l", col="chocolate3")
  #lines(Survival ~ time, data=dung160, type="l", col="chocolate4")
  #lines(Survival ~ time, data=dung200, type="l", col="black")
  points(Survival ~ Day, data = control_df, pch=21, bg="turquoise1")
  points(Survival ~ Day, dung_40, pch=21, bg="chocolate1", col="chocolate1", typ="b")
  points(Survival ~ Day, dung_80, pch=21, bg="chocolate2", col="chocolate2", typ="b")
  points(Survival ~ Day, dung_120, pch=21, bg="chocolate3", col="chocolate3", typ="b")
  #points(Survival ~ Day, dung_160, pch=21, bg="chocolate4", col="chocolate4", typ="b")
  #points(Survival ~ Day, dung_200, pch=21, bg="black", col="black", typ="b")
  
  
  par(mfrow=c(1,1))
}

params.t = DEB_parameter_trans(pars)
params.t["rho"] = 1e-5 # smaller number favors dung (ratio of consumption rates)
params.t["yEF3"] = .2
params.t["yED"] = 2.5e-1
params.t["kR2"] = 0.1
params.t["mR2"] = 1e-3
params.t["kkM"] = 0#2e-2
# params.t["d02"] = 80
params.t["d0M"] = 60
# params.t["kk2"] = 1e-3
# # #params.t["d02"] = 1

DEB_dung_vis(params.t)
