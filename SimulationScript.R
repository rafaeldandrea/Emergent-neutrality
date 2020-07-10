## Libraries
library(dplyr)    ## for data wrangling
library(ggplot2)  ## for plotting
library(tidyr)    ## for function crossing()

## Load functions
dummy <-
  list.files(
    path = '~/R_Functions/',
    pattern = "[.][R]$",
    full.names = TRUE,
    ignore.case = TRUE,
    recursive = TRUE
  ) %>%
  lapply(source)

## Set working directory
setwd('~/EmergentNeutrality/')

## Resume simulation (TRUE) or start from scratch (FALSE)
resume_simulation = FALSE

## Keep track of species extinction times
keep_track_extinctions = TRUE

## Save parameter data along with abundances
save_parameters = TRUE

## Set folder to save data files in
datadir = '~/EmergentNeutrality/Data/20200122/'

## Indicator variable; 1 if running on High Power Computing cluster, 0 if running on my PC.
hpcc = as.numeric((nodename <- as.character(Sys.info()['nodename'])) != 'RAF-PC')

## Read scenario from either input on command line (if running on hpcc)
## or chosen scenario listed in the line below (if running on PC)
simind = ifelse(hpcc, as.numeric(commandArgs(TRUE)[1]), 991)

## Start simulation from scratch
if (!resume_simulation) {
  ## Define parameter table
  {
    scenarios <- {
      crossing(
        limit = c('Generalists'),
        ## equalize resource budget of all consumers?
        parametrization = c('FixedAbuns'),
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        NR = 50,
        ## number of resources
        NS = 50,
        ## number of consumers
        Rstar = 100,
        ## expected resource abun
        Nstar = seq(500, 1000, l = 5),
        ## expected consumer abun
        mu = seq(.05, 2, l = 5),
        ## mean consumption rate in C matrix
        b = seq(0, 1, l = 11),
        ## semirange of C matrix
        run = 1                                              ## controls seed of random number generator
      ) %>%
        mutate(b = b * mu) %>%
        filter(b > 0) %>%
        arrange(Nstar, mu, b) %>%
        mutate(scenario = 1:n())
    }
    
    higherNstar <- {
      crossing(
        limit = c('Generalists'),
        ## equalize resource budget of all consumers?
        parametrization = c('FixedAbuns'),
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        NR = 50,
        ## number of resources
        NS = 50,
        ## number of consumers
        Rstar = 100,
        ## expected resource abun
        Nstar = c(2000, 5000, 10000),
        mu = seq(.05, 2, l = 5),
        b = c(seq(.1, 1, l = 5), seq(1.5, 5, l = 5)),
        run = 1
      ) %>%
        mutate(b = mu * b * exp(3.9677479) * Nstar ^ -0.7906098) %>%
        arrange(Nstar, mu, b) %>%
        mutate(scenario = nrow(scenarios) + 1:n())
    }
    
    specialists <- {
      tibble(
        limit = 'Specialists',
        ## equalize resource budget of all consumers?
        parametrization = 'FixedAbuns',
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        Rstar = 100,
        ## expected resource abun
        Nstar = 100,
        NS = rep(c(25, 10 * (5:10)), each = 10),
        NR = NS,
        mu = as.numeric(sapply(unique(NS), function(n)
          seq(1.1, 1.5335 + 0.1653 * n, l = 10))),
        b = 1,
        run = 1
      ) %>%
        arrange(NR, NS, mu, b) %>%
        mutate(scenario = nrow(rbind(scenarios, higherNstar)) + 1:n())
    }
    
    neutral <- {
      tibble(
        limit = 'Neutral',
        ## equalize resource budget of all consumers?
        parametrization = 'FixedAbuns',
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        Rstar = 100,
        ## expected resource abun
        Nstar = 500,
        NS = 50,
        NR = NS,
        mu = .05,
        b = 0,
        run = 1:10
      ) %>%
        mutate(scenario = nrow(rbind(
          scenarios, higherNstar, specialists
        )) + 1:n())
    }
    
    sad_vs_extime <- {
      crossing(
        limit = 'Generalists',
        ## equalize resource budget of all consumers?
        parametrization = 'FixedAbuns',
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        Rstar = 100,
        ## expected resource abun
        Nstar = 500,
        NS = 50,
        NR = 50,
        mu = .05,
        b = seq(.1, 1, l = 5),
        run = 1:10
      ) %>%
        mutate(b = mu * b,
               scenario = nrow(rbind(
                 scenarios, higherNstar, specialists, neutral
               )) + 1:n())
    }
    
    specialists_sad_vs_extime <- {
      crossing(
        limit = 'Specialists',
        ## equalize resource budget of all consumers?
        parametrization = 'FixedAbuns',
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        Rstar = 100,
        ## expected resource abun
        Nstar = 100,
        NS = 50,
        NR = 50,
        mu = as.numeric(seq(1.1, 1.5335 + 0.1653 * 50, l = 5)),
        b = 1,
        run = 1:10
      ) %>%
        mutate(scenario = nrow(
          rbind(
            scenarios,
            higherNstar,
            specialists,
            neutral,
            sad_vs_extime
          )
        ) + 1:n())
    }
    
    specialists_higher_NS <- {
      tibble(
        limit = 'Specialists',
        ## equalize resource budget of all consumers?
        parametrization = 'FixedAbuns',
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        Rstar = 100,
        ## expected resource abun
        Nstar = 100,
        NS = rep(c(500, 1000, 2000), each = 10),
        NR = NS,
        mu = as.numeric(sapply(unique(NS), function(n)
          seq(1.1, 1.5335 + 0.1653 * n, l = 10))),
        b = 1,
        run = 1
      ) %>%
        arrange(NS, mu) %>%
        mutate(scenario = nrow(
          rbind(
            scenarios,
            higherNstar,
            specialists,
            neutral,
            sad_vs_extime,
            specialists_sad_vs_extime
          )
        ) + 1:n())
    }
    
    specialists_sad_vs_extime_lower_Cd <- {
      crossing(
        limit = 'Specialists',
        ## equalize resource budget of all consumers?
        parametrization = 'FixedAbuns',
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        Rstar = 100,
        ## expected resource abun
        Nstar = 100,
        NS = 50,
        NR = 50,
        mu = 1 + seq(0, .08, l = 5),
        b = 1,
        run = 1:10
      ) %>%
        mutate(scenario = nrow(
          rbind(
            scenarios,
            higherNstar,
            specialists,
            neutral,
            sad_vs_extime,
            specialists_sad_vs_extime,
            specialists_higher_NS
          )
        ) + 1:n())
    }
    
    full_niche <- {
      tibble(
        limit = 'Specialists',
        ## equalize resource budget of all consumers?
        parametrization = 'FixedAbuns',
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        Rstar = 100,
        ## expected resource abun
        Nstar = 100,
        NS = 50,
        NR = NS,
        mu = 1,
        b = 0,
        run = 1
      ) %>%
        mutate(scenario = nrow(
          rbind(
            scenarios,
            higherNstar,
            specialists,
            neutral,
            sad_vs_extime,
            specialists_sad_vs_extime,
            specialists_higher_NS,
            specialists_sad_vs_extime_lower_Cd
          )
        ) + 1:n())
    }
    
    specialists_speciation <- {
      crossing(
        limit = 'Specialists',
        ## equalize resource budget of all consumers?
        parametrization = 'FixedAbuns',
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = FALSE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 1e-3,
        ## speciation rate
        Rstar = 100,
        ## expected resource abun
        Nstar = 100,
        NS = 50,
        NR = 50,
        mu = as.numeric(seq(1.1, 1.5335 + 0.1653 * 50, l = 5)),
        b = 1,
        run = 1:10
      ) %>%
        mutate(scenario = nrow(
          rbind(
            scenarios,
            higherNstar,
            specialists,
            neutral,
            sad_vs_extime,
            specialists_sad_vs_extime,
            specialists_higher_NS,
            specialists_sad_vs_extime_lower_Cd,
            full_niche
          )
        ) + 1:n())
    }
    
    specialists_speciation2 <- {
      tibble(
        limit = 'Specialists',
        ## equalize resource budget of all consumers?
        parametrization = 'FixedAbuns',
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = FALSE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 1e-3,
        ## speciation rate
        Rstar = 100,
        ## expected resource abun
        Nstar = 100,
        NS = rep(10 * (2:8), each = 10),
        NR = NS,
        mu = as.numeric(sapply(unique(NS), function(n)
          seq(1.1, 1.5335 + 0.1653 * n, l = 10))),
        b = 1,
        run = 1
      ) %>%
        mutate(scenario = nrow(
          rbind(
            scenarios,
            higherNstar,
            specialists,
            neutral,
            sad_vs_extime,
            specialists_sad_vs_extime,
            specialists_higher_NS,
            specialists_sad_vs_extime_lower_Cd,
            full_niche,
            specialists_speciation
          )
        ) + 1:n())
    }
    
    specialists_logisticresources <- {
      tibble(
        limit = 'Specialists',
        ## equalize resource budget of all consumers?
        parametrization = 'FixedAbuns_LogisticResources',
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        Rstar = 100,
        ## expected resource abun
        Nstar = 100,
        NS = rep(10 * (2:8), each = 10),
        NR = NS,
        mu = as.numeric(sapply(unique(NS), function(n)
          seq(1.1, 1.5335 + 0.1653 * n, l = 10))),
        b = 1,
        run = 1
      ) %>%
        mutate(scenario = nrow(
          rbind(
            scenarios,
            higherNstar,
            specialists,
            neutral,
            sad_vs_extime,
            specialists_sad_vs_extime,
            specialists_higher_NS,
            specialists_sad_vs_extime_lower_Cd,
            full_niche,
            specialists_speciation,
            specialists_speciation2
          )
        ) + 1:n())
    }
    
    generalists_gaussian <- {
      crossing(
        limit = c('Generalists_Gaussian'),
        ## equalize resource budget of all consumers?
        parametrization = c('FixedAbuns'),
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        NR = 50,
        ## number of resources
        NS = 50,
        ## number of consumers
        Rstar = 100,
        ## expected resource abun
        Nstar = exp(seq(4, 7, l = 7)),
        ## expected consumer abun
        mu = 2,
        ## mean consumption rate in C matrix
        b = seq(.01, .6, l = 10),
        ## desired CV of C matrix
        run = 1                                ## controls seed of random number generator
      ) %>%
        arrange(Nstar, b) %>%
        mutate(scenario = nrow(
          rbind(
            scenarios,
            higherNstar,
            specialists,
            neutral,
            sad_vs_extime,
            specialists_sad_vs_extime,
            specialists_higher_NS,
            specialists_sad_vs_extime_lower_Cd,
            full_niche,
            specialists_speciation,
            specialists_speciation2,
            specialists_logisticresources
          )
        ) + 1:n())
    }
    
    generalists_uniform <- {
      crossing(
        limit = c('Generalists_Uniform'),
        ## equalize resource budget of all consumers?
        parametrization = c('FixedAbuns'),
        ## fix rates eta, rho, or equil abuns N*, R*?
        bounce = TRUE,
        ## do consumers immediately bounce from extinction?
        epsilon = 1,
        ## consumer conversion efficiency
        nu = 0,
        ## speciation rate
        NR = 50,
        ## number of resources
        NS = 50,
        ## number of consumers
        Rstar = 100,
        ## expected resource abun
        Nstar = exp(seq(4, 7, l = 7)),
        ## expected consumer abun
        mu = 2,
        ## mean consumption rate in C matrix
        b = seq(.01, .6, l = 10),
        ## desired CV of C matrix
        run = 1                                ## controls seed of random number generator
      ) %>%
        arrange(Nstar, b) %>%
        mutate(scenario = nrow(
          rbind(
            scenarios,
            higherNstar,
            specialists,
            neutral,
            sad_vs_extime,
            specialists_sad_vs_extime,
            specialists_higher_NS,
            specialists_sad_vs_extime_lower_Cd,
            full_niche,
            specialists_speciation,
            specialists_speciation2,
            specialists_logisticresources,
            generalists_gaussian
          )
        ) + 1:n())
    }
    
    scenarios <- {
      scenarios %>%
        union(higherNstar) %>%
        union(specialists) %>%
        union(neutral) %>%
        union(sad_vs_extime) %>%
        union(specialists_sad_vs_extime) %>%
        union(specialists_higher_NS) %>%
        union(specialists_sad_vs_extime_lower_Cd) %>%
        union(full_niche) %>%
        union(specialists_speciation) %>%
        union(specialists_speciation2) %>%
        union(specialists_logisticresources) %>%
        union(generalists_gaussian) %>%
        union(generalists_uniform) %>%
        arrange(scenario)
    }
  }
  
  scen = scenarios[simind, ]
  
  ## Assign each column of scen to a variable matching its name
  list2env(as.list(scen), envir = globalenv())
  
  ## Initial conditions
  resources = rep(round(Rstar), NR)			## initial resource abundances
  consumers = rep(round(Nstar), NS)			## initial consumer abundances
  RN = outer(resources, consumers)				## R*N (abundance products)
  simtime = 0                           ## number of simulation steps
  phystime0 = 0                         ## total clock time (in seconds) since beginning of simulation
  
  
  ## Set seed to ensure reproducibility of results
  set.seed(run)
  
  ## Consumption matrix, resources-by-consumers
  if (limit == 'Neutral')
    C = matrix(mu, NR, NS)
  if (limit == 'Specialists') {
    C = matrix(b, NR, NS)
    diag(C) = mu
  }
  if (limit == 'Generalists')
    C = matrix(runif(NR * NS, min = mu - b, max = mu + b), NR, NS)
  if (limit == 'Generalists_Gaussian')
    C = matrix(pmax(0, rnorm(
      NR * NS, mean = mu, sd = mu * b
    )), NR, NS)
  if (limit == 'Generalists_Uniform')
    C = matrix(pmax(0, runif(
      NR * NS,
      min = mu - mu * b * sqrt(3),
      max = mu + mu * b * sqrt(3)
    )), NR, NS)
  if (limit == 'ConstrainedBudget') {
    C = matrix(runif(NR * NS, min = mu - b, max = mu + b), NR, NS)
    C = t(t(C) / colSums(C) * NS * mu)
  }
  
  ## Resource supply and consumer mortality
  if (parametrization == 'FixedAbuns') {
    eta = epsilon * Rstar * colSums(C)			## consumer mortality
    rho = Rstar * Nstar * rowSums(C)				## resource supply rate
  }
  if (parametrization == 'FixedRates') {
    eta = epsilon * mu * NR * Rstar						## consumer mortality
    rho = eta * NS / NR * Nstar / epsilon				## resource supply rate
  }
  if (limit == 'Specialists') {
    rho = rep((mu + b * (NS - 1)) * Nstar * Rstar, NR)
    eta = rep(epsilon * (mu + b * (NR - 1)) * Rstar, NS)
  }
  if (parametrization == 'FixedAbuns_LogisticResources') {
    intrinsic_growth_rate = 2 * Nstar * rowSums(C)  ## resource intrinsic growth rate
    carrying_capacity = 2 * Rstar                 ## resource carrying capacity
    rho = rep(NA, NR)
  }
  
}

## If resuming simulation from a previous run
if (resume_simulation) {
  datafilename = paste0(datadir, formatC(simind, width = 3, flag = "0"), '.RData')
  dat = get(load(datafilename))
  dat = dat[[length(dat)]]
  list2env(as.list(dat), envir = globalenv())
  list2env(as.list(parms), envir = globalenv())
  RN = outer(resources, consumers)
  phystime0 = phystime
  scen = parms[, 1:13]
  
  ## Set seed to ensure reproducibility of results
  set.seed(run)
}

## If keeping track of extinctions
if (keep_track_extinctions) {
  tempo = current_tempo = 0
  biotime = 0
  extinctiontimes =
    tibble(
      species = seq_along(consumers),
      abundance = consumers,
      clocktime = tempo,
      biotime = biotime,
      simtime = simtime
    )
}

extinction_event = 0                  ## identifies extinction events
speciation_counter = 0                ## tracks speciation events
extinction_counter = 0                ## tracks extinction events

## Dynamics
t0 = as.numeric(Sys.time())
Sys.sleep(1)
while (as.numeric(Sys.time()) - t0 < 3600 * 24 * 7) {
  simtime = simtime + 1
  
  if (keep_track_extinctions)
    biotime = biotime + rexp(1, rate = sum(NR * mean(rho), sum(C * RN), sum(eta *
                                                                              consumers)))
  
  event = sample(
    c('supply', 'consume', 'die'),
    size = 1,
    prob = c(
      ifelse(
        parametrization == 'FixedAbuns_LogisticResources',
        sum(
          intrinsic_growth_rate * resources * pmax(0, 1 - resources / carrying_capacity)
        ),
        NR * mean(rho)
      ),
      sum(C * RN),
      sum(eta * consumers)
    )
  )
  
  if (event == 'supply') {
    if (parametrization != 'FixedAbuns_LogisticResources') {
      probabilities = rho
    } else{
      probabilities = intrinsic_growth_rate * resources * pmax(0, 1 - resources /
                                                                 carrying_capacity)
    }
    i = sample(seq_along(resources), size = 1, prob = probabilities)
    resources[i] = resources[i] + 1
    RN[i, ] = RN[i, ] + consumers
  }
  
  if (event == 'die') {
    j = sample(seq_along(consumers),
               size = 1,
               prob = eta * consumers)
    if (consumers[j] == 1) {
      extinction_event = 1
      extinction_counter = extinction_counter + 1
      # writeLines(paste('extinction',j,sum(consumers>0)-1))
      # print(consumers)
    }
    if (consumers[j] > 1 | nu > 0) {
      consumers[j] = consumers[j] - 1
      RN[, j] = RN[, j] - resources
    }
    # else if(bounce==FALSE & limit=='Generalists' & length(consumers)>1){
    #   consumers=consumers[-j]
    #   RN=RN[,-j]; C=C[,-j]
    # }
    if (extinction_event == 1 & keep_track_extinctions) {
      x = extinctiontimes %>% filter(species == j) %>% tail(1)
      if (x$abundance > 0) {
        extinctiontimes =
          rbind(
            extinctiontimes,
            tibble(
              species = j,
              abundance = 0,
              clocktime = tempo,
              biotime = biotime,
              simtime = simtime
            )
          )
      }
      extinction_event = 0
    }
  }
  
  if (event == 'consume') {
    pair = sample(seq_along(RN), size = 1, prob = as.numeric(C * RN))
    i = 1 + (pair - 1) %% length(resources)		## which resource got eaten
    j = 1 + (pair - 1) %/% length(resources)	## which consumer ate
    
    ## Resource depletion
    if (resources[i] > 1) {
      resources[i] = resources[i] - 1
      RN[i, ] = RN[i, ] - consumers
    }
    
    ## Consumer reproduction / Speciation
    if (runif(1) <= epsilon) {
      if (runif(1) <= nu & any(consumers == 0)) {
        speciation_counter = speciation_counter + 1
        j = if (sum(consumers == 0) > 1)
          sample(which(consumers == 0), size = 1)
        else
          which(consumers == 0)
        # writeLines(paste('speciation',j))
      }
      consumers[j] = consumers[j] + 1
      RN[, j] = RN[, j] + resources
    }
    
  }
  
  ## Update simulation clock time
  tempo = round(as.numeric(Sys.time()) - t0)
  
  ## Every 30 min (1800 sec) update extinctiontimes
  ## by adding the current state of the community
  if (keep_track_extinctions) {
    if (tempo - current_tempo >= 1800) {
      extinctiontimes =
        rbind(
          extinctiontimes,
          tibble(
            species = seq_along(consumers),
            abundance = consumers,
            clocktime = tempo,
            biotime = biotime,
            simtime = simtime
          )
        )
      current_tempo = tempo
    }
  }
  
  ## Save data (saves current state every hour if running on hpcc, rewriting previous save)
  if (tempo %% 1800 == 0 & hpcc == 1) {
    if (keep_track_extinctions) {
      save(extinctiontimes,
           file = paste0(
             datadir,
             'ExtinctionTimes',
             formatC(simind, width = 3, flag = "0"),
             '.RData'
           ))
    }
    if (save_parameters) {
      datafilename = paste0(datadir, formatC(simind, width = 3, flag = "0"), '.RData')
      com = list(
        list(
          parms = data.frame(
            scen,
            simtime = simtime,
            phystime = phystime0 + tempo,
            node = nodename
          ),
          rho = rho,
          eta = eta,
          C = C,
          resources = resources,
          consumers = consumers
        )
      )
      if (file.exists(datafilename)) {
        load(datafilename)
        data = c(data, com)
      } else
        data = com
      save(data, file = datafilename)
    }
    Sys.sleep(1)
  }
  
  ## Plot data (plots current SADs every 10 secs if running on PC)
  if (tempo %% 10 == 0 & hpcc == 0) {
    trimmed_resources = resources[resources > 0]
    trimmed_consumers = consumers[consumers > 0]
    wr = plyr::count(trimmed_resources) %>% mutate(cum = (sum(freq) - cumsum(freq) +
                                                            freq) / sum(freq))
    wc = plyr::count(trimmed_consumers) %>% mutate(cum = (sum(freq) - cumsum(freq) +
                                                            freq) / sum(freq))
    
    fit.r = fitSAD(data = trimmed_resources, dbn = 'Pois')
    fit.c = fitSAD(data = trimmed_consumers, dbn = 'LS')
    
    pr <-
      wr %>%
      ggplot(aes(x = x, y = cum)) +
      geom_point() +
      scale_x_log10() +
      geom_line(fit.r$dtf,
                mapping = aes(x = x, y = cum),
                color = red) +
      coord_cartesian(xlim = range(wr$x)) +
      xlab('Abundance') +
      ylab('Cumulative No. Species') +
      ggtitle('Resources') +
      theme(legend.position = 'none') +
      annotate(
        'text',
        x = min(trimmed_resources),
        y = 0,
        hjust = 0,
        vjust = 0,
        label = paste('Pois.pval =', round(fit.r$p.value, 3))
      ) +
      annotate(
        'text',
        x = Inf,
        y = Inf,
        hjust = 1,
        vjust = 1,
        label = paste('Res.mean = ', (round(
          mean(trimmed_resources)
        )))
      )
    
    pc <-
      wc %>%
      ggplot(aes(x = x, y = cum)) +
      geom_point() +
      scale_x_log10() +
      geom_line(fit.c$dtf,
                mapping = aes(x = x, y = cum),
                color = blue) +
      coord_cartesian(xlim = range(wc$x)) +
      xlab('Abundance') +
      ylab('Cumulative No. Species') +
      ggtitle('Consumers') +
      theme(legend.position = 'none') +
      annotate(
        'text',
        min(trimmed_consumers),
        0,
        hjust = 0,
        vjust = 0,
        label = paste('LS.pval =', round(fit.c$p.value, 3))
      ) +
      annotate(
        'text',
        Inf,
        Inf,
        hjust = 1,
        vjust = 1,
        label = paste('Cons.mean = ', (round(
          mean(trimmed_consumers)
        )))
      )
    
    gridExtra::grid.arrange(pr,
                            pc,
                            nrow = 1,
                            top = paste0(limit, '\nT = ', round(simtime / 1e3), 'k'))
    
    Sys.sleep(1)
  }
}
