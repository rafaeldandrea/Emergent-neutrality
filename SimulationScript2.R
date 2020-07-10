## Libraries
library(tidyverse)    ## for data wrangling

## Load custom functions
dummy <- {
  list.files(
    path = '~/R_Functions/',
    pattern = "[.][R]$",
    full.names = TRUE,
    ignore.case = TRUE,
    recursive = TRUE
  ) %>%
  lapply(source)
}

## Read scenario from either manual input or command line from Bash file
simind <- make_index(986)

## Set working directory
setwd('~/EmergentNeutrality/')

## Resume simulation (TRUE) or start from scratch (FALSE)
resume_simulation = FALSE

## Keep track of species extinction times
keep_track_extinctions = TRUE

## Set folder to save data files in
datadir = '~/EmergentNeutrality/Data/20200416/'

## Name of data file
datafilename <- paste0(datadir, formatC(simind, width = 4, flag = "0"), '.RData')

## Indicator variable; 1 if running on High Power Computing cluster, 0 if running on my PC.
hpcc = as.numeric((nodename <- as.character(Sys.info()['nodename'])) != 'RAF-PC')

## Parameters of simulation
parms <- 
  tibble(
    algorithm = 'Tau-leaping',
    max_years = 1e6,
    tau = 5e-2,                ## average fraction of the resource-consumer community that turns over every year
    census_interval = 1e2
  )

list2env(parms, envir = environment())


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
        arrange(Nstar, mu, b)
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
        arrange(Nstar, mu, b)
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
        arrange(NR, NS, mu, b) 
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
      ) 
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
        mutate(b = mu * b)
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
      )
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
        arrange(NS, mu)
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
      ) 
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
      ) 
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
      )
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
      ) 
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
      ) 
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
      ) 
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
      )
    }
    
    benchmarking <- {
      crossing(
        limit = 'Benchmarking',
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
        Nstar = c(100,10),
        NS = c(25,50,100,200,500,1000,2000),
        mu = 1,
        b = 1,
        run = 1
      ) %>%
        mutate(NR = NS) %>%
        arrange(desc(Nstar))
    }
    
    specialists2 <- {
      crossing(
        limit = 'Specialists',
        parametrization = 'FixedAbuns',
        bounce = TRUE,
        epsilon = 1,
        nu = 0,
        Nstar = c(10,50,100),
        NS = c(25, 50, 100, 200, 400, 800, 1600),
        b = 1,
        run = 1
      ) %>%
      mutate(
        Rstar = Nstar,
        NR = NS
      ) %>%
      left_join(
        tibble(NS = c(25, 50, 100, 200, 400, 800, 1600)) %>%
          #ddply(.(NS), function(.) mu = seq(1.1, 1.5335 + 0.1653 * .$NS, l = 25)) %>% 
          ddply(.(NS), function(.) mu = pmax(1, qnorm(1:30/31, mean = 5/2 + 1/40 * .$NS, sd = 5))) %>% 
          pivot_longer(cols = -NS, values_to = 'mu', names_to = 'bla') %>%
          select(-bla),
          by = 'NS'
      ) %>%
      unique %>%
      arrange(NS, Nstar, mu, b) 
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
        union(benchmarking) %>%
        union(specialists2) %>%
        rowid_to_column(var = 'scenario')
    } 
  }
  
  ## subset scenarios according to simind
  scen <- filter(scenarios, scenario == simind)
  
  ## Assign each column of scen to a variable matching its name
  list2env(scen, envir = environment())
  
  ## Analytic prediction for variance of log-series (aka logarithmic) distribution
  var_logseries <- function(nstar) {
    p = uniroot(function(x) - x / (1 - x) / log(1 - x) - nstar, c(1e-10, 1 - 1e-10))$root
    - (p ^ 2 + p * log(1 - p)) / (1 - p) ^ 2 / log(1 - p) ^ 2
  } 
  varlogseries <- var_logseries(Nstar)
  
  ## Set seed to ensure reproducibility of results
  set.seed(run)
  
  ## Consumption matrix [resources-by-consumers]
  {
    if (limit %in% c('Neutral', 'Benchmarking'))
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
  }
  
  ## Resource supply and consumer mortality
  {
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
    if(limit == 'Benchmarking'){
      rho = (mu + b * (NS - 1)) * Nstar * Rstar
      eta = epsilon * (mu + b * (NR - 1)) * Rstar
    }
  }
  
  ## Initial conditions
  {
    resources = rep(round(Rstar), NR)			## initial resource abundances
    
    consumers = rep(round(Nstar), NS)			## initial consumer abundances
    
    RN = outer(resources, consumers)				## R*N (abundance products)
    
    simtime = 0                           ## number of simulation steps (Gillespie)
    
    phystime0 = 0                         ## total clock time (in seconds) since beginning of simulation (Gillespie)
    
    extinction_event = 0                  ## identifies extinction events
    
    speciation_counter = 0                ## tracks speciation events
    
    extinction_counter = 0                ## tracks extinction events
    
    next_census = census_interval         ## how many sim steps to skip before censusing
    
    years = 0                             ## how many years since start of simulation (tau-leaping)
    
    tot_events = 0                        ## how many supply + consume + die events since start of simulation (tau-leaping)
    
    t0 = as.numeric(Sys.time())
    
    benchmark <-
      tibble(
        simulation_hours = 0, 
        community_years = 0,
        num_events = 0,
        sad_variance = 0, 
        logseries_fit_statistic = NA
      )
    
    data <-
      list(
        parms = cbind(scen, parms),
        rho = rho,
        eta = eta,
        C = C,
        benchmark = benchmark
      )
    
    resources_df <- 
      tibble(
        resource_id = 1:NR,
        resources = resources
      ) %>%
      pivot_wider(
        names_from = resource_id, 
        values_from = resources, 
        names_prefix = 'resource_'
      )
    
    consumers_df <- 
      tibble(
        consumer_id = 1:NR,
        consumers = consumers
      ) %>%
      pivot_wider(
        names_from = consumer_id, 
        values_from = consumers, 
        names_prefix = 'species_'
      )
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


## Dynamics
if(algorithm == 'Gillespie'){ 
  while (as.numeric(Sys.time()) - t0 < 3600 * 24) {
  
    simtime = simtime + 1
    
    if (keep_track_extinctions)
      biotime = biotime + rexp(1, rate = sum(NR * mean(rho), sum(C * RN), sum(eta * consumers)))
    
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
      if(limit == 'Benchmarking'){ 
        i = seq(NR)[runif(1,1,NR+1)]
      }else{
        if (parametrization != 'FixedAbuns_LogisticResources') {
          probabilities = rho
        } else{
          probabilities = intrinsic_growth_rate * resources * pmax(0, 1 - resources / carrying_capacity)
        }
        i = sample(seq_along(resources), size = 1, prob = probabilities)
      }
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
      }
      
      if (consumers[j] > 1 | nu > 0) {
        consumers[j] = consumers[j] - 1
        RN[, j] = RN[, j] - resources
      }
      
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
      if (keep_track_extinctions & limit != 'Benchmarking') {
        save(extinctiontimes,
             file = paste0(
               datadir,
               'ExtinctionTimes',
               formatC(simind, width = 3, flag = "0"),
               '.RData'
             ))
      }
      if (save_parameters & limit != 'Benchmarking') {
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
      
      if(limit == 'Benchmarking') save(benchmark, file = paste0(datadir, 'benchmark_', formatC(simind, width = 3, flag = "0"), '.RData'))
      
      Sys.sleep(1)
    }
    
    ## Plot data (plots current SADs every 10 secs if running on PC)
    if (tempo %% 10 == 0) {
      
      trimmed_resources = resources[resources > 0]
      trimmed_consumers = consumers[consumers > 0]
      wr = plyr::count(trimmed_resources) %>% mutate(cum = (sum(freq) - cumsum(freq) + freq) / sum(freq))
      wc = plyr::count(trimmed_consumers) %>% mutate(cum = (sum(freq) - cumsum(freq) + freq) / sum(freq))
      
      fit.r = fitSAD(data = trimmed_resources, dbn = 'Pois')
      fit.c = fitSAD(data = trimmed_consumers, dbn = 'LS')
      
      foo <- c(tempo, var(consumers)/varlogseries, as.numeric(fit.c$statistic), simtime)
      
      benchmark <-
        benchmark %>%
        rbind(foo)
      
      print(foo, quote = FALSE)
      
      if(!hpcc){
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
              label = paste('Cons.mean = ', (round(mean(trimmed_consumers))))
            )
        
        gridExtra::grid.arrange(pr,
                                pc,
                                nrow = 1,
                                top = paste0(limit, '\nT = ', round(simtime / 1e3), 'k'))
      }
      
      Sys.sleep(1)
    }
  }
}

if(algorithm == 'Tau-leaping'){
  
  while(years < max_years){
  
    num_events = rpois(1, lambda = tau * sum(resources, consumers))
    
    tot_events = tot_events + num_events
    
    years = years + 1
    
    events = sample(
      c('supply', 'consume', 'die'), 
      size = num_events, 
      replace = TRUE, 
      prob = c(NR * mean(rho), sum(C * RN), sum(eta * consumers))
    )
    
    num_supply = sum(events == 'supply')
    supplied = seq(NR)[runif(num_supply, 1, NR + 1)]
    
    num_consume = sum(events == 'consume')
    pair = sample(seq_along(RN), size = num_consume, replace = TRUE, prob = as.numeric(C * RN))
    pair_i = 1 + (pair - 1) %% NR		## which resource got eaten
    pair_j = 1 + (pair - 1) %/% NR	## which consumer ate
    
    num_die = sum(events == 'die')
    dead = sample(seq(NS), size = num_die, replace = TRUE, prob = consumers)
    
    delta_resources = sapply(seq(NR), function(i) sum(supplied == i) - sum(pair_i == i))
    delta_consumers = sapply(seq(NS), function(j) sum(pair_j == j) - sum(dead == j))
    
    resources = pmax(1, resources + delta_resources)
    consumers = pmax(1, consumers + delta_consumers)
    
    RN = outer(resources, consumers)
    
    ## Census and Save data 
    if (years == next_census) {
      
        next_census = next_census + census_interval
        
        tempo = round(as.numeric(Sys.time()) - t0)
        
        trimmed_consumers = consumers[consumers > 0]
        wc = plyr::count(trimmed_consumers) %>% mutate(cum = (sum(freq) - cumsum(freq) + freq) / sum(freq))
        fit.c = fitSAD(data = trimmed_consumers, dbn = 'LS')
        
        benchmark <- 
          benchmark %>% 
          rbind(
            c(
              simulation_hours = round(tempo / 3600, 3), 
              community_years = years,
              num_events = tot_events,
              sad_variance = var(consumers) / varlogseries, 
              logseries_fit_statistic = as.numeric(fit.c$statistic) 
            )
          )
        
        resources_df <- rbind(resources_df, resources)
        
        consumers_df <- rbind(consumers_df, consumers)
        
        data$benchmark <- benchmark
        
        data$resources <- resources_df
        
        data$consumers <- consumers_df
        
        if(hpcc) save(data, file = datafilename)
        
    }
    
    if(!hpcc & round(as.numeric(Sys.time()) - t0) %% 10 == 0){
      
      trimmed_resources = resources[resources > 0]
      trimmed_consumers = consumers[consumers > 0]
      wr = plyr::count(trimmed_resources) %>% mutate(cum = (sum(freq) - cumsum(freq) + freq) / sum(freq))
      wc = plyr::count(trimmed_consumers) %>% mutate(cum = (sum(freq) - cumsum(freq) + freq) / sum(freq))
      
      fit.r = fitSAD(data = trimmed_resources, dbn = 'Pois', method = 'CVM')
      fit.c = fitSAD(data = trimmed_consumers, dbn = 'LS', method = 'CVM')
      
      current <-
        c(
          simulation_hours = round(tempo / 3600, 3), 
          community_years = years,
          num_events = tot_events,
          sad_variance = var(consumers) / varlogseries, 
          logseries_fit_statistic = as.numeric(fit.c$statistic) 
        )
      
      print(current, quote = FALSE)
      
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
          label = paste('Cons.mean = ', (round(mean(trimmed_consumers))))
        )
      
      gridExtra::grid.arrange(
        pr,
        pc,
        nrow = 1,
        top = paste0(limit, '\nT = ', round(simtime / 1e3), 'k')
      )
      
      Sys.sleep(1)
      
    }
    
  } 
  
}