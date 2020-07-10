## Analyze the extinction times of species

do.pdf = FALSE

library(gridExtra)
## Load my functions ------
dummy <-
  list.files(
    path = '~/R_Functions/',
    pattern = "[.][R]$",
    full.names = TRUE,
    ignore.case = TRUE,
    recursive = TRUE
  ) %>%
  lapply(source)


## set work dir
## 20190307: res-cons model, GENERALISTS scenarios, 500 species, 5 mu values, 10 b values, 1 run
## 20190514: res-cons model, NEUTRAL scenarios, 500 species, 10 runs
## 20190507: res-cons model, SPECIALISTS scenarios, 50 species, 10 Cd values, 1 run
## 20190520: res-cons model, GENERALISTS scenarios, N*=500, NS=50, 10 runs
## 20190530: res-cons model, SPECIALISTS scenarios, 50 species, 5 Cd values, 10 runs, lower Cd values

datadir = '20190530'
setwd(paste0('~/emergentneutrality/data/', datadir, '/'))

## set data treatment.
## simple: collect one extinction per species, counting from the beginning of the simulation.
## all: collect all exticntions
## correction: collect all extinctions by abundance group excluding groups where some
## species remain throughout the observation window. This removes biasing the estimated
## extinction times towards lower values.
sample_all = TRUE
sample_correction = FALSE
sample_simple = FALSE


## Define parameter table
if (datadir == '20190307') {
  ## declare scenarios to analyze
  index_dtf = tibble(index = 601:650)
  
  ## cutoffs for data reading
  burnin_cutoff = 0
  maxtime_cutoff = 1e9
  abundance_cutoff = 1e9
  
  Nstars = c(100, 1800, 250, 50, 30, 500, 750, 1000, 850, 950, 1100)
  scenarios <-
    expand.grid(
      limit = c('Generalists', 'ConstrainedBudget'),
      ## equalize resource budget of all consumers?
      parametrization = c('FixedRates', 'FixedAbuns'),
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
      Nstar = Nstars,
      ## expected consumer abun
      mu = seq(.05, 2, l = 5),
      ## mean consumption rate in C matrix
      b = seq(0, 1, l = 11),
      ## semirange of C matrix
      run = 1                                              ## controls seed of random number generator
    ) %>%
    mutate(limit = as.character(limit), b = b * mu) %>%
    arrange(parametrization, match(Nstar, Nstars), desc(limit), mu, b, run) %>%
    filter(b > 0) %>%
    filter(limit == 'Generalists' |
             (limit == 'ConstrainedBudget' & Nstar == 100)) %>%
    filter(parametrization == 'FixedAbuns' |
             !(Nstar > 250 & Nstar < 1800)) %>%
    mutate(index = seq(n())) %>%
    as_tibble
}
if (datadir == '20190514') {
  ## declare scenarios to analyze
  index_dtf = tibble(index = 471:480)
  
  ## cutoffs for data reading
  burnin_cutoff = 5 * 3600
  maxtime_cutoff = 1e9
  abundance_cutoff = 1e9
  
  scenarios <-
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
      NR = 50,
      ## number of resources
      NS = 50,
      ## number of consumers
      Rstar = 100,
      ## expected resource abun
      Nstar = 500,
      ## expected consumer abun
      mu = .05,
      ## mean consumption rate in C matrix
      b = 0,
      ## semirange of C matrix
      run = 1:10                                              ## controls seed of random number generator
    ) %>%
    mutate(index = 471:480)
}
if (datadir == '20190520') {
  ## declare scenarios to analyze
  index_dtf = tibble(index = 481:530)
  
  ## cutoffs for data reading
  burnin_cutoff = 10 * 3600
  maxtime_cutoff = 1e9
  abundance_cutoff = 1e9
  
  scenarios <-
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
      NR = NS,
      mu = .05,
      b = mu * seq(.1, 1, l = 5),
      run = 1:10
    ) %>%
    mutate(index = index_dtf$index)
}
if (datadir == '20190530') {
  ## declare scenarios to analyze
  index_dtf = tibble(index = c(531:580))
  
  ## cutoffs for data reading
  burnin_cutoff = 10 * 3600
  maxtime_cutoff = 1e9
  abundance_cutoff = 1e9
  
  scenarios <-
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
      NR = NS,
      # mu=c(seq(1.1,1.5335+0.1653*NS,l=5),1+seq(0,.08,l=5)),
      mu = seq(1.1, 1.5335 + 0.1653 * NS, l = 5),
      b = 1,
      run = 1:10
    ) %>%
    mutate(index = index_dtf$index)
}

if (TRUE) {
  ## analyze extinction data ---------
  
  if (sample_all) {
    dtf <-
      index_dtf %>%
      ddply(.(index), function(.) {
        char = .$index
        print(char)
        dat = get(load(
          paste0(
            '~/emergentneutrality/data/',
            datadir,
            '/ExtinctionTimes',
            char,
            '.Rdata'
          )
        ))
        pardat = get(load(
          paste0(
            '~/emergentneutrality/data/',
            datadir,
            '/',
            char,
            '.Rdata'
          )
        ))[[1]]
        pardtf <-
          tibble(
            species = seq_along(pardat$rho),
            rho = pardat$rho,
            eta = pardat$eta
          )
        res <-
          dat %>%
          filter(clocktime >= burnin_cutoff) %>%
          filter(clocktime <= maxtime_cutoff) %>%
          filter(abundance < abundance_cutoff) %>%
          ddply(.(species), function(df) {
            if (any(df$abundance == 0)) {
              y = df %>%
                filter(abundance == 0) %>%
                mutate(slot = seq(n()))
              
              x = df %>%
                filter(abundance != 0 & simtime < max(y$simtime)) %>%
                mutate(slot = 1 + findInterval(simtime, y$simtime)) %>%
                inner_join(y, by = 'slot') %>%
                mutate(
                  delta_biotime = biotime.y - biotime.x,
                  delta_simtime = simtime.y - simtime.x,
                  delta_clocktime = clocktime.y - clocktime.x
                )
            }
          }) %>%
          as_tibble %>%
          transmute(
            species,
            n0 = abundance.x,
            biotime = delta_biotime,
            clocktime = delta_clocktime,
            simtime = delta_simtime
          ) %>%
          inner_join(pardtf, by = 'species')
        return(res)
      }) %>%
      inner_join(scenarios, by = 'index') %>%
      mutate(bmu = if (datadir != '20190530')
        round(b / mu, 1)
        else
          mu) %>%
      as_tibble
  }
  
  if (sample_correction) {
    dtf <-
      index_dtf %>%
      ddply(.(index), function(.) {
        char = .$index
        dat = get(load(
          paste0(
            '~/emergentneutrality/data/',
            datadir,
            '/ExtinctionTimes',
            char,
            '.Rdata'
          )
        )) %>%
          mutate(window = clocktime %/% 1800,
                 group = findInterval(abundance, vec = exp(seq(
                   log(1), log(max(1 + abundance))
                 )))) %>%
          group_by(window, species) %>%
          mutate(window_permanent = ifelse(n() > 1, 0, 1))
        
        pardat = get(load(
          paste0(
            '~/emergentneutrality/data/',
            datadir,
            '/',
            char,
            '.Rdata'
          )
        ))[[1]]
        pardtf = tibble(species = 1:50,
                        rho = pardat$rho,
                        eta = pardat$eta)
        
        unbiased_dat <-
          dat %>%
          ddply(.(window), function(df) {
            ddply(subset(df, window_permanent == 1), .(species), function(.) {
              y = subset(dat,
                         window > unique(.$window) & species == unique(.$species))
              return(c(extinct = 1 * any(y$window_permanent == 0)))
            })
          }) %>%
          left_join(dat, by = c('window', 'species')) %>%
          group_by(window, group) %>%
          filter(all(extinct == 1)) %>%
          ungroup %>%
          full_join(dat %>%
                      filter(abundance == 0) %>%
                      mutate(extinct = 1)) %>%
          as_tibble
        
        res <-
          unbiased_dat %>%
          ddply(.(species), function(df) {
            if (any(df$abundance == 0)) {
              y = df %>%
                filter(abundance == 0) %>%
                mutate(slot = seq(n()))
              
              x = df %>%
                filter(abundance != 0 & simtime < max(y$simtime)) %>%
                mutate(slot = 1 + findInterval(simtime, y$simtime)) %>%
                inner_join(y, by = 'slot') %>%
                mutate(
                  delta_biotime = biotime.y - biotime.x,
                  delta_simtime = simtime.y - simtime.x,
                  delta_clocktime = clocktime.y - clocktime.x
                )
            }
          }) %>%
          as_tibble %>%
          transmute(
            species,
            n0 = abundance.x,
            biotime = delta_biotime,
            clocktime = delta_clocktime,
            simtime = delta_simtime
          ) %>%
          inner_join(pardtf, by = 'species')
        return(res)
      }) %>%
      inner_join(scenarios, by = 'index') %>%
      mutate(bmu = round(b / mu, 1)) %>%
      as_tibble
  }
  
  if (sample_simple) {
    dtf <-
      index_dtf %>%
      ddply(.(index), function(.) {
        char = .$index
        print(char)
        dat = get(load(
          paste0(
            '~/emergentneutrality/data/',
            datadir,
            '/ExtinctionTimes',
            char,
            '.Rdata'
          )
        ))
        pardat = get(load(
          paste0(
            '~/emergentneutrality/data/',
            datadir,
            '/',
            char,
            '.Rdata'
          )
        ))[[1]]
        pardtf <-
          tibble(species = 1:50,
                 rho = pardat$rho,
                 eta = pardat$eta)
        res <-
          dat %>%
          filter(clocktime >= burnin_cutoff) %>%
          filter(clocktime <= maxtime_cutoff) %>%
          filter(abundance < abundance_cutoff) %>%
          ddply(.(species), function(df) {
            if (any(df$abundance == 0)) {
              y = df %>%
                filter(abundance == 0) %>%
                head(1)
              
              x = df %>%
                filter(abundance != 0) %>%
                head(1) %>%
                mutate(
                  delta_biotime = y$biotime - biotime,
                  delta_simtime = y$simtime - simtime,
                  delta_clocktime = y$clocktime - clocktime
                )
            }
          }) %>%
          as_tibble %>%
          transmute(
            species,
            n0 = abundance,
            biotime = delta_biotime,
            clocktime = delta_clocktime,
            simtime = delta_simtime
          ) %>%
          inner_join(pardtf, by = 'species')
        return(res)
      }) %>%
      inner_join(scenarios, by = 'index') %>%
      mutate(bmu = round(b / mu, 5)) %>%
      as_tibble
  }
}

dtf %>%
  filter(index == index_dtf$index[1]) %>%
  ggplot(aes(x = n0 / eta, y = clocktime / 3600)) +
  geom_point() +
  stat_smooth(method = 'lm') +
  #scale_x_log10() +
  #scale_y_log10() +
  #facet_wrap(~round(b/mu,1),ncol=5) +
  theme(aspect.ratio = 1) +
  xlab(expression(paste(N[0], ' / ', eta))) +
  ylab('Extinction Times') +
  ggtitle(' N* = 500  Same random number seeds')

dtf %>%
  mutate(predictor = n0 / eta) %>%
  ddply(.(bmu), function(.) {
    coef(summary(lm(biotime ~ predictor - 1, data = .)))
  }) %>%
  ggplot(aes(bmu, Estimate)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab(expression(paste('b / ', mu))) +
  ylab(expression(paste(
    'slope of lm(extinction time ~ ', N[0], ' / ', eta, ')'
  )))


v = uniroot(function(v)
  - 1 / log(v) * (1 - v) / v - unique(dtf$Nstar),
  interval = c(1e-10, 1 - 1e-10))$root
foo = tibble(n0 = seq(max(dtf$n0)),
             H = cumsum(1 / (1:max(dtf$n0))))

dtf_new <-
  dtf %>%
  mutate(eta = if (!datadir %in% c('20190507', '20190530'))
  {
    epsilon * Rstar * NR * mu
  } else{
    epsilon * Rstar * (mu + b * (NR - 1))
  },
  B = pbeta(1 - v, 1 + n0, 1e-20) * beta(1 + n0, 1e-20)) %>%
  left_join(foo, by = 'n0') %>%
  mutate(T = 1 / eta * 1 / v * ((1 / (1 - v)) ^ n0 * B + H + log(v)))

dtf_new %>%
  filter(mu == .05 & bmu == .1) %>%
  ggplot(aes(T, biotime)) +
  geom_abline(intercept = 0,
              slope = 1,
              color = red) +
  geom_point() +
  stat_smooth(method = 'lm') +
  scale_x_log10() +
  scale_y_log10() +
  #facet_wrap(mu~bmu,ncol=10) +
  theme(aspect.ratio = 1) +
  xlab('Predicted T') +
  ylab('Observed T') +
  ggtitle(' N* = 500  Same random number seeds')

dtf_new %>%
  ggplot(aes(n0 / eta, T)) +
  geom_point() +
  stat_smooth(method = 'lm') +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(mu ~ bmu, ncol = 10) +
  theme(aspect.ratio = 1) +
  xlab('n0 / eta') +
  ylab('Predicted T') +
  ggtitle(' N* = 500  Same random number seeds')

dtf_new %>%
  ddply(.(bmu), function(.) {
    coef(summary(lm(biotime ~ T - 1, data = .)))
  }) %>%
  ggplot(aes(bmu, Estimate)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab(expression(paste('b / ', mu))) +
  ylab('slope of lm(Observed T ~ Predicted T)')


base = exp(1)

dtf_groups <-
  dtf_new %>%
  group_by(mu, bmu) %>%
  mutate(group = findInterval(n0, vec = base ^ (seq(
    log(min(n0, base)), log(max(1 + n0, base)), l = 15
  )))) %>%
  ungroup %>%
  group_by(mu, bmu, group) %>%
  summarize_at(c('n0', 'T', 'biotime', 'eta', 'clocktime'),
               funs(n(), mean, se = sd(.) / sqrt(n()))) %>%
  ungroup



p1 <-
  dtf_groups %>%
  filter(bmu >= 1.1 | bmu == 1) %>%
  mutate(bmu = round(bmu, 1)) %>%
  ggplot(aes(n0_mean, biotime_mean)) +
  geom_errorbar(
    aes(ymin = biotime_mean - biotime_se,
        ymax = biotime_mean + biotime_se),
    width = 0,
    color = red
  ) +
  geom_errorbarh(aes(xmin = n0_mean - n0_se, xmax = n0_mean + n0_se), color =
                   red) +
  geom_line(aes(n0_mean, T_mean), color = blue) +
  geom_point() +
  facet_wrap( ~ bmu, ncol = 10, scale = 'free') +
  theme_bw() +
  xlab('Initial abundance') +
  ylab('Observed Extinction Time') +
  ggtitle('Specialists')


gridExtra::grid.arrange(p1 + scale_x_log10() + scale_y_log10())


p2 <-
  dtf_groups %>%
  filter(mu >= 1.1 | bmu == 1) %>%
  mutate(bmu = round(bmu, 1)) %>%
  ggplot(aes(
    T_mean,
    biotime_mean,
    group = as.factor(mu),
    color = as.factor(mu)
  )) +
  geom_errorbar(
    aes(ymin = biotime_mean - biotime_se, ymax = biotime_mean + biotime_se),
    width = 0,
    color = 'grey'
  ) +
  geom_errorbarh(aes(xmin = T_mean - T_se, xmax = T_mean + T_se), color =
                   'grey') +
  geom_smooth(method = 'lm', color = blue, se = FALSE) +
  geom_abline(slope = 1,
              intercept = 0,
              color = red) +
  geom_point() +
  facet_wrap( ~ bmu, ncol = 6) +
  theme_bw() +
  
  theme(legend.position = 'none') +
  xlab('Predicted Extinction Time') +
  ylab('Observed Extinction Time') +
  ggtitle(' N* = 500  Same random number seeds')


gridExtra::grid.arrange(p2)

thexlab = if (datadir != '20190530')
  expression(paste('b /', mu))
else
  expression(C[d] / C[o])

res <-
  ddply(dtf_groups, .(mu, bmu), function(.) {
    coefs = coef(summary(lm(biotime_mean ~ T_mean - 1, data = .)))
    tibble(slope = coefs[1],
           se = coefs[2])
  }) %>%
  as_tibble


res %>%
  ggplot(aes(bmu, slope)) +
  geom_point() +
  stat_smooth(method = 'lm') +
  facet_wrap( ~ mu, scale = 'free') +
  theme_bw()


res2 <-
  ddply(dtf_groups, .(bmu), function(.) {
    coefs = coef(summary(lm(biotime_mean ~ T_mean - 1, data = .)))
    tibble(slope = coefs[1],
           se = coefs[2])
  }) %>%
  as_tibble

res2 %>%
  ggplot(aes(bmu, slope)) +
  stat_smooth(method = 'lm') +
  geom_errorbar(aes(ymax = slope + 1.96 * se, ymin = slope - 1.96 * se),
                color = red,
                width = .02) +
  geom_point() +
  ylab('Slope between Observed and Predicted ET') +
  xlab(thexlab) +
  annotate('text',
           x = .25,
           y = 2,
           label = paste('R2 =', round(summary(
             lm(slope ~ bmu, data = res2)
           )$adj.r.squared, 2))) +
  theme_bw()

## Testing for normality of residuals of linear regression
## between observed and predicted extinction times
dtf_groups %>%
  group_by(bmu) %>%
  summarize(p.value = shapiro.test(summary(lm(
    biotime_mean ~ T_mean
  ))$residuals)$p.value)

dtf_groups %>%
  group_by(bmu) %>%
  summarize(
    intercept = coef(summary(lm(
      biotime_mean ~ T_mean
    )))[1, 1],
    intercept_se = coef(summary(lm(
      biotime_mean ~ T_mean
    )))[1, 2],
    slope = coef(summary(lm(
      biotime_mean ~ T_mean
    )))[2, 1],
    slope_se = coef(summary(lm(
      biotime_mean ~ T_mean
    )))[2, 2]
  )

p3 <-
  dtf_groups %>%
  group_by(bmu) %>%
  summarize(Chi2 = sum((biotime_mean - T_mean) ^ 2 / T_mean),
            G = 2 * sum(biotime_mean * log(biotime_mean / T_mean))) %>%
  gather(key = 'statistic', value = 'value', -bmu) %>%
  ggplot(aes(bmu, value, group = statistic, color = statistic)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  xlab(thexlab) +
  ylab('Statistic') +
  ggtitle('Observed vs Expected Extinction Times') +
  theme(
    legend.justification = c(0, 1),
    legend.position = c(.05, .95),
    legend.background = element_rect(color = blue),
    aspect.ratio = 1
  )

p4 <-
  dtf_groups %>%
  filter(bmu >= 1.1 | bmu == 1) %>%
  group_by(bmu) %>%
  summarize(chi2 = sum((biotime_mean - T_mean) ^ 2 / T_mean)) %>%
  ggplot(aes(bmu, chi2)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  xlab(thexlab) +
  ylab(expression(chi ^ 2)) +
  ggtitle('Observed vs Expected Extinction Times')

pAD <-
  dtf_groups %>%
  mutate(bmu = if (datadir %in% c('20190507', '20190530'))
    round(bmu, 1)
    else
      round(bmu / sqrt(3), 2)) %>%
  filter(bmu %in% range(bmu)) %>%
  ggplot(aes(
    n0_mean,
    biotime_mean * eta_mean,
    group = factor(bmu),
    color = factor(bmu)
  )) +
  labs(color = ifelse(datadir %in% c('20190507', '20190530'), expression(C[d] /
                                                                           C[o]), 'CV')) +
  geom_errorbar(aes(
    ymin = (biotime_mean - biotime_se) * eta_mean,
    ymax = (biotime_mean + biotime_se) * eta_mean,
    color = factor(bmu)
  ),
  width = .02) +
  geom_errorbarh(aes(
    xmin = n0_mean - n0_se,
    xmax = n0_mean + n0_se,
    color = factor(bmu)
  )) +
  geom_line(aes(n0_mean, T_mean * eta_mean, color = factor(bmu))) +
  geom_point() +
  theme_bw() +
  theme(
    legend.justification = c(1, 0),
    legend.position = c(.95, .05),
    legend.background = element_rect(color = blue)
  ) +
  xlab('Initial abundance') +
  ylab('Species extinction time (generations)') +
  ggtitle(ifelse(datadir %in% c('20190530', '20190507'), 'D', 'A'))

pBE <-
  dtf_groups %>%
  #filter(mu>=1.1 | bmu==1) %>%
  mutate(bmu = if (datadir %in% c('20190507', '20190530'))
    round(bmu, 1)
    else
      round(bmu / sqrt(3), 2)) %>%
  ggplot(aes(
    T_mean * eta_mean,
    biotime_mean * eta_mean,
    color = factor(bmu),
    group = factor(bmu)
  )) +
  labs(color = ifelse(datadir %in% c('20190530', '20190507'), expression(C[d] /
                                                                           C[o]), 'CV')) +
  geom_errorbar(aes(
    ymin = (biotime_mean - biotime_se) * eta_mean,
    ymax = (biotime_mean + biotime_se) * eta_mean,
    color = as.factor(bmu)
  ),
  width = 0) +
  geom_errorbarh(aes(
    xmin = (T_mean - T_se) * eta_mean,
    xmax = (T_mean + T_se) * eta_mean
  ),
  color = 'grey') +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', se = FALSE) +
  #geom_line() +
  geom_point() +
  theme_bw() +
  theme(
    legend.justification = c(1, 0),
    legend.position = c(.98, .02),
    legend.background = element_rect(color = blue)
  ) +
  xlab('Predicted extinction time (generations)') +
  ylab('Observed extinction time (generations)') +
  ggtitle(ifelse(datadir %in% c('20190530', '20190507'), 'E', 'B'))

pCF <-
  res2 %>%
  # filter(bmu>=1.1 | bmu==1) %>%
  ggplot(aes(bmu, slope)) +
  geom_line() +
  geom_errorbar(
    aes(ymax = slope + se, ymin = slope - se),
    color = red,
    width = ifelse(datadir %in% c('20190530', '20190507'), .2, .02)
  ) +
  geom_point() +
  ylab('Slope: observed vs predicted extinction time') +
  xlab(ifelse(
    datadir %in% c('20190530', '20190507'),
    expression(C[d] / C[o]),
    'Consumption matrix CV'
  )) +
  theme_bw() +
  ggtitle(ifelse(datadir %in% c('20190530', '20190507'), 'F', 'C'))


if (do.pdf) {
  plotname = 'Fig4_ET_results.pdf'
  if (file.exists(plotname))
    unlink(plotname)
  ggsave(
    plotname,
    plot = gridExtra::grid.arrange(
      pA + scale_x_log10() + scale_y_log10(),
      pB,
      pC,
      pD + scale_x_log10() + scale_y_log10(),
      pE,
      pF,
      ncol = 3,
      respect = TRUE
    ),
    width = 27,
    height = 18,
    units = 'cm',
    path = '~/emergentneutrality/figures'
  )
}
