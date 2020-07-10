redo_ls_analysis = FALSE
do_resources = FALSE
do.pdf = FALSE

phystime_cutoff = 48 * 1800

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

## load data ------
#datafolder='20190520'
datafolder = '20190307'
setwd(paste0('~/emergentneutrality/data/', datafolder, '/'))
lf <-
  setdiff(list.files(pattern = '*.RData'),
          c(
            list.files(pattern = 'ExtinctionTimes*'),
            'LogSeriesFits.RData'
          ))
dat = do.call('c', lapply(lf, function(l)
  get(load(l))))

metadata <-
  dat %>%
  ldply(function(l)
    l$parms) %>%
  as_tibble %>%
  mutate(index = 1:n())

raw_consumers <-
  dat %>%
  ldply(function(l)
    c(consumers = l$consumers)) %>%
  mutate(index = 1:n()) %>%
  as_tibble

consumers <-
  raw_consumers %>%
  gather(key = 'key', value = 'x', -index) %>%
  group_by(index) %>%
  count(x) %>%
  mutate(cum = (sum(n) - cumsum(n) + n) / sum(n)) %>%
  rename(freq = n) %>%
  inner_join(metadata, by = 'index') %>%
  filter(limit == 'Generalists') %>%
  ungroup


if (redo_ls_analysis) {
  results <-
    consumers %>%
    ddply(.(index), function(.) {
      print(unique(.$index))
      mod = fitSAD(., dbn = 'LS', method = 'CVM')
      return(data.frame(
        p.value = mod$p.value,
        statistic = mod$statistic
      ))
    }) %>%
    mutate(LS = 1 * (p.value >= .05)) %>%
    inner_join(metadata) %>%
    as_tibble
}

if (!redo_ls_analysis) {
  results = get(load('LogSeriesFits.Rdata')) %>%
    as_tibble
}

results = results %>% filter(phystime > phystime_cutoff)

statistic_over_time <-
  function(theNstar, thecutoff) {
    results %>%
      filter(Nstar == theNstar & phystime >= thecutoff) %>%
      mutate(bmu = b / mu) %>%
      ggplot(aes(
        phystime / 1800,
        statistic,
        color = (LS == 1),
        group = (LS == 1)
      )) +
      geom_point() +
      theme(legend.position = 'none') +
      facet_wrap(mu ~ bmu, scale = 'free', ncol = 10)
  }

phystime_cutoff <-
  tibble(
    Nstar = c(50, 100, 250, 500, 750, 1000, 1100),
    phystime_cutoff = 1800 * c(0, 0, 0, 20, 70, 70, 70)
  )

results_summary <-
  results %>%
  inner_join(phystime_cutoff, by = 'Nstar') %>%
  filter(phystime >= phystime_cutoff) %>%
  group_by(Nstar, mu, b) %>%
  summarize_at('LS', funs(mean, se = sd(.) / sqrt(n()))) %>%
  ungroup() %>%
  mutate(bmu = b / mu)

results_summary2 <-
  results %>%
  mutate(bmu = round(b / mu, 1)) %>%
  inner_join(phystime_cutoff, by = 'Nstar') %>%
  filter(phystime >= phystime_cutoff) %>%
  group_by(Nstar, bmu) %>%
  summarize_at('LS', funs(mean, se = sd(.) / sqrt(n()))) %>%
  ungroup()

logit <-
  results %>%
  mutate(bmu = round(b / mu, 1)) %>%
  inner_join(phystime_cutoff, by = 'Nstar') %>%
  filter(phystime >= phystime_cutoff) %>%
  dlply(.(Nstar), function(.) {
    glm(LS ~ bmu, data = ., family = 'binomial')
  })

bmus = seq(.1, 1, l = 1000)
predicted <-
  logit %>%
  ldply(function(mod)
    predict(mod, data.frame(bmu = bmus), type = 'response')) %>%
  as_tibble %>%
  gather(key = 'bmu', value = 'ProbLS', -Nstar) %>%
  arrange(Nstar) %>%
  plyr::mutate(bmu = bmus)

significance_cutoff = .5
bthresh <-
  logit %>%
  ldply(function(mod) {
    intercept = as.numeric(coefficients(mod)[1])
    slope = as.numeric(coefficients(mod)[2])
    s_intercept = as.numeric(coefficients(summary(mod))[1, 2])
    s_slope = as.numeric(coefficients(summary(mod))[2, 2])
    b.threshold = (log(significance_cutoff / (1 - significance_cutoff)) -
                     intercept) / slope
    c(
      intercept = intercept,
      slope = slope,
      b.threshold = b.threshold,
      s_b.threshold = abs(1 / slope) * sqrt((s_intercept) ^ 2 + (b.threshold *
                                                                   s_slope) ^ 2)
    )
  })

model <-
  lm(logb ~ logN,
     data = bthresh %>%
       mutate(logN = log(Nstar),
              logb = log(b.threshold)))
if (datafolder != '20190520') {
  model_slope = round(coef(model)[2], 2)
  model_se = round(coef(summary(model))[2, 2], 2)
}

sad_data <-
  consumers %>%
  filter(Nstar == 500 & run == 1) %>%
  mutate(cv = round(b / sqrt(3) / mu, 2)) %>%
  filter(phystime == max(phystime)) %>%
  filter(mu == min(mu)) %>%
  filter(cv %in% range(cv))

sad_fit <-
  sad_data %>%
  filter(cv == max(cv)) %>%
  ddply(.(cv), function(.)
    fitSAD(., dbn = 'LS', method = 'CVM')$dtf) %>%
  ungroup %>%
  as_tibble %>%
  filter(x <= max(sad_data$x))

sad_fit_stats <-
  sad_data %>%
  ddply(.(cv), function(.) {
    foo = fitSAD(., dbn = 'LS', method = 'CVM')
    return(tibble(
      p.value = foo$p.value,
      statistic = foo$statistic
    ))
  }) %>%
  ungroup %>%
  as_tibble

## Plots
pA <-
  sad_data %>%
  ggplot(aes(x, cum, group = cv, color = factor(cv))) +
  labs(color = 'CV') +
  geom_line(aes(x, cum), data = sad_fit, color = 'black') +
  geom_point() +
  theme_bw() +
  theme(
    legend.justification = c(0, 0),
    legend.position = c(.05, .05),
    legend.background = element_rect(color = blue)
  ) +
  scale_x_log10() +
  scale_y_log10() +
  xlab('Abundance') +
  ylab('Cumulative frequency') +
  ggtitle(label = 'A')



p1 <-
  results_summary %>%
  ggplot(aes(bmu, mean, color = as.factor(mu), group = as.factor(mu))) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                color = red,
                width = .02) +
  geom_point() +
  facet_wrap( ~ Nstar, scale = 'free') +
  theme_bw() +
  ggtitle(label = 'B')

p1b <-
  results_summary2 %>%
  ggplot(aes(bmu, mean)) +
  geom_line(aes(bmu, ProbLS), data = predicted, color = blue) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                color = red,
                width = .02) +
  geom_point() +
  facet_wrap( ~ Nstar, scale = 'free', ncol = 4) +
  theme_bw() +
  xlab(expression(paste('b / ', mu))) +
  ylab('Probability of log-series SAD') +
  ggtitle(label = 'B')

thenstar = 500
pB <-
  results_summary2 %>%
  filter(Nstar == thenstar) %>%
  ggplot(aes(bmu / sqrt(3), mean)) +
  geom_line(
    aes(bmu / sqrt(3), ProbLS),
    data = predicted %>% filter(Nstar == thenstar),
    color = blue
  ) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                color = red,
                width = .02) +
  geom_point() +
  theme_bw() +
  xlab('Consumption matrix CV') +
  ylab('Probability of logseries SAD') +
  ggtitle(label = 'B')


pC <-
  bthresh %>%
  ggplot(aes(log(Nstar), log(b.threshold))) +
  stat_smooth(method = 'lm') +
  geom_errorbar(aes(
    ymin = log(b.threshold - s_b.threshold),
    ymax = log(b.threshold + s_b.threshold)
  ),
  color = red,
  width = .1) +
  geom_point() +
  theme_bw() +
  annotate('text',
           x = 5.9,
           y = .7,
           label = as.expression(bquote(R ^ 2 == .(
             round(summary(model)$adj.r.squared, 2)
           )))) +
  annotate(
    'text',
    x = 5.9,
    y = .4,
    label = paste('slope =', model_slope, '\u00B1', model_se)
  ) +
  xlab('log(Mean abundance)') +
  ylab(expression(paste('log(', CV ^ threshold, ')'))) +
  ggtitle(label = 'C')

gridExtra::grid.arrange(pA, pB, pC, ncol = 3, respect = TRUE)

results %>%
  filter(Nstar == min(Nstar)) %>%
  mutate(bmu = round(b / mu, 1)) %>%
  ggplot(aes(phystime / 3600, statistic, group = LS, color = LS)) +
  geom_point() +
  facet_wrap(mu ~ bmu, ncol = 10, scale = 'free') +
  theme_bw() +
  theme(legend.position = 'none')

if (do.pdf) {
  plotname = 'Generalist_SAD_results.pdf'
  if (file.exists(plotname))
    unlink(plotname)
  ggsave(
    plotname,
    plot = gridExtra::grid.arrange(pA, pB, pC, ncol = 3, respect = TRUE),
    width = 21,
    height = 7,
    units = 'cm',
    path = '~/emergentneutrality/figures'
  )
  
  plotname = 'Fig2_SAD_results.pdf'
  if (file.exists(plotname))
    unlink(plotname)
  ggsave(
    plotname,
    plot = gridExtra::grid.arrange(pA, pB, pC,
                                   pD, pE, pF,
                                   ncol = 3,
                                   respect = TRUE),
    width = 24,
    height = 16,
    units = 'cm',
    path = '~/emergentneutrality/figures'
  )
}


if (do_resources) {
  raw_resources <-
    dat %>%
    ldply(function(l)
      c(resources = l$resources)) %>%
    mutate(index = 1:n()) %>%
    as_tibble
  
  resources <-
    raw_resources %>%
    gather(key = 'key', value = 'x', -index) %>%
    group_by(index) %>%
    count(x) %>%
    mutate(cum = (sum(n) - cumsum(n) + n) / sum(n)) %>%
    rename(freq = n) %>%
    inner_join(metadata, by = 'index') %>%
    filter(limit == 'Generalists') %>%
    ungroup
  
  sad_data <-
    resources %>%
    filter(Nstar == 500 & mu == min(mu)) %>%
    mutate(cv = round(b / sqrt(3) / mu, 2)) %>%
    filter(phystime == max(phystime))
  
  sad_fit <-
    sad_data %>%
    ddply(.(mu, cv), function(.) {
      fitSAD(., dbn = 'Pois', method = 'CVM')$dtf %>%
        filter(x >= min(.$x) & x <= max(.$x))
    }) %>%
    ungroup %>%
    as_tibble
  
  sad_fit_stats <-
    sad_data %>%
    ddply(.(mu, cv), function(.) {
      foo = fitSAD(., dbn = 'Pois', method = 'CVM')
      return(tibble(
        p.value = foo$p.value,
        statistic = foo$statistic
      ))
    }) %>%
    ungroup %>%
    as_tibble
  
  pR <-
    sad_data %>%
    ggplot(aes(x, cum)) +
    labs(color = 'CV') +
    geom_line(aes(x, cum), data = sad_fit, color = red) +
    geom_point() +
    theme_bw() +
    xlab('Abundance') +
    ylab('Cumulative frequency') +
    xlim(range(sad_data$x)) +
    facet_wrap( ~ cv, ncol = 5, scale = 'free')
  
  plotname = 'Generalist_SAD_results_Resources.pdf'
  if (file.exists(plotname))
    unlink(plotname)
  ggsave(
    plotname,
    plot = pR,
    width = 26,
    height = 10,
    units = 'cm',
    path = '~/emergentneutrality/figures'
  )
  
}
