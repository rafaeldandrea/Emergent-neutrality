do_resources=FALSE
do.pdf=FALSE

## Load my functions ------
dummy <-
  list.files(path='~/R_Functions/',pattern = "[.][R]$",full.names=TRUE,ignore.case=TRUE,recursive=TRUE) %>%
  lapply(source)

## specify model scenario --------
datafolder='20190507'

##
redo_ls_analysis=TRUE

fitting_test='KS'

foo <-
  tibble(
    thelimit='Specialists',
    theNstar=100,
    theNS=c(25,10*(5:10),500,1000,2000),
    phystime_cutoff=c(rep(1800,7),rep(90*3600,3)),
    theparametrization='FixedAbuns',
    datafolder=datafolder
  )

if(redo_ls_analysis){
  dtf=NULL
  for(thecase in 2:7){
  
    setwd(paste0('~/emergentneutrality/data/',datafolder,'/'))
    list2env(foo[thecase,],envir=globalenv())
  
    if(thecase<=7){ 
      lf=paste0(401:410+10*(thecase-1),'.RData')
      dat=do.call('c',lapply(lf,function(l) get(load(l))))
      
      if(thecase==2){
        setwd(paste0('~/emergentneutrality/data/',datafolder,'/20190507_2/'))
        lf2=paste0(411:420,'.RData')
        dat2=do.call('c',lapply(lf,function(l) get(load(l))))
        dat=c(dat,dat2)
      }
      
      metadata <-
        dat %>%
        ldply(function(l) l$parms) %>% 
        as_tibble %>%
        mutate(index=1:n()) %>% 
        arrange(b,mu) %>%
        filter(limit==thelimit) %>%
        filter(Nstar==theNstar) %>%
        filter(NS==theNS) %>%
        filter(phystime>phystime_cutoff)
      
      raw_consumers <-
        dat %>%
        ldply(function(l) c(consumers=l$consumers)) %>% 
        mutate(index=1:n()) %>%
        as_tibble
    }
    
    if(thecase %in% 8:10){ 
        lf=paste0(581:590+10*(thecase-8),'.RData')
        metadata<-
          ddply(tibble(name=lf),.(name),function(.){
            foo=get(load(.$name))
            ldply(foo,function(v){
              v$parms
            })
          }) %>%
          as_tibble %>%
          select(-name) %>%
          mutate(index=1:n()) %>% 
          arrange(b,mu) %>%
          filter(limit==thelimit) %>%
          filter(Nstar==theNstar) %>%
          filter(NS==theNS) %>%
          filter(phystime>phystime_cutoff)
        
        raw_consumers<-
          ddply(tibble(name=lf),.(name),function(.){
            ldply(get(load(.$name)),function(v) consumers=v$consumers)
          }) %>%
          mutate(index=1:n()) %>%
          select(-name) %>%
          as_tibble
    }
    
    consumers <- 
      raw_consumers %>%
      gather(key='key',value='x',-index) %>%
      group_by(index) %>%
      count(x) %>%
      mutate(cum=(sum(n)-cumsum(n)+n)/sum(n)) %>%
      rename(freq=n) %>%
      inner_join(metadata,by='index') %>%
      ungroup()
    
    results <-
      consumers %>%
      ddply(.(index),function(.){ 
        mod=fitSAD(.,dbn='LS',method=fitting_test)
        return(data.frame(p.value=mod$p.value,statistic=mod$statistic))  
      }) %>%
      mutate(LS=1*(p.value>=.05)) %>%  
      inner_join(metadata) %>%  
      as_tibble
    
    dtf=rbind(dtf,results)
  }
  results=dtf
}
if(!redo_ls_analysis){ 
  setwd(paste0('~/emergentneutrality/data/',datafolder,'/'))
  results=get(load('LogSeriesFits.RData'))
}

if(fitting_test=='KS') results=filter(results,NS>30) 

thecase=2
list2env(foo[thecase,],envir=globalenv())

lf=paste0(401:410+10*(thecase-1),'.RData')
dat=do.call('c',lapply(lf,function(l) get(load(l))))

if(thecase==2){
  setwd(paste0('~/emergentneutrality/data/',datafolder,'/20190507_2/'))
  lf2=paste0(411:420,'.RData')
  dat2=do.call('c',lapply(lf,function(l) get(load(l))))
  dat=c(dat,dat2)
}

metadata <-
  dat %>%
  ldply(function(l) l$parms) %>% 
  as_tibble %>%
  mutate(index=1:n()) %>% 
  arrange(b,mu) %>%
  filter(limit==thelimit) %>%
  filter(Nstar==theNstar) %>%
  filter(NS==theNS) %>%
  filter(phystime>phystime_cutoff)

raw_consumers <-
  dat %>%
  ldply(function(l) c(consumers=l$consumers)) %>% 
  mutate(index=1:n()) %>%
  as_tibble

consumers <- 
  raw_consumers %>%
  gather(key='key',value='x',-index) %>%
  group_by(index) %>%
  count(x) %>%
  mutate(cum=(sum(n)-cumsum(n)+n)/sum(n)) %>%
  rename(freq=n) %>%
  inner_join(metadata,by='index') %>%
  ungroup


logit <-
  results %>%
  filter(phystime>phystime_cutoff) %>%
  dlply(.(NS),function(.){
    glm(LS~mu,data=.,family='binomial')
  })

sad_data <-
  consumers %>%
  filter(node=='rosalind') %>%
  filter(phystime==max(phystime)) %>%
  filter(mu %in% range(mu))

sad_fit <-
  sad_data %>%
  filter(mu==13) %>%
  ddply(.(mu),function(.)  fitSAD(.,dbn='LS',method=fitting_test)$dtf) %>%
  ungroup %>%
  as_tibble %>%
  filter(x<=max(sad_data$x))

sad_fit_stats <-
  sad_data %>%
  ddply(.(mu),function(.){ 
    foo=fitSAD(.,dbn='LS',method=fitting_test)
    return(tibble(p.value=foo$p.value,statistic=foo$statistic))
  }) %>%
  ungroup %>%
  as_tibble

significance_cutoff=.5 
Cd_thresh <-
  logit %>%
  ldply(function(mod){
    coefs=coefficients(summary(mod))
    intercept=as.numeric(coefs[1,1])
    slope=as.numeric(coefs[2,1])
    s_intercept=as.numeric(coefs[1,2])
    s_slope=as.numeric(coefs[2,2])
    cd.threshold=(log(significance_cutoff/(1-significance_cutoff))-intercept)/slope
    c(
      intercept=intercept,
      slope=slope,
      cd.threshold=cd.threshold,
      s_cd.threshold=abs(1/slope)*sqrt((s_intercept)^2+(cd.threshold*s_slope)^2)
    )
  }) %>%
  mutate(logc=log(cd.threshold),logn=log(NS))

model_log=lm(logc~logn,data=Cd_thresh)
model_log_slope=round(coef(model_log)[2],2)
model_log_se=round(coef(summary(model_log))[2,2],2)

model=lm(cd.threshold~NS,data=Cd_thresh)
model_slope=round(coef(model)[2],2)
model_se=round(coef(summary(model))[2,2],4)


themus=seq(1,max(results$mu),l=100)
df_logit_fit <-
  tibble(
    NS=rep(sort(unique(dtf$NS)),each=length(themus)),
    mus=rep(themus,length(logit)),
    predLS=unlist(llply(logit,function(l) as.numeric(predict(l,data.frame(mu=themus),type='response'))))
  )

results_summary <-
  results %>% 
  filter(phystime>phystime_cutoff) %>%
  group_by(NS,mu) %>% 
  summarize_at('LS',funs(n=n(),mean,se=sd(.)/sqrt(n()))) %>%
  ungroup

mod2=summary(lm(cd.threshold~NS,data=Cd_thresh))

## Plots
pD <- {
  sad_data %>%
  ggplot(aes(x,cum,group=mu,color=factor(mu))) +
  labs(color=expression(paste(C[d]/C[o]))) +
  geom_line(aes(x,cum),data=sad_fit,color='black') +
  geom_point() +
  theme_bw() +
  theme(
    legend.justification=c(0,0),
    legend.position=c(.05,.05),
    legend.background=element_rect(color=blue)
  ) +
  scale_x_log10() +
  scale_y_log10() +
  xlab('Abundance') +
  ylab('Cumulative frequency') +
  ggtitle(label='D')
}
  
p1 <- {
  results_summary %>%
  ggplot(aes(mu,mean,group=NS,color=factor(NS))) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.05) +
  xlim(c(1,10)) +
  geom_line(aes(mus,predLS),data=df_logit_fit) +
  xlab(expression(paste(C[d])))+
  ylab('Prob Logseries (logit)') +
  ylim(c(0,1)) +
  ggtitle(
    bquote(
      paste(
        .(thelimit),
        ' ',
        N[s],
        ' = ',
        .(theNS),
        '  Burn-in cutoff = ',
        .(phystime_cutoff/1800),
        '* 30 min'
      )
    )
  ) +
  theme(plot.title=element_text(size=15)) +
  theme_bw() +
  theme(
    legend.justification=c(1,1),
    legend.position=c(.95,.7),
    legend.background=element_rect(color=blue)
  )
}

pE <- {
  results_summary %>%
  filter(NS==50) %>%
  ggplot(aes(mu,mean)) +
  geom_line(aes(mus,predLS),data=df_logit_fit %>% filter(NS==50),color=blue) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.4,color=red) +
  geom_point() +
  xlim(c(1,10)) +
  xlab(expression(paste(C[d]/C[o])))+
  ylab('Probability of logseries SAD') +
  ylim(c(0,1)) +
  ggtitle('E') +
  theme_bw() +
  theme(
    legend.justification=c(1,1),
    legend.position=c(.95,.7),
    legend.background=element_rect(color=blue)
  )
}

pF<- {
  Cd_thresh %>% 
  ggplot(aes(NS,cd.threshold)) + 
  geom_smooth(method='lm') +
  geom_errorbar(
    aes(
      ymin=cd.threshold-s_cd.threshold,
      ymax=cd.threshold+s_cd.threshold
    ),
    color=red,
    width=2
  ) +
  geom_point() + 
  ylab(expression(paste('(',C[d]/C[o],')'^threshold))) +
  xlab('Number of species') +
  theme_bw() +
  annotate(
    'text',
    x=ifelse(fitting_test=='KS',60,40),
    y=ifelse(fitting_test=='KS',5.75,5.25),
    label=as.expression(bquote(R^2 == .(round(summary(model)$adj.r.squared,2))))
  ) +
  #annotate('text',x=55,y=5.5,label=paste('slope =',model_slope,'\u00B1',model_se)) +
  ggtitle('F') 
}

p2b<- {
  Cd_thresh %>% 
  ggplot(aes(logn,logc)) + 
  geom_smooth(method='lm') +
  geom_errorbar(
    aes(
      ymin=log(cd.threshold-s_cd.threshold),
      ymax=log(cd.threshold+s_cd.threshold)
    ),
    color=red,
    width=.04
  ) +
  geom_point() + 
  ylab(expression(paste('log(',C[d],'*)'))) +
  xlab(expression(paste('log(',N[S],')'))) +
  theme_bw() +
  annotate('text',x=4.2,y=1.1,label=paste('R2 = ',round(summary(model_log)$adj.r.squared,2))) +
  annotate('text',x=4.2,y=1.05,label=paste('slope =',model_log_slope,'\u00B1',model_log_se)) +
  ggtitle('C') 
}

gridExtra::grid.arrange(pD,pE,pF,ncol=3,respect=TRUE)

if(do.pdf){
  plotname='Specialist_SAD_results.pdf'
  if(file.exists(plotname)) unlink(plotname)
  ggsave(
    plotname,
    plot=gridExtra::grid.arrange(pD,pE,pF,ncol=3,respect=TRUE),
    width=21,
    height=7,
    units='cm',
    path='~/emergentneutrality/figures/'
  )
}

newcd <-
  Cd_thresh %>%
  mutate(
    log_ns=log(NS-min(NS)),
    log_cd=log(cd.threshold-min(cd.threshold))
  ) %>%
  filter(is.finite(log_ns)) 

newcd %>%
  ggplot(aes(log_ns,log_cd)) +
  geom_point()
coef(summary(lm(log_cd~log_ns,data=newcd)))


if(do_resources){
  raw_resources <-
    dat %>%
    ldply(function(l) c(resources=l$resources)) %>% 
    mutate(index=1:n()) %>%
    as_tibble
  
  resources <- 
    raw_resources %>%
    gather(key='key',value='x',-index) %>%
    group_by(index) %>%
    count(x) %>%
    mutate(cum=(sum(n)-cumsum(n)+n)/sum(n)) %>%
    rename(freq=n) %>%
    inner_join(metadata,by='index') %>%
    filter(limit=='Specialists') %>%
    ungroup
  
  sad_data <-
    resources %>%
    mutate(mu=round(mu,1)) %>%
    filter(phystime==max(phystime) & node=='wilson')
  
  sad_fit <-
    sad_data %>%
    ddply(.(mu),function(.){ 
      fitSAD(.,dbn='Pois',method=fitting_test)$dtf %>%
        filter(x>=min(.$x) & x<=max(.$x))  
    }) %>%
    ungroup %>%
    as_tibble
  
  sad_fit_stats <-
    sad_data %>%
    ddply(.(mu),function(.){ 
      foo=fitSAD(.,dbn='Pois',method=fitting_test)
      return(tibble(p.value=foo$p.value,statistic=foo$statistic))
    }) %>%
    ungroup %>%
    as_tibble
  
  pR <-
    sad_data %>%
    ggplot(aes(x,cum)) +
    geom_line(aes(x,cum),data=sad_fit,color=red) +
    geom_point() +
    theme_bw() +
    xlab('Abundance') +
    ylab('Cumulative frequency') +
    xlim(range(sad_data$x)) +
    facet_wrap(~mu,ncol=5,scale='free')
  
  plotname='Specialist_SAD_results_Resources.pdf'
  if(file.exists(plotname)) unlink(plotname)
  ggsave(
    plotname,
    plot=pR,
    width=26,
    height=10,
    units='cm',
    path='~/emergentneutrality/figures'
  )
  
}


if(FALSE){
  res_sum <-
    results_summary_cvm %>%
    mutate(`Goodness-of-fit test`='Cramer-von Mises') %>%
    union(
      results_summary_ks %>%
        mutate(`Goodness-of-fit test`='Kolmogorov-Smirnov')
    ) %>%
    filter(NS==50)
  
  df_logitfit <-
    df_logit_fit_cvm %>%
    mutate(`Goodness-of-fit test`='Cramer-von Mises') %>%
    union(
      df_logit_fit_ks %>%
        mutate(`Goodness-of-fit test`='Kolmogorov-Smirnov')
    ) %>%
    filter(NS==50)
  
  cd_thresh <-
    cd_thresh_cvm %>%
    mutate(test='Cramer-von Mises') %>%
    union(
      cd_thresh_ks %>%
        mutate(test='Kolmogorov-Smirnov')
    )
  
  pEall <-
    res_sum %>%
    ggplot(aes(color=`Goodness-of-fit test`,group=`Goodness-of-fit test`)) +
    geom_line(aes(mus,predLS),data=df_logitfit) +
    geom_point(aes(mu,mean)) +
    xlim(c(1,10)) +
    xlab(expression(paste(C[d]/C[o])))+
    ylab('Probability of logseries SAD') +
    ylim(c(0,1)) +
    ggtitle('A') +
    theme_bw() +
    theme(
      legend.position=c(7.5,8.5)/10,
      legend.box.background = element_rect(colour = "black"),
      legend.background = element_blank()
    )
  
  pFall <-
    cd_thresh %>%
    mutate(`Goodness-of-fit test`=test) %>%
    ggplot(aes(NS,cd.threshold,group=`Goodness-of-fit test`,color=`Goodness-of-fit test`)) +
    geom_smooth(aes(fill=`Goodness-of-fit test`),method='lm',se=TRUE) +
    geom_point() +
    theme_bw() +
    ylab(expression(paste('(',C[d]/C[o],')'^threshold))) +
    xlab('Number of species') +
    theme(
      legend.position=c(2.5,8.5)/10,
      legend.box.background = element_rect(colour = "black"),
      legend.background = element_blank()
    ) +
    ggtitle('B')
    
  gridExtra::grid.arrange(pEall,pFall,nrow=1,respect=TRUE)
  
  plotname='Specialist_SAD_results_CVM_vs_KS.pdf'
  if(file.exists(plotname)) unlink(plotname)
  ggsave(
    plotname,
    plot=gridExtra::grid.arrange(pEall,pFall,nrow=1,respect=TRUE),
    width=24,
    height=12,
    units='cm',
    path='~/emergentneutrality/figures/'
  )
}    
  
