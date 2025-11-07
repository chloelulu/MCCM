
summary_plot <- function(
  main.dirs,
  dirs,
  figure.dir,
  
  alpha.summary = T,
  alpha.plot = T,
  beta.summary = T,
  beta.plot = T,
  daa.file = 'DAA_P_R2.RData',
  DAA.summary = T,
  DAA.plot = T,
  volcano.plot =T,
  taxa.levels = c("Phylum","Class","Order","Family","Genus","Species"),
  volcano.top =NULL,
  volcano.level ='Genus',
  
  ## parameter for DAA Heatmap
  heat.p = T,
  bar.p = F,
  taxa.level.heatmap='Species',
  q.cut, p.cut, abund.cut, prev.cut,
  q.cut1, q.cut2, q.cut3, q.cut4,# stars on the plot
  plot.method = c('each.top','all.top'),
  top.n = 5,
  
  ## barplot paras
  taxa.level.bar = 'Species',
  
  heatmap.width,
  heatmap.height,
  fd = "/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/mforge_clean/Figure/",
  wd = '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/mforge_clean/',
  rd = '/research/bsi/projects/staff_analysis/m216453/2023_01_09_Oncobiome/mforge_clean/Result/'
){
  
  if(length(main.dirs)>1){len <- (main.dirs)}else{len <- dirs}
  
  setwd(fd)
  if(!dir.exists(figure.dir)){dir.create(figure.dir)}
  
  if(grepl('\\_func', rd)){alpha.summary=F; alpha.plot =F} # functional data does not have alpha diversity result, in case disable here
  
  # Alpha 
  if(alpha.summary==T){
    cat('----- Alpha Diveristy Summary -----\n')
    alpha.measures <- c("Observed","Chao1","Shannon","InvSimpson")
    setwd(wd)
    if(length(main.dirs)>length(dirs)){
      len <- (main.dirs)
    }else{
      len <- (dirs)
    }
    
    pval_coef_unadj.All <- pval_coef_adj.All <- 
      array(NA, c(length(alpha.measures), 3, length(len)), dimnames = list(alpha.measures, c('P','R2','coef'),len))
    for(main.dir in main.dirs){
      for(dir in dirs){
        setwd(paste0(rd,main.dir, '/Alpha/',dir))
        cat('[',main.dir,'-',dir,']\n')
        if(length(main.dirs)>length(dirs)) {dir.name = main.dir}else{dir.name = dir}
        fit2 <- fit1 <- NULL
        load('Alpha.RData')
        pval_coef_unadj.All[,,dir.name] <- fit1$res
        pval_coef_adj.All[,,dir.name] <- fit2$res
        }
      }
    
    setwd(fd)
    setwd(figure.dir)
    save(pval_coef_unadj.All, pval_coef_adj.All, file= 'Alpha_P_R2.RData')
  }

  if(alpha.plot==T){
    cat('----- Alpha Diveristy Plot -----\n')
    setwd(fd)
    setwd(figure.dir)
    load('Alpha_P_R2.RData')
    for(alpha.measure in alpha.measures){
      pval_coef_unadj <- t(pval_coef_unadj.All[alpha.measure,,])
      pval_coef_unadj[pval_coef_unadj[,'R2']<0,'R2'] <- 0
      pval_coef_unadj <- as.data.frame(pval_coef_unadj) %>% 
        mutate(direction = ifelse(is.na(coef),'none',ifelse(coef>0, 'positive','negative'))) %>%
        dplyr::select(-coef) %>%
        tibble::rownames_to_column('var') %>% melt()
      
      sub.r <- pval_coef_unadj[pval_coef_unadj$variable=='R2',]
      sub.p <- pval_coef_unadj[pval_coef_unadj$variable=='P',c('var','value')] %>% 
        mutate(P.col = ifelse(value<0.05,'*','')) %>% dplyr::select(-value)
      sub.rp <- inner_join(sub.p, sub.r) %>% mutate(value = value *100)
      ord <- (sub.rp[order(sub.rp$value),'var'])
      sub.rp <- within(sub.rp, var <- factor(var, levels=ord))
      p1 <- ggplot(sub.rp, aes(x = reorder(var,-value),y= value, fill = direction)) + 
        geom_bar(position="dodge", stat="identity",color = 'black') + 
        geom_text(aes(label = P.col), vjust = 0, colour = "black", size = 6) + 
        scale_fill_manual(values = c('positive'=brewer.pal(9, 'Set1')[5],'negative'=brewer.pal(9, 'Set1')[3], 'none'= brewer.pal(9, 'Pastel1')[9])) + 
        scale_y_continuous(limits = c(0, max(sub.rp$value)*1.2))+
        theme_classic() + 
        labs(y = ('Alpha-diversity variance explained(R2),%'), fill = '') +
        theme(panel.grid = element_blank(), 
              # legend.position = 'none',
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust =1, color = 'black'),
              axis.text.y = element_text(color = 'black')) + 
        ggtitle(alpha.measure)
      print(p1)
      setwd(fd)
      setwd(figure.dir)
      ggsave(file = paste0('AlphaR2_AllDodge_',alpha.measure,'_unAdj.pdf'), width = 5,height = 5.5)
      
      
      
      pval_coef_adj <- t(pval_coef_adj.All[alpha.measure,,])
      pval_coef_adj[pval_coef_adj[,'R2']<0,'R2'] <- 0
      pval_coef_adj <- as.data.frame(pval_coef_adj) %>% 
        mutate(direction = ifelse(is.na(coef),'none',ifelse(coef>0, 'positive','negative'))) %>%
        dplyr::select(-coef) %>%
        tibble::rownames_to_column('var') %>% melt()
      
      sub.r <- pval_coef_adj[pval_coef_adj$variable=='R2',]
      sub.p <- pval_coef_adj[pval_coef_adj$variable=='P',c('var','value')] %>% 
        mutate(P.col = ifelse(value<0.05,'*','')) %>% dplyr::select(-value)
      sub.rp <- inner_join(sub.p, sub.r) %>% mutate(value = value *100)
      ord <- (sub.rp[order(sub.rp$value),'var'])
      sub.rp <- within(sub.rp, var <- factor(var, levels=ord))
      p2 <- ggplot(sub.rp, aes(x = reorder(var,-value),y= value, fill = direction)) + 
        geom_bar(position="dodge", stat="identity",color = 'black') + 
        geom_text(aes(label = P.col), vjust = 0, colour = "black", size = 6) + 
        scale_fill_manual(values = c('positive'=brewer.pal(9, 'Set1')[5],'negative'=brewer.pal(9, 'Set1')[3], 'none'= brewer.pal(9, 'Pastel1')[9])) + 
        scale_y_continuous(limits = c(0, max(sub.rp$value)*1.2))+
        theme_classic() + 
        labs(y = ('Alpha-diversity variance explained(R2),%'), fill = '') +
        theme(panel.grid = element_blank(), 
              # legend.position = 'none',
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust =1, color = 'black'),
              axis.text.y = element_text(color = 'black')) + 
        ggtitle(alpha.measure)
      print(p2)
      setwd(fd)
      setwd(figure.dir)
      ggsave(file = paste0('AlphaR2_AllDodge_',alpha.measure,'_Adj.pdf'), width = 5,height = 5.5)
      
    }
  }
  
  # Beta 
  if(beta.summary==T){
    cat('----- Beta Diveristy Summary-----\n')
    R2.adj <- R2.unadj <- P.adj <- P.unadj <- NULL
    for(main.dir in main.dirs){
      for(dir in dirs){
        setwd(rd)
        setwd(main.dir)
        getwd()
        setwd('Beta')
        cat('[',main.dir, '-',dir,']\n')
        setwd(dir)
        load('R2_pvalue.RData')
        
        R2.adj <- cbind(R2.adj, r2.adj.mat)
        P.adj <- cbind(P.adj, pv.adj.mat)
        
        R2.unadj <- cbind(R2.unadj, r2.unadj.mat)
        P.unadj <- cbind(P.unadj, pv.unadj.mat)
      }
    }
    if(length(main.dirs)>1){name = main.dirs}else{name = dirs}
    colnames(P.adj) <- colnames(P.unadj) <- colnames(R2.adj) <- colnames(R2.unadj) <- name
    
    setwd(fd)
    setwd(figure.dir)
    save(R2.unadj, R2.adj, P.unadj, P.adj, file= 'Beta_P_R2.RData')
  }

  if(beta.plot ==T){
    cat('----- Beta Diveristy Plot -----\n')
    setwd(fd)
    setwd(figure.dir)
    load('Beta_P_R2.RData')
    dist.names <- rownames(R2.unadj)
    # plot.list <- list()
    # for(dist.name in dist.names){
    #   R2.adj.sub <- as.data.frame(R2.adj[dist.name,]) %>%
    #     tibble::rownames_to_column(dist.name)
    #   colnames(R2.adj.sub)[2] <- 'value'
    #   R2.m <- melt(R2.adj.sub)
    #   R2.m$value <- R2.m$value * 100
    #   R2.m <- R2.m[order(R2.m$value),]
    #   plot.list[[dist.name]] <- ggplot(R2.m, aes(x = variable,y= (value),
    #                                              fill = reorder(!!as.name(dist.name), value))) +
    #     geom_bar(position="stack", stat="identity") +
    #     scale_fill_manual(values = c(brewer.pal(8, 'Pastel1'),brewer.pal(8, 'Pastel2')))+
    #     theme_classic() +
    #     labs(x = 'Taxa', y = 'Beta-diversity variance explained (R2),%', fill = '') +
    #     theme(panel.grid = element_blank(),
    #           legend.position = 'right',
    #           axis.text.x = element_blank()) +
    #     ggtitle(dist.name)
    # }
    # ggarrange(plot.list[[1]], plot.list[[2]],plot.list[[3]], plot.list[[4]],
    #           common.legend = T, legend = 'right')
    # setwd(fd)
    # setwd(figure.dir)
    # ggsave(file = 'BetaR2_AllStack_Adj.pdf', width = 5,height = 7)
    # 
    # 
    # 
    # plot.list <- list()
    # for(dist.name in dist.names){
    #   R2.adj.sub <- as.data.frame(R2.unadj[dist.name,]) %>%
    #     tibble::rownames_to_column(dist.name)
    #   colnames(R2.adj.sub)[2] <- 'value'
    #   R2.m <- melt(R2.adj.sub)
    #   R2.m$value <- R2.m$value * 100
    #   R2.m <- R2.m[order(R2.m$value),]
    #   plot.list[[dist.name]] <- ggplot(R2.m, aes(x = variable,y= (value),
    #                                              fill = reorder(!!as.name(dist.name), value))) +
    #     geom_bar(position="stack", stat="identity") +
    #     scale_fill_manual(values = c(brewer.pal(8, 'Pastel1'),brewer.pal(8, 'Pastel2')))+
    #     theme_classic() +
    #     labs(x = 'Taxa', y = 'Beta-diversity variance explained (R2),%', fill = '') +
    #     theme(panel.grid = element_blank(),
    #           legend.position = 'right',
    #           axis.text.x = element_blank()) +
    #     ggtitle(dist.name)
    # }
    # ggarrange(plot.list[[1]], plot.list[[2]],plot.list[[3]], plot.list[[4]],
    #           common.legend = T, legend = 'right')
    # setwd(fd)
    # setwd(figure.dir)
    # ggsave(file = 'BetaR2_AllStack_unAdj.pdf', width = 5,height = 7)

    
    setwd(fd)
    setwd(figure.dir)
    load('Beta_P_R2.RData')
    dist.names <- rownames(R2.unadj)
    setwd(fd)
    setwd(figure.dir)
    pdf(file = paste0('BetaR2_AllDodge_Adj.pdf'), width = 5,height = 5)
    for(dist.name in dist.names){
      P.adj.sub <- as.data.frame(P.adj[dist.name,]) %>% 
        tibble::rownames_to_column(dist.name)
      colnames(P.adj.sub) <- c('variable','P')
      P.adj.sub$color <- ifelse(P.adj.sub$P<0.05, 'sig','nosig')
      
      R2.adj.sub <- as.data.frame(R2.adj[dist.name,]) %>% 
        tibble::rownames_to_column(dist.name)
      colnames(R2.adj.sub) <- c('variable','R2')
      R2_P <- inner_join(R2.adj.sub,P.adj.sub)
      R2_P$R2[R2_P$R2<0] <- 0
      ord <- R2_P[order(-R2_P$R2),'variable']
      R2_P$sig <- ifelse(R2_P$P<0.05,'*','')
      R2_P$R2 <- R2_P$R2 *100
      R2_P <- within(R2_P, variable <- factor(variable, levels=ord))
      p1 <- ggplot(R2_P, aes(x = variable,y= (R2))) + 
        geom_bar(position="dodge", stat="identity",color='black', aes(fill =color)) + 
        theme_classic() + 
        labs(y = 'Beta-diversity variance explained (R2),%', fill = '') + 
        geom_text(aes(label = sig), vjust = 0, colour = "black", size = 6) + 
        scale_fill_manual(values = c('sig'=brewer.pal(9, 'Paired')[1],'nosig'= brewer.pal(9, 'Pastel1')[9])) + 
        scale_y_continuous(limits = c(0, max(R2_P$R2)*1.2))+
        theme(panel.grid = element_blank(), 
              legend.position = 'none',
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust =1)) + 
        ggtitle(dist.name)
      print(p1)
    }
    dev.off()
    
    setwd(fd)
    setwd(figure.dir)
    pdf(file = paste0('BetaR2_AllDodge_unAdj.pdf'), width = 5,height =5)
    for(dist.name in dist.names){
      P.unadj.sub <- as.data.frame(P.unadj[dist.name,]) %>% 
        tibble::rownames_to_column(dist.name)
      colnames(P.unadj.sub) <- c('variable','P')
      P.unadj.sub$color <- ifelse(P.unadj.sub$P<0.05, 'sig','nosig')
      
      R2.unadj.sub <- as.data.frame(R2.unadj[dist.name,]) %>% 
        tibble::rownames_to_column(dist.name)
      colnames(R2.unadj.sub) <- c('variable','R2')
      R2_P <- inner_join(R2.unadj.sub,P.unadj.sub)
      R2_P$R2[R2_P$R2<0] <- 0
      ord <- R2_P[order(-R2_P$R2),'variable']
      R2_P$sig <- ifelse(R2_P$P<0.05,'*','')
      R2_P$R2 <- R2_P$R2 *100
      R2_P <- within(R2_P, variable <- factor(variable, levels=ord))
      p1 <- ggplot(R2_P, aes(x = variable,y= (R2))) + 
        geom_bar(position="dodge", stat="identity",color='black', aes(fill =color)) + 
        theme_classic() + 
        labs(y = 'Beta-diversity variance explained (R2),%', fill = '') + 
        geom_text(aes(label = sig), vjust = 0, colour = "black", size = 6) + 
        scale_fill_manual(values = c('sig'=brewer.pal(9, 'Paired')[1],'nosig'= brewer.pal(9, 'Pastel1')[9])) + 
        scale_y_continuous(limits = c(0, max(R2_P$R2)*1.2))+
        theme(panel.grid = element_blank(), 
              legend.position = 'none',
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust =1)) + 
        ggtitle(dist.name)
      print(p1)
      
    }
    dev.off()
}
    
  # DAA
  if(DAA.summary==T){
    cat('----- DAA Summary -----\n')
    R2.All <- P.All <- Q.All <- list()
    for(taxa.level in taxa.levels){
      cat('----',taxa.level,'----\n')
      R2 <- P <- Q <- list()
      dd <- c()
      for(main.dir in main.dirs){
        for(dir in dirs){
          cat('[',main.dir,'-',dir,']\n')
          setwd(rd)
          setwd(main.dir)
          setwd('DAA')
          setwd(dir)
          file <- list.files(pattern = 'ZicoSeq.Rdata$')
          if(length(file.exists(file))==1){
            load(file)
            p <- diff.obj$pv.list[[taxa.level]]
            q <- diff.obj$qv.list[[taxa.level]]
            r2 <- diff.obj$R2.list[[taxa.level]]
            coef <- diff.obj$coef.list[[taxa.level]]
            idx <- colnames(coef)[grep(dir,colnames(coef))]
            if(length(idx)==1){
              coef <- coef[,idx, drop =F]
            }else{ # for multi-level categorical data, coefficient does not have meaning
              coef <- coef[,idx[1], drop =F]
              coef[,1] <- 1
            }
            r2 <- r2[,1] * sign(coef)
            
            if(length(main.dirs)==1){name <- dir}else{name = main.dir}
            colnames(r2) <- colnames(q) <- colnames(p) <- name
            r2 <- as.data.frame(r2) %>% tibble::rownames_to_column(taxa.level)
            p <- as.data.frame(p) %>% tibble::rownames_to_column(taxa.level)
            q <- as.data.frame(q) %>% tibble::rownames_to_column(taxa.level)
            
            R2[[name]] <- r2
            P[[name]] <- p
            Q[[name]] <- q
            gc()
          }else{
            dd <- c(dd, dir)
            cat(dir,'no sig! \n')
          }
        }
      }
      P.All[[taxa.level]] <- Reduce(function(...) merge(..., by=taxa.level, all=TRUE), P)
      Q.All[[taxa.level]] <- Reduce(function(...) merge(..., by=taxa.level, all=TRUE), Q)
      R2.All[[taxa.level]] <- Reduce(function(...) merge(..., by=taxa.level, all=TRUE), R2)
    }
    
    setwd(fd)
    setwd(figure.dir)
    save(P.All, Q.All, R2.All, file= 'DAA_P_R2.RData')
  }
  
 if(DAA.plot==T |volcano.plot==T){
   setwd(fd)
   setwd(figure.dir)
   load(daa.file)
 }
  
  if(volcano.plot==T){
    if(length(dirs) >= length(main.dirs)){
      drs <- dirs
    }else{
        drs <- main.dirs
    }
    setwd(fd)
    setwd(figure.dir)

    pdf(file = paste0('DAA_volcano_',volcano.level,'.pdf'), width = 8,height = 7)
    for(dir in drs){
      if(dir %in% colnames(Q.All[[volcano.level]])){
        R2 <- R2.All[[volcano.level]][,c(volcano.level,dir), drop =F] %>% na.omit()
        colnames(R2)[2] <- 'R2'
        Q <- Q.All[[volcano.level]][,c(volcano.level,dir), drop =F] %>% na.omit()
        colnames(Q)[2] <- 'Q'
        R2_P <- inner_join(R2,Q) %>% mutate(`-log10P` = -log10(Q))
        R2_P$sig <- 'nosig'
        R2_P[(R2_P$Q<q.cut),'sig'] <- 'sig'
        R2_P <- R2_P %>% arrange(Q, desc(abs(R2)), !!as.name(volcano.level))
        R2_P$label <- ''
        
        if(sum(R2_P$Q <q.cut) >volcano.top){
          volcano.top.n <- volcano.top
        }else{
          volcano.top.n  <- sum(R2_P$Q<q.cut)
        }
        if(volcano.top.n >0){
          R2_P[1:volcano.top.n,'label'] <- R2_P[1:volcano.top.n,volcano.level]
          R2_P$label <- gsub('.*g__|.*s__|.*f__|.*o__|.*c__|.*p__','',R2_P$label )
        }
        text.size <- 16
        volcano.p <- ggplot(R2_P, aes(x = R2, y = `-log10P`)) + 
          geom_point(aes(color = sig), size = 0.8) + 
          geom_vline(aes(xintercept = 0), color = "gray", linetype = "dashed") + 
          geom_hline(aes(yintercept = -log10(q.cut)), color = "gray", 
                     linetype = "dashed") + 
          theme_bw()+
          scale_colour_manual(values = c('sig'='#c91f37','nosig'='#549688')) + #c('sig'='#f73668','nosig'='#f0ead8')
          scale_y_continuous(limits = c(0,max((R2_P$`-log10P`)) * 1.1)) + 
          scale_x_continuous(limits = c(min(R2_P$R2)*1.1,max(R2_P$R2) * 1.1)) + 
          ggrepel::geom_text_repel(aes(label = label),max.overlaps = Inf, color = "black") + 
          theme(axis.text = element_text(color = "black", size = text.size), 
                axis.title = element_text(color = "black", size = text.size), 
                panel.grid = element_blank(),
                legend.position = 'none',
                legend.text = element_text(color = "black", size = text.size), 
                legend.title = element_text(color = "black", size = text.size))+
          labs(x = bquote(R^2), y = '-log10(FDR adjusted P)') + 
          ggtitle(paste(dir, ': ',sum(R2_P$Q <q.cut),' ',volcano.level,' with FDR <', q.cut))
        print(volcano.p)
      }

    }
    dev.off()
  }
  
  
  if(DAA.plot==T){
    cat('----- DAA Plot -----\n')
    ## color scheme: R2 * effect direction (multi-class variable has no direction, all positive)
    ## Significance(Q<0.1): "+" positive coefficient; "-" negative coefficient; "." no direction
    ## load abundance >0.01 or prevelance > 10% samples
    if(grepl('\\_func',rd) | grepl('\\_func',figure.dir)){
      ## [to be revised] for a better coding strategy
      setwd(rd)
      if(grepl('PanCancer_func|CancerOnly_func',figure.dir)){
        setwd(main.dirs)
        load('data.obj.wk2.RData')
      }else if(grepl('subCancerX-Ex_func',figure.dir)){
        load('base_of_tongue/data.obj.wk2.RData')
      }
      data.obj <- data.obj3
    }else{
      ## for otu data
      if(length(grep('Ex-|CancerOnly',c(main.dirs, dirs)))>1){
        setwd(rd)
        load('CancerOnly/data.obj.wk.RData')
      }else{
        setwd(wd)
        load('Data/data.obj.raw.core.RData')
      }
    }

      Q.sig <- as.matrix(Q.All[[taxa.level.heatmap]] %>% column_to_rownames(taxa.level.heatmap))
      R2.sig <- as.matrix(R2.All[[taxa.level.heatmap]] %>% column_to_rownames(taxa.level.heatmap))
      P.sig <- as.matrix(P.All[[taxa.level.heatmap]] %>% column_to_rownames(taxa.level.heatmap))
      
      if(heat.p ==T){
        if(length(grep('^Control',c(dirs, main.dirs)))>0){
          idx <- apply(Q.sig, 1, function(x) (sum(x<q.cut,na.rm = T)>0 & sum(is.na(x))==0)) # To aviod some taxa exist in this but not in others
          Q.sig <- Q.sig[idx,, drop =F]
          P.sig <- P.sig[idx,, drop =F]
          R2.sig <- R2.sig[idx,, drop =F]
        }
        if(plot.method =='each.top'){ # the smalled q value
          tt <- apply(Q.sig, 2, function(x){
            x.ord <- order(x, decreasing = F)
            if(sum(x[x.ord]<q.cut, na.rm = T) >= top.n){
              xx <- x.ord[1:top.n]
            }else{
              xx <- which(x < q.cut)
            }
            return(xx)
          }) %>% unlist() %>% unname %>% as.vector  %>%unique
          sig.taxa <- rownames(Q.sig)[1:nrow(Q.sig) %in% as.vector(tt)]
        }
        
        if(plot.method =='each.top.R2'){
          # extract each variable's taxa any q< q.cut
          tt <- apply(Q.sig, 2, function(x){
            x.ord <- order(x, decreasing = F)
            xx <- which(x < q.cut)
            return(xx)
          }) %>% unlist() %>% unname %>% as.vector  %>%unique
          sig.taxa <- rownames(Q.sig)[1:nrow(Q.sig) %in% as.vector(tt)]
          
          # then extract each variable's taxa top.n R2
          R2.sig0 <- R2.sig[sig.taxa,,drop =F]
          rr <- apply(R2.sig0, 2, function(x){
            x.ord <- order(abs(x), decreasing = T)
            return((x.ord[1:top.n]))
            }) %>% as.vector() %>% unique
          sig.taxa <- rownames(R2.sig0[rr,,drop =F])
        }
          
        
        if(plot.method =='all.top'){
          # all taxa tested by ZicoSe & with q < q.cut
          if(grepl('\\_func',rd)){
            comm <- data.obj$abund.list[[taxa.level.heatmap]]
          }else{
            comm <- data.obj.rff$abund.list[[taxa.level.heatmap]]
          }
          prop <- t(t(comm) /colSums(comm))
          filter.ind <- names(which(rowMeans(prop != 0) >= prev & rowMaxs(prop, useNames = T) >= minp))
          sig.taxa <- names(which(apply(Q.sig, 1, function(x) sum(x<q.cut)) >0));length(sig.taxa)
          sig.taxa <- intersect(filter.ind,sig.taxa);length(sig.taxa)
        }
        
        cat(length(sig.taxa),' ',taxa.level.heatmap, 'will be shown on heatmap!\n')
        Q.sig <- Q.sig[sig.taxa,, drop =F]
        R2.sig <- R2.sig[sig.taxa,, drop =F]
        P.sig <- P.sig[sig.taxa,, drop =F]
        
        rownames(R2.sig) <- gsub('^p__|.*g__|.*f__|.*o__|.*c__|.*s__','',rownames(R2.sig))
        rownames(Q.sig) <- gsub('^p__|.*g__|.*f__|.*o__|.*c__|.*s__','',rownames(Q.sig))
        rownames(P.sig) <- gsub('^p__|.*g__|.*f__|.*o__|.*c__|.*s__','',rownames(P.sig))
        if(ncol(R2.sig)>1){
          R2.sig.scale <- t(apply(R2.sig,1, function(x) scale(x)))
          colnames(R2.sig.scale) <- colnames(R2.sig)
        }
        
        
        idd <- names(which(apply(Q.sig, 2, function(x) sum(x<q.cut))>0))
        R2.sig <- R2.sig[,idd, drop =F]
        Q.sig <- Q.sig[,idd, drop =F]
        
        if('CancerOnly' %in% main.dirs){
          tech.vars <- c("Bristol_score","Sample_season")
          demo.vars <- c("Sex", "Age", "BMI", "Urban")
          clin.vars <- c("Abx_day_365","PPI_day_365","Cancer_class","Metastasis","GI_nonGI","Charlson_score", "Elix_score")
          
          variables <- c(tech.vars, demo.vars, clin.vars)
          variables <- variables[variables %in% colnames(Q.sig)]
          Q.sig <- Q.sig[,variables]
          R2.sig <- R2.sig[,variables]
          P.sig <- P.sig[,variables]
          types <- c(rep('Technical',length(tech.vars)),rep('Demographic',length(demo.vars)),rep('Clinical',length(clin.vars)))
          column_ha = HeatmapAnnotation(Factors = types, annotation_height = 1,border = TRUE, simple_anno_size = unit(0.3, "cm"), 
                                        col  = list(Factors = c("Technical" =brewer.pal(8,'Dark2')[1],
                                                                "Demographic" = brewer.pal(8,'Dark2')[6],
                                                                "Clinical"=brewer.pal(8,'Dark2')[3])),
                                        annotation_legend_param = list(
                                          Factors = list(
                                            title = "Factors",
                                            at = c("Demographic", "Clinical","Technical"),
                                            labels = c("Demographic", "Clinical","Technical"))
                                        )
          )
          
          heatmap.obj <- Heatmap(R2.sig, name = "R2 x effect direction", 
                                 heatmap_legend_param = list(direction = "vertical"),
                                 col = colorRamp2(c(min(R2.sig), 0, max(R2.sig)), 
                                                  c(brewer.pal(9, 'Blues')[7],'white', 
                                                    brewer.pal(9, 'YlOrRd')[4])),#
                                 top_annotation = column_ha,
                                 column_gap = unit(1, "mm"), 
                                 border = TRUE,
                                 cell_fun = function(j, i, x, y, width, height, fill) {
                                   if(Q.sig[i,j] < q.cut1) {
                                     grid.text('***', x = x, y = y,r = unit(1/30,'cm'))
                                   }else if(Q.sig[i,j] < q.cut2){
                                     grid.text('**', x = x, y = y,r = unit(1/30,'cm'))
                                   }else if(Q.sig[i,j] < q.cut3){
                                     grid.text('*', x = x, y = y,r = unit(1/30,'cm'))
                                   }
                                 },
                                 column_split = factor(c(rep('Technical',sum(tech.vars %in% variables)),
                                                         rep('Demographic',sum(demo.vars %in% variables)),
                                                         rep('Clinical',sum(clin.vars %in% variables)))),
                                 rect_gp = gpar(col= "white"),
                                 show_column_names = T, show_row_names = T,
                                 show_column_dend = F, show_row_dend = T,
                                 show_heatmap_legend = T,
                                 # row_order = c(1:nrow(R2.sig)),
                                 row_names_gp = gpar(fontface = 'italic'),
                                 row_names_max_width = max_text_width(rownames(R2.sig), gp = gpar(fontsize = 18))
          ) 
        }else if(length(grep('^Ex-',main.dirs))>0){
          heatmap.obj  <- Heatmap(R2.sig, name = "R2 x effect direction", 
                                  heatmap_legend_param = list(direction = "vertical"),
                                  col = colorRamp2(c(min(R2.sig), 0, max(R2.sig)), 
                                                   c(brewer.pal(9, 'Blues')[7],'white', 
                                                     brewer.pal(9, 'YlOrRd')[4])),#brewer.pal(9, 'Blues')[5] brewer.pal(9, 'YlOrRd')[4]
                                  column_gap = unit(1, "mm"), 
                                  border = TRUE,
                                  cell_fun = function(j, i, x, y, w, h, fill) {
                                    if(Q.sig[i,j] < q.cut1) {
                                      grid.text('***', x = x, y = y, r = unit(1/30,'cm'))
                                    }else if(Q.sig[i,j] < q.cut2){
                                      grid.text('**', x = x, y = y, r = unit(1/30,'cm'))
                                    }else if(Q.sig[i,j] < q.cut3){
                                      grid.text('*', x = x, y = y,r = unit(1/30,'cm'))
                                    }else if(Q.sig[i,j] < q.cut4){
                                      grid.circle(x = x, y = y, r = unit(1/40,'cm'), gp = gpar(fill = 'black', col = 'black'))
                                    }
                                  },
                                  rect_gp = gpar(col= "white"),
                                  show_column_names = T, show_row_names = T,
                                  show_column_dend = F,
                                  show_row_dend = T,
                                  row_names_gp = gpar(fontface = 'italic'),
                                  row_names_max_width = max_text_width(rownames(R2.sig), gp = gpar(fontsize = 18))
          ) 

        }else{
          heatmap.obj  <- Heatmap(R2.sig, name = "R2 x effect direction", 
                                  heatmap_legend_param = list(direction = "vertical"),
                                  col = colorRamp2(c(min(R2.sig), 0, max(R2.sig)), 
                                                   c(brewer.pal(9, 'Blues')[7],'white', 
                                                     brewer.pal(9, 'YlOrRd')[4])),#brewer.pal(9, 'Blues')[5] brewer.pal(9, 'YlOrRd')[4]
                                  column_gap = unit(1, "mm"), 
                                  border = TRUE,
                                  cell_fun = function(j, i, x, y, w, h, fill) {
                                    if(Q.sig[i,j] < q.cut1) {
                                      grid.text('***', x = x, y = y, r = unit(1/30,'cm'))
                                    }else if(Q.sig[i,j] < q.cut2){
                                      grid.text('**', x = x, y = y, r = unit(1/30,'cm'))
                                    }else if(Q.sig[i,j] < q.cut3){
                                      grid.text('*', x = x, y = y,r = unit(1/30,'cm'))
                                    }
                                  },
                                  rect_gp = gpar(col= "white"),
                                  show_column_names = T, show_row_names = T,
                                  show_column_dend = F,
                                  show_row_dend = T,
                                  row_names_gp = gpar(fontface = 'italic'),
                                  row_names_max_width = max_text_width(rownames(R2.sig), gp = gpar(fontsize = 18))
          ) 
        }
        return(heatmap.obj)
      }
      
      if(bar.p == T){ # not finished
        Q.sig <- as.matrix(Q.All[[taxa.level.bar]] %>% column_to_rownames(taxa.level.bar))
        R2.sig <- as.matrix(R2.All[[taxa.level.bar]] %>% column_to_rownames(taxa.level.bar))
        P.sig <- as.matrix(P.All[[taxa.level.bar]] %>% column_to_rownames(taxa.level.bar))
        
        for(var.bar in colnames(Q.sig)){
          res <- cbind.data.frame(as.data.frame(Q.sig[,var.bar,drop =F]) %>% dplyr::rename(Q=1),
                                  as.data.frame(P.sig[,var.bar,drop =F]) %>% dplyr::rename(P=1),
                                  as.data.frame(R2.sig[,var.bar,drop =F])%>% dplyr::rename(R2=1)) %>% dplyr::filter(Q < q.cut)
          if(nrow(res)>0){
              if(grepl('\\_func',rd) | grepl('\\_func',figure.dir)){
              comm <- data.obj$abund.list[[taxa.level.bar]]
              meta <- data.obj$meta.dat[,var.bar, drop =F]%>% na.omit()
            }else{
              comm <- data.obj.rff$abund.list[[taxa.level.bar]]
              meta <- data.obj.rff$meta.dat[,var.bar, drop =F] %>% na.omit()
            }
            prop <- t(t(comm) /colSums(comm))
            setwd(fd)
            setwd(figure.dir)
            pdf(paste0('barplot_',var.bar, '_',taxa.level.bar,'.pdf'), width = 6, height = 5)
            if(length(unique(meta[,var.bar]))>2){
              for(i in rownames(res)){
                df.bar <- merge(t(prop[i,, drop =F]), na.omit(meta), by = 0)
                bar <- ggplot(df.bar, aes(x = !!as.name(var.bar), y = !!as.name(i))) + geom_point() + theme_classic() + 
                  scale_y_continuous(trans = sqrt_trans(),
                                     # limits = c(0,0.05),
                                     breaks = trans_breaks("sqrt", function(x) x^2),
                                     labels = trans_format("sqrt", math_format(.x^2)))
                print(bar)
                
              }
            }else{
              for(i in rownames(res)){
                df.bar <- merge(t(prop[i,, drop =F]), na.omit(meta), by = 0)
                bar <- ggplot(df.bar, aes(x = !!as.name(var.bar), y = !!as.name(i))) + 
                  geom_violin(aes(fill = !!as.name(var.bar)), trim = F) +
                  geom_jitter(width = 0.1, size = 0.5) + 
                  # geom_boxplot(width= 0.2) +
                  theme_classic() + 
                  scale_fill_brewer(palette = 'Dark2') + 
                  scale_y_continuous(trans = sqrt_trans(),
                                     # limits = c(0,0.05),
                                     breaks = trans_breaks("sqrt", function(x) x^2),
                                     labels = trans_format("sqrt", math_format(.x^2)))
                print(bar)
              }
            }
            
            dev.off()
          }
          
        }
      }
  }
}









# tree nodes selected by 3 criteria: qcut, minp, prev
generate_tree <- function(data.obj, dir, Q.All, R2.All, taxa.level, 
                          filter.mode='or', q.cut, minp, prev, # filter taxa show on the tree by abundance, prevealcance amd minp
                          subset.taxa = NULL, subset.level =NULL,subset.group=F, # subset specific level of taxa; For example, only want to retain Firmicute and Bacteriodetes
                          tip=T,
                          branch.level# branch color by group
                          ){
  
  Q.sig <- as.matrix(Q.All[[taxa.level]][,c(taxa.level,dir), drop =F] %>% column_to_rownames(taxa.level))
  R2.sig <- as.matrix(R2.All[[taxa.level]][,c(taxa.level,dir), drop =F] %>% column_to_rownames(taxa.level))
  rownames(Q.sig) <- gsub('s__|.*g__','',rownames(Q.sig))
  rownames(Q.sig) <- gsub('\\ ','_',rownames(Q.sig))
  rownames(R2.sig) <- gsub('s__|.*g__','',rownames(R2.sig))
  rownames(R2.sig) <- gsub('\\ ','_',rownames(R2.sig))
  
  abund <- (data.obj$abund.list[[taxa.level]])
  rownames(abund) <- gsub('\\ ','_',rownames(abund))
  rownames(abund) <- gsub('s__|.*g__','',rownames(abund))
  prop <- t(t(abund) /colSums(abund))
  
  if(filter.mode == 'or'){
    abund.prev1 <- names(which(rowSums(prop > minp) > 0))
    abund.prev2 <- names(which(rowSums(prop > 0) > prev * ncol(prop)))
    abund.prev <- unique(c(abund.prev1, abund.prev2))
    hypen = '|'
  }
  
  if(filter.mode == 'both'){
    abund.prev <- names(which(rowSums(prop > minp) > prev * ncol(prop)))
    hypen = '&'
  }
  
  if(subset.group ==T){
    subset.nodes <- data.obj$otu.name[,taxa.level][which(data.obj$otu.name[,subset.level] %in% subset.taxa)]
    subset.nodes <- gsub('^s\\_\\_','',names(subset.nodes))
  }
  
  sig.taxa <- names(which(Q.sig[,dir] <q.cut))
  sig.nodes <- intersect(abund.prev,sig.taxa) 
  
  if(subset.group ==T){
    sig.nodes <- intersect(subset.nodes,sig.nodes)
    }# species for final analysis
  
  if(length(sig.nodes)>0){
    sig.coefs <- R2.sig[sig.nodes,dir]
    direction <- ifelse(sig.coefs>0,'pos','neg')
    direction <- lapply(split(direction, direction), function(x) names(x))
    
    otu.name <- (data.obj$otu.name)
    rownames(otu.name) <- gsub('s\\_\\_','',rownames(otu.name))
    taxa.split <- otu.name[sig.nodes,branch.level]
    taxa.split <- gsub('f\\_\\_|p\\_\\_|o\\_\\_|c\\_\\_|g\\_\\_','',taxa.split)
    taxa.split <- lapply(split(taxa.split, taxa.split), function(x) names(x))
    
    data.obj$tree$tip.label <- gsub('^s__','',data.obj$tree$tip.label)
    tips <- data.obj$tree$tip.label[!(data.obj$tree$tip.label %in% sig.nodes)]
    tree.sub <- drop.tip(data.obj$tree, tips)
    tree.sub <- groupOTU(tree.sub, direction, group_name = 'group')
    tree.sub <- groupOTU(tree.sub, taxa.split, group_name = 'taxa.split')
    
    cols <- brewer.pal(8,'Set1')[1:3]
    names(cols) <- c('pos','neg')
    tree.sub$node.label <- gsub('g__|f__|o__|p__|c__|s__','',tree.sub$node.label)
    
    cols.shape <- c(24,25)# increase; decrease 
    names(cols.shape) <- c('pos','neg')
    cols.fill <- c('red','forestgreen')
    names(cols.fill) <- c('pos','neg')
    
    ord <- gsub('f\\_\\_|p\\_\\_|o\\_\\_|c\\_\\_|g\\_\\_','',unique(otu.name[get_taxa_name(ggtree(tree.sub)),branch.level]))
    
    p <- ggplot(tree.sub, aes(x, y)) +
      geom_tree() +
      theme_tree() +
      geom_tiplab(align=T,aes(color =taxa.split), fontface='italic',size =5) +
      lims(x=c(0,3.5)) +
      scale_shape_manual(values=cols.shape)+
      scale_color_manual(values = c(brewer.pal(8,'Dark2'),brewer.pal(8,'Set2'),
                                    brewer.pal(12,'Paired')[-11],brewer.pal(9,'Set1')[-6],
                                    brewer.pal(8,'Pastel2')[-c(6:7)],brewer.pal(8,'Accent')[-4],brewer.pal(8,'Set3')[-4]),
                         breaks = ord) +
      scale_fill_manual(values = cols.fill) +
      theme(legend.position="right", legend.text = element_text(size = 16, color = 'black', face = 'italic'),legend.title = element_text(size = 16, color = 'black')) +
      labs(color = branch.level, shape = 'direction',
           caption = paste0('[We found ',length(sig.taxa), ' taxa with FDR<',q.cut,']')) +
      ggtitle(paste0(dir,'(FDR<',q.cut,' & minp>',minp,hypen,' prev>',prev,')'))
    if(tip ==T){
      p <- p + geom_tippoint(size = 2, aes(shape = group, fill =group), show.legend = F)
    }
      
  }else{
    cat('No sig taxa\n')
  }

  return(list(plot=p, sig.taxa= sig.taxa))
}






apply_func <- function(i, j, cutoff = 0.05) {
  x <- tree.tmp0[i, j]
  y <- tree.tmp1[i, j]
  if (x < cutoff) {
    if (y > 0) {
      return("Positive")
    } else {
      return("Negative")
    }
  } 
  return('No')
}
