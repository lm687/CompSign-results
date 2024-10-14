
give_plot_betas_all_ct <- function(nucleotide1, names_trinucleotide, prepare_for_tikz=T){
  
  nucleotide_colours_logR <- c('C$>$A/T$>$G'= '#3cb371', 'C$>$G/T$>$G'= '#90ee90', 'C$>$T/T$>$G'= '#66cdaa',
                               'T$>$A/T$>$G'= '#cd5c5c', 'T$>$C/T$>$G'= '#f4a460')
  nucleotide_colours <- c('C>A' = '#3cb371', 'C>G'= '#90ee90', 'C>T'= '#66cdaa',
                          'T>A'= '#cd5c5c', 'T>C'= '#f4a460', 'T>G'='red')
  nucleotide_colours_dollar <- c('C$>$A' = '#3cb371', 'C$>$G'= '#90ee90', 'C$>$T'= '#66cdaa',
                                 'T$>$A'= '#cd5c5c', 'T$>$C'= '#f4a460', 'T$>$G'='red')
  
  nucleotide_colours_logR <- c('C$>$A/T$>$G'= '#a53606', 'C$>$G/T$>$G'= '#b32db5', 'C$>$T/T$>$G'= '#881a58',
                               'T$>$A/T$>$G'= '#0e288e', 'T$>$C/T$>$G'= '#164c64')
  nucleotide_colours <- c('C>A' = '#a53606', 'C>G'= '#b32db5', 'C>T'= '#881a58',
                          'T>A'= '#0e288e', 'T>C'= '#164c64', 'T>G'='red')
  nucleotide_colours_dollar <- c('C$>$A' = '#a53606', 'C$>$G'= '#b32db5', 'C$>$T'= '#881a58',
                                 'T$>$A'= '#0e288e', 'T$>$C'= '#164c64', 'T$>$G'='red')
  
  nucleotide_colours_logR <- c('C$>$A/T$>$G'= '#377eb8', 'C$>$G/T$>$G'= '#ff7f00', 'C$>$T/T$>$G'= '#984ea3',
                               'T$>$A/T$>$G'= '#f781bf', 'T$>$C/T$>$G'= '#a65628')
  nucleotide_colours <- c('C>A' = '#377eb8', 'C>G'= '#ff7f00', 'C>T'= '#984ea3',
                          'T>A'= '#f781bf', 'T>C'= '#a65628', 'T>G'='red')
  nucleotide_colours_dollar <- c('C$>$A' = '#377eb8', 'C$>$G'= '#ff7f00', 'C$>$T'= '#984ea3',
                                 'T$>$A'= '#f781bf', 'T$>$C'= '#a65628', 'T$>$G'='red')
  
  betas_nucleotides <- lapply(nucleotide1, function(i) plot_betas(i, return_df = T))
  betas_nucleotides <- lapply(betas_nucleotides, function(i){
    i$LogR <- names_trinucleotide[i$LogR]
    # rownames(i) <- make.names(i$LogR, unique = T)
    i
  })
  
  betas_nucleotides_slopes <- do.call('cbind', lapply(betas_nucleotides, function(i) i%>% filter(type_beta == 'Slope' ) %>% select(Estimate)))
  colnames(betas_nucleotides_slopes) <- names(nucleotide1)
  rownames(betas_nucleotides_slopes) <- names_trinucleotide

  
  if(prepare_for_tikz){
    titley=("$\\hat{\\beta}_1$")
    rownames(betas_nucleotides_slopes) <- gsub(">", "$>$", rownames(betas_nucleotides_slopes))
  }else{
    names(nucleotide_colours_logR) <- gsub("\\$>\\$", ">", names(nucleotide_colours_logR))
    titley=latex2exp::TeX(r"(Estimated $\beta_1$)")
  }
  
  df_to_plot = melt(as(betas_nucleotides_slopes, 'matrix'))
  df_to_plot$Var2 = factor(Var2,levels=names(sort(colMeans(betas_nucleotides_slopes))))
  ggplot(df_to_plot,
         aes(x=x,
             col=Var1, shape=Var1, y=value))+geom_point()+
    geom_hline(yintercept = 0, lty='dashed')+theme_bw()+geom_line(aes(group=Var1))+
    theme_bw()+theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
    theme(axis.title.x = element_blank(),
          legend.position = "bottom",
          # legend.position = c(.5,-.5),
          # plot.margin = unit(c(0,0,4,0), 'lines'),
          legend.title = element_blank())+
    # labs(y=("$\\widehat{\\betab}_1$"))+
    labs(y=titley)+
    guides(col=guide_legend(nrow=1,byrow=TRUE))+
    scale_color_manual(values = nucleotide_colours_logR)+
    # geom_segment(aes(x = 5, y = 0, xend = 5, yend = Inf),
    #              arrow = arrow(length = unit(0.2, "cm")), col='black')+
    # geom_segment(aes(x = 19, y = 0, xend = 18, yend = -Inf),
    #              arrow = arrow(length = unit(0.2, "cm")), col='black')+
    annotate("text", x = 5.5, y = -0.6, label="More clonal than baseline")+
    annotate("text", x = 18, y = 0.65, label="More subclonal than baseline")+
    ylim(c(-0.7, 0.8))
}
