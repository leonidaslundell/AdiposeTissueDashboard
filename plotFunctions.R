library(ggpubr)
library(ggplot2)
library(ggiraph)

####################################################
#functions

boxPlot <- function(highlight = NULL,
                    datBox, 
                    datRes){
  if(is.null(highlight)){
    g <- ggplot() +
      annotate(geom = "text", label = "Type in you favorite gene in the bar above,\nor click on the plots below to explore a gene", x = .6,y = .5, size = 18) +
      theme_void()
    
    return(g)
  }
  if(length(highlight)>1){
    g <- ggplot() +
      annotate(geom = "text", label = "Select a single gene", x = .6,y = .5, size = 18) +
      theme_void()
    
    return(g)
  }
  g <- ggplot(datBox[Symbol == highlight, ],
              aes(x = Group, y = value, fill = Diagnosis)) +
    geom_boxplot() +
    scale_x_discrete(labels = rep(c("PRE", "POST", "REC"), 2)) +
    scale_y_log10(labels = scales::comma) +
    xlab("") +
    ylab("Expression (counts per million)") +
    theme_bw() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 18))
  
  if(nrow(datRes[Symbol == highlight & value != "ns", ])>0){
    sigs <- nrow(datRes[Symbol == highlight & value != "ns", ])
    positions <- datBox[Symbol == highlight, ]$value |> ceiling() |> sort(decreasing = T) |> log10()
    g <- g +
      ggpubr::stat_pvalue_manual(datRes[Symbol == highlight & value != "ns", ],
                                 label = "value",
                                 y.position = if(sigs==1){positions[1]}else{positions[1:sigs]},
                                 label.size = 6, 
                                 bracket.size = 1, bracket.nudge.y = 0.5)
    g
  }
  g <- g + ggtitle(highlight)
  return(g)
}

exp2logFC <- function(aveExp,
                      logFC,
                      sig,
                      Symbol,
                      title = "", 
                      xZoom = c(NA, NA),
                      yZoom = c(NA, NA),
                      highlight = NULL){
  
  y <- data.table(Symbol,
                  aveExp,
                  logFC,
                  sig,
                  stringsAsFactors = F)
  y <- y[!is.na(y$aveExp),]
  
  g <- ggplot(mapping = aes(x = aveExp, y = logFC, tooltip = Symbol, data_id = Symbol)) +
    geom_point(data = y[y$sig == "*",], color = "black", size = 1.4) +
    geom_point_interactive(data = y[y$sig == "ns",], color = "gray", alpha = .5, size = 0.5) +
    geom_point_interactive(data = y[y$sig == "*",], color = "red", alpha = .9, size = 1) +
    geom_hline(yintercept = 0, color = "gray") +
    scale_x_log10() +
    xlab("Expression (counts per million)") + 
    ylab("logFC") +
    # coord_cartesian(xlim = xZoom, ylim = yZoom) +
    theme_bw() +
    theme(text = element_text(size = 14),
          plot.margin = margin(0,0,0,0))
  g
  if(!is.null(highlight)){
    g <- g + geom_point_interactive(data = y[Symbol %in% highlight,], size = 6, color = "blue")
  }
  
  return(g)
}

girafize <- function(g){
  girafe(ggobj = g, 
         options = list(
           opts_sizing(width = .7),
           opts_zoom(max = 5),
           opts_toolbar(saveaspng = F),
           opts_selection(
             type = "single")
         ))
}

