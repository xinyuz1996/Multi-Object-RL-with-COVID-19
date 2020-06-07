# Environment --------------------------------------------------
rm(list = ls())

## Path
main_path = "./"
data_path = paste(main_path, "data.rds", sep = "")
out_path = paste(main_path, "output/", sep = "")

## Code
miceadds::source.all(paste(main_path, "codes/", sep = ""))

## Preparation    --------
data = readRDS(data_path)
name = names(data)
finalPlot <- function(t0 = 12,
           t_end = 120,
           alpha = .1,
           title.size = 12,
           MAP.only = F) {
    theta = getSIRMAP(data = data, t0 = 12)
    cat("R0:", theta[13:15], "\n")
    plot_one_city <- function(l) {
      cat(l, "\n")
      result = plotNewValidation(
        data = data,
        l = l,
        t0 = t0,
        theta = theta,
        t_end = (t_end + 9),
        alpha = alpha,
        date_break = "30 day"
      )
      plots_i = result$gg + theme(axis.title = element_blank()) +
        theme(
          plot.title = element_text(size = 10, margin = ggplot2::margin(0, 0, .5, 0)),
          plot.margin = unit(c(.25, 0, 0, 0), "lines")
        )
      return(plots_i)
    }
    ggs = lapply(c(1:6), plot_one_city)
    
    g = ggarrange(
      ggs[[1]],
      ggs[[2]],
      ggs[[3]],
      ggs[[4]],
      ggs[[5]],
      ggs[[6]],
      ncol = 3,
      nrow = 2,
      legend = "none",
      label.x = 0,
      label.y = 0,
      common.legend = T,
      align = "v"
    )
    yvar =  c("count")
    xvar = c("date")
    g1 = annotate_figure(
      g,
      left = text_grob(
        yvar,
        just = "centre",
        size = title.size,
        rot = 90,
        hjust = .5,
        vjust = 1.5
      ),
      bottom = text_grob(
        xvar,
        just = "centre",
        size = title.size,
        rot = 0,
        hjust = .48,
        vjust = .5
      )
    )
    return(g1)
  }




# GET RESULT FOR TABLE 1 -------------------------------------------
writeTable1 <- function(path) {
  sink(path)
  cat("---- Result for Table 1: ----", "\n")
  for (t0 in c(12,  61)) {
    a = getSIRMAP(
      data = data,
      t0 = t0,
      brief = FALSE,
    ) # t_begin = t_begin,
    cat("\n", "t0 = ", t0, "\n")
    print(round(a[1:dim(a)[1], c(1, 2)], 3))
    cat("------------------------------", "\n")
  }
  sink()
}
writeTable1(path = paste(out_path, "Table1", ".txt", sep = ""))


#  FIGURE 2 IN PAPER  ---------------------------------------------------------
gg_T = finalPlot(
  t0 = 12,
  t_end = 120,
  alpha = .01
)
ggsave(
  filename = paste(out_path, "Figure2.png", sep = ""),
  plot = gg_T,
  width = 12,
  height = 3,
  units = c("in"),
  dpi = 500
)


### GET COST described in C.1  --------
printCostinC1 <- function(path) {
  sink(path)
  dat.all = do.call(rbind, data)
  dat.all = dat.all[which(dat.all$date <= "2020-03-15"), 1:8]
  cat(paste("C_j", ":", sep = ""), " mean ", " sd", "\n", sep = " ")
  for (i in c(2, 3)) {
    C = round(c(mean(unlist(dat.all[dat.all$action == i, 8]), na.rm = T),
                sd(unlist(dat.all[dat.all$action == i, 8]), na.rm = T)), 3)
    cat(paste("C_", i, ":", sep = ""), C, "\n", sep = " ")
  }
  sink()
}
printCostinC1(path = paste(out_path, "Cost_in_C1.txt", sep = ""))
