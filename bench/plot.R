library(ggpubr)
library(ggplot2)

plot <- function(path, x, xlab, outfile) {
    data = lapply(Sys.glob(path), function(x) {
      res = tryCatch({
        return(read.table(x, sep = "\t", header = T))
      }, error = function(e) {
        print(e)
        return(NULL)
      })
      return(res)
    })
    data = do.call(rbind, data)
    data$memory = data$memory / 1024 / 1024

    ncol = length(unique(data$event))

    p1 = ggboxplot(
      data, x = x, y = "time",
      color = "software",
      xlab = xlab,
      ylab="Time (s)",
      facet.by = "event",
      palette = c("#00AFBB", "#E7B800", "#FC4E07"),
      ncol = ncol, add = "jitter",
      scales="free_y", size=.3
      # position=position_dodge(width = .3),
    ) + ylim(0, NA) + grids(linetype = "dashed")


    p2 = ggboxplot(
      data,
      x = x,
      y = "memory",
      color = "software",
      xlab = xlab,
      ylab = "Max Memory (MB)",
      facet.by = "event",
      palette = c("#00AFBB", "#E7B800", "#FC4E07"),
      ncol = ncol, add = "jitter",
      scales="free_y", size=.3
      # position=position_dodge(width = .3)
    ) + ylim(0, NA) + grids(linetype = "dashed")

    ggsave(
      cowplot::plot_grid(p1, p2, ncol = 1),
      filename = outfile,
      width = 4 * ncol, height = 8
    )
}

plot(path="./benchmark_files/*/stats.txt", x="num_of_files", xlab="Number of files", outfile="./n_files.pdf")
plot(path="./benchmark_threads/*/stats.txt", x="n_jobs", xlab="Number of processes", outfile="./n_process.pdf")