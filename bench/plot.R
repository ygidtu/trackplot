library(ggpubr)
library(ggplot2)

plot <- function(path, x, xlab, outfile, header = TRUE, log10 = FALSE) {
    data = lapply(Sys.glob(path), function(x) {
      res = tryCatch({
        return(read.table(x, sep = "\t", header = header))
      }, error = function(e) {
        print(e)
        return(NULL)
      })
      return(res)
    })
    data = do.call(rbind, data)

    if (!header) {
        colnames(data) <- c("event", "time", "memory", "software", "num_of_files", "n_jobs")
    }

    data$memory = data$memory / 1024 / 1024
    data$memory_per_process = data$memory / data$n_jobs
    data[data$software == "ggsashimi", "memory_per_process"] = data[data$software == "ggsashimi", "memory"]

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

    p3 = ggboxplot(
      data,
      x = x,
      y = "memory_per_process",
      color = "software",
      xlab = xlab,
      ylab = "Max Memory Per process (MB)",
      facet.by = "event",
      palette = c("#00AFBB", "#E7B800", "#FC4E07"),
      ncol = ncol, add = "jitter",
      scales="free_y", size=.3
      # position=position_dodge(width = .3)
    ) + ylim(0, NA) + grids(linetype = "dashed")

    if (log10) {
        p1 = p1 + scale_y_log10()
        p2 = p2 + scale_y_log10()
        p3 = p3 + scale_y_log10()
    }

    ggsave(
      cowplot::plot_grid(p1, p2, p3, ncol = 1),
      filename = outfile,
      width = 4 * ncol, height = 12
    )
}

plot(path="./benchmark_files/*/stats.txt", x="num_of_files", xlab="Number of files", outfile="./n_files.pdf")
plot(path="./benchmark_threads/*/stats.txt", x="n_jobs", xlab="Number of processes", outfile="./n_process.pdf", log10 = TRUE)