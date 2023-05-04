library(ggpubr)
library(ggplot2)
library(reshape2)

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

    data = data[data$software != "miso_preprocess", ]

    ncol = 2

    p1 = ggboxplot(
      data, x = x, y = "time",
      fill = "software",
      xlab = xlab,
      ylab="Time (s)",
      facet.by = ifelse(x == "num_of_files", "n_jobs", "num_of_files" ),
      palette = c(
        "ggsashimi"="#00AFBB",
        "trackplot"="#E7B800",
        "miso"="#FC4E07",
        "miso_preprocess"="lightgrey"
      ),
      # add = "jitter",
      scales="free_y", size=.3,
    ) +
      ylim(0, NA) +
      grids(linetype = "dashed")


    p2 = ggboxplot(
      data, x = x, y = "memory",
      fill = "software",
      xlab = xlab,
      ylab = "Max Memory (MB)",
      facet.by = ifelse(x == "num_of_files", "n_jobs", "num_of_files" ),
      palette = c(
        "ggsashimi"="#00AFBB",
        "trackplot"="#E7B800",
        "miso"="#FC4E07",
        "miso_preprocess"="lightgrey"
      ),
      # add = "jitter",
      scales="free_y", size=.3,
    ) +
      ylim(0, NA) +
      grids(linetype = "dashed")


    p3 = ggboxplot(
      data, x = x, y = "memory_per_process",
      fill = "software",
      xlab = xlab,
      ylab = "Max Memory Per process (MB)",
      facet.by = ifelse(x == "num_of_files", "n_jobs", "num_of_files" ),
      palette = c(
        "ggsashimi"="#00AFBB",
        "trackplot"="#E7B800",
        "miso"="#FC4E07",
        "miso_preprocess"="lightgrey"
      ),
      # add = "jitter",
      scales="free_y", size=.3,
    ) +
      ylim(0, NA) +
      grids(linetype = "dashed")

    if (log10) {
        p1 = p1 + scale_y_log10()
        p2 = p2 + scale_y_log10()
        p3 = p3 + scale_y_log10()
    }

    ggsave(
      cowplot::plot_grid(p1, p2, p3, ncol = 1),
      filename = paste0(outfile, ".pdf"),
      width = 4 * ncol, height = 12,
      limitsize = FALSE
    )

    ## scaled
    data = data[data$software != "miso_preprocess", ]
    scale_usage = function(data, value.var = "time") {
      temp = dcast(
        data,
        event+num_of_files+n_jobs~software,
        value.var = value.var,
        fun.aggregate = mean,
        fill = 0
      )

      for (i in c("ggsashimi", "miso")) {
        temp[, i] = log10(temp[, i] / temp[, "trackplot"] + 1)
      }
      temp = melt(temp[, colnames(temp) != "trackplot"], id.vars = c("event", "num_of_files", "n_jobs"))
      return(temp)
    }

    temp = scale_usage(data)
    p1 = ggboxplot(
      temp, x = x, y = "value",
      fill = "variable",
      xlab = xlab,
      ylab="log10(Time (s) / trackplot)",
      facet.by = ifelse(x == "num_of_files", "n_jobs", "num_of_files" ),
      palette = c(
        "ggsashimi"="#00AFBB",
        "trackplot"="#E7B800",
        "miso"="#FC4E07",
        "miso_preprocess"="lightgrey"
      ),
      scales="free_y", size=.3
    ) +
      grids(linetype = "dashed")

    temp = scale_usage(data, value.var = "memory")
    p2 = ggboxplot(
      temp, x = x, y = "value",
      fill = "variable",
      xlab = xlab,
      ylab = "lgo10(Max Memory (MB) / trackplot)",
      facet.by = ifelse(x == "num_of_files", "n_jobs", "num_of_files" ),
      palette = c(
        "ggsashimi"="#00AFBB",
        "trackplot"="#E7B800",
        "miso"="#FC4E07",
        "miso_preprocess"="lightgrey"
      ),
      # add = "jitter",
      scales="free_y", size=.3,
    ) +
      grids(linetype = "dashed")

    temp = scale_usage(data, value.var = "memory_per_process")
    p3 = ggboxplot(
      temp, x = x, y = "value",
      fill = "variable",
      xlab = xlab,
      ylab = "log10(Max Memory Per process (MB) / trackplot)",
      facet.by = ifelse(x == "num_of_files", "n_jobs", "num_of_files" ),
      palette = c(
        "ggsashimi"="#00AFBB",
        "trackplot"="#E7B800",
        "miso"="#FC4E07",
        "miso_preprocess"="lightgrey"
      ),
      # add = "jitter",
      scales="free_y", size=.3,
    ) +
      grids(linetype = "dashed")

    ggsave(
      cowplot::plot_grid(p1, p2, p3, ncol = 1),
      filename = paste0(outfile, "_scaled.pdf"),
      width = 4 * ncol, height = 12,
      limitsize = FALSE
    )
}

plot(
    path="./benchmark_*/stats.txt",
    x="num_of_files",
    xlab="Number of files",
    outfile="./n_files",
    log10 = TRUE
)
plot(
    path="./benchmark_*/stats.txt",
    x="n_jobs",
    xlab="Number of processes",
    outfile="./n_process",
    log10 = TRUE
)