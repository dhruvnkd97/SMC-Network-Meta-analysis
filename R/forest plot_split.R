library(dplyr)
library(forestplot)
library(grid)
library(svglite)


# ---------- 1) Dataset ----------
forest_split_plot <-  read.csv("data", "forest plot_split.csv")

forest_split_plot <- forest_split_plot %>%
  mutate(
    label_source = tools::toTitleCase(as.character(source)),
    label_ci     = sprintf("%.2f (%.2f, %.2f)", est, lo, hi),
    label_p      = ifelse(source == "network", fmt_p(p_incoherence), "")
  )


# ---------- 2) Label columns ----------
fmt_p <- function(p) ifelse(is.na(p), "",
                            ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p))) #Labelling p values



# Order rows: Direct → Indirect → Network
src_order <- c("direct","indirect","network")
forest_split_plot <- forest_split_plot %>% 
  mutate(source = factor(source, levels = src_order))  


# ---------- 3) Insert bold comparison headers and mark Network as summary (diamond) ----------
rows <- list()
for (cmp in unique(forest_split_plot$comparison)) {
  g <- forest_split_plot %>% filter(comparison == cmp) %>% arrange(source)
  # Header row (bold, no estimate)
  header <- tibble(
    comparison = cmp,
    source = factor(NA, levels = src_order),
    est = NA_real_, lo = NA_real_, hi = NA_real_,
    label_source = "", label_ci = "", label_p = "",
    is_summary = TRUE,  # bold header
    is_network = FALSE,
    boxsize = NA_real_
  )
  # Body rows (network row drawn as diamond/summary)
  g <- g %>%
    mutate(is_summary = (source == "network"),
           is_network = (source == "network"))
  rows[[length(rows)+1]] <- bind_rows(header, g)
}

tbl <- bind_rows(rows)



# ---------- 4) Compose label matrix (left & right columns) ----------
labeltext <- cbind(
  Comparison = ifelse(is.na(tbl$source), tbl$comparison, ""), # show cmp only on header rows
  Source     = tbl$label_source,
  `IRR [95% CrI]` = tbl$label_ci,
  `p (incoherence)`     = tbl$label_p
)

# ---------- 5) Vectors for forestplot ----------
mean  <- tbl$est
lower <- tbl$lo
upper <- tbl$hi
is.sum <- tbl$is_summary


# ---------- 6) Axis settings (log ratio) ----------
xmin <- min(lower[!is.na(lower)], na.rm = TRUE)
xmax <- max(upper[!is.na(upper)], na.rm = TRUE)
clip <- c(max(0.05, xmin/1.5), xmax*1.5)
ticks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
ticks <- ticks[ticks >= clip[1] & ticks <= clip[2]]


tbl$is_summary <- is.na(tbl$source)
row_type <- ifelse(is.na(tbl$source), "header", as.character(tbl$source))


#Draw functions
draw_direct <- function(mean, lower, upper, size, ...){
  fpDrawNormalCI(mean, lower, upper, size,
                 gp = gpar(fill = "grey75", col = "grey30"), ...)
}

draw_indirect <- function(mean, lower, upper, size, ...){
  fpDrawCircleCI(mean, lower, upper, size,
                 gp = gpar(fill = NA, col = "#1f78b4", lwd = 1.3), ...)
}

draw_network <- function(mean, lower, upper, size, ...){
  fpDrawDiamondCI(mean, lower, upper, size,
                  gp = gpar(fill = "black", col = "black"), ...)
}

draw_header <- function(mean, lower, upper, size, ...){
  NULL  # no symbol for comparison header rows
}

# ---------- 7) Draw ----------
forestplot(
  labeltext = labeltext,
  mean  = mean, lower = lower, upper = upper,
  is.summary = tbl$is_summary,
  xlog = TRUE, zero = 1, #clip = clip, 
  xticks = ticks,
  col = fpColors(box = "black", line = "black", summary = "black"),
  graph.pos = 3,
  lwd.ci = 1, lwd.zero = 1, vertices = TRUE,  
  txt_gp = fpTxtGp(
    label   = gpar(cex = 0.9),
    ticks   = gpar(cex = 0.9),
    xlab    = gpar(cex = 0.85),
    summary = gpar(fontface = "bold")         # bold comparison headers & diamonds
  ),
  lineheight = unit(0.6, "cm"),               # tighten vertical spacing
  colgap     = unit(6, "mm"),                 # gap between plot and right columns
  xlab = "Incidence Rate Ratio",
  title = "Incoherence in network estimates\n Antimalarials vs placebo/no SMC",
)
























