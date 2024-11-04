library(RColorBrewer)
library(viridis)

# palettes
ylrd <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "OrRd"))(n = 100)
rdbu <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(n = 100))
rdbu2   <- colorRampPalette(c("navy", "white", "red"))(n = 100)
rdbu3   <- colorRampPalette(c("navy", "#2E2E97", "#4848A4", "white",
                              "#FF4E4E", "#FF3434", "red"))(n = 100)

# from https://github.com/gabrielakinker/CCLE_heterogeneity/blob/master/custom_magma.R
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))


# white -> purple
purples <- colorRampPalette(brewer.pal(9, "Purples"))(n = 100)

# gray -> red
gryrd   <- grDevices::colorRampPalette(c("gray83", "#E09797", "red"))(n = 200)

palette_tumor_normal_mm <- c("Tumor" = "red4", "Normal" = "seagreen", "Doublet" = "cornflowerblue")
palette_malig_normal    <- c("Malignant" = "black", "Likely normal" = "gray70", "Normal" = "gray40")

palette_groups <- c("NB-FOXR2"         = "red",
                    "DIPG-H3K27M-FOXR2"= "turquoise3", 
                    "DIPG-H3K27M"      = "darkgreen",
                    "HGG-IDH"          = "darkolivegreen4", 
                    "EP-PFA"           = "darkolivegreen2",
                    "HGG-H3.3G34R/V"   = "turquoise4",
                    "ETMR"             = "darksalmon",
                    "MB-SHH"           = "hotpink4",
                    "MB-WNT"           = "darkred",
                    "EC-NB"            = "gray70",
                    "EC-NB (MYCN-amp)" = "gray40",
                    "Unknown"          = "gray90")

palette_nb_mycn <- c("Amp"    = "gray40",
                     "NonAmp" = "gray70")

palette_nb_epi <- c("MES"    = "#420A68",
                    "MYCN"   = "#2FB47C",
                    "MNA-HR" = "#b22222",
                    "MNA-LR" = "#00afaf",
                    "Unknown" = "gray90")

palette_class <- c("Radial glia" = "#ffcc00",
                   "Neural stem cell" = "#ffcc00",
                   "Proliferating progenitors" = "#e8a805",
                   "Neuronal progenitors" = "#ffbda3",
                   "Excitatory neurons" = "darkred",
                   "GE neural precursors" = "#6885e3",
                   "GE/Inhibitory neurons" = "#234edb",
                   "Striatal neurons" = "#3c214a",
                   "Other neurons" = "#e3bfcd",
                   "Neural crest lineages" = "#b580d1",
                   "Oligodendrocyte precursors" = "#b4e04e",
                   "Oligodendrocytes" = "#57a607",
                   "Astrocytes" = "#00a385",
                   "Ependymal" = "#8ee5cf",
                   "Glial progenitors" = "#d5d98b",
                   "Immune" = "#aca2b2",
                   "Meninges" = "#6e5335",
                   "Mesenchyme" = "#ffe5af",
                   "Mural" = "gray40",
                   "Muscle" = "#c7a7b4",
                   "Unlabelled" = "gray70", 
                   "Other" = "gray90")


#' Apply a clean theme to a ggplot2 object
#'
#' @references https://github.com/sjessa/ggmin
#'
#' @importFrom ggplot2 theme_light theme
theme_min <- function(base_size = 11, base_family = "",
                      border_colour = "grey90",
                      border_size = 1) {
  
  theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = border_colour, size = border_size),
      axis.ticks = element_line(colour = border_colour),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black", size = rel(1.2)),
      strip.text.y = element_text(colour = "black", size = rel(1.2)),
      title = element_text(size = rel(1.2), face = "bold"),
      axis.text = element_text(colour = "black", size = rel(1.2)),
      axis.title = element_text(colour = "black", size = rel(1.2), face = "plain"),
      legend.title = element_text(colour = "black", size = rel(1.2), face = "plain"),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)
    )
}


#' Rotate the x axis labels in a ggplot
#'
#' @param angle Integer, value in degrees to rotate labels. Default: 90.
#'
#' @return A theme element to rotate labels
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(colour = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + rotate_x()
rotate_x <- function(angle = 90) {
  
  theme(axis.text.x = element_text(angle = angle, hjust = 1))
  
}


#' Remove the legend in a ggplot
#'
#' @return A theme element to hide legend
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(colour = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + no_legend()
no_legend <- function() {
  
  theme(legend.position = "none")
  
}



#' modify the colour legend - by default, just makes the keys bigger and opaque,
#' in one column. should be used with geom_*(aes(colour = ___))
#' 
#' @param n numeric, size of keys
#' @param ncol numeric, number of columns
#' @param title set to NULL to hide legend title
mod_legend_col <- function(n = 3, ncol = 1, title = waiver()) {
    
    guides(colour = guide_legend(title = title,
                                 ncol = ncol,
                                 override.aes = list(size = n, alpha = 1)))
}


#' modify the fill legend - by default, just makes the keys bigger and opaque,
#' in one column. should be used with geom_*(aes(fill = ___))
#' 
#' @param n numeric, size of keys
#' @param ncol numeric, number of columns
#' @param title set to NULL to hide legend title
mod_legend_fill <- function(n = 3, ncol = 1, title = waiver()) {
    
    guides(fill = guide_legend(title = title,
                                 ncol = ncol,
                                 override.aes = list(size = n, alpha = 1)))
}



set_qual_pal <- function(n) {
  
  pal_qual <- c("#646BA8",
                "#ef893b",
                "#e2445e",
                "#5E7A41",
                "#FFFC63",
                "#5DAD3B",
                "#6E3688",
                "#A2ACD3",
                "#f2c7d8",
                "#62babd",
                "#519674",
                "#f0c992",
                "#BEBEBE",
                "#2e3082",
                "#61cfe8")
  
  if (n <= length(pal_qual)) { pal_qual_ramped <- head(pal_qual, n)}
  else {
    
    pal_qual_ramped <- colorRampPalette(pal_qual)(n = n)
    
  }
  
  return(pal_qual_ramped)
  
}


#' Remove axis ticks and tick labels from a ggplot
#'
#' @return A theme element to remove ticks
#' @export
#'
#' @author Selin Jessa
no_ticks <- function() {
  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
  
}

color_legend_ncol <- function(n) guides(color = guide_legend(ncol = n))
fill_legend_ncol <- function(n) guides(fill = guide_legend(ncol = n))


#' Adapted from Samantha Worme
kable_format_numeric <- function(df, colour = TRUE) {
  
  require(kableExtra)
  
  if (colour) {
    
    df %>%
      mutate_if(is.numeric, function(x) {
        cell_spec(x, bold = T,
                  color = spec_color(x, end = 0.9),
                  font_size = spec_font_size(x))
      }) %>%
      kable(escape = F, align = "c") %>%
      kable_styling(c("striped", "condensed"), full_width = F)
    
  } else {
    
    df %>%
      kable(escape = F, align = "c") %>%
      kable_styling(c("striped", "condensed"), full_width = F)
    
  }
  
}

