glimpse(contrasts_shape)

shape_results <- 



library(dplyr)
library(tidyr)
library(tidybayes)
library(purrr)



dominant_patterns |>  print(n=30)

library(dplyr)
library(ggplot2)
library(ggh4x)




library(ggplot2)
library(ggh4x)
library(ggtext)
library(grid)

fig_4c <- ggplot(plot_df,
       aes(x = Predictor, y = shape)) +
  # centered label
  geom_label(
    aes(label = glue::glue("{prob_lab}"), colour = direction),
    hjust=.5,
    fill = "white",
  ) +
  ggh4x::facet_nested(
    response   ~ facet +class ,
    scales = "free_x",
    space  = "free_x",
    strip = strip_nested(
      text_y = elem_list_text(face = "bold", size = c(10, 9)),
      by_layer_y = TRUE,
      text_x = elem_list_text(face = "bold", size = c(10, 9)),
      by_layer_x = TRUE
    )
  ) +
  
  scale_colour_manual(
    values = direction_colors
  ) +
  
#  scale_x_continuous(limits = c(0.84, 1.12), expand = c(0,0)) +
  theme_bw(base_family = "sans", base_size = 12) +
  theme(
    strip.text.x =  element_text(face="bold",margin=margin(t=5)),
    strip.text.y = element_text(size=10,face="bold",angle=270, margin=margin(l=5)),
    axis.title.x       = element_blank(),
    axis.text.y        = element_text(size = 10, hjust = 1),
    axis.text.x        = ggtext::element_markdown(size = 8, margin=margin(t=3),angle=30,vjust=1),
    axis.title.y       = element_blank(),
    axis.line.x        = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.x       = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.length  = unit(0.2, "lines"),
    plot.tag           = element_text(face = "bold", size = 14),
    panel.spacing.x = unit(2,"lines"),
    plot.margin = margin(0,0,0,0),
    legend.position="bottom",
    legend.title = element_text(face="bold"),
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect_round(
      fill = NA, 
      color = "black", 
      radius = unit(0.5, "cm") # Controla o arredondamento
    )
  ) +
  guides(
    x = guide_axis(cap = TRUE)
  ) + labs(
           colour="Dominant direction",
           tag="c")


fig_4c
remotes::install_github("teunbrand/elementalist")
library(elementalist)
library(ggplot2)

dev.off()
ggplot(mpg, aes(cty, hwy)) +
  geom_point() +
  facet_wrap(~class) +
  theme_minimal() + 
  theme(
    # Remove as bordas individuais e o fundo do strip
    strip.background = element_blank(),
    panel.border = element_blank(),
    
    # Aplica a borda arredondada ao fundo de cada painel
    panel.background = element_rect_round(
      fill = "grey95", 
      color = "black", 
      radius = unit(0.5, "cm") # Controla o arredondamento
    ),
    
    # Adiciona espaçamento para que as caixas não encostem uma na outra
    panel.spacing = unit(1, "lines")
  )
