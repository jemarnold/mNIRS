library(hexSticker)
library(sysfonts)
library(showtext)

font_families_google()[which(grepl("Merriweather", font_families_google()))]
font_add_google("Merriweather Sans", "merriweather_sans")
font_add_google("Merriweather", "merriweather")

showtext_auto()

hexSticker::sticker(
    ## imported image
    subplot = r"(C:\R-Projects\mnirs\man\figures\mnirs-hex-transparent-dark3.png)",
    s_x = 1, # centred horizontally
    s_y = 0.9, # slightly below centre to leave room for text
    s_width = 0.84, # fill most of the hex width
    s_height = 0.7, # constrain height proportionally

    ## package name parameters
    package = "mnirs",
    p_x = 1,
    p_y = 1.4,
    # p_size = 12, ## for .svg
    p_size = 35, ## for .png
    p_color = "#004992",
    p_family = "merriweather_sans",
    p_fontface = 2, # "plain", "bold", "italic", "bolditalic"

    ## hex geometry
    h_fill = "#bcced7", ## light blue-grey background
    # h_color = "#fb653a", ## orange of muscle icon
    # h_color = "#1a284f", ## dark navy
    h_color = "#004992", ## blue
    h_size = 1.5,
    white_around_sticker = FALSE,

    ## spotlight / white lines echo — leave off unless you want glow
    # spotlight = TRUE,
    # l_x = 1,
    # l_y = 0.35,

    ## ggsave
    ## `*mnirs-hex*.png` will be gitignored
    filename = "man/figures/mnirs-hex-temp.png",
    dpi = 300
)
