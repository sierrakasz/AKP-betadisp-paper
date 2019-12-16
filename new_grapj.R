library(grid)
library(pBrackets)

sam_numb <- as.numeric(1:10)
beta_disp <- c(rep('Full communities', 5), rep('Random forest indicators', 5))
body_site <- c(rep(c('Nose', 'Mouth', 'Nose', 'Ears', 'Mouth'), 2))
of_Deaths <- c(rep(c('MOD', 'MOD', 'COD', 'COD', 'COD'), 2))
mcfadden <- c(0.277, 0.250, 0.361, 0.304, 0.298,
              0.310, 0.255, 0.377, 0.359, 0.291)
percent_correct <- c(61.0, 55.3, 62.8, 62.9, 60.6,
                     61.0, 54.1, 59.9, 41.3, 54.7)
n <- c(172, 170, 172, 167, 170,
       172, 170, 172, 167, 170)

df_for_graph <- as.data.frame(cbind(sam_numb, beta_disp, body_site, of_Deaths, mcfadden,
                                    percent_correct, n))
df_for_graph$sam_numb <- factor(df_for_graph$sam_numb, levels = c(1,2,3,4,5
                                                                     ,6,7,8,9,10))
df_for_graph$percent_correct <- as.numeric(as.character(df_for_graph$percent_correct))
df_for_graph$mcfadden <- as.numeric(as.character(df_for_graph$mcfadden))
df_for_graph$of_Deaths <- factor(df_for_graph$of_Deaths, levels = c('MOD', 'COD'))


grob <- grobTree(textGrob("Full Communities", x=0.07,  y=0.02, 1, hjust=0,
                          gp=gpar(col="black", fontsize=14)))
grob2 <- grobTree(textGrob("RF Indicators", x=0.62,  y=0.02, hjust=0,
                          gp=gpar(col="black", fontsize=14)))

theme_set(theme_classic(base_size = 20))
newlabs <- c("n=172", "n=170", "n=172", "n=167", "n=170",
                    "n=172", "n=170", "n=172", "n=167", "n=170")

p <- ggplot(df_for_graph, aes(sam_numb, percent_correct, fill = body_site)) + geom_col() +
  scale_fill_manual("Body Sites", values = c("Nose" = "#37F196", "Mouth" = "#F1375C", "Ears" = "#F1B337")) +
  facet_wrap(~of_Deaths, scales = 'free_x', strip.position = "bottom") + 
  ylab('Percent Correct') + xlab("") +
  annotation_custom(grob) + annotation_custom(grob2) + 
  scale_x_discrete(labels= newlabs) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
  
df_for_graph$sam_numb <- factor(df_for_graph$sam_numb, levels = c(1,2,6,7,3,4,5,8,9,10))
p2 <- ggplot(df_for_graph, aes(sam_numb, mcfadden)) + geom_line(group = 10) +
  geom_point(aes(color = body_site), size = 4) +
  scale_color_manual("Body Sites", 
                     values = c("Nose" = "#37F196", "Mouth" = "#F1375C", 
                                "Ears" = "#F1B337")) +
  ylab('McFadden R-squared') + theme(axis.line.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank(),
                                     axis.title.x=element_blank()
                                    )

p2

tiff("rf_full_grpah.TIF", width = 3500, height = 3500, res=300)
ggarrange(p2, p, nrow = 2)
dev.off()
