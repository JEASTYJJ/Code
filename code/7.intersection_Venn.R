inter<-intersect(modProbes,cygene)
write.table(inter,"intersect_mod_DEG_gene.txt",row.names = F,sep = '\t')

library(ggvenn)

a <- list("Modle"=modProbes,"DEG"=cygene)

opar <- par(family = "Roboto Condensed")
ggvenn(a,fill_color=c("turquoise", "orangered2"),fill_alpha = .7,
       stroke_linetype = "solid",
       stroke_color = "white",set_name_size = 5,
       text_size=4) 
ggsave("intersect_DEG.pdf",width = 6,height = 6)

b<-list("Modle&DEG"=inter,"PPI>700"=unique(ppi1$protein1))
ggvenn(b,fill_color=c("orangered4", "cyan2"),fill_alpha = .7,
       stroke_linetype = "solid",
       stroke_color = "white",set_name_size = 5,
       text_size=4) 

ggsave("DEG_intersect_PPI>700.pdf",width = 6,height = 6)

