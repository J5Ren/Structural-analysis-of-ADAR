library(ggplot2)
library(ggnewscale)
library(ggforce)
point_normal=function(pdb){
  pdb=read.table(paste0('./',pdb))
  pdb=pdb[pdb$V4%in%'CA',]
  deaminase=pdb[pdb$V9%in%c(324:724)&pdb$V4%in%"CA",]
  RBM=pdb[pdb$V9%in%c(1:323)&pdb$V4%in%"CA",]
  deaminase$label='deaminase'
  RBM$label='RBM'
  df=rbind(deaminase,RBM)
  return(df)
}
files=list.files('./')
df=point_normal(files[1])
for(i in 2:length(files)){
  tmp=point_normal(files[i])
  df=rbind(df,tmp)
}
library(ggpointdensity)
p=ggplot() +
  geom_pointdensity(aes(df_DRBM[,11], df_DRBM[,12]), adjust=100, size=1.5) +
  scale_color_gradientn(colors = colorRampPalette(c("#F7F2EB","#D2B48C"))(100)) +
  new_scale_color() +  # 确保在每个需要独立颜色映射的图层之前调用
  geom_pointdensity(aes(df_deam[,11], df_deam[,12]), adjust=100, size=1.5) +
  scale_color_gradientn(colors = colorRampPalette(c("#DFF7F6","#48D1CC"))(100)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  ) +
  geom_mark_hull(aes(df[,11],df[,12]),alpha=0.1, expand = unit(0.1, "inches")) +
  coord_cartesian(clip = "off") +
  coord_fixed(ratio = 1)
ggsave('heatmap.png',p)
