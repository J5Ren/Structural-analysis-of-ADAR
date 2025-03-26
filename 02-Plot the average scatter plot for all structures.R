files=list.files('./') #./ is directory where the transformed pdb file in
pdb=read.table(paste0('./',files[1]))
pdb=pdb[pdb$V4%in%'CA',]
pdb$id=paste0(pdb$V3,pdb$V9)
for(i in 2:length(files)){
  tmp=read.table(paste0('./',files[i]))
  tmp=tmp[tmp$V4%in%'CA',]
  tmp$id=paste0(tmp$V3,tmp$V9)
  for(i in unique(tmp$id)){
    pdb[pdb$id%in%i,11]=pdb[pdb$id%in%i,11]+tmp[tmp$id%in%i,11]
    pdb[pdb$id%in%i,12]=pdb[pdb$id%in%i,12]+tmp[tmp$id%in%i,12]
    pdb[pdb$id%in%i,13]=pdb[pdb$id%in%i,13]+tmp[tmp$id%in%i,13]
  }
}
pdb[,11]=pdb[,11]/length(files)
pdb[,12]=pdb[,12]/length(files)
pdb[,13]=pdb[,13]/length(files)
deaminase_new=pdb[pdb$V9%in%c(324:724)&pdb$V4%in%"CA",]
RBM=pdb[pdb$V9%in%c(1:323)&pdb$V4%in%"CA",]
deaminase_new$label='deaminase_new'
RBM$label='RBM_new'
df=rbind(deaminase_new,RBM)
pdf("Point.pdf")
ggplot(df, aes(x = df[, 11], y = df[, 12], color = df[, 'label'])) +
  geom_point()+
  scale_color_manual(values = c("#48D1CC",
                                '#D2B48C',
                                '#D2B48C',
                                '#D2B48C'))+
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 0))+
  xlim(-max(abs(df[,11])),max(abs(df[,11])))+
  ylim(-max(abs(df[,12])),max(abs(df[,12])))+
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  xlab('')+ylab('')
dev.off()