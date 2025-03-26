###################Obtain the pdb coordinates after coordinate transformation###############################
library(optparse)
parser=OptionParser()
parser=add_option(parser,c("-n","--name"),help = 'file name')
name=parse_args(parser)$name
# pdb preprocess ---------------------------------------------------------
#name is input_file: grep ATOM pdb_file > input_file
pdb=read.table(name)
pdb_pro=pdb[pdb$V4%in%"CA",]
deaminase=pdb_pro[pdb_pro$V9%in%c(324:724),]
RBM1=pdb_pro[pdb_pro$V9%in%c(1:69),]
RBM2=pdb_pro[pdb_pro$V9%in%c(112:180),]
RBM3=pdb_pro[pdb_pro$V9%in%c(224:292),]
dsRNA=pdb[pdb$V4%in%"P",]

# define x-axis -----------------------------------------------------------
dsRNA_point=dsRNA[,11:13]
centroid <- colMeans(dsRNA_point)  # Calculate the centroid
X <- sweep(dsRNA_point, 2, centroid)
svd_res <- svd(X)
x_axis <- svd_res$v[,1]  #x-axis

# define y-axis -----------------------------------------------------------
P0 <- centroid  
deaminase_center=colMeans(deaminase[,11:13]) #deaminase center point
v <- deaminase_center - P0
proj_length <- sum(v * x_axis)
foot_point <- P0 + proj_length * x_axis #foot point coordinates
if(abs(sum((foot_point-deaminase_center)*x_axis))<0.00000001){
  print('foot point right!')
}else{
  print('foot point wrong!')
}
y_axis=(foot_point-deaminase_center)/sqrt(sum((foot_point-deaminase_center)**2)) #y axis

# define z-xais -----------------------------------------------------------
library(pracma)
z_axis <- cross(x_axis, y_axis)
z_axis <- z_axis / sqrt(sum(z_axis^2))



# transformer matrix ------------------------------------------------------
transformation_matrix <- rbind(x_axis, y_axis, z_axis)
new_coordinates <- t(apply(pdb[,11:13], 1, function(point) {
  transformed_point <- point - foot_point 
  transformed_point %*% t(transformation_matrix)
}))
new_coordinates[,1:2]=-new_coordinates[,1:2]
pdb_new=pdb
pdb_new[,11:13]=new_coordinates
pdb_new_dsRNA=pdb_new[pdb_new$V4%in%"P",]
# if(max(pdb_new_dsRNA$V11)<abs(min(pdb_new_dsRNA$V11))){
#   pdb_new$V11=-pdb_new$V11
# }
if(mean(pdb_new[pdb_new$V6%in%"LEU"&pdb_new$V9%in%353,'V11'])<0){
  pdb_new$V11=-pdb_new$V11
}
if(mean(pdb_new[pdb_new$V6%in%"LEU"&pdb_new$V9%in%353,'V12'])<0){
  pdb_new$V12=-pdb_new$V12
}
write.table(pdb_new,paste0('./',name),row.names = F,col.names = F,quote = F)