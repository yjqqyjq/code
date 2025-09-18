library("gadgetry")
library("rhdf5")
S100=h5read("/Users/24756376/data/Flamingo/L1000N0900/S_compare.hdf5","S_r100")
Sub=h5read("/Users/24756376/data/Flamingo/L1000N0900/S_compare.hdf5","S_cen_unbound")
for (i in 1:100){
  file=paste0("/Users/24756376/data/Flamingo/L1000N0900/halos/-",as.character(i),"_1.0_r100.hdf5")
  sn=readsnapshot(file)
  sn$PartType1$lum=1.5
  
  plot4(sn,length.unit='mpc',kde=4,width=8,npixel=300,center=c(0,0,0), ,arrows="True",scale = "True",
        pdffile=paste0("/Users/24756376/plot/Flamingo/L1000N0900/visualization/r100/",as.character(i),'_',as.character(round(Sub[i+1],2)),"_",as.character(round(S100[i+1],2)),".pdf"))
}
#  plot4(sn,length.unit='mpc',kde=4,width=7,npixel=300,center=c(0,0,0), ,arrows="True",scale = "True")


