library("gadgetry")
library("rhdf5")
S100=h5read("/Users/24756376/data/Flamingo/L1000N0900/S_compare.hdf5","S_r100")
Sub=h5read("/Users/24756376/data/Flamingo/L1000N0900/S_compare.hdf5","S_cen_unbound")
r200=h5read("/Users/24756376/data/Flamingo/L1000N0900/halos_ranked.hdf5","r200")
r200=r200[r200!=0]
for (i in 1:50){
  file=paste0("/Users/24756376/data/Flamingo/L1000N0900/halos/-",as.character(i),"_5.0_r200.hdf5")
  sn=readsnapshot(file)
  sn$PartType1$lum=0
#  sn$PartType0$lum=0
#  sn$PartType1$col="white"
  sn$PartType0$value=log(sn$PartType0$Temperatures)
  sn$PartType0$col=rev(rainbow(100,end=2/3))
  plot4(sn,length.unit='mpc',kde=4,width=10*r200[i+1],npixel=300,center=c(0,0,0), ,arrows="True",scale = "True",
        pdffile=paste0("/Users/24756376/plot/Flamingo/L1000N0900/visualization/5_r200_gas_T/",as.character(i),'_',as.character(round(Sub[i+1],2)),"_",as.character(round(S100[i+1],2)),".pdf"))
}
#  plot4(sn,length.unit='mpc',kde=4,width=7,npixel=300,center=c(0,0,0), ,arrows="True",scale = "True")


