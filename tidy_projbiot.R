library(tidyverse)
library(readxl)
# library(magrittr)
library(egg)
library(ggtern)


### We define the projection matrix
projectionMatrix <- matrix(c(3,0,2,0,
                              1,0,1,0,
                              1,1,0,0,
                              1,0,1,3),
                            byrow=F,nrow=4)
rownames(projectionMatrix)<-c("Al","Ca","NaK","FM")
colnames(projectionMatrix)<-c("ms1","fsp","CaAl","bio")

invProjMat <- solve(projectionMatrix)

# This function will do the job
# It can probably be generalized to work with all projection matrices
# And/or wirtten in a nicer way... 
projBiotTidy <- function(Al,Ca,NaK,FM){
  x<- cbind(Al,Ca,NaK,FM)

  prd<-apply(x,MARGIN=1,FUN=function(z){
    return(invProjMat%*%z)
  })
  prd<-t(prd)
  colnames(prd)<-c("ms1","fsp","CaAl","bio")
  
  data.frame(prd)
}

##### START HERE #####
# This line has to be adapted to your file, of course
# Change the path, and perhaps also the loading options to fit your file
mydata <- read_excel("U:\\Recherche\\@Databases\\all_experiments.xlsx") 

# The following has to be adjusted based on the actual colnmaes in your file (and on the way Fe is defined)
# There are no safeguards, so it WILL fail if you are not careful
# You have been warned....
mydata %>%
  # Correct for Fe, which in the case of Sazava comes as FeO and Fe2O3
  # Adapt to your case... 
  replace_na(list(FeO=0,Fe2O3=0)) %>%
  mutate(FeOt = FeO + Fe2O3/1.111) %>%
  # We first calculate millications
  mutate(Si.m = SiO2 / 60.078 * 1000,
         Ti.m = TiO2 / 79.865 * 1000,
         Al.m = Al2O3 / 101.961 * 2 * 1000,
         Fe.m = FeOt / 71.839 * 1000,
         Mn.m = MnO / 70.937 * 1000,
         Mg.m = MgO / 40.299 * 1000,
         Ca.m = CaO / 56.077 * 1000,
         Na.m = Na2O / 61.979 * 2 * 1000,
         K.m = K2O / 94.195 * 2 * 1000,
         Al = Al.m,
         Ca = Ca.m,
         NaK = Na.m + K.m,
         FM = Fe.m + Mg.m) %>% 
  # Actuel coordinate mapping step
  mutate(projBiotTidy(Al,Ca,NaK,FM)) %>%  
  # Ploting coordinates 
  mutate(ss = ms1 + fsp + CaAl,
         aa = CaAl / ss,
         bb = ms1 / ss,
         cc = fsp / ss,
         xcoord = sqrt(3)/2*cc,
         ycoord = (bb-aa)/2 ) %>%
  # Angle from feldspar
  mutate(ang=atan(-ycoord/(xcoord-sqrt(3)/2))*180/pi) %>%
  # Store safely !
  {.} -> mydataProjected


write_csv(mydataProjected,"U:\\Recherche\\@Databases\\all_experiments_withangle.csv")

# Plotting with ggplot (why not ?)

mydataProjected %>%
  # Basic plot definition
  ggplot(aes(x=xcoord,y=ycoord))+
  # We plot points
  # Color and Symbol are taken from the file (i.e. the columns must exist, or else !)
  geom_point(aes(color=`Src_type(Jensen)`))+
  # We use the col and pch values found in the file, no mapping is done
  #scale_color_identity()+
  # Draw the decorations
  annotate("path",x=c(-2,1),y=c(0,0),color="black")+ # Horiz line
  annotate("path",x=c(0,0,sqrt(3)/2,0),y=c(-1/2,1/2,0,-1/2),color="black")+ # Triangle
  annotate("path",x=c(-2,sqrt(3)/2),y=c(1/6*(1+2/(sqrt(3)/2)),0),col="grey",lty="dashed" )+ # Feldspar line
  # Text
  annotate("text",x=0,y=-0.55,label="Ca+Al",adj=0.5)+
  annotate("text",x=0,y=0.55,label="3 Al+2(Na+K)",adj=0.5)+
  annotate("text",x=sqrt(3)/2+.03,y=0+0.03,label="Al+(Na+K)",adj=0)+
  # We fix the plot boundaries, and critically set aspect ratio to 1
  coord_fixed(ratio=1,xlim=c(-2,1),ylim=c(-1/2-0.1,1/2+0.1))+
  # Niceties, to remove the "true" axes and such
  theme_article()+theme(axis.line = element_blank(),
                         panel.border = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank())


# Or, even better...
mydataProjected %>%
  # Basic plot definition
  ggtern(aes(x=CaAl,y=ms1,z=fsp))+
  geom_point(aes(color=`Src_type(Jensen)`))+
  theme_article()+
  theme_nogrid()+theme_ticksoutside()+theme_ticklength(major=unit(5,"mm"))+
  annotate(geom="line",x=c(0.5,0),y=c(0.5,0),z=c(0,1))+
  theme_rotate(degrees = 30)
  


  