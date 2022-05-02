#####################################################################
#                 Projection from biotite                           #
#                                                                   #
#####################################################################
# Plugin version, to be used in .../library/GCDkit/Diagrams/Plugins

projbiocoords<-function(where=WR,add=FALSE){
ee<-millications(where)
plane<-matrix(c(3,0,2,0,
				1,0,1,0,
				1,1,0,0),
				byrow=F,nrow=4)
rownames(plane)<-c("Al","Ca","NaK","FM")
colnames(plane)<-c("ms1","fsp","CaAl")
bio<-c(1,0,1,3)
names(bio)<-c("Al","Ca","NaK","FM")
aa<-cbind(plane,bio)

#Extract and calculate
ox<-cbind(ee[,"Al2O3"],ee[,"CaO"],ee[,"Na2O"]+ee[,"K2O"],ee[,"FeOt"]+ee[,"MgO"])
colnames(ox)<-c("ms1","fsp","CaAl","bio")

iaa<-solve(aa)
ox.p<-t(apply(ox,MARGIN=1,FUN=function(z){
return(iaa%*%z)
}))
colnames(ox.p)<-c("ms1","fsp","CaAl","bio")
results<<-ox.p

if(add){
addResults()
}

return(results)
}

.rotatedTernaryCoordinates<-function(where=WR,alab,blab,clab){
aa <- where[, alab]
bb <- where[, blab]
cc <- where[, clab]

ss<-aa+bb+cc
aa<-aa/ss
bb<-bb/ss
cc<-cc/ss

xx<-sqrt(3)/2*cc
yy<-(bb-aa)/2

return(cbind(aa,bb,cc,xx,yy))
}


projbioplot<-function(mins=FALSE,addWR=FALSE,ticks=FALSE,xmin=-2,xmax=1,ymin=-0.5,ymax=0.5){

### Rotated ternary !
coords<-.rotatedTernaryCoordinates(projbiocoords(where=WR,add=addWR),alab="CaAl",blab="ms1",clab="fsp")
x.data<<-coords[,"xx"]
y.data<<-coords[,"yy"]

temp1<-list(
        lines0=list("lines",x=c(0,0,sqrt(3)/2,0),y=c(-1/2,1/2,0,-1/2),col="black"), # triangle
		lines1=list("lines",x=c(xmin,1),y=c(0,0),col="black"), # horiz line
		lines2=list("lines",x=c(xmin,sqrt(3)/2),y=c(1/6*(1-xmin/(sqrt(3)/2)),0),col="grey",lty="dashed"), #feldspar line
        A=list("text",x=0,y=-0.55,text="Ca+Al",adj=0.5),
        B=list("text",x=0,y=0.55,text="3 Al+2(Na+K)",adj=0.5),
        C=list("text",x=sqrt(3)/2+.03,y=0+0.03,text="Al+(Na+K)",adj=0),
        GCDkit=list("NULL",plot.type="ternary",plot.position=666,plot.name="Projection from bt")
        )

## Ticks
if(ticks){
t1<-seq(0.1,0.9,0.1)

c1<-0
a1<-t1
b1<-1-t1

c2<-0.03
a2a<-t1-c2
b2a<-1-t1

a2b<-t1
b2b<-1-t1-c2

t.coords<-cbind(a2a,b2a,c2,a1,b1,c1,a2b,b2b,c2)

# Left hand side
tl1<-.rotatedTernaryCoordinates(where=t.coords,alab="a2a",blab="b2a",clab="c2")
tl2<-.rotatedTernaryCoordinates(where=t.coords,alab="a1",blab="b1",clab="c1")
tl3<-.rotatedTernaryCoordinates(where=t.coords,alab="a2b",blab="b2b",clab="c2")
tl<-cbind(tl1[,"xx"],tl1[,"yy"],tl2[,"xx"],tl2[,"yy"],tl3[,"xx"],tl3[,"yy"])


ee<-apply(tl,FUN=function(z) {
	paste("=list(\"lines\",x=c(",z[1],",",z[3],",",z[5],"),y=c(",z[2],",",z[4],",",z[6],"),col=\"black\")",
	sep="")
},MARGIN=1)
ee<-paste("ticks",seq(1:length(ee)),ee,sep="")
eee<-paste(ee,collapse=",")
maketemp2<-paste("temp2<-list(",eee,")",sep="")
eval(parse(text=maketemp2))

# top (right)
tt1<-.rotatedTernaryCoordinates(where=t.coords,alab="c2",blab="b2a",clab="a2a")
tt2<-.rotatedTernaryCoordinates(where=t.coords,alab="c1",blab="b1",clab="a1")
tt3<-.rotatedTernaryCoordinates(where=t.coords,alab="c2",blab="b2b",clab="a2b")
tt<-cbind(tt1[,"xx"],tt1[,"yy"],tt2[,"xx"],tt2[,"yy"],tt3[,"xx"],tt3[,"yy"])

ee<-apply(tt,FUN=function(z) {
	paste("=list(\"lines\",x=c(",z[1],",",z[3],",",z[5],"),y=c(",z[2],",",z[4],",",z[6],"),col=\"black\")",
	sep="")
},MARGIN=1)
ee<-paste("ticks",seq(1:length(ee)),ee,sep="")
eee<-paste(ee,collapse=",")
maketemp3<-paste("temp3<-list(",eee,")",sep="")
eval(parse(text=maketemp3))


# Bottom (right)
tb1<-.rotatedTernaryCoordinates(where=t.coords,alab="a2a",blab="c2",clab="b2a")
tb2<-.rotatedTernaryCoordinates(where=t.coords,alab="a1",blab="c1",clab="b1")
tb3<-.rotatedTernaryCoordinates(where=t.coords,alab="a2b",blab="c2",clab="b2b")
tb<-cbind(tb1[,"xx"],tb1[,"yy"],tb2[,"xx"],tb2[,"yy"],tb3[,"xx"],tb3[,"yy"])


ee<-apply(tb,FUN=function(z) {
	paste("=list(\"lines\",x=c(",z[1],",",z[3],",",z[5],"),y=c(",z[2],",",z[4],",",z[6],"),col=\"black\")",
	sep="")
},MARGIN=1)
ee<-paste("ticks",seq(1:length(ee)),ee,sep="")
eee<-paste(ee,collapse=",")
maketemp4<-paste("temp4<-list(",eee,")",sep="")
eval(parse(text=maketemp4))
}else{
temp2<-NULL
temp3<-NULL
temp4<-NULL}

## Ideal minerals
if(mins){
im<-c(1,0,0,0,0,	#q
	  3,1,0,1,0,	#fsp
	  2,2,1,0,0,	#an
	  2.5,1.5,0.5,0.5,0,	#an50
	  3,3,2,0,0,	#cz
	  3,2,2,0,1,	#Ep
	  1,1,0,0,0,	#sill
	  2,0,0,0,2,	#opx
	  2,0,1,0,1,	#cpx
	  1,0,0,0,2,	#olv
	  3,2,3,0,0,	#grs
	  3,2,0,0,3,	#Gt
	  5,4,0,1,2,	#NaCrd
	  5,4,0,0,2,	#Crd
	  6,2,0,2,6,	#bio
	  6,6,0,2,0,	#ms
	  7,2,2,0,4,	#MgHbl
	  7,1,2,1,5,	#Edn
	  6,3,2,1,4)	#Pgs
id.min<-matrix(im,byrow=T,ncol=5)
colnames(id.min)<-c("Si","Al","Ca","NaK","FM")
rownames(id.min)<-c("q","fsp","an","an50","cz","Ep","sill","opx","cpx","olv","grs-Gt","Gt","NaCrd","Crd","bio","ms","MgHbl","Edn","Pgs")

aa<-matrix(c(3,0,2,0,	#ms1
			 1,0,1,0,	#fsp
			 1,1,0,0,	#CaAl
			 1,0,1,3),	#bio
			byrow=F,nrow=4)
rownames(aa)<-c("Al","Ca","NaK","FM")
colnames(aa)<-c("ms1","fsp","CaAl","bio")

mins.p<-t(apply(id.min,MARGIN=1,FUN=function(z){
return(solve(aa)%*%z[2:5])
}))
colnames(mins.p)<-c("ms1","fsp","CaAl","bio")

neg.ph<-apply(mins.p,FUN=sum,MARGIN=1)<0
rownames(mins.p)[neg.ph]<-paste("-",rownames(mins.p)[neg.ph],sep="")

symb<-rep(16,nrow(mins.p))
names(symb)<-rownames(mins.p)
symb[neg.ph]<-1

imcoords<-.rotatedTernaryCoordinates(where=round(mins.p,2),alab="CaAl",blab="ms1",clab="fsp")

temp5<-list(
        pts=list("points",x=imcoords[,"xx"],y=imcoords[,"yy"],pch=symb,cex=0.8,col=plt.col[3]), # points
        labs=list("text",x=imcoords[,"xx"],y=imcoords[,"yy"]+0.08,text=rownames(mins.p),adj=0,cex=0.6,col=plt.col[3])
        )
}else{
temp5<-NULL
}

temp<-c(temp1,temp2,temp3,temp4,temp5)

sheet<<-list(demo=list(
				fun="plot",
				call=list(
					xlim=c(xmin-0.1,xmax+0.1),
					ylim=c(ymin-0.1,ymax+0.1),
					main=annotate("Projected from biotite, onto A-C-NK plane"),
					xlab="",
					ylab="",
					bg="transparent",
					fg="black",
					asp=1,
					axes=FALSE),
				template=temp))
}


.projbioGUI<-function(){

## Ticks
plt.tcks<-winDialog(type="yesno",message="Draw ticks?")
if(plt.tcks=="YES"){ticks<-TRUE}else{ticks<-FALSE}

## Ideal mins
plt.im<-winDialog(type="yesno",message="Plot ideal minerals?")
if(plt.im=="YES"){mins<-TRUE}else{mins<-FALSE}

## Ideal mins
add<-winDialog(type="yesno",message="Add new coords to WR?")
if(add=="YES"){addWR<-TRUE}else{addWR<-FALSE}

## Limits
def.lim<-"-2,1,-0.5,0.5"
plt.lim<-winDialogString("Plot limits\n(xmin,xmax,ymin,ymax)",def.lim)
lims<-as.numeric(unlist(strsplit(plt.lim,",")))
if(length(lims)!=4){lims<-c(-2,1,-0.5,0.5)}

ee<-paste("plotDiagram(\"projbioplot\",mins=",mins,
		  ",addWR=",addWR,
		  ",ticks=",ticks,
		  ",xmin=",lims[1],
		  ",xmax=",lims[2],
		  ",ymin=",lims[3],
		  ",ymax=",lims[4],")",sep="")

cat("GCDkit->",ee,"\n")
.save2hist(ee)
eval(parse(text=ee))  
}

### Plugin
if(.Platform$OS.type=="windows"){winMenuAddItem("Plugins","ACNK proj from bio",".projbioGUI()")}
