mat1<-c(1:11)
yields1<-c(0.1,0.15,0.25,0.5,1,1.5,2.5,3,3.4,4.1,4.2)
dat1<-data.frame(cbind(mat1,yields1))
piece_lin(mat1,yields1)
af1<-af
ticks<-c("1 mo","3 mo", "6 mo", "1 yr", "2 yr", "3 yr", "5 yr", "7 yr", "10 yr", "20 yr", "30 yr")

mat2<-c(1,2,3,3.5,4,5,6,6.5,7,8,9,9.5,10,10.5,11)
yields2<-c(0.17,0.36,0.65,0.95,1.35,1.85,2.8,3.4,3.9,4.2,4.9,5,5.2,5.7,5.9)
dat2<-data.frame(cbind(mat2,yields2))
piece_lin(mat2,yields2)
af2<-af

missing<-c(3.5,6.5,9.5,10.5)
missingy<-af1(missing)
dat_m<-data.frame(cbind(missing,missingy))

yplot<-ggplot(dat2,aes(x=mat2, y=yields2)) + ggtitle("Yield Curve")+
  geom_point(size=3, col='red')+
  stat_function(fun = af1, size=0.8, aes(colour="Government Treasury"))+
  geom_point(data=dat1,aes(x=mat1, y=yields1),size=3, col='blue')+
  geom_point(data=dat_m,aes(x=missing, y=missingy),size=3, col='black',shape=17)+
  ylim(0,6.5)+
  stat_function(fun = af2, size=0.8, aes(colour="Corporate bond"))+
  xlab("Time to maturity")+ylab("Yield [%]")+
  scale_x_continuous(breaks=c(1:11),labels=ticks)+
  geom_segment(aes(x=mat2,y=af1(mat2),xend=mat2,yend=af2(mat2)),linetype="dotted",size=1.3)+
  geom_segment(aes(x=3.5,y=af1(3.5),xend=3.5,yend=af2(3.5)),size=1.3)+
  geom_segment(aes(x=6.5,y=af1(6.5),xend=6.5,yend=af2(6.5)),size=1.3)+
  geom_segment(aes(x=9.5,y=af1(9.5),xend=9.5,yend=af2(9.5)),size=1.3)+
  geom_segment(aes(x=10.5,y=af1(10.5),xend=10.5,yend=af2(10.5)),size=1.3)+
  scale_colour_manual("", values = c("red","blue"))+ theme(legend.position="top")

print(yplot)

ggsave(filename="yields.pdf", plot=yplot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
