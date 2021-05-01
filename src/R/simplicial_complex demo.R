library(ggplot2)
library(ggforce)
library(gridExtra)

points <- data.frame(x=c(4,3,5,4,4),y=c(4,3,3,2,.3))
path1 <- data.frame(x=c(4,3,4,5,4),y=c(4,3,2,3,4))
path2 <- data.frame(x=c(4,4),y=c(2,.3))

e0 <- ggplot(data=points) + geom_point(aes(x=x,y=y), size=2) + xlim(1,7) + ylim(-1,6) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ggtitle('\u03f5 = 0')

e.5 <- ggplot(data=points) + geom_circle(aes(x0=x,y0=y,r=.5), fill='blue', alpha=.5) + geom_point(aes(x=x,y=y), size=2) + xlim(1,7) + ylim(-1,7) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ggtitle('\u03f5 = 0.5')

e.75 <- ggplot(data=points) + geom_circle(aes(x0=x,y0=y,r=.75), fill='blue', alpha=.5) + geom_point(aes(x=x,y=y), size=2) + xlim(1,7) + ylim(-1,7) + 
  geom_path(data=path1, aes(x, y), size=2)  +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ggtitle('\u03f5 = 0.75')
e.85 <- ggplot(data=points) + 
  geom_circle(aes(x0=x,y0=y,r=.85), fill='blue', alpha=.5) + 
  geom_point(aes(x=x,y=y)) + xlim(1,7) + ylim(-1,7) +
  geom_path(data=path1, aes(x, y), size=2)  +
  geom_path(data=path2, aes(x, y), size=2)  +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ggtitle('\u03f5 = 0.85')

e1 <- ggplot(data=points) + 
  geom_circle(aes(x0=x,y0=y,r=1), fill='blue', alpha=.5) + 
  geom_point(aes(x=x,y=y)) + xlim(1,7) + ylim(-1,7) +
  geom_path(data=path1, aes(x, y), size=2)  +
  geom_path(data=path2, aes(x, y), size=2)  +
  geom_polygon(data=path1, aes(x, y), size=2, fill ='black', alpha=.8) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ggtitle('\u03f5 = 1')

pddf <- data.frame(Birth=c(.75,0,0,0),Death=c(1,.75,.8,1), Feature=as.factor(c('Hole', 'Connected Component', 'Connected Component','Connected Component')))
pd <- ggplot(data=pddf) + geom_point(aes(x=Birth,y=Death, shape=Feature, color=Feature), size=3) + 
  xlim(c(0,1.2)) + ylim(c(0,1.2)) + 
  scale_color_manual(values=c('black','red')) + 
  geom_abline(intercept = 0, slope=1) + theme(legend.position="bottom", legend.title = element_blank()) + ggtitle('Persistence Diagram')
pd


gridExtra::grid.arrange(e0,e.5,e.75,e.85,e1,pd, ncol=3)



