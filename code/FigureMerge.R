library(ggplot2)
library(png)
library(ggpubr)
library(patchwork)

img1<-readPNG("3A.png")
p1<-ggplot()+background_image(img1)+theme_void()

img2<-readPNG("3B.png")
p2<-ggplot()+background_image(img2)+theme_void()

img3<-readPNG("3C.png")
p3<-ggplot()+background_image(img3)+theme_void()

img4<-readPNG("3D.png")
p4<-ggplot()+background_image(img4)+theme_void()

img5<-readPNG("3E.png")
p5<-ggplot()+background_image(img5)+theme_void()

img6<-readPNG("3F.png")
p6<-ggplot()+background_image(img6)+theme_void()

p123456<-p1+p2+p3+p4+p5+p6+plot_layout(ncol=2)

ggsave(filename="all1.pdf",
       p123456,
       width=16,
       heigh=24,
       dpi = 1000)
