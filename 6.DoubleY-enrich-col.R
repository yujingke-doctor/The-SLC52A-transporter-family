library(tidyverse)

#读入数据
data1 <- read.csv('XXXXXXXXXX.csv',header = T, check.names=F)

#初步绘制柱状图
p1 <- ggplot(data1)+
  # geom_col(aes(x = reorder(GO_terms,-`-Log10(P value)`),
  geom_col(aes(x = reorder(GO_terms,-Log10P),
               y = Gene_count),
           color='black',   #柱子边框颜色
           width=0.6,       #柱子宽度
           fill='#A59ACA')+ #柱子填充色 文献的颜色：#d5a478
  labs(x=NULL,y='Gene count',
       title = 'Enriched GO Terms in vitamin B2 deficiency group (expanded gene families)')+
  theme_test(base_size = 15)+ #主题基本大小
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text = element_text(color = 'black',face = 'bold',size=13),
        plot.margin = margin(1,0.5,0.5,2.5,'cm'),
        panel.border = element_rect(size = 1),
        axis.title = element_text(face = 'bold',size = 18),
        plot.title = element_text(face = 'bold',
                                  size=15,hjust = 0.5))
#双y轴绘制
p2 <- p1+
  scale_y_continuous(expand = c(0,0),limits = c(0,75),
                     sec.axis = sec_axis(~./2.5,
                                         # name = '-Log10(P value)',
                                         breaks = seq(0,30,5)))+
  geom_line(aes(x= reorder(GO_terms,-Log10P),
                y=Log10P*2.5,
                group=1),
            linetype=3,cex=0.6)+
  geom_point(aes(x= reorder(GO_terms,-Log10P),
                 y=Log10P*2.5),
             color = "black", fill = '#589c47', shape = 21, size=3.5)

p3 <- p2+
  geom_text(aes(x= reorder(GO_terms,-Log10P),
                     y=Gene_count,
                     label=Gene_count),
                 vjust=-0.5,size=3.5,fontface='bold')

#添加图例
#先定义一个用来画框的数据框：
df <- data.frame(a=c(2.5,2.5,18.5,18.5), 
                 b=c(75,65,65,75))
p4 <- p3+
  geom_line(data = df,aes(a,b),cex=0.5)+
  geom_rect(aes(xmin=4,xmax=5.4,ymin=68,ymax=72),
            fill='#A59ACA',color='black')+ #文献中颜色 #d5a478
  annotate(geom='text',x=7.5,y=70,label='Gene Count',
           fontface='bold',size=5)+
  annotate('segment',x=10.9,xend = 12.6,y=70,yend = 70,
           linetype=3,cex=0.5)+
  annotate(geom='point', x=11.7,y=70,
           color = "black", fill = '#589c47', shape = 21, size = 3)+
  annotate('text',x=15.4,y=70,label='-Log10(P value)',
           fontface='bold',size=5)
p4
ggsave("FigS4I-双轴富集分析-enrih-bar+自加次轴坐标.png", width = 9, height = 7, dpi = 600)
ggsave("FigS4I-双轴富集分析-enrih-bar+自加次轴坐标.pdf", width = 9, height = 7)
ggsave("FigS4I-双轴富集分析-enrih-bar+自加次轴坐标.tiff", width = 9, height = 7, dpi = 600)


#####加-Log10(P value)名字#####
p2 <- p1+
  scale_y_continuous(expand = c(0,0),limits = c(0,75),
                     sec.axis = sec_axis(~./2.5,
                                         name = '-Log10(P value)',
                                         breaks = seq(0,30,5)))+
  geom_line(aes(x= reorder(GO_terms,-Log10P),
                y=Log10P*2.5,
                group=1),
            linetype=3,cex=0.6)+
  geom_point(aes(x= reorder(GO_terms,-Log10P),
                 y=Log10P*2.5),
             color = "black", fill = '#589c47', shape = 21, size=3.5)

p3 <- p2+
  geom_text(aes(x= reorder(GO_terms,-Log10P),
                y=Gene_count,
                label=Gene_count),
            vjust=-0.5,size=3.5,fontface='bold')

#添加图例
#先定义一个用来画框的数据框：
df <- data.frame(a=c(2.5,2.5,18.5,18.5), 
                 b=c(75,65,65,75))
p4 <- p3+
  geom_line(data = df,aes(a,b),cex=0.5)+
  geom_rect(aes(xmin=4,xmax=5.4,ymin=68,ymax=72),
            fill='#A59ACA',color='black')+ #文献中颜色 #d5a478
  annotate(geom='text',x=7.5,y=70,label='Gene Count',
           fontface='bold',size=5)+
  annotate('segment',x=10.9,xend = 12.6,y=70,yend = 70,
           linetype=3,cex=0.5)+
  annotate(geom='point', x=11.7,y=70,
           color = "black", fill = '#589c47', shape = 21, size = 3)+
  annotate('text',x=15.4,y=70,label='-Log10(P value)',
           fontface='bold',size=5)
p4
ggsave("FigS4I-双轴富集分析-enrih-bar.png", width = 9, height = 7, dpi = 600)
ggsave("FigS4I-双轴富集分析-enrih-bar.pdf", width = 9, height = 7)
ggsave("FigS4I-双轴富集分析-enrih-bar.tiff", width = 9, height = 7, dpi = 600)
