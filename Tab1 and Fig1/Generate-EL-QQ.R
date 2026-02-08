#==============================================#
#            Model Generate QQ-Plot            #
#==============================================#
rm(list = ls()) 
# EL-Chi QQ图
Model_EL_QQ <- function(EL, m, a=0.95, qtime=200){
  EL = EL[,1]
  m = m
  nsim = length(EL)
  qtime = qtime                     # 分位数个数
  
  
  as = (1:qtime-0.5)/qtime
  qas <- qchisq(as,m)               # 卡方分布的分位数点
  
  x <- seq(0,max(EL)+5,2) 
  c1 = sort(EL)[ceiling(as*nsim)]   # mel的分位数点

  plot(x,x,xaxs = 'i', yaxs = 'i',xaxt="n",yaxt="n",
       xlim =c(0, max(c1)),
       ylim =c(0, max(c1)),
       ann = F, type = 'l',lwd=1.5)
  axis(1,seq(0,max(c1)+5,2),mgp=c(4,0,  0),tcl=0.5,font=2,lwd=1.5,las=1) 
  axis(2,seq(0,max(c1)+5,2),mgp=c(4,0.2,0),tcl=0.5,font=2,lwd=1.5,las=1) 
  axis(3,seq(0,max(c1)+5,2),tck=0.01,labels=F,tcl=0,lwd=1.5)
  axis(4,seq(0,max(c1)+5,2),         labels=F,tcl=0,lwd=1.5)
  # 根据modelname显示对应的标题
  if(modelname == 'Model1'){
    title(main = paste0('Model5.1',' n = ',n,' m = ',m))
  } else if(modelname == 'Model2'){
    title(main = paste0('Model5.2',' n = ',n,' m = ',m))
  } else {
    title(main = paste0(modelname,' n = ',n,' m = ',m))
  }
  title(xlab='Chi-Square quantile',line=1,col.lab = 1, font.lab =2, cex.lab = 1)
  title(ylab ='EL quantile',       line=2,col.lab = 1, font.lab =2, cex.lab = 1)
  points(qchisq(as,m),c1,        pch=20,cex=0.7,col='blue1')
  points(qchisq(a,m),qchisq(a,m),pch=5 ,cex=1,  col='blue1')
  
}

# EEL-Chi QQ图
Model_EEL_QQ <- function(EEL, m, a=0.95, qtime=200){
  EEL = EEL[,1]
  m = m
  nsim = length(EEL)
  qtime = qtime                     # 分位数个数
  
  as = (1:qtime-0.5)/qtime
  qas <- qchisq(as,m)               # 卡方分布的分位数点
  
  x <- seq(0,max(EEL)+5,2) 
  c2 = sort(EEL)[ceiling(as*nsim)]
  plot(x,x,xaxs = 'i', yaxs = 'i',xaxt="n",yaxt="n",
       xlim =c(0, max(c2)),
       ylim =c(0, max(c2)),
       ann = F, type = 'l',lwd=1.5)
  axis(1,seq(0,max(c2)+5,2),mgp=c(4,0,  0),tcl=0.5,font=2,lwd=1.5,las=1) 
  axis(2,seq(0,max(c2)+5,2),mgp=c(4,0.2,0),tcl=0.5,font=2,lwd=1.5,las=1) 
  axis(3,seq(0,max(c2)+5,2),tck=0.01,labels=F,tcl=0,lwd=1.5)
  axis(4,seq(0,max(c2)+5,2),         labels=F,tcl=0,lwd=1.5)
  # 根据modelname显示对应的标题
  if(modelname == 'Model1'){
    title(main = paste0('Model5.1',' n = ',n,' m = ',m))
  } else if(modelname == 'Model2'){
    title(main = paste0('Model5.2',' n = ',n,' m = ',m))
  } else {
    title(main = paste0(modelname,' n = ',n,' m = ',m))
  }
  title(xlab='Chi-Square quantile',line=1,col.lab = 1, font.lab =2, cex.lab = 1)
  title(ylab ='EEL quantile',       line=2,col.lab = 1, font.lab =2, cex.lab = 1)
  points(qchisq(as,m),c2,         pch=20,cex=0.7,col='aquamarine4')
  points(qchisq(a, m),qchisq(a,m),pch=5 ,cex=1,  col='aquamarine4')
}


#===============================================
#                   赋值画图            
#=============================================== 
# 选择模型
# modelname = c('Model2')
# # 图EL
# par(mfrow = c(2, 3))
# ns = c(400,2000,4000)
# ms = c(1,5)
# for(m in ms){
#   for(n in ns){
#     EL_path <- file.path("data", paste0(modelname, '-', 'EL', '-', n, '-', m,'.csv'))
#     EL = read.csv( EL_path)
#     Model_EL_QQ(EL=EL,m=m,a=0.95,qtime=100)
#   }
# }
# # 图EEL
# par(mfrow = c(2, 3))
# ns = c(400,2000,4000)
# ms = c(1,5)
# for(m in ms){
#   for(n in ns){
#     EEL_path <- file.path("data", paste0(modelname, '-', 'EEL', '-', n, '-', m,'.csv'))
#     EEL = read.csv( EEL_path)
#     Model_EEL_QQ(EEL=EEL,m=m,a=0.95,qtime=100)
#   }
# }

ns = c(2000,4000,8000)
ms = c(1,5,10)

# 保存图片
for(modelname in c('Model1','Model2')){
  # 保存图EL
  for(m in ms){
    # setEPS()
    # postscript(paste0(modelname,'-QQ-EL-',m,'.eps'), width=8.9, height=3)
    par(mfrow = c(1, 3),mar = c(4, 4, 2, 1),oma = c(0, 0, 0, 0))
    for(n in ns){
      EL_path <- file.path("data", paste0(modelname, '-', 'EL', '-', n, '-', m,'.csv'))
      EL = read.csv(EL_path)
      Model_EL_QQ(EL=EL,m=m,a=0.95,qtime=100)
    }
    # dev.off()
  }
  
  # 保存图EEL
  for(m in ms){
    # setEPS()
    # postscript(paste0(modelname,'-QQ-EEL-',m,'.eps'), width=8.9, height=3)
    par(mfrow = c(1, 3),mar = c(4, 4, 2, 1),oma = c(0, 0, 0, 0))
    for(n in ns){
      EEL_path <- file.path("data", paste0(modelname, '-', 'EEL', '-', n, '-', m,'.csv'))
      EEL = read.csv(EEL_path)
      Model_EEL_QQ(EEL=EEL,m=m,a=0.95,qtime=100)
    }
    # dev.off()
  }
}
