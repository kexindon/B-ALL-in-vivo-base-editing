pdf(file='ABE_d15_input.pdf',width=4.5,height=4.5);
gstable=read.table('ABE_d15_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Ctnnb1","Vhl","Rad21","Smarca4","Inpp4a","Dis3","Cbfb","Foxa1","Socs1","Recql")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='d15_rep1,d15_rep2,d15_rep3_vs_input_rep1,input_rep2,input_rep3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1546.4278087155587,839.4168529945387,1015.79797644797,924.3424027264667,897.3898410881176,1315.4850026310758),c(382.2907661099069,508.56377142400834,452.66138639736516,432.91036231588583,520.0782033578863,608.0498792649291),c(1859.6595977217405,697.7396670744848,646.3516785166814,417.1545259057756,431.3589804321293,504.7518302778369),c(426.6858228194444,540.5025590013615,322.81715353219386,207.07670710430588,143.78632681070977,242.5939029242317),c(1973.1136315350032,841.8736828081812,675.046536608432,730.7706982593983,693.4376044771818,806.8203674673641),c(298.43343676966924,389.8169970979516,210.90720697436663,92.28418468778848,110.13420776990534,55.561829379420814),c(900.2330943878452,507.7448281527941,289.1006952743869,154.55725240393843,189.6755800481703,215.20426872310878),c(838.5732934023763,497.09856562700975,355.816240337707,51.76917677607647,145.8258491768191,292.67780546342794))
targetgene="Ctnnb1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(83.85732934023764,271.0702227718949,106.17097493947709,141.80252769099206,84.64017819353836,106.42829289579197),c(416.8202546617694,203.09793126111765,114.77943236700226,171.06336673833962,151.94441627514718,136.16561002843974),c(224.4416755871066,382.446507657024,238.88469361382343,366.8859049782811,302.8690713672397,420.23524474294334),c(3768.647036231856,1040.876897713228,1489.2631349618543,518.4420456850556,665.9040525347054,866.2950017326597),c(1087.6788893836706,531.4941830180055,789.8259689754343,573.9626120825869,458.8925323746056,525.0984156843854),c(1292.389428655427,462.70294823601404,465.5740725386529,184.56836937557696,269.21695232643526,287.19987862320335),c(342.8284934792068,108.1005118002723,118.36628962847108,266.34866312329194,201.91271424482647,239.46365901553193),c(295.96704473025045,293.1816910946779,104.73623203488955,146.30419523673785,179.47796821762353,199.55304917960996),c(273.7695163754817,140.03929937762547,44.47703004221337,90.78362883920656,7.138328281382754,30.5198781098227),c(3117.519537825305,895.1049954371032,1547.3702225976492,841.8118310544609,1283.879329465841,1132.3657339721397),c(3267.969452229849,1052.3421035102265,1488.5457635095604,936.3468495151222,984.0695416477654,1124.5401242003902),c(207.17693131117534,281.7164852976793,150.64800498169046,260.34643972896424,364.0547423505204,266.8532932166549))
targetgene="Vhl"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2343.0724374478164,981.0940389145925,1077.4919213452335,742.7751450480537,1060.5516303768663,744.2154892933689),c(1531.6294564790462,610.9316803257814,530.8548746973854,307.61394895929493,296.75050426891164,461.7109765332152))
targetgene="Rad21"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(912.5650545849389,532.3131262892197,429.7054999239647,238.58837992452632,156.0234610073659,402.23634226791967),c(717.7200834708574,493.0038492709388,568.8755616689549,516.9414898364737,453.7937264593322,387.3676837015958),c(594.4004814999197,740.3247171776225,535.1591034111481,785.5409867326385,748.5047083621345,795.864513786915),c(1003.8215600434329,737.0489440927656,756.8268821699211,997.1193613826903,683.239992646635,729.346830727045),c(641.2619302488761,475.80604057544093,695.8503087249512,469.673980606143,359.97569761830175,350.58731777437356),c(986.5568157675016,455.3324587950863,274.0358947762179,240.83921369739923,168.26059520402205,226.16012240355795),c(1736.3399957508027,774.7203345686181,601.874648474468,460.6706455146515,533.3350987375971,875.6857334587589))
targetgene="Smarca4"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3531.8734004476555,965.534116761523,1322.1155865774072,931.8451819693764,933.0814824950314,1160.5379291504375),c(3127.38510598298,972.0856629312365,1567.4566232618745,1302.4824765691124,1432.7644621918241,1380.4375637365958),c(2140.8282902154783,869.7177540294635,1034.449634207608,955.1037976223963,1271.6421952691849,1202.0136609407093),c(3667.524962615687,1247.2506020592023,1684.3881699857582,1458.5402848216327,1381.7764030390902,1511.907807901986))
targetgene="Inpp4a"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(83.85732934023764,132.6688099366978,71.73714522937641,18.75694810727408,25.494029576366977,24.25939029242317),c(912.5650545849389,601.9233043424254,383.07635552487,554.4553860510218,338.56071277415344,299.7208542580024),c(3087.92283335228,1293.930368518411,1059.5576350378897,1139.6721669979731,916.7653035661565,933.5952457697047),c(2076.702097190591,868.0798674870351,546.6370466478483,902.5843429220289,883.1131845253522,1083.8469533872933),c(170.18105071989402,452.8756289814438,142.75691900645904,386.3931310098461,263.0983852281072,460.1458545788653))
targetgene="Dis3"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2643.9722662569043,917.2164637598862,1190.8366108076484,1108.9107721020437,1055.4528244615929,945.3336604273287),c(1220.8640595122831,578.9928927484282,785.5217402616717,567.9603886882592,700.5759327585646,577.5300011551064),c(419.28664670118815,415.2042385055913,419.662299591852,317.3675619750775,180.4977294006782,203.46585406548465))
targetgene="Cbfb"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1978.0464156138407,717.3943055836253,781.9348830002028,535.6984379437478,566.9872177784016,811.5157333304138),c(2185.223346925016,937.6900455402408,1005.0374046635635,961.106021016724,974.8916910002732,949.2464653132034),c(1213.4648833940269,573.2602898499289,388.09795569092637,672.9992980889941,545.5722329342533,576.7474401779315),c(720.1864755102762,633.0431486485643,356.53361179000075,441.16341948308644,481.32727840180854,609.6150012192791),c(2079.1684892300095,774.7203345686181,692.2634514634824,778.7884854140199,792.3544392334857,882.7287822533334),c(453.81613525305073,353.7834931645275,154.23486224315928,407.4009128899931,370.1733094488485,448.4074399212412),c(180.04661887756905,356.24032297817007,104.73623203488955,286.60616707914795,349.77808578775495,293.46036644060285),c(3642.8610422214997,871.3556405718919,847.2156851589353,839.5609972815879,813.769424077634,981.331465377376),c(1450.2385191782273,959.8015138630237,413.9233279735019,489.18120663770804,315.10620556389586,430.40853744621757),c(1827.5965012092968,836.1410799096819,652.0906501350315,431.4098064673039,455.83324882544156,418.67012278859346),c(3519.541440250562,1097.3839834270066,1170.750210143423,874.0737817989723,1110.5199283465456,1389.8282954626952),c(2007.6431200868658,877.9071867416053,779.0653971910278,682.7529111047766,745.4454248129704,949.2464653132034))
targetgene="Foxa1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(631.396362091201,504.4690550679374,190.10343485784747,180.81697975412214,113.19349131906938,219.8996345861584))
targetgene="Socs1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1642.6170982528902,802.5644057899004,971.3209464057566,627.2323447072453,800.5125286979231,874.9031724815841),c(483.4128397260758,484.814416558797,313.4913246523749,558.2067756724767,444.6158758118401,431.1910984233925),c(2774.6910443460984,845.149455893038,1251.0958128003247,728.5198644865253,990.1881087460935,938.2906116327542),c(1770.8694843026653,773.901391297404,1298.442328651713,1040.635480991566,1165.5870322314981,1202.7962219178844))
targetgene="Recql"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Erbb3","Atrx","Runx1","Ret","Rara","Gli1","Gata3","Tsc1","Keap1","Abl1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='d15_rep1,d15_rep2,d15_rep3_vs_input_rep1,input_rep2,input_rep3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(683.1905949189949,625.6726592076367,474.8999014184718,584.4665030226604,1107.4606447973815,1094.8028070677424),c(974.2248555704078,723.1269084821246,524.3985316267416,616.7284537671718,781.1370662198842,857.6868309837354),c(1627.8187460163776,623.2158293939941,713.0672235800015,579.2145575526237,895.3503187220083,948.4639043360285),c(414.3538626223507,259.60501697489633,139.17006174499022,519.9426015336376,475.20871130348047,558.7485377029079),c(78.92454526140013,294.81957763710625,105.45360348718332,608.4753965999712,405.8649508557623,610.397562196454),c(330.49653328211303,442.22936645565943,216.64617859271675,583.7162250983695,660.8052466194321,616.6580500138535),c(175.11383479873152,515.1153175937218,190.10343485784747,486.1800949405442,591.4614861717139,503.969269300662),c(2330.7404772507225,972.0856629312365,642.0474498029189,776.537651641147,768.8999320232281,885.8590261620332),c(1543.96141667614,881.182959826462,525.1159030790353,415.65397005719365,515.9991586256676,676.9152452563239),c(192.3785790746628,459.42717515115726,326.4040107936627,468.1734247575611,446.65539817794945,367.02109829504735),c(1622.88596193754,662.5251064122749,699.43716598642,596.4709498113158,579.2243519750577,749.6934161335935),c(574.6693451845697,517.5721474073644,318.51292481843126,476.4264819247617,735.2478129824236,864.7298797783099),c(693.0561630766699,443.8672529980878,402.44538473680166,650.4909603602651,548.6315164834174,399.1060983592199),c(355.16045367630056,452.8756289814438,418.22755668726444,405.1500791171202,592.4812473547686,835.775123622837),c(919.9642307031952,456.97034533751474,220.95040730647935,491.43204041058095,506.8213079781755,415.53987887989365),c(947.0945431368016,773.0824480261898,349.3598972670631,934.8462936665403,681.2004702805257,981.331465377376),c(421.75303874060694,308.7416132477474,114.77943236700226,711.2634722278332,688.3387985619084,671.4373184160994))
targetgene="Erbb3"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1048.2166167529704,736.2300008215515,628.4173922093373,961.106021016724,856.5993937659305,849.8612212119859),c(599.3332655787572,741.1436604488366,541.6154464817919,889.8296182090825,772.9789767554468,806.0378064901893),c(619.0644018941073,495.4606790845814,796.999683498372,855.3168336916982,574.1255460597844,622.918537831253),c(1001.3551680040141,665.8008794971316,599.7225341175867,633.234568101573,586.3626802564405,665.1768305986999),c(1371.3139739168273,701.0154401593416,741.0447102194584,858.317945388862,972.8521686341638,849.078660234811),c(483.4128397260758,257.96713043246797,121.23577543764613,365.38534912969914,434.4182639812933,658.9163427813003),c(98.65568157675015,191.6327254641191,100.43200332112697,440.4131415587955,524.1572480901051,379.54207392984637),c(545.0726407115446,788.6423701792593,973.4730607626378,1076.6488213575324,676.1016643652523,939.0731726099292),c(601.799657618176,469.2544944057275,440.46607170837115,447.16564287741414,466.03086065598836,394.4107324961703))
targetgene="Atrx"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(411.8874705829319,176.89174658226375,149.21326207710294,128.29752505375473,128.48990906488956,140.07841491431444),c(209.64332335059407,270.25127950068077,175.03863435967844,138.80141599382821,129.50967024794426,129.9051222110402),c(357.6268457157193,429.1262741162325,284.0790951083306,329.3720087637329,481.32727840180854,444.4946350353665),c(17.264744275931278,149.8666186321957,75.32400249084523,42.015563760293944,38.7509249560778,70.43048794574469),c(342.8284934792068,506.10694161036577,201.58137809454772,536.4487158680388,750.5442307282439,698.8269526172223),c(397.0891183464194,362.7918691478836,260.4058371826364,291.10783462489377,233.52531091952153,256.68000051338066),c(78.92454526140013,162.15076770040844,73.17188813396393,269.3497748204558,287.57265362141953,316.9371957558511))
targetgene="Runx1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(771.98070833807,603.5611908848537,519.3769314606852,504.9370430478183,637.3507394091745,715.2607331378961),c(246.6392039418754,312.01738633260413,180.77760597802856,510.939266442146,688.3387985619084,519.6204888441608),c(789.2454526140012,361.15398260545516,390.25007004780764,435.1611960887587,471.12966657126174,597.0940255844799),c(271.3031243360629,232.57988902482828,190.10343485784747,922.0915689535939,994.2671534783121,1125.322685177565))
targetgene="Ret"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(81.39093730081888,240.7693217369701,72.45451668167017,185.31864729986793,247.801967482287,318.502317710201),c(2846.216413489242,919.6732935735288,1260.4216416801435,927.3435144236306,898.4096022711723,1083.0643924101184),c(340.36210143978803,180.98646293833468,130.56160431746505,382.64174138839127,352.837369336919,264.50561028513005),c(1561.2261609520713,585.5444389181416,632.7216209231,622.7306771614996,554.7500835817455,775.5179283803666))
targetgene="Rara"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(369.9588059128131,452.8756289814438,293.4049239881495,266.34866312329194,374.25235418106723,223.81243947203313),c(520.408720317357,259.60501697489633,152.80011933857176,240.83921369739923,262.07862404505255,337.2837811623996),c(78.92454526140013,134.3066964791262,120.51840398535236,371.38757252402684,337.5409515910988,452.3202448071159))
targetgene="Gli1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(135.65156216803146,323.4825921296027,149.9306335293967,383.39201931268224,419.12184623547313,544.662440113759),c(1245.5279799064708,615.8453399530664,462.70458672947785,763.0326490039097,653.6669183380493,719.1735380237708))
targetgene="Gata3"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(315.6981810456005,366.88658550395445,123.38788979452742,282.85477745769316,253.9205345806151,309.89414696127665),c(328.0301412426943,118.7467743260567,65.28080215873253,338.37534385522446,296.75050426891164,244.1590248785816))
targetgene="Tsc1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1640.1507062134715,576.5360629347856,535.8764748634418,407.4009128899931,360.9954588013564,668.3070745073996),c(2246.8831479104847,644.5083544455629,811.3471125442472,593.4698381141519,650.6076347888853,700.3920745715722),c(288.5678686119942,236.6746053808992,115.49680381929602,541.7006613380755,314.08644438084116,438.23414721796695),c(135.65156216803146,250.59664099154034,77.47611684772652,396.1467440256286,283.4936088892008,462.4935375103901),c(165.24826664105652,282.5354285688935,197.27714938078512,225.08337728728898,208.03128134315455,187.81463452198585),c(214.5761074294316,361.97292587666936,135.5832044835214,399.89813364708345,439.5170698965667,719.9560990009458))
targetgene="Keap1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(59.1934089460501,298.9142939931772,289.8180667266807,390.8947985555919,257.9995793128338,139.2958539371395),c(4.932784078837508,176.89174658226375,104.01886058259579,236.33754615165344,201.91271424482647,76.69097576314422))
targetgene="Abl1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("ABE_d15_input_summary.Rnw");
library(tools);

texi2dvi("ABE_d15_input_summary.tex",pdf=TRUE);

