pdf(file='ABE_meninges_input.pdf',width=4.5,height=4.5);
gstable=read.table('ABE_meninges_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Ctnnb1","Traf7","Ctcf","Smarca4","Dis3","Arid1a","Stat6","Rhoa","Bcor","Rara")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='meninges1,meninges2,meninges3,meninges4,meninges5_vs_input_rep1,input_rep2,input_rep3 neg.'


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
targetmat=list(c(1030.9960755007621,550.6842472375045,793.0445214946421,1457.5638573780507,340.8315726382892,2494.6746409661478,2684.1335619527363,2094.568488772317),c(254.87143812219796,333.6340658873076,353.3976645925983,230.51251611849005,44.3873675994051,119.55309780987928,403.7505357931312,82.31442813019194),c(1239.8262215750792,457.7394913622964,504.61378097928855,165.08319931629853,26.949473185353096,808.9759618468498,255.1703386212589,258.0668557595207),c(284.46941158155,354.58692992853946,252.02686064448372,98.64727764022717,63.41052514200729,46.49287137050861,101.74513501986907,163.5164991234894),c(1315.4654870823122,552.2960060099069,527.0161685921315,70.46234117159085,50.72842011360583,55.79144564461033,888.2511787448886,1206.907493530517),c(198.9641549211997,255.73239188785573,164.65754895439602,57.37647781115254,68.16631452765783,106.93360415216979,33.91504500662302,74.52792817193055),c(600.1811284813049,333.0968129631734,225.7040551993932,24.15851697311686,55.48420949925637,33.209193836077574,69.44509215641857,94.5503566360313),c(88.79392037805607,192.33654684002596,40.88435739343847,10.066048738798692,28.534736313903277,10.626942027544825,557.9832404661073,5.561785684472429),c(559.0728320099827,326.1125249494295,277.78960639925316,36.23777545967529,20.608420671152366,24.57480343869741,83.98011144497129,1.1123571368944858))
targetgene="Ctnnb1"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(162.78885402643613,44.05473977900036,54.88584965146534,14.092468234318169,6.341052514200729,2.656735506886206,3.2300042863450495,0.0),c(2290.554279382076,609.782068892261,1105.5578286938019,435.8599103899834,286.932626267583,130.84422371414564,1227.4016288111188,123.47164219528793))
targetgene="Traf7"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(21.37631416508757,10.207805558548863,16.801790709632247,2.0132097477597384,0.0,61.10491665838274,0.8075010715862624,0.0),c(866.5628896154731,539.4019358306873,663.1106733401527,1770.61797315469,943.2315614873584,1078.6346157957996,641.9633519110786,1209.132207804306),c(62.484610636409826,23.10187573776848,15.121611638669023,2.0132097477597384,12.682105028401457,0.0,12.112516073793936,40.04485692820149),c(841.8979117326797,423.8925571418449,511.8945569534625,627.1148364271585,529.4778849357608,2341.912349320191,542.6407201059683,205.78607032547987),c(2820.0291379327064,1027.7648438686304,1138.0412907324244,1423.339291666135,334.49052012408845,2026.4250078774537,834.9561080201953,813.1330670698691),c(942.202155122706,503.4059899136992,441.32703597300707,6807.668762049555,462.89683353665316,130.84422371414564,937.5087441116507,235.81971302163097),c(1460.1666906613664,546.3862238444312,592.5431523596973,7.046234117159084,3.1705262571003643,23.246435685254305,74.29009858593614,229.14557020026407),c(715.2843586010072,470.09630861738185,271.62894980572133,295.9418329206815,488.2610435934561,197.92679526302237,155.04020574456237,204.6737131885854))
targetgene="Ctcf"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(608.4027877755693,349.21440068719795,335.47575450232387,64.42271192831163,1.5852631285501821,393.19685501915853,96.90012859035149,80.08971385640298),c(478.500570926191,323.4262603287587,444.12733442461246,193.26813578493488,245.71578492527823,101.62013313839739,1897.6275182277166,24.471857011678686),c(396.2839779835465,485.67664341727226,417.80452897952193,157.0303603252596,3531.966250409806,834.8791330389903,427.1680668691328,308.12292691977257),c(669.2430665531263,483.52763172073566,590.8629732887341,193.26813578493488,456.55578102245244,329.43520285388956,629.8508358372846,520.5831400666193),c(427.52628330175145,312.14394892194156,543.2578996114427,65.4293168021915,33.290525699553825,9.298574274101721,110.62764680731794,73.41557103503607),c(657.7327435411561,298.7126258185878,213.94280170265063,47.310429072353855,653.128408962675,2.656735506886206,50.872567509934534,15.5729999165228),c(1157.6096286324346,508.24126623090655,469.89008017938187,592.890270715243,2214.6125905846043,455.63013943098434,1075.5914273529015,716.3579961600489))
targetgene="Smarca4"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(55.907283200998265,87.03497370973241,56.005969032107494,49.32363882011359,1.5852631285501821,1.328367753443103,1.6150021431725248,243.60621297989238),c(608.4027877755693,394.88089923860076,299.071874631454,382.5098520743503,58.654735756356736,21.9180679318112,188.14774967959914,1271.4242074703973),c(2058.7034872838185,848.8596201319581,827.2081626042277,103.68030200962653,775.1936698610391,844.8418911898136,2879.5488212766118,898.7845666107445),c(1384.5274251541334,569.4880995821998,426.7654840246591,187.22850654165566,1020.9094547863173,449.6524845404904,258.40034290760394,367.07785517518033),c(113.45889826084942,297.10086704618533,111.45187837389392,84.55480940590901,118.89473464126367,1071.9927770285842,37.14504929296807,76.75264244571952))
targetgene="Dis3"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(1075.3930356897902,553.9077647823094,455.32852823103394,306.00788165948023,291.6884156532335,1490.4286193631617,623.3908272645946,8466.15016890393),c(394.6396461246936,484.60213756900396,462.0492445148868,504.3090418138145,204.4989435829735,455.63013943098434,692.0284183494268,380.42614081791413),c(513.0315399621018,246.5990921775752,212.2626226316874,1.0066048738798692,11.096841899851276,13.947861411152582,11.305015002207673,32.25835696994009))
targetgene="Arid1a"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(310.77872132319624,217.58743427433103,192.1004737801287,95.62746301858758,85.60420894170984,168.7027046872741,105.78264037780038,24.471857011678686),c(555.7841682922768,307.84592552886835,371.3195746828727,185.21529679389593,431.19157096564953,7.306022643937067,71.86759537117736,23.359499874784202),c(832.0319205795624,557.6685352512484,552.7789143469009,48.31703394623372,39.63157821375455,61.76910053510429,364.18298328540436,1074.5369942400732),c(925.7588365341771,284.2067968669657,296.83163587016975,1917.5822847411507,71.3368407847582,4991.34183356246,176.0352336058052,52.28078543404083))
targetgene="Stat6"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(439.0366063137217,416.9082691281009,357.31808242484584,190.2483211632953,85.60420894170984,27.231538945583612,505.49567081300023,102.3368565942927),c(1042.5063985127324,505.01774868610164,1069.153948822932,188.23511141553553,190.23157542602186,112.2470751659422,1093.3564509277992,692.9984962852647),c(483.43356650274967,207.37962871578216,152.3362357673324,210.38041864089266,415.33893968014775,130.84422371414564,566.8657522535563,338.1565696159237),c(116.7475619785552,218.66194012259933,127.13354970288401,3.0198146216396076,104.62736648431202,57.11981339805343,338.34294899464396,153.50528489143903),c(266.3817611341682,305.69691383233175,194.90077223173407,325.13337426319777,450.2147285082517,220.5090470715551,108.20514359255915,72.30321389814158),c(802.4339471202104,499.6452194447602,444.6873941149335,131.86523847826285,90.35999832736039,123.53820107020859,528.9132018890018,385.9879265023866),c(947.1351506992647,473.85707908632094,614.9455399725402,463.03824198473984,534.2336743214114,866.7599591216248,1656.9921988950105,463.85292608500055),c(703.774035589037,628.5859212369563,893.2952060621145,1152.5625805924503,664.2252508625263,248.40476989386028,1090.1264466414543,521.6954972035138))
targetgene="Rhoa"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(2617.776319293801,541.0136946030898,836.1691176493649,517.3949051742528,253.64210056802915,1959.3424363285772,448.16309473037563,401.5609264189094),c(2847.9827795332058,809.102903746031,963.86272704257,703.6168068420286,2948.589419103339,206.56118566040254,2336.9081011706435,265.8533557177821),c(723.5060178952716,372.3162764249664,224.58393581875106,269.77010619980496,331.31999386698806,74.38859419281377,72.67509644276362,131.25814215354933),c(1343.4191286828113,469.02180276911355,664.7908524111159,89.58783377530835,3533.551513538356,616.3626375975998,453.0081011598932,1955.523846660506),c(1178.9859427975223,457.7394913622964,754.4004028624879,317.0805352721588,550.0863056069132,2796.8783048744535,1546.3645520876926,310.3476411935615),c(1516.0739738623647,445.9199270313451,631.7473306821726,879.7726597710057,358.26946705234116,55.127261767888776,163.115216460425,289.2128555925663))
targetgene="Bcor"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(54.26295134214537,157.95235969544032,56.56602872242857,94.62085814470771,240.9599955396277,39.18684872657154,227.715302187326,484.9877116859958),c(1897.5589651162352,603.3350338026512,984.0248758941286,1027.7435762313464,204.4989435829735,340.72632875815594,299.58289755850336,2458.3092725368138),c(226.91779652169885,118.73289623364731,101.93086363843564,28.184936468636337,85.60420894170984,27.895722822305164,102.55263609145533,35.595428380623545),c(1040.8620666538795,384.13584075591774,493.9726468631881,208.36720889313293,391.559992751895,180.65801446826202,897.9411916039238,62.2919996660912))
targetgene="Rara"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetgenelist=c("Daxx","Ep300","Keap1","Syk","Pik3cg","Abl1","Tbl1xr1","Ncoa3","Arid5b","Raf1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='meninges1,meninges2,meninges3,meninges4,meninges5_vs_input_rep1,input_rep2,input_rep3 pos.'


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
targetmat=list(c(564.0058275865413,434.10036270039376,398.76249950860534,872.7264256538466,745.0736704185856,295.5618251410904,836.5711101633678,305.8982126459836),c(106.88157082543786,113.89761991643995,53.7657302708232,3573.4473022735356,1228.5789246263912,838.2000524225981,65.40758679848726,180.2018561769067))
targetgene="Daxx"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(891.2278674982665,402.4024401764789,405.4832157924583,161.05677982077907,692.7599871764296,731.9306321471498,406.98054007947627,152.39292775454456),c(154.56719473217169,232.63051615008726,288.99080020567465,794.2112454912168,169.6231547548695,212.5388405508965,276.972867554088,493.8865687811517),c(392.99531426584076,423.8925571418449,514.6948554050679,189.24171628941542,199.74315419732295,3030.0068456037184,231.7528075452573,1716.3670622281916),c(13.154654870823121,65.54485674436638,16.801790709632247,95.62746301858758,7.926315642750911,110.9187074124991,176.84273467739146,259.1792128964152),c(236.78378767481618,224.03446936394084,150.0959970060481,54.356663189512936,423.2652553228986,62.433284411825845,177.65023574897774,171.3029990817508),c(88.79392037805607,95.63102049587883,71.68764036109759,33.21796083803568,174.37894414052005,905.2826239714748,68.63759108483231,2444.9609868940797),c(1415.7697304723383,639.8682326437735,999.1464875327977,1347.843926125145,229.8631536397764,2138.672083043396,1255.664166316638,107.89864227876512),c(430.81494701945724,275.07349715668516,301.3121133927383,33.21796083803568,206.08420671152368,423.7493133483499,1174.9140591580117,40.04485692820149),c(154.56719473217169,284.7440497910999,81.20865509655587,33.21796083803568,359.85473018089135,52.470526261002576,369.83549078650816,261.40392717020416),c(2129.409757214493,774.1814636773112,1306.0591978287468,3330.855527668487,321.808415095687,221.8374148249982,1293.6167166811924,874.3127095990658),c(859.9855621800615,474.93158493458924,713.5160454690495,2403.7724388251277,144.25894469806659,493.4886204041128,934.2787398253056,453.8417118529502))
targetgene="Ep300"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(1093.480686137172,378.2260585904421,418.364588669843,128.84542385662326,114.13894525561311,124.86656882365169,205.10527218291065,116.79749937392101),c(1497.986323414983,422.8180512935766,633.4275097531357,333.1862132542367,149.01473408371712,1178.926381180754,323.00042863450494,109.01099941565961),c(192.38682748578813,155.26609507476957,90.16961014169307,1696.1292124875795,294.8589419103339,2452.1668728559684,1464.80694385748,33.370714106834576),c(90.43825223690895,164.3993947850501,60.48644655467609,433.8467006422236,38.04631508520437,261.02426355156973,378.71800257395705,474.97649745394546),c(110.17023454314364,185.352258826282,154.01641483829562,5524.247547852722,193.40210168312223,468.9138169654154,169.5752250331151,81.20207099329747),c(143.05687172020143,237.46579246729462,105.85128147068316,999.5586397627101,225.10736425412586,812.2968812304575,664.573381915494,503.89778301320206))
targetgene="Keap1"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(113.45889826084942,47.815510247939415,58.2462077933918,1173.7012829439275,177.5494703976204,158.07576265972926,230.94530647367105,2110.1414886888397))
targetgene="Syk"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(667.5987346942734,494.8099431275528,320.3541428636549,433.8467006422236,185.4757860403713,320.13662857978784,320.5779254197462,441.60578334711084),c(475.21190720848523,455.5904796657598,189.8602350188444,56.369872937272675,1.5852631285501821,27.231538945583612,233.36780968842982,35.595428380623545),c(47.68562390673382,147.2073012127573,52.64561089018105,282.85596956024324,280.59157375338225,296.226009017812,297.9678954153308,212.4602131468468),c(698.8410400124783,224.03446936394084,193.22059316077085,85.56141427978888,107.79789274141238,37.85848097312844,1245.1666523860165,24.471857011678686))
targetgene="Pik3cg"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(39.46396461246936,196.09731730896502,226.26411488971428,77.50857528874992,117.30947151271347,55.79144564461033,624.1983283361808,275.86456994983246),c(3.2886637177057803,140.7602661231475,31.363342657980198,139.9180774693018,63.41052514200729,1973.2902977397296,134.85267895490583,72.30321389814158),c(3.2886637177057803,116.04663161297655,81.20865509655587,107.706721505146,118.89473464126367,25.903171192140512,146.15769395711348,17.797714190311773))
targetgene="Abl1"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(529.4748585506306,341.1556068251857,344.9967692377822,934.1293229605186,2355.701009025571,5094.954518331022,279.39537076884676,1708.5805622699302))
targetgene="Tbl1xr1"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(287.7580752992558,219.73644597086763,120.97289310935219,430.826886020584,596.0589363348685,120.88146556332238,226.10030004415347,1082.3234941983346),c(108.52590268429074,189.6502822193552,117.61253496742574,535.5137929040905,416.9242028086979,461.6077943214783,550.7157308218309,2361.5342016269933))
targetgene="Ncoa3"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(21.37631416508757,43.517486854866206,36.40387987086987,37.244380333555156,139.50315531241603,820.2670877511162,152.61770252980358,126.80871360597138),c(167.72184960299478,180.51698250907464,85.68913261912446,145.95770671258103,39.63157821375455,43.17195198690085,2697.053579098116,563.9650684055043),c(82.21659294264451,117.12113746124486,100.81074425779349,261.717267208766,374.122098337843,175.3445434544896,96.90012859035149,56.730213981618775))
targetgene="Arid5b"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
targetmat=list(c(2267.5336333581354,873.5732546421291,835.0489982687227,637.1808851659572,386.80420336624445,83.02298459019394,2360.325632246645,1409.3564924453135),c(443.96960189028033,455.0532267416257,294.5913971088854,357.3447302273536,169.6231547548695,197.92679526302237,341.572953280989,233.59499874784203),c(70.70626993067428,74.1409035305128,29.1231038966959,25.16512184699673,76.09263017040874,9.962758150823273,452.20060008830694,1120.1436368527472),c(304.20139388778466,334.7085717355759,208.90226448976097,2759.1039593047212,2209.8568011989537,48.485423000673265,539.4107158196233,2135.725702837413))
targetgene="Raf1"
collabel=c("input_rep1","input_rep2","input_rep3","meninges1","meninges2","meninges3","meninges4","meninges5")

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
Sweave("ABE_meninges_input_summary.Rnw");
library(tools);

texi2dvi("ABE_meninges_input_summary.tex",pdf=TRUE);

