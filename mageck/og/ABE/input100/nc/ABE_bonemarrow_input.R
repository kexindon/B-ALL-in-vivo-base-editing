pdf(file='ABE_bonemarrow_input.pdf',width=4.5,height=4.5);
gstable=read.table('ABE_bonemarrow_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Smarca4","Ctnnb1","Vhl","Ezh2","Cbfb","Rhoa","Crebbp","Arid1a","Cd79b","Cbl")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='bonemarrow1,bonemarrow2,bonemarrow3,bonemarrow4,bonemarrow5_vs_input_rep1,input_rep2,input_rep3 neg.'


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
targetmat=list(c(1204.6790851839587,542.1870809927696,471.32110991754274,154.28670042687736,90.96003677604577,408.85570918352147,74.22508761394782,317.7401516799366),c(947.4638210500865,502.14865039638045,623.9693491896685,136.67128974821816,445.174783163187,163.96450987910586,3682.7062700766423,49.57647756708231),c(784.6693500792812,754.0571095653288,586.9875592629163,1344.23806316803,281.542970973475,294.85463364525907,550.0269312931006,2038.2692708829977),c(1325.1469937023546,750.7205736822963,830.1231568664567,1146.8239779761593,758.9628465387523,1561.533250737064,466.2858068055697,594.9177308049877),c(846.5312490481872,484.6318370104602,763.2411963606285,383.28703924944733,202.1334150578795,351.85517141439027,451.0601478078368,573.5097064010204),c(1302.3557677664419,463.7784877415075,300.5753989791341,170.07982724222703,120.7987789988756,53.48198605498732,121.8052719818631,221.9674109253458),c(2292.146151268938,789.0907363371692,660.1642925222343,397.2578822014874,1221.4633509926148,394.0777919841171,1436.9215679110412,479.9904418994787))
targetgene="Smarca4"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(2041.4426659738976,854.9873200270597,1114.1747773676805,808.4866073550149,776.2885678294277,2605.1356891521455,943.9908578594391,3217.964089354252),c(504.6628600094962,517.9971958407845,496.50020093150164,764.751794635585,480.78876581624195,876.1193768218318,323.54525370182387,228.72783968449338),c(2454.940622239743,710.6821430859072,708.9487813617797,569.159993307024,408.11699040257577,139.33464788009854,373.0286454444557,194.92569588875543),c(563.2688695589861,550.5284207003506,354.08096738379675,68.03193089689081,127.05528946495284,117.51962725240635,650.8969221530809,75.49145447714807),c(2604.7115355328838,857.489721939334,740.4226451292283,2167.910369383958,1254.670983466409,976.0462455035185,864.0561481213414,1571.7996865018142),c(393.9626197493487,397.0477700808589,231.3328986907472,185.26552610314016,219.94040638440697,216.7427855912644,243.6105439637262,104.7866457667876),c(1188.3996380868782,517.1630618700264,317.0991774570446,185.87295405757666,232.9346973524135,318.0770749586088,683.2514475232633,267.0369359863297),c(1107.0024026014755,506.31932025017096,390.2759107163626,54.66851589928726,148.71244107829705,27.44470337032244,312.1260094535242,109.293598272886))
targetgene="Ctnnb1"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(110.70024026014757,276.09834432093345,116.45329593955981,108.12217588970147,56.308594194695004,163.96450987910586,34.257732744898995,329.0075329451826),c(550.2453118813218,206.86522474801055,125.8954550697944,58.31308362590641,184.8076937672041,49.259723998014636,91.35395398639733,69.85776384452507),c(296.28593716686555,389.540564344036,262.0199158640096,564.3005696715318,186.2515038747604,180.8535581069966,348.2869495731398,322.247104186035),c(4974.999032867808,1060.1842768335541,1633.4935295305822,616.5393737530729,654.5272487588479,467.26366763831027,1124.795558457517,1276.5942973523695),c(1435.8472339625023,541.3529470220114,866.3181001990226,458.60810559957645,221.38421649196326,452.4857504389059,443.44731830897035,230.98131593754258),c(1706.0860557740389,471.2856934783305,510.6634396268535,267.87572790650756,398.972859721386,254.74314410401854,142.7405531037458,479.9904418994787),c(452.5686292988386,110.10568414007012,129.82968804072547,65.6022190791447,64.00891476832851,116.81591690957757,215.06243334297704,161.1235520930175),c(390.70673032993255,298.61996153140234,114.87960275118738,224.14091518707775,283.4680511168834,206.1871304488327,116.09564985771327,777.4493073019726),c(361.40372555518763,142.6369089996363,48.784488839545325,10.933703179857451,38.501602868167524,31.666965427295125,5.709622124149833,27.04171503659035),c(4115.444226141956,911.708430038611,1697.2281036596655,838.8580050768411,1121.8404535712311,769.1554047118572,947.7972726088723,602.8048976906599),c(4314.053480726338,1071.8621524241676,1632.706682936396,677.2821691967255,1406.2710447598188,1218.122603436619,2112.5601859354383,1780.2462399088647),c(273.4947112309528,286.94208594078884,165.23778477910514,244.18603768348308,253.62930889405357,446.85606769627566,692.7674843968464,1176.3146040916802))
targetgene="Vhl"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(1601.8975943527234,960.9223343133393,1140.1407149758254,716.1575582806631,1037.6181972971149,745.9329633985074,1374.115724545393,644.49420837207),c(2809.8325689560984,762.3984492729098,809.6651454176151,509.63205377224455,915.3756081906829,300.4843163878893,397.7703413157717,246.75564970888695),c(970.2550469859992,588.0644493844654,592.4954854222199,659.6667585180662,161.2254620104515,309.63255084466346,253.12658083730926,395.4850824101339))
targetgene="Ezh2"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(3490.313457614064,934.2300472490799,1306.1653463491168,1130.4234232063732,666.0777296192981,834.6004665949337,1741.434747865699,1002.7969326068921),c(1611.665262610972,589.7327173259816,861.5970206339053,427.0218519688771,556.8294314808728,1123.1217071547337,972.5389684801881,410.13267805495366),c(553.5012013007378,422.9059231743603,460.3052575989357,136.06386179378163,186.2515038747604,177.3350063928527,165.57904160034516,633.2268271068241))
targetgene="Cbfb"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(869.3224749840999,647.2879613082911,502.00812709080515,995.5744173214646,376.35316803633754,286.41010953131365,287.3843135822083,196.05243401528003),c(2064.2338919098106,784.0859325126206,1502.0901483014843,165.82783156117134,162.1880020821557,227.2984407336961,167.48224897506176,682.8033046739064),c(957.2314893083349,321.97571271262933,214.02227361865047,253.29745700003096,502.9271874654383,769.1554047118572,228.3848849659933,2077.7051053113587),c(231.16814877854344,339.49252609854955,178.6141768802708,437.95555514873456,137.6432302536989,294.1509233024303,199.83677434524415,335.76796170433016),c(527.454085945409,474.6222293613629,273.8226147768028,261.1940204077058,744.0434754273374,209.00197182014782,285.4811062074916,381.9642248918387),c(1588.8740366750592,775.7445928050396,624.7561957838547,551.5445826283648,421.1112813705823,534.1161502070444,593.8007009115826,424.78027369977343),c(1875.3923055836763,735.7061622086504,863.957560416464,1410.4477102016112,480.30749578038984,923.9716801341889,591.897493536866,1106.4568402471552),c(1393.5206715100928,975.9367457869852,1255.0203177270128,1014.4046839089968,1156.01062611673,643.1912533455054,1039.1512265952695,2408.96611450959))
targetgene="Rhoa"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(739.0868982074558,357.8434734552279,261.2330692698234,153.07184451800433,107.32321799501698,41.51891022689805,180.80470059807803,152.1096470808207),c(1722.3655028711194,1095.2179036053944,886.7761116478642,296.42484176502427,526.5094192221909,302.59544741637563,595.7039082862992,539.7075626052824),c(1748.4126182264483,613.9226024779667,585.4138660745439,174.3318229232827,96.25400717041882,109.07510313846099,234.09450709014314,140.84226581557473),c(654.4337733026371,283.6055500577564,291.1332398488995,900.2082284749301,381.1658683948585,3861.258651101519,6954.319747214497,1846.723789373816),c(289.7741583280333,331.15118639096846,392.63645049892125,182.22838633095753,138.60577032540309,80.92668942530976,220.77205546712688,267.0369359863297),c(84.65312490481872,387.03816243176163,369.817899267521,274.5574354053093,165.07562229726827,249.11346136138832,123.70847935657972,343.6551285900024),c(159.53858155138914,547.1918848173182,605.8718775233855,1036.2720902687117,241.5975579977512,547.4866467207912,738.444461390045,464.21610812813435),c(3021.465381218145,1111.0664490497986,926.118441357175,480.47551195929134,1060.7191590180153,745.2292530556786,2040.238305696207,998.2899801007937),c(1321.8911042829386,744.8816358869896,787.6334407804012,386.32417902162996,247.85406846382844,567.1905363199971,216.96564071769365,340.27491421042856),c(2142.375237975797,644.7855593960168,734.1278723757385,406.9767294724718,858.1044739242836,693.8583980291776,923.0555767375563,623.0861839681027),c(2578.664420177555,779.9152626588301,872.6128729525124,427.6292799233137,967.352772062709,582.6721638622303,765.0893646360776,600.5514214376107),c(1947.0218728108307,982.6098175530501,1123.616936497915,857.0808437099369,3403.541693546009,1256.1229619493733,1396.9542130419925,1006.1771469864659),c(3841.9495149110035,974.268477845469,1351.0156022177312,883.2002457507075,692.0663115553112,641.7838326598478,1086.7314109631848,390.9781299040355),c(2432.14939630383,864.9969276761569,832.4836966490154,975.5292948250593,1789.8432633339378,342.7069369576161,462.4793920561365,1064.767529565745),c(973.5109364054154,633.1076838054032,519.3187521629019,450.10411423746507,648.2707382927707,189.29808222094195,411.09279293878797,319.9936279329858),c(849.7871384676033,477.1246312736372,319.45971723960326,323.15167176023135,590.5183339905194,887.3787423070922,319.73883895239067,647.8744227516438),c(605.5954320113955,576.386573793852,508.30289984429487,844.9322846212063,365.2839572117394,234.3355441619839,336.86770532484013,319.9936279329858))
targetgene="Crebbp"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(2129.3516802981326,859.9921238516083,639.7062810733928,601.3536748921598,1426.0031162297546,439.81896426798784,744.1540835141949,1410.6761344087965),c(781.4134606598651,752.3888416238125,649.1484402036273,245.40089359235614,366.72776731929565,938.0458869907644,154.15979735204547,1030.965385770007),c(1015.8374988578247,382.8674925779711,298.21485919657545,31.586253630699304,42.83303319083637,14.777917199404392,72.32188023923122,103.65990764026301))
targetgene="Arid1a"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(856.2989173064356,361.18000933826033,259.65937608145094,221.10377541489513,220.42167642025908,299.07689570223175,121.8052719818631,299.71234165554307),c(3799.622952458594,1328.7754154176646,1310.099579320048,916.0013552902798,717.5736234554722,1427.124575256767,3144.0985830318414,1878.2724569165048),c(2620.9909826299645,601.4105929165952,746.717417882718,710.0832787362979,272.39884029228523,284.2989785028273,234.09450709014314,710.9717578370213),c(5193.143623968687,998.4583629974541,1554.0220235177744,664.5261821535585,500.0395672503257,641.0801223170191,1924.1426558384937,838.2931661343008),c(6199.213454568264,1427.2032239671212,1776.6996096724733,2393.2661404799087,1773.4800821149665,848.6746734515093,740.3476687647617,2848.3939838541837))
targetgene="Cd79b"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(1002.8139411801602,536.3481431974628,553.9400023070954,260.58659245326925,142.93720064807192,471.48592969528295,203.64318909467738,164.5037664725913),c(1745.156728807032,641.4490235129844,681.4091505652622,622.0062253430017,991.8975438911658,1645.2747815336888,1195.2142313220318,377.4572723857403),c(1920.9747574555017,809.1099516353638,1024.4742656304518,938.4761896044313,492.3392466766922,463.0414055813376,635.671263155348,246.75564970888695),c(433.03329278234196,343.6631959523401,74.75042644769042,61.350223398089035,144.3810107556282,220.9650476482371,102.77319823469699,128.44814642380416),c(963.743268147167,199.35801901118757,159.72985861980163,178.58381860433838,226.6781868863363,67.55619291156293,205.54639646939398,421.4000593201996),c(937.6961527918381,548.0260187880763,735.701565564111,371.7459081151533,367.69030739099986,453.8931711245634,612.8327746587487,881.1092149422356))
targetgene="Cbl"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetgenelist=c("Keap1","Ret","Inppl1","Erbb3","Tsc1","Lyn","Rtel1","Ep300","Lck","Atrx")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='bonemarrow1,bonemarrow2,bonemarrow3,bonemarrow4,bonemarrow5_vs_input_rep1,input_rep2,input_rep3 pos.'


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
targetmat=list(c(2165.16646391171,587.2303154137073,587.7744058571026,143.3529972470199,277.21154065080617,627.0059154604435,327.35166845125707,798.8573317059399),c(2966.1152610880713,656.4634349866302,889.9234980246091,931.1870541511929,503.88972753714245,781.8221908827752,1010.6031159745204,1238.285201050533),c(380.93906207168425,241.06471754909293,126.68230166398061,761.1072269089659,1145.9039553638358,1754.3498846721498,1383.631761418976,885.616167448334),c(179.07391806788576,255.24499505198074,84.97943217211122,445.2446906019729,412.4484207252446,675.5619291156293,588.0910787874328,1514.3360420490596),c(218.144591100879,287.77621991154695,216.38281340120912,283.6688547218572,371.5404676778166,738.8958599702196,119.90206460714649,252.38934034150992),c(283.2623794892011,368.6872150750833,148.71400630119462,781.7597773598078,1226.7573213869878,1665.6823814757236,1828.9822871026631,1541.37775708565))
targetgene="Keap1"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(1019.0933882772408,614.7567364487248,569.6769341908196,909.319647791478,460.0941542746019,396.8926333554322,422.5120371870876,348.16208109610074),c(325.58894194161047,317.80504285883876,198.28534173492616,869.8368307531039,740.1933151405207,771.9702460831722,740.3476687647617,910.4044062318751),c(1041.8846142131536,367.8530811043252,428.04454723730095,359.5973490264229,586.6681737037027,933.8236249337917,688.9610696474132,1331.8044655520748),c(358.1478361357715,236.8940476953024,208.51434745934696,1460.2568024654063,810.4587403749264,669.2285360301703,2502.7176977523436,1797.1473118067338))
targetgene="Ret"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(159.53858155138914,301.12236344367665,155.79562564887055,343.8042222110732,720.4612436705847,828.2670735094747,2622.61976235949,446.1882981037408),c(980.0227152442475,904.201224301788,697.9329290431726,736.8101087315049,602.5500848868218,387.74439889865806,658.5097516519473,1723.9093335826349))
targetgene="Inppl1"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(901.881369178261,637.2783536591937,520.8924453512743,1971.103712146524,438.9182726971098,1217.4188930937903,1604.403816886103,283.93800788419867),c(1286.0763206693614,736.5402961794085,575.1848603501231,556.404006263857,507.73988782395924,1130.8625209258503,521.4788206723514,537.4540863522332),c(2148.887016814629,634.7759517469194,782.1255146210976,462.2526733261956,1763.3734113620726,367.33679895662345,1526.372314522722,313.2331991738382),c(546.9894224619056,264.42046873031995,152.6482392721257,328.61852335016005,593.8872242414841,605.1908948327513,888.7978439926574,363.93641486744514),c(104.18846142131535,300.28822947291854,115.6664493453736,512.6691935444271,224.7531067429279,509.4862882080371,721.3155950175956,139.71552768905013),c(436.289182201758,450.4323442093778,237.6276714442369,934.2241939233755,388.3849189326399,1126.6402588688777,3166.937071528441,2045.0296996421453),c(231.16814877854344,524.6702676068493,208.51434745934696,379.03504356839164,429.2928719800679,695.2658187148352,388.25430444218864,1122.2311740184996),c(3076.815501348219,990.117023289873,704.2277017966624,688.8233003310195,1750.3791203940661,2448.2082827013273,471.99542892971954,1096.3161971084337),c(2038.1867765544816,897.5281525357232,575.9717069443093,464.07495718950514,466.3506647406791,774.7850874544873,371.1254380697391,401.11877304275686),c(253.95937471445617,467.94915759529806,358.0152003547278,168.86497133335396,193.9518244483939,666.4136946588552,173.1918710992116,486.7508706586263),c(2142.375237975797,674.8143823433086,767.1754293315596,362.6344887986055,1392.7954837559603,698.0806600861503,816.4759637534261,1189.8354616099755),c(758.6222347239524,527.1726695191236,349.35988781867945,521.780612860975,221.86548652781536,424.33733672575465,393.96392656633844,233.23479219059178),c(914.9049268559254,452.100612150894,441.4209393384666,487.76464741252965,800.8333396578845,377.89245409905516,1358.8900655476602,734.6332584940378),c(468.84807639591907,461.2760858292332,458.7315644105633,1556.2304192663773,905.750207473641,765.6368529977133,285.4811062074916,369.5701055000681),c(1214.4467534422072,465.4467556830237,242.3487510093542,428.2367078777502,216.57151613344232,1438.3839407420273,348.2869495731398,519.4262763278397),c(1250.2615370557842,787.422468395653,383.19429136868666,574.6268448969528,1575.6780973797559,1021.7874177873894,1442.6311900351911,806.7444985916121),c(556.7570907201539,314.46850697580635,125.8954550697944,607.427954436525,930.7762493379499,555.9311708347366,1023.9255675975367,726.7460916083656))
targetgene="Erbb3"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(416.7538456852614,373.69201889963193,135.33761420002898,734.9878248681954,362.3963369966268,364.5219575853083,108.48282035884682,705.3380672043983),c(433.03329278234196,120.94942575992552,71.60304007094555,146.39013701920254,146.78736093488868,496.81950203711904,382.5446823180388,678.296352167808))
targetgene="Tsc1"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(644.6661050443887,367.0189471335671,207.72750086516075,1198.455354103264,357.10236660225377,338.4846749006434,4910.275026768856,1453.4921832167313))
targetgene="Lyn"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(439.54507162117415,220.21136828014025,172.31940412678108,1186.3067950145335,844.628912920425,538.338412264017,1941.2715222109432,325.6273185656088),c(1253.5174264752004,474.6222293613629,359.58889354310026,431.88127560436936,558.2732415884291,752.9700668267951,1023.9255675975367,829.2792611221041))
targetgene="Rtel1"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(1764.6920653235288,624.7663440978222,569.6769341908196,1172.94338001693,662.2275693324814,259.66911650382,869.7657702454912,315.4866754268874),c(306.05360542511386,361.18000933826033,406.0128426000869,1691.6868531057223,2448.2206723796025,304.706578444862,5266.174805840862,2576.8500953617554),c(778.157571240449,658.1317029281464,723.1120200571315,736.8101087315049,250.2604186430889,744.5255427128498,502.44674692518527,788.7166885672185),c(468.84807639591907,347.83386580613063,210.8748872419056,133.63414997603553,347.958235921064,162.5570891934483,190.3207374716611,148.72943270124694),c(175.81802864846966,148.47584679494304,100.71636405583551,352.91564152762106,332.0763247379449,474.3007710665981,350.1901569478564,372.95031987964194),c(2803.3207901172664,993.4535591729054,1403.7343240282075,846.1471405300795,638.6453375757288,1210.3817896655025,1809.950213355497,452.94872686288835),c(853.0430278870194,427.0765930281508,423.32346767218365,1535.5778688155353,355.1772864588454,197.74260633488734,137.03093097959598,158.8700758399683),c(306.05360542511386,442.0910045017967,114.09275615700116,349.27107380100193,204.05849520128788,146.37175130838637,445.35052568368695,267.0369359863297),c(4216.376798143856,1201.9870518624323,1834.9262576422532,8058.1392435549415,3607.600188747297,2227.946945395919,1147.6340469541165,4775.116180211246),c(1702.8301663546229,737.3744301501666,1002.4425609932379,1293.8215429497984,746.4498256065979,551.0051984349352,1471.1793006559403,1245.0456298096808))
targetgene="Ep300"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(459.08040813767076,236.8940476953024,101.50321065002173,731.3432571415763,604.4751650302301,478.52303312357077,1545.4043882698882,327.880794818658))
targetgene="Lck"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
targetmat=list(c(1383.7530032518446,749.8864397115382,689.2776165071243,525.4251805875942,1710.9149774541943,408.1519988406927,550.0269312931006,280.5577935046249),c(791.1811289181135,754.8912435360869,594.0691786105923,710.0832787362979,1415.4151754410086,558.042301863223,1387.4381761684094,772.9423547958742),c(817.2282442734423,504.65105230865476,874.1865661408848,883.2002457507075,210.31500566736509,244.18748896158684,485.3178805527358,617.4524933354796),c(1321.8911042829386,678.150918226341,657.8037527396757,527.8548924053403,447.09986330659535,654.4506188307658,253.12658083730926,437.174393091544),c(1810.2745171953543,714.0186789689396,812.81253179436,1474.8350733718828,424.48017162154696,795.192687396522,1326.5355401774777,1103.0766258675815),c(638.1543262055566,262.7522007888037,132.97707441747033,253.9048849544675,194.91436452009808,181.5572684498254,175.0950784739282,343.6551285900024),c(130.2355767766442,195.18734915739705,110.1585231860701,544.862875129563,218.0153262409986,339.1883852434722,827.8952080017258,1167.3006990794834),c(719.5515616909591,803.2710138400571,1067.7508283106936,877.1259662063422,1301.8354469799144,916.2308663630722,1549.2108030193212,839.4199042608254),c(794.4370183375296,477.95876524439535,483.123808830336,179.79867451321144,134.27434000273425,168.18677193607854,312.1260094535242,963.3610981785313))
targetgene="Atrx"
collabel=c("input_rep1","input_rep2","input_rep3","bonemarrow1","bonemarrow2","bonemarrow3","bonemarrow4","bonemarrow5")

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
Sweave("ABE_bonemarrow_input_summary.Rnw");
library(tools);

texi2dvi("ABE_bonemarrow_input_summary.tex",pdf=TRUE);

