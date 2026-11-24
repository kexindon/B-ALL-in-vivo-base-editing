pdf(file='ABE_spleen_input.pdf',width=4.5,height=4.5);
gstable=read.table('ABE_spleen_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Ctnnb1","Smarca4","Crebbp","Epha7","Cbfb","Vhl","Epha3","Foxa1","Max","Ankrd11")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='spleen1,spleen2,spleen3,spleen4_vs_input_rep1,input_rep2,input_rep3 neg.'


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
targetmat=list(c(1423.8818441783983,624.6934912282044,853.1252857204411,798.5532344399991,921.2828785556869,1109.7587512617283,969.2550918768274),c(351.9963091669087,378.4728371245999,380.17094300112876,623.8270566479993,819.9061911470496,485.4465200051304,555.4057084615122),c(1712.2917233022526,519.2574190501757,542.8431373122298,488.6872785119994,450.5630551494994,345.41387000365046,192.04775308815275),c(392.8732999088723,402.2416626445023,271.1203238518351,38.22135139199995,126.27622466689918,28.006530000295985,49.5897953983493),c(1816.7551440872708,626.521862422043,566.9427216546152,1679.009364719998,937.8825700611949,535.6248862556607,525.6518312225027),c(274.7842155431997,290.10156275573195,177.13194491653226,121.48929549599985,120.34776341493207,322.0750950034038,293.0306091720641),c(828.8945344898173,377.86338005998704,242.8033122495323,226.5980118239997,464.1985160290237,246.2240762526022,210.0804059602798),c(122.63097222589077,218.18562913141187,43.98174142485325,363.1028382239996,889.2691877950647,672.1567200071037,1708.5938596340352),c(772.1209362370901,369.94043822001953,298.83484584557823,87.3630888959999,59.87745864486769,36.17510125038231,40.573468962285794))
targetgene="Ctnnb1"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(840.2492541403627,396.14709199837347,360.8912755272205,184.28151563999978,128.64760916768603,46.67755000049331,29.753877239009583),c(660.8446836617447,366.89315289695514,477.7742595877894,334.4368246799996,348.59352161566534,161.0375475017019,945.8126431430621),c(547.2974871562903,550.9491864100456,449.4572479854866,520.0833885839994,507.47628316838353,442.2697862546741,681.6342785664013),c(924.274179554399,548.5113581515941,635.6265370304134,610.1765740079993,637.309584586463,969.7261012602485,484.1767296166105),c(590.445421828363,354.09455454008463,584.4149203028445,222.50286703199973,197.41775969050434,364.0848900038478,375.98081238384833),c(908.3775720436354,338.8581279247625,230.15103046978,147.42521251199983,422.10644114005737,297.56938125314485,120.81877424325103),c(1598.7445267967983,576.5463831237867,505.4887815815325,382.2135139199995,951.5180309407192,423.5987662544768,665.404890981487))
targetgene="Smarca4"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(515.504272134763,261.4570807189265,200.02655004179832,55.96697882399993,33.19938301101575,145.8673437515416,89.26163171702875),c(1201.3293390277076,800.2171258367144,679.005788846707,363.1028382239996,550.1612041825466,301.0701975031818,293.0306091720641),c(1219.4968904685804,448.56039955508135,448.25226876836734,121.48929549599985,131.61183979366956,308.07183000325585,18.03265287212702),c(456.45972995192676,207.21540196837998,222.9211551670644,330.3416798879996,252.55244933379836,653.4857000069063,114.50734573800658),c(202.11400977970885,241.95445465131428,300.6423146712571,217.04267397599975,206.9032976936517,236.88856625250355,232.62122205043855),c(59.044542182836295,282.7880779803774,283.17011602302773,343.9921625279996,158.88276155271822,214.71673000226923,203.76897745503533),c(111.27625257534532,399.80383438605077,463.9169985909178,607.4464774799993,253.73814158419177,714.1665150075477,1074.7461111787704),c(2107.435967141234,811.7968100643592,709.1302692746887,303.04071460799963,526.4473591746782,851.8652875090029,866.4689705057033),c(922.0032356242899,544.2451586993039,603.0920981681932,303.04071460799963,742.2433487462806,413.0963175043658,699.6669314385283),c(1494.2811060117801,471.110310945758,562.1228047861381,438.1804927439995,513.9975905455474,457.43999000483444,477.865301111366),c(1798.587592646398,569.842355413045,668.1609758926336,592.4309465759993,1044.5948725966025,635.9816187567213,670.8146868431252),c(1358.0244702052348,717.9404221139754,860.3551610231567,834.044489303999,1706.8039944413274,1127.2628325119133,1069.3363153171322),c(2679.713837528724,711.8458514678465,1034.4746578968907,443.64068579999946,308.8728312274858,469.10937750495776,564.4220348975757),c(1696.395115791489,632.006976003559,637.4340058560923,520.0833885839994,602.9245093250538,723.5020250076462,576.1432592644583),c(679.0122351026174,462.5779120411777,397.64314164935814,992.3900879279988,538.8971278038092,381.5889712540328,476.96366846775965),c(592.7163657584721,348.60944095856865,244.6107810752112,139.23492292799983,558.4610499353006,338.4122375035765,167.7036717107813),c(422.39557100029043,421.13483164750164,389.2082871295232,369.92807954399956,246.03114195663454,329.0767275034778,442.7016280107183))
targetgene="Crebbp"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(1612.3701903774527,432.10505881053354,555.4954190919821,313.9611007199996,352.1505983668456,221.7183625023432,582.4546877697028),c(1696.395115791489,543.635701634691,531.9983243581564,327.6115833599996,322.5082921070101,656.9865162569433,348.03020043205146),c(1821.297031947489,600.3152086436891,863.9700986745145,1491.9977525519982,846.5842667809015,2403.8938250254055,608.6020344342869),c(2309.5499769209428,774.6199291229734,963.3808840868539,349.4523555839996,326.6582149833871,557.796722505895,414.75101605892144),c(1700.937003651707,499.1453359179506,575.9800657830096,225.23296355999972,196.82491356530764,177.37469000187457,95.5730602222732))
targetgene="Epha7"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(2434.4518930769427,682.5919123664281,1000.1327502089916,384.94361044799956,1067.7158714792743,301.0701975031818,574.3399939772456),c(1124.1172454039988,430.88614468130777,659.7261213727987,401.3241896159995,626.0455082077256,451.6052962547728,303.85020089534027),c(386.06046811854503,308.99473175873135,352.4564210073856,156.9805503599998,103.1552257842275,243.89019875257753,146.06448826422886))
targetgene="Cbfb"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(86.29586934414536,73.74430481815875,29.521990819422044,15.015530903999982,81.81276527714594,25.67265250027132,13.524489654095264),c(77.21209362370901,201.73028838686403,89.16846206682575,58.69707535199993,58.098920269277556,70.01632500073997,115.40897838161293),c(383.7895241884359,151.1453520239948,96.39833736954137,77.8077510479999,135.76176267004652,145.8673437515416,86.55673378620969),c(206.65589763992705,284.61644917421603,200.62903965035795,391.76885176799954,443.4489016471389,264.89509625279953,306.55509882615934),c(3470.002325206687,774.6199291229734,1250.7684273697992,390.40380350399954,432.1848252684014,550.7950900058211,528.3567291533217),c(1001.4862731781079,395.5376349337606,663.3410590241565,207.48733612799975,321.3225998566167,255.55958625270085,561.7171369667567),c(1189.9746193771623,344.3432415062785,391.01575595520217,139.23492292799983,164.2183766794886,354.74938000374914,235.32611998125762),c(315.6612062851633,80.44833252890047,99.41078541233954,125.58444028799985,106.11945641021104,32.674285000345314,586.0612183441282),c(272.5132716130906,218.18562913141187,87.9634828497065,218.40772223999974,294.64452422276474,122.52856875129494,155.9824473438987),c(252.0747762421088,104.21715804880287,37.35435573069728,1.3650482639999983,0.0,0.0,34.26204045704134),c(0.0,0.0,3.0124480427981677,0.0,0.0,0.0,0.0),c(2870.4731276578877,666.1365716218803,1299.5700856631295,948.7085434799989,1147.1572522556335,730.5036575077203,520.2420353608645),c(3009.000707394542,783.1523280275537,1250.1659377612395,780.8076070079991,1187.4707887690097,1096.9224250115928,1263.1873336924978),c(190.75929012916342,209.6532302268315,126.52281779752303,571.9552226159993,518.1475134219244,596.305701256302,174.91673285963208))
targetgene="Vhl"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(1625.9958539581073,344.95269857089136,833.2431286379732,521.4484368479993,403.13536513376266,1256.7930337632824,978.2714183128908),c(1625.9958539581073,727.0822780831686,602.4896085596336,376.75332086399953,370.5288282479436,639.4824350067584,596.8808100674044),c(390.60235597876317,282.7880779803774,404.27052734351406,230.69315661599973,257.8880644605687,225.2191787523802,419.2591792769532),c(699.4507304735993,471.110310945758,595.2597332569179,264.81936321599966,265.595064088126,371.0865225039218,409.34122019728335))
targetgene="Epha3"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(1821.297031947489,533.8843886008849,656.7136733300006,106.47376459199987,215.20314344640565,200.71346500212124,175.81836550323845),c(2012.0563220766524,697.8283389817502,844.0879415920466,393.13390003199953,290.4946013463878,549.6281512558087,176.7199981468448),c(1117.3044136136714,426.6199452290176,325.94687823076174,1265.3997407279985,531.7829743014487,436.43509250461244,1778.9212058353305),c(663.1156275918538,471.110310945758,299.43733545413784,342.6271142639996,471.31266953138424,632.4808025066843,366.0628533041785),c(1914.4057330819614,576.5463831237867,581.4024722600464,281.19994238399966,406.6924418849429,217.05060750229387,276.80122158714977),c(417.85368314007223,263.28545191276515,129.53526584032122,438.1804927439995,182.5966065605866,214.71673000226923,244.34244641732113),c(165.77890689796345,265.1138231066038,87.9634828497065,223.86791529599972,218.7602201975859,241.55632125255286,259.6702013586291),c(45.41887860218176,104.82661511341576,48.19916868477068,96.91842674399989,208.08898994404512,204.21428125215823,71.22897884490173),c(3354.1841847711235,648.4623167481068,711.5402277089272,330.3416798879996,608.2601244518243,705.9979437574613,302.04693560812757),c(1335.315030904144,714.283679726298,347.63650413890855,277.10479759199967,324.28683048260024,373.4204000039465,177.62163079045115),c(1682.7694522108345,622.2556629697528,547.6630541807069,585.6057052559993,313.0227541038628,730.5036575077203,528.3567291533217),c(3240.636988265669,816.6724665812623,983.2630411693219,877.7260337519989,511.0333599195638,665.1550875070296,825.8955015434175),c(1848.548359108798,653.3379732650099,654.303714895762,343.9921625279996,678.2159672250359,662.821210007005,326.39101698549905))
targetgene="Foxa1"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(2643.378734646979,627.1313194866559,859.752671414597,364.4678864879996,243.066911330651,135.3648950014306,292.1289765284577))
targetgene="Max"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(1971.1793313346886,607.0192363544307,795.2862832987163,382.2135139199995,627.2312004581189,758.5101875080163,533.7665250149598),c(1407.9852366676348,440.6374577151139,599.4771605168354,375.38827259999954,235.9527578282905,453.9391737547974,415.6526487025278),c(1619.1830221677799,586.9071532222057,818.1808884239823,335.80187294399957,962.18926119426,666.322026257042,600.4873406418297))
targetgene="Ankrd11"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetgenelist=c("Syk","Ercc2","Ep300","Arid5b","Lck","Raf1","Smad4","Card11","Abl1","Nfe2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='spleen1,spleen2,spleen3,spleen4_vs_input_rep1,input_rep2,input_rep3 pos.'


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
targetmat=list(c(156.6951311775271,54.24167875054653,62.65891929020189,719.3804351279991,642.0523535880367,506.4514175053524,1298.3510067931454))
targetgene="Syk"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(0.0,0.0,30.72697003654131,25.93591701599997,101.96953353383408,8.168571250086329,0.0),c(944.7126749253807,538.760045117788,500.0663751044958,520.0833885839994,663.3948140951182,457.43999000483444,511.225708924801),c(724.4311137047991,523.5236185024659,399.450610475037,981.4697018159989,1058.2303334761268,1017.5705900107541,768.191012352611),c(797.1013194682901,402.2416626445023,153.0323605741469,515.9882437919994,740.4648103706905,527.4563150055744,258.76856871502275),c(2822.7833051255966,660.6514580403644,827.8207221609365,1004.6755223039988,477.24113078335137,772.5134525081643,1291.1379456442946),c(2270.9439301090883,730.1295634062329,677.198320021028,711.1901455439992,547.7898196817598,436.43509250461244,548.1926473126614),c(331.5578137959269,441.2469147797268,186.16928904492676,305.7708111359996,903.4974947997857,260.2273412527502,553.6024431742995),c(1326.2312551837076,461.9684549765648,603.0920981681932,167.90093647199978,344.44359873928835,199.5465262521089,248.85060963535287),c(302.03554270450877,338.8581279247625,206.05144612739466,316.69119724799964,1091.4297164871427,567.1322325059937,472.45550524972793),c(815.2688709091627,298.6339616603123,222.31866555850476,313.9611007199996,354.52198286763246,131.8640787513936,287.62081331042594),c(2872.7440715879966,982.4447881559663,1038.0895955482486,1194.4172309999985,1389.6313174610877,1301.136706263751,944.0093778558495),c(451.91784209170856,355.92292573392325,189.78422669628455,293.48537675999967,138.72599329603008,168.03918000177592,107.29428458915577),c(865.2296373715627,553.99647173311,450.0597375940462,1296.7958507999983,1827.7446039814563,1427.166091265083,1108.1065189922053),c(1328.5021991138167,510.1155630809825,599.4771605168354,300.31061807999964,391.87128875502515,214.71673000226923,731.2240739647507),c(104.46342078501806,124.93869824564086,58.44149203028445,784.902751799999,661.0234295943313,1229.9534425129987,136.146529184559),c(1051.447039640508,682.5919123664281,602.4896085596336,1381.4288431679984,1577.5635391484448,801.6869212584726,1294.74447621872),c(656.3027958015265,450.38877074892,284.375095240147,109.20386111999987,133.98322429445642,541.4595800057224,156.88407998750506),c(2073.3718081895977,769.7442726060704,1051.3443669365604,887.2813715999989,1077.2014094824217,1253.2922175132453,900.7310109627447),c(2761.4678190126515,505.23990656407943,850.7153272862025,335.80187294399957,584.5462794439558,140.03265000147994,633.8477484552648))
targetgene="Ercc2"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(1230.8516101191258,456.48334139504885,436.20247659717467,868.170695903999,476.64828465815464,1107.4248737617038,426.472240425804),c(213.4687294302543,263.894908977378,310.8846380167709,600.6212361599993,323.6939843574035,383.92284875405744,344.42366985762607),c(542.7555992960721,480.8616239795641,553.6879502663032,633.3823944959993,393.05698100541855,970.8930400102608,635.6510137424774),c(18.167551440872707,74.35376188277164,18.074688256789006,185.64656390399978,246.03114195663454,340.7461150036012,190.24448780094005),c(327.0159259357087,254.1435959435719,161.4672150939818,113.29900591199986,254.9238338345852,306.9048912532435,57.70448919080646),c(122.63097222589077,108.48335750109305,77.11866989563309,173.3611295279998,168.36829955586558,217.05060750229387,209.17877331667341),c(1955.2827238239252,725.8633639539428,1074.8414616703863,615.6367670639993,640.8666613376432,1653.5522087674753,610.4052997214997),c(594.9873096885811,312.04201708179573,324.13940940508286,185.64656390399978,119.16207116453866,404.92774625427944,240.7359158428957),c(213.4687294302543,323.0122442448276,87.36099324114686,363.1028382239996,145.2473006731939,275.39754500291053,162.29387584914318),c(2940.8723894912696,878.2276301071633,1405.0057671610653,1348.6676848319983,816.9419605210661,2338.5452550247146,1482.284066088841),c(1187.7036754470532,538.760045117788,767.571761304973,713.9202420719992,489.09805328728555,319.7412175033792,609.5036670778933))
targetgene="Ep300"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(29.522271091418148,49.36602223364346,39.16182455637618,223.86791529599972,153.54714642594783,231.0538725024419,91.9665296478478),c(231.636280871127,204.77757370992845,92.18091010962392,425.89505836799947,391.87128875502515,506.4514175053524,1288.4330477134756),c(113.54719650545442,132.86164008560834,108.44812954073403,485.9571819839994,155.91853092673466,94.52203875099895,189.3428551573337))
targetgene="Arid5b"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(320.20309414538144,173.08580635005856,77.72115950419273,337.16692120799956,265.595064088126,95.68897750101128,196.5559163061845),c(163.50796296785435,106.6549863072544,26.509542776623874,440.91058927199947,524.6688207990882,508.78529500537707,560.8155043231503))
targetgene="Lck"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(3131.631679620433,990.9771870605466,898.3120063624136,936.4231091039989,771.885655006116,865.8685525091508,815.0759098201413),c(613.1548611294538,516.2101337271113,316.9095341023672,327.6115833599996,385.942827503058,564.798355005969,425.57060778219767),c(97.65058899469079,84.10507491657776,31.329459645100943,487.3222302479994,508.661975418777,668.6559037570667,231.7195894068322),c(420.12462707018136,379.69175125382566,224.7286239927433,429.99020315999945,1872.8009094964061,458.60692875484676,337.2106087087753))
targetgene="Raf1"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(93.10870113447262,102.9982439195771,13.254771388311937,96.91842674399989,129.24045529288273,43.17673375045631,49.5897953983493),c(38.606046811854505,116.40629934106052,37.35435573069728,102.37861979999988,352.74344449204233,133.03101750140593,256.96530342781006),c(910.6485159737445,603.9719510313663,476.5692803706701,1102.9589973119987,1284.69755330127,1047.9109975110748,1286.6297824262629),c(213.4687294302543,142.61295311941447,46.99418946765142,151.52035730399982,298.79444709914173,119.02775250125794,146.9661209078352),c(372.4348045378905,248.04902529744308,131.3427346660001,1097.4988042559987,740.4648103706905,875.2040625092495,367.8661185913912),c(828.8945344898173,333.9824714078595,357.2763378758627,309.8659559279996,452.93443965028627,333.74448250352714,336.3089760651689),c(0.0,80.44833252890047,22.894605125266075,199.29704654399976,330.80813785976403,243.89019875257753,75.73714206293349),c(946.9836188554898,524.7425326316917,414.51285068902786,937.7881573679989,329.62244560937063,609.1420275064377,491.38979076546127),c(788.0175437478537,421.74428871211455,392.82322478088105,526.9086299039993,730.3864262423464,484.2795812551181,414.75101605892144),c(1623.7249100279982,475.98596746266105,380.17094300112876,328.9766316239996,460.64143927784346,310.4057075032805,183.03142665208924),c(277.0551594733088,290.71101982034486,123.51036975472488,144.6951159839998,252.55244933379836,78.1848962508263,62.21265240883822),c(692.6378986832719,314.47984534024727,162.67219431110104,304.4057628719996,290.4946013463878,856.5330425090522,325.48938434189273),c(560.9231507369449,451.6076848781458,235.5734369468167,296.21547328799966,217.57452794719248,326.74285000345316,201.96571216782263),c(545.0265432261812,316.3082165340859,371.7360884812939,174.7261777919998,378.82867400069756,261.3942800027625,444.504893297931),c(177.1336265485089,210.2626872914444,148.21244370566984,99.64852327199988,368.1574437471568,390.9244812541315,118.11387631243198))
targetgene="Smad4"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(281.597047333527,343.7337844416656,187.97675787060567,345.3572107919996,291.0874474715845,245.05713750258988,445.4065259415374),c(36.335102881745414,191.97897535305793,90.97593089250466,1501.5530903999982,567.9465879384479,544.9603962557594,221.80163032716234),c(940.1707870651626,485.7372804964672,428.3701116858994,171.9960812639998,239.50983457947075,506.4514175053524,357.046526868115),c(483.7110571132358,352.87564041085886,298.23235623701856,226.5980118239997,362.8218286203864,334.9114212535395,246.14571170453382),c(2161.938621463852,449.16985661969426,601.2846293425142,416.3397205199995,805.0850380171319,665.1550875070296,871.8787663673414))
targetgene="Card11"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(54.50265432261812,222.45182858370202,243.40580185809193,84.6329923679999,106.71230253540776,186.71020000197325,210.0804059602798),c(4.541887860218177,159.67775092857516,33.73941807933948,88.72813715999989,146.4329929235873,319.7412175033792,18.03265287212702),c(4.541887860218177,131.64272595638258,87.36099324114686,102.37861979999988,93.66968778108014,96.85591625102361,79.34367263735889))
targetgene="Abl1"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(174.86268261839982,180.39929112541316,87.9634828497065,308.5009076639996,1023.8452582147178,1120.2612000118395,1071.139580604345),c(340.64158951636324,213.9194296791217,111.4605775835322,360.37274169599954,495.61936066444935,406.09468500429176,640.1591769605092),c(461.00161781214496,291.3204768849577,219.90870712426624,286.66013543999964,313.0227541038628,66.51550875070296,67.62244827047633))
targetgene="Nfe2"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
Sweave("ABE_spleen_input_summary.Rnw");
library(tools);

texi2dvi("ABE_spleen_input_summary.tex",pdf=TRUE);

