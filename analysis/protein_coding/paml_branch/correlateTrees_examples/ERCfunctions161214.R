if(!require(ape)){
  install.packages("ape")
}
if (!require(compiler)){
  install.packages("compiler")
}
if(!require(phytools)){
  install.packages("phytools")
}

simpleAUCmat=function(lab, value){
  value=t(apply(rbind(value),1,rank))
  posn=sum(lab>0)
  negn=sum(lab<=0)
  if(posn<2||negn<2){
    auc=rep(NA, nrow(value))
    pp=rep(NA, nrow(value))
  }
  else{
    stat=apply(value[,lab>0, drop=F],1,sum)-posn*(posn+1)/2
    
    auc=stat/(posn*negn)
    mu=posn*negn/2
    sd=sqrt((posn*negn*(posn+negn+1))/12)
    stattest=apply(cbind(stat, posn*negn-stat),1,max)
    pp=(2*pnorm(stattest, mu, sd, lower.tail = F))
  }
  return(list(auc=auc, pp=pp))
}


# reads trees from a 2 column file
# this will take a while for a whole genome so max.read is useful for testing
readTrees=function(file, max.read=NA, rearrange=F, computePaths=F){  
  tmp=scan(file, sep="\t", what="character")
  trees=vector(mode = "list", length = min(length(tmp)/2,max.read, na.rm = T))
  treenames=character()
  maxsp=0; # maximum number of species
  show(length(trees))
  for ( i in 1:min(length(tmp),max.read*2, na.rm = T)){  
    if (i %% 2==1){
      treenames=c(treenames, tmp[i])
    }
    else{
      trees[[i/2]]=unroot(read.tree(text=tmp[i]))
      #check if it has more species
      if(length(trees[[i/2]]$tip.label)>maxsp){
        maxsp=length(trees[[i/2]]$tip.label)
        allnames=trees[[i/2]]$tip.label
      }
    }
    
  }
  names(trees)=treenames
  treeObj=vector(mode = "list")
  treeObj$trees=trees
  treeObj$numTrees=length(trees)
  treeObj$maxSp=maxsp
  
  message(paste("max is ", maxsp))
  
  report=matrix(nrow=treeObj$numTrees, ncol=maxsp)
  colnames(report)=allnames
  
  rownames(report)=treenames
  for ( i in 1:nrow(report)){
    ii=match(allnames, trees[[i]]$tip.label)
    report[i,]=1-is.na(ii)
  }
  treeObj$report=report
    
  ii=which(rowSums(report)==maxsp)
  
  #Create a master tree with no edge lengths
  master=trees[[ii[1]]]
  master$edge.length[]=1
  treeObj$masterTree=master
  
  
  if(rearrange){
    show(treeObj$masterTree$tip)
    treeObj$masterTree=rotateConstr(treeObj$masterTree, sort(treeObj$masterTree$tip.label))
    #this gets the absolute alphabetically constrained order when all branches are present

    tiporder=treeTraverse(treeObj$masterTree)
    show(tiporder)
    #treeObj$masterTree=CanonicalForm(treeObj$masterTree)
    show(treeObj$masterTree$tip)
    for ( i in 1:treeObj$numTrees){
      treeObj$trees[[i]]=rotateConstr(treeObj$trees[[i]], tiporder)
    }  
  }
  
  
  if (computePaths){
    ap=allPaths(master)
    #  show(length(ap))
    paths=matrix(nrow=treeObj$numTrees, ncol=length(ap$dist))
    for( i in 1:treeObj$numTrees){
      paths[i,]=allPathMasterRelative(treeObj$trees[[i]], master, ap)
    }
    colnames(paths)=namePaths(ap$nodeId)
    treeObj$paths=paths
  }
  
  treeObj
}


treePlot=function(tree, vals=NULL,rank=F, nlevels=8, type="c", col=NULL){
  if(is.null(vals)){
    vals=tree$edge.length
  }
  vals=as.numeric(vals)
  layout(matrix(c(1,2), ncol=1),heights=c(10,2))
  if(is.null(col)){
    col=colorpanel(nlevels, "blue", "red")
  }
  plot.phylo(tree, use.edge.length = F,type=type,edge.color=col[cut(vals, nlevels)], edge.width=3.5, lab4ut="axial", cex=0.6)
  
  
  min.raw <- min(vals, na.rm = TRUE)
  max.raw <- max(vals, na.rm = TRUE)
  z <- seq(min.raw, max.raw, length = length(col))
  
  par(mai=c(1,0.5,0,0.5))
  image(z = matrix(z, ncol = 1), col = col, breaks = seq(min.raw, max.raw, length.out=nlevels+1), 
        xaxt = "n", yaxt = "n")
  #par(usr = c(0, 1, 0, 1))
  lv <- pretty(seq(min.raw, max.raw, length.out=nlevels+1))
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  xv <- scale01(as.numeric(lv), min.raw, max.raw)
  axis(1, at = xv, labels = lv)
}


getRoot = function(phy) phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]

treeTraverse=function(tree, node=NULL){
  if(is.null(node)){
    rt=getRoot(tree)
    ic=getChildren(tree,rt)
    return(c(treeTraverse(tree, ic[1]), treeTraverse(tree, ic[2])))
    
  }
  else{
    if (node<=length(tree$tip)){
      return(tree$tip[node])
    }
    else{
      ic=getChildren(tree,node)
      return(c(treeTraverse(tree, ic[1]), treeTraverse(tree, ic[2])))
      
    }
  }
}

#mathces node if descendands in tr1 are a subset of descendands in tr2
matchNodesInject=function (tr1, tr2){  
  desc.tr1 <- lapply(1:tr1$Nnode + length(tr1$tip), function(x) extract.clade(tr1, 
                                                                              x)$tip.label)
  names(desc.tr1) <- 1:tr1$Nnode + length(tr1$tip)
  desc.tr2 <- lapply(1:tr2$Nnode + length(tr2$tip), function(x) extract.clade(tr2, 
                                                                              x)$tip.label)
  names(desc.tr2) <- 1:tr2$Nnode + length(tr2$tip)
  Nodes <- matrix(NA, length(desc.tr1), 2, dimnames = list(NULL, 
                                                           c("tr1", "tr2")))
  for (i in 1:length(desc.tr1)) {
    Nodes[i, 1] <- as.numeric(names(desc.tr1)[i])
    for (j in 1:length(desc.tr2)) if (all(desc.tr1[[i]] %in% 
                                            desc.tr2[[j]])) 
      Nodes[i, 2] <- as.numeric(names(desc.tr2)[j])
  }
  
  iim=match(tr1$tip.label, tr2$tip.label)
  Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
  Nodes
}
matchNodesInject_c=cmpfun(matchNodesInject)
getChildren=function(tree, nodeN){
  tree$edge[tree$edge[,1]==nodeN,2]
}
getAncestors=function(tree, nodeN){
  if(is.character(nodeN)){
    nodeN=which(tree$tip.label==nodeN)
  }
  im=which(tree$edge[,2]==nodeN)
  if(length(im)==0){
    return()
  }
  else{
    anc=tree$edge[im,1]
    return(c(anc, getAncestors(tree, anc)))
  }
  
}

#computes all paths from the tips to ancestors
allPaths=function(tree){
  dd=dist.nodes(tree)
  allD=double()
  nn=matrix(nrow=0, ncol=2)
  nA=length(tree$tip.label)+tree$Nnode
  for ( i in 1:nA){
    ia=getAncestors(tree,i)
    if(length(ia)>0){
      allD=c(allD, dd[i, ia])
      nn=rbind(nn,cbind(rep(i, length(ia)), ia))   
    }
  }
  return(list(dist=allD, nodeId=nn))
}

#computes paths relative to a master tree, returns a vector that is the same length as the number of 
#paths in the master tree though with missing values
allPathMasterRelative=function(tree, masterTree, masterTreePaths=NULL){
  if(! is.list(masterTreePaths)){
    masterTreePaths=allPaths(masterTree)
  }
  nnMaster=masterTreePaths$nodeId[,1]*1000+masterTreePaths$nodeId[,2]
  
  treePaths=allPaths(tree)
  
  map=matchAllNodes_c(tree,masterTree)
  #remap the nodes
  treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
  treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]
  
  nnTree=treePaths$nodeId[,1]*1000+treePaths$nodeId[,2]
  ii=match(nnTree, nnMaster)
  #show(nnTree)
  vals=double(length(masterTreePaths$dist))
  vals[]=NA
  vals[ii]=treePaths$dist
  vals  
}

#Maps edges to paths in a master tree
edgeIndexRelativeMaster=function(tree, masterTree){
  map=matchAllNodes(tree,masterTree)
  newedge=tree$edge
  newedge[,1]=map[newedge[,1],2]
  newedge[,2]=map[newedge[,2],2]
  newedge
}

#a hack, turns a 2-tuble into a single value for indexing
namePaths=function(nodeMat, invert=F, mult=1000){
  if(invert){
    nodeMat=nodeMat[,c(2,1)]
  }
  return(nodeMat[,1]*mult+nodeMat[,2])
}
printTipDist=function(tree, node){
  if(is.character(node)){
    node=match(node, tree$tip.label)
  }
  ii=which(tree$edge[,2]==node)
  print(tree$edge.length[ii])
}
printLeaveEdges=function(tree){
  nA=tree$Nnode+length(tree$tip)
  nT=length(tree$tip.label)
  ii=match(1:nT, tree$edge[,2])
  tmp=as.data.frame(tree$edge)
  tmp[ii,2]=tree$tip.label[tmp[ii,2]]
  show(tmp[ii,])
  tmp
}
printClade=function(tree,node){
  nT=length(tree$tip.label)
  clade=character()
  if(node<=nT){
    clade=c(clade, tree$tip.label[node])
  }
  else{
    ic=getChildren(tree,node)
    for ( i in 1:length(ic)){
      clade=c(clade, printClade(tree,ic[i]))
    }
  }
  paste(clade, collapse=",")
}

#puts tree it tip sorted order
CanonicalForm=function(tree){
  par(mfrow=c(1,2))
  oo=order(tree$tip.label)
  tree$tip.label=tree$tip.label[oo]
  ii=match(1:length(oo), tree$edge[,2])
  tree$edge[ii,2]=order(oo)
  rotateConstr(tree, sort(tree$tip.label))
  
}



rescaleTree=function(tree){
  tree$edge.length=tree$edge.length/sqrt(sum(tree$edge.length^2))
  tree
}

treeSum=function(tree){
  sum(tree$edge.length^2)
}

distToVec=function(dist){
  vec=as.vector(dist)
  names=character()
  for(i in 1:nrow(dist)){
    for(j in 1:ncol(dist)){
      names=c(names, paste(rownames(dist)[i], colnames(dist)[j], sep="_"))
    }
  }
  names(vec)=names
  vec
}
#drop.tip wrapper
pruneTree=function(tree, tip.names){
  keep=intersect(tree$tip.label, tip.names)
  torm=setdiff(tree$tip.label, keep)
  tree=drop.tip(tree, torm)
  tree
}

#linear fit with intercept
residLN=function(x,y, plot=F){
  if(plot){
    plot(y,x);abline(a=0,b=1)
  }
  return(as.vector(resid(x, cbind(rep(1,length(y)), y))))
}
#linear fit no intersept
residLN0=function(x,y, plot=F){
  if(plot){
    plot(y,x);abline(a=0,b=1)
  }
  return(as.vector(resid(x, cbind(y))))
}


#linear fit in sqroot space
residSQ=function(x,y, plot=F){
  x=sqrt(x);y=sqrt(y)
  if(plot){
    plot(y,x);abline(a=0,b=1)
  }
  return(as.vector(resid(x, cbind(rep(1,length(y)), y))))
}
#linear in (1/5) space
residOptPower=function(x,y, plot=F){
  x=(x)^0.2;y=(y)^0.2
  if(plot){
    plot(y,x);abline(a=0,b=1)
  }
  return(as.vector(resid(x, cbind(rep(1,length(y)), y))))
}

#loess fit
residLO=function(x,y, plot=F){
  x=as.vector(sqrt(x));y=as.vector(sqrt(y))
  
  fit=loess(x~y, span = 0.99, family = "s")
  if(plot){
    plot(fit$x, fit$y)
    lines(sort(fit$x), fit$fitted[order(fit$x)])
  }
  return(as.vector(fit$residuals))
}

matchAllNodes=function(tree1, tree2){
  map=matchNodesInject_c(tree1,tree2)  
  map=map[order(map[,1]),]
  map
}
matchAllNodes_c=cmpfun(matchAllNodes)

checkOrder=function(tree1, tree2, plot=F){
  both=intersect(tree1$tip.label, tree2$tip.label)
  # show(c(length(tree1$tip), length(tree2$tip), length(both)))
  tree1=unroot(pruneTree(tree1, both))
  tree2=unroot(pruneTree(tree2, both))
  # show(tree1)
  #show(tree2)
  tmpe1=as.data.frame(tree1$edge)
  tmpe1[match(1:length(both),tmpe1[,2]),2]=tree1$tip
  tmpe2=as.data.frame(tree2$edge)
  tmpe2[match(1:length(both),tmpe2[,2]),2]=tree2$tip
  #  show(cbind(tmpe1, tmpe2))
  #   show(cbind(tree1$tip, tree2$tip))
  map=matchNodes(tree1,tree2, method = "descendant")  
  n=length(both)
  im=match(tree1$tip, tree2$tip)
  map=rbind(map, cbind(1:n, im))
  map=map[order(map[,1]),]
  #show(n)
  #show(map)
  edge1remap=tree1$edge
  #show(nrow(edge1remap))
  #show(nrow(tree2$edge))
  
  edge1remap[,1]=map[edge1remap[,1],2]
  edge1remap[,2]=map[edge1remap[,2],2]
  if(plot){
    tmpd1=tree1$edge.length
    tmpd2=tree2$edge.length
    nn=character(length(tree1$edge))
    iim1=match(1:length(tree1$tip.label), tree1$edge[,2])
    nn[iim1]=tree1$tip.label 
    par(mfrow=c(1,3))
    plot(tree1, use.edge.length = T); plot(tree2, use.edge.length = T)
    plotWtext(tmpd1, tmpd2, nn[])
  }
  return(all(edge1remap[,1]==tree2$edge[,1]) &&all(edge1remap[,2]==tree2$edge[,2]))
  
  #   #show(tmpd1)
  
  #  return(cbind(tmpd1, tmpd2))
}

checkOrderAll=function(treeObj){
  for (i in 1:(treeObj$numTrees-1)){
    for(j in (i+1):treeObj$numTrees){
      if(! checkOrder(treeObj$trees[[i]], treeObj$trees[[j]])){
        return(c(i,j))
      }
    }
  }
}



correlateTrees=function(tree1, tree2, mastertree, residfun=residLN, plot=F, cutoff=0.00001, Tree1Bin=F){
  both=intersect(tree1$tip.label, tree2$tip.label)
  if(length(both)<10){
    return(0)
  }
  iibad1=which(tree1$edge.length<cutoff)
  iibad2=which(tree2$edge.length<cutoff)
  show(c(length(iibad1), length(iibad2)))
  show(tree1$edge.length)
  if (!Tree1Bin){
    tree1=rescaleTree(tree1)
  }
  tree2=rescaleTree(tree2)
  tree1$edge.length[iibad1]=mastertree$edge.length[iibad1]
  tree2$edge.length[iibad2]=mastertree$edge.length[iibad2]
  if(!Tree1Bin){
    e1=residfun(t(tree1$edge.length), mastertree$edge.length, plot=F)
  }
  else{
    e1=tree1$edge.length
  }
  e2=residfun(t(tree2$edge.length), mastertree$edge.length, plot=F)
  cc=cor((e1), (e2))
  nn=character(length(e1))
  iim=match(1:length(tree1$tip.label), tree1$edge[,2])
  
  nn[iim]=tree1$tip.label
  
  if(plot){   
    plotWtext(e1, e2, nn)
    title(paste("R=", round(cc,2)))
  }
  return(list(l1=tree1$edge.length, l2=tree2$edge.length,e1=e1, e2=e2, cor=cc, names=nn, tree1=tree1, tree2=tree2))  
}


getNV=function(name1, name2, treeObj, residfun=residLN, plot=T){
  report=treeObj[["report"]]
  both=names(which(colSums(report[c(name1,name2),])==2))
  show(length(both))
  mastertree=pruneTree(treeObj[["master"]], both)
  allbranch=matrix(nrow=0, ncol=length(mastertree$edge.length))
  for ( i in 1:(length(treeObj)-3)){
    if(sum(is.na(match(both, treeObj[[i]]$tip.label)))==0){
      tmptree=(pruneTree(treeObj[[i]], both, mastertree))
      
      show(c(length(tmptree$edge.length), ncol(allbranch)))
      allbranch=rbind(allbranch, tmptree$edge.length)
    }
  }
  nv=projection(t(allbranch), method="AVE", returnNV = T)
  mastertree$edge.length=nv  
  # par(mfrow=c(1,3))
  res=correlateTrees(treeObj[[name1]], treeObj[[name2]], mastertree, residfun=residfun, plot=plot)
  res$nv=nv 
  res$master=mastertree
  return(res)
}

scaleDist=function(x){
  x/sqrt(sum(x^2)) 
}
scaleDist_c=compiler::cmpfun(scaleDist)

getProjection=function(treeObj, tree1, tree2, maxT=treeObj$numTrees){
  both=intersect(tree1$tip.label, tree2$tip.label)
  tree1=unroot(pruneTree(tree1, both))
  tree2=unroot(pruneTree(tree2, both))
  allreport=treeObj$report[1:maxT,both]

  ss=rowSums(allreport)
  iiboth=which(ss==length(both))
  torm=setdiff(treeObj$masterTree$tip.label, both)
  allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
  for ( k in 1:length(iiboth)){        
    tmptree=rescaleTree(unroot(drop.tip(treeObj$trees[[iiboth[k]]], torm)))         
    allbranch[k, ]=tmptree$edge.length           
  }
  allbranch
}

getProjectionPaths=function(treeObj, tree1, tree2, maxT=treeObj$numTrees){
  both=intersect(tree1$tip.label, tree2$tip.label)
  tree1=unroot(pruneTree(tree1, both))
  tree2=unroot(pruneTree(tree2, both))
  allreport=treeObj$report[1:maxT,both]
  ss=rowSums(allreport)
  iiboth=which(ss==length(both))
  allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
  ee=edgeIndexRelativeMaster(tree1, treeObj$masterTree)
  ii= match(namePaths(ee,T), colnames(treeObj$paths)) 
  allbranch=treeObj$paths[iiboth,ii]
  allbranch=scaleMat_c(allbranch)
  allbranch
}


######
correlateTreesAll=function(treeObj,  usePaths=T, maxDo=NULL, species.threshold=10, species.list=NULL){
  maxT=treeObj$numTrees
  if (is.null(maxDo)){
    maxDo=maxT*(maxT-1)
  }
  corout=matrix(nrow=maxT, ncol=maxT)
  maxn=treeObj$report[, species.list]%*%t(treeObj$report[, species.list])
  colnames(corout)=rownames(corout)=names(treeObj$trees)

  done=0
  todo=length(maxn[upper.tri(maxn)]>= species.threshold)
  corout[maxn<species.threshold]=0
  diag(corout)=1
  corout[lower.tri(corout)]=0
  for (i in 1:(maxT-1)){
    for(j in (i+1):maxT){
      #show(c(i,j))
      if (is.na(corout[i,j])){
        tree1=treeObj$trees[[i]]
        tree2=treeObj$trees[[j]]
        if(!is.null(species.list)){
          tree1=pruneTree(tree1, species.list)
        }
        both=intersect(tree1$tip.label, tree2$tip.label)
        tree1=unroot(pruneTree(tree1, both))
        tree2=unroot(pruneTree(tree2, both))
        allreport=treeObj$report[,both] 
        ss=rowSums(allreport)
        iiboth=which(ss==length(both))
        torm=setdiff(treeObj$masterTree$tip.label, both)
        allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
        if(length(both)<species.threshold){
          next
        }
       
        if(! usePaths){  
          torm=setdiff(treeObj$masterTree$tip.label, both)
          allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
          for ( k in 1:length(iiboth)){        
            tmptree=rescaleTree(unroot(drop.tip(treeObj$trees[[iiboth[k]]], torm)))         
            allbranch[k, ]=tmptree$edge.length           
          }
        }
        else{
          ee=edgeIndexRelativeMaster(tree1, treeObj$masterTree)
          ii= match(namePaths(ee,T), colnames(treeObj$paths))
          allbranch=treeObj$paths[iiboth,ii]
          allbranch=scaleMat_c(allbranch)
        }
        nb=length(both)
        proj=t(projection(t(allbranch), method="AVE", returnNV = F))
        #i1=match(i, iiboth)
        #j1=match(j,iiboth)
        #corout[i,j]=cor(proj[i1, ], proj[j1,])
        #  tmpcor=cor(t(proj)) 
        ai=which(maxn[iiboth, iiboth]==nb, arr.ind = T)
        for (m in 1:nrow(ai)){
          k=sort(ai[m,])[1]
          l=sort(ai[m,])[2]
     
          tmpcor=cor(proj[k,], proj[l,])
          if (is.na(tmpcor)){
            tmpcor=0
          }
          corout[iiboth[k], iiboth[l]]=tmpcor
        }
        done=done+nrow(ai)
        message(paste("Done with",done, "out of", todo))
        message(paste(sum(is.na(corout))," left"), appendLF = T)
        if(done>=maxDo){
          message("Finished. Exiting because done>=maxDo")
		  diag(corout)=1
          return(corout)
        }
      }
    }
  }
  message("Finished. Exiting due to completed traversal of matrix")
  diag(corout)=1
  corout
}


removeEmpties=function(table){
	#create vector of column/row indices to keep
	keep = apply(derc, 1, sum) != "1" | apply(derc, 2, sum) != "1"  	#the strange "1" is due to an unexplained phenomenon of the first column being summed to a nonnumeric 1. ???
	table[keep,keep]
}


correlateTreesChar=function(treeObj,  char, usePaths=F, maxDo=NULL){
  maxT=treeObj$numTrees
  if (is.null(maxDo)){
    maxDo=maxT
  }
  corout=matrix(nrow=maxT, ncol=1)
  rownames(corout)=names(treeObj$trees)
  charReport=as.vector(as.numeric(names(char) %in% colnames(treeObj$report)))

  maxn=treeObj$report%*%(charReport)
  
  done=0
  todo=length(maxn>=10)
  corout[maxn<10]=0

  for (i in 1:maxT){

      if (is.na(corout[i,1])){
        tree1=treeObj$trees[[i]]

        both=tree1$tip.label
        charTree=edgeVars(tree1, char)
       
        allreport=treeObj$report[,both] 
        ss=rowSums(allreport)
        iiboth=which(ss==length(both))
        torm=setdiff(treeObj$masterTree$tip.label, both)
        allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
        if(length(both)<10){
          next
        }
        
        if(! usePaths){  
          torm=setdiff(treeObj$masterTree$tip.label, both)
          allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
          for ( k in 1:length(iiboth)){        
            tmptree=rescaleTree(unroot(drop.tip(treeObj$trees[[iiboth[k]]], torm)))         
            allbranch[k, ]=tmptree$edge.length           
          }
        }
        else{
         
          ii= match(namePaths(edgeIndexRelativeMaster(tree1, treeObj$masterTree),T), colnames(treeObj$paths)) 
          ii2=match(namePaths(edgeIndexRelativeMaster(charTree, treeObj$masterTree),T), colnames(treeObj$paths)) 
       stopifnot(all(ii=ii2))
          allbranch=treeObj$paths[iiboth,ii]
          allbranch=scaleMat_c(allbranch)
        }
       # message("done")
        nb=length(both)
        proj=t(projection(t(allbranch), method="AVE", returnNV = F))
        # i1=match(i, iiboth)
        #j1=match(j,iiboth)
        #corout[i,j]=cor(proj[i1, ], proj[j1,])
        #  tmpcor=cor(t(proj)) 
        ai=which(maxn[iiboth, 1]==nb)
    
    

        corout[iiboth[ai]]=cor(t(proj[ai, ,drop=F]), cbind(charTree$edge.length), method="s")

        done=done+length(ai)
      #  message(paste("Done with",done, "out of", todo))
      #  message(paste(sum(is.na(corout))," left"), appendLF = T)
        if(done>=maxDo){
          message("DONE")
          return(corout)
        }
        #generate the projection
      }
    }

  corout
  }
  
  
#assuming the binTree is already in canonical form
correlateTreesBinary=function(treeObj, binTree, usePaths=T, maxDo=NULL, species.list=NULL, useSQ=F){
     print(Sys.time())
     maxT=treeObj$numTrees
     if (is.null(maxDo)){
          maxDo=maxT
     }
     corout=matrix(nrow=maxT, ncol=1)
     pout=matrix(nrow=maxT, ncol=1)
     rownames(corout)=rownames(pout)=names(treeObj$trees)
     binReport=as.vector(as.numeric(binTree$tip.label %in% colnames(treeObj$report)))
     names(binReport)=colnames(treeObj$report)
     if(is.null(species.list)){
          maxn=treeObj$report%*%(binReport)
     }else{
          maxn=treeObj$report[, species.list]%*%(binReport[species.list])
     }
     done=0
     todo=length(maxn>=10)
     corout[maxn<10]=0
     for (i in 1:maxT){
          if(i==18979){print('lim2')}
          if (is.na(corout[i,1])){
               tree1=treeObj$trees[[i]]
               if(! is.null(species.list)){
                    tree1=pruneTree(tree1, species.list)
               } 
               both=tree1$tip.label
               bothIndex=match(both, colnames(treeObj$report))
               allreport=treeObj$report[,bothIndex] 
               ss=rowSums(allreport)
               iiboth=which(ss==length(both))
               allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
               if(length(both)<10){
                    next
               }
               if(! usePaths){
                    binTreeUse=unroot(pruneTree(binTree,tree1$tip.label))
                    torm=setdiff(treeObj$masterTree$tip.label, both)
                    allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
                    for ( k in 1:length(iiboth)){        
                         tmptree=rescaleTree(unroot(drop.tip(treeObj$trees[[iiboth[k]]], torm)))         
                         allbranch[k, ]=tmptree$edge.length           
                    }
               }
               else{
                    binTreeUse=unroot(pruneTree(binTree,tree1$tip.label))
                    ii= match(namePaths(edgeIndexRelativeMaster(tree1, treeObj$masterTree),T), colnames(treeObj$paths)) 
                    ii2=match(namePaths(edgeIndexRelativeMaster(binTreeUse, treeObj$masterTree),T), colnames(treeObj$paths)) 
                    stopifnot(all(ii=ii2))
                    allbranch=treeObj$paths[iiboth,ii]
                    allbranch=scaleMat_c(allbranch)
               }
               nb=length(both)
               if(!useSQ){
                    proj=t(projection(t(allbranch), method="AVE", returnNV = F))
               }
               else{
                    proj=projectionSQ(allbranch)
               }
               ai=which(maxn[iiboth, 1]==nb)
               tmp=simpleAUCmat(binTreeUse$edge.length, (proj[ai, ,drop=F]))
               corout[iiboth[ai]]=tmp$auc
               pout[iiboth[ai]]=tmp$pp
               if(i==18979){
                    print(length(iiboth))
                    print(iiboth[ai])
                    print(tmp$auc)
                    print(tmp$pp)
               }
               done=done+length(ai)
               if(done>=maxDo){
                    message("DONE")
                    print(Sys.time())
                    return(list(r=corout, p=pout))
               }
          }
     }  
     print(Sys.time())
     return(list(r=corout, p=pout))
}

# NATHAN #
#pseudo is table of genes x species with 0=lesion and 1=none
readPseudo=function(pseudoFile){
	pseudoTable=read.table(pseudoFile, header=TRUE, row.names=1, sep="\t")
	print("pseudoFile table read with dimensions:")
	show(dim(pseudoTable))
	pseudoTable
}


### ZELIA ###  
#assuming the binTree is already in canonical form
correlateTreesBinaryPseudo=function(treeObj,  binTreeObj, pseudoObj, usePaths=F, maxDo=NULL, species.list=NULL, useSQ=F){
  maxT=treeObj$numTrees
  if (is.null(maxDo)){    
    maxDo=maxT
  }
  corout=matrix(nrow=maxT, ncol=1)
  pout=matrix(nrow=maxT, ncol=1)
  rownames(corout)=rownames(pout)=names(treeObj$trees)

  binTree=binTreeObj$trees[[1]]
  binTreePruned=pruneTree(binTree,species.list)
  show(binTreePruned$tip.label )
  binReport=as.vector(as.numeric(binTree$tip.label %in% colnames(treeObj$report)))
  show((binReport))
  names(binReport)=colnames(treeObj$report)
  maxn=treeObj$report[, species.list]%*%(binReport[species.list])
  
  done=0
  todo=length(maxn>=10)
  corout[maxn<10]=0
  
  for (i in 1:maxT){
    if (is.na(corout[i,1])){
      tree1=treeObj$trees[[i]]
      if(! is.null(species.list)){
        tree1=pruneTree(tree1, species.list)
      } 
      
      #both is vector of species in tree1
      both=tree1$tip.label     
      #bothIndex are those species columns in treeObj that are in both
      bothIndex=match(both, colnames(treeObj$report))
      #report for just those both species columns
      allreport=treeObj$report[,bothIndex] 
      #total number of species for each GENE in truncated report
      ss=rowSums(allreport)
      #row(gene) indices of ss(genes) for which minimum number of species are present
      iiboth=which(ss==length(both))
      #  torm=setdiff(treeObj$masterTree$tip.label, both)
      binTreeUse=unroot(pruneTree(binTree,tree1$tip.label))
      allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
      if(length(both)<10){
        next
      }
      
      if(! usePaths){  
        torm=setdiff(treeObj$masterTree$tip.label, both)
        allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
        for ( k in 1:length(iiboth)){        
          tmptree=rescaleTree(unroot(drop.tip(treeObj$trees[[iiboth[k]]], torm)))         
          allbranch[k, ]=tmptree$edge.length           
        }
      }
      else{
        ii= match(namePaths(edgeIndexRelativeMaster(tree1, treeObj$masterTree),T), colnames(treeObj$paths)) 
        ii2=match(namePaths(edgeIndexRelativeMaster(binTreeUse, treeObj$masterTree),T), colnames(treeObj$paths)) 
        stopifnot(all(ii=ii2))
        allbranch=treeObj$paths[iiboth,ii]    
        allbranch=scaleMat_c(allbranch)
      }
      # message("done")
      nb=length(both)
      if(!useSQ){
      	proj=t(projection(t(allbranch), method="AVE", returnNV = F))
      }
      else{
      	proj=projectionSQ(allbranch)
      }
      # i1=match(i, iiboth)
      # j1=match(j, iiboth)
      # corout[i,j]=cor(proj[i1, ], proj[j1,])
      # tmpcor=cor(t(proj)) 
      ai=which(maxn[iiboth, 1]==nb)
      # show(iiboth[1])
      
	  #ai is a vector of indeces for genes to compute with the current projection
	  ainames = names(ai)
	  #pcheck are species differences between the binTree and tree in question, for asking if it is a pseudogene
	  pcheck=setdiff(binTreePruned$tip.label,tree1$tip.label)
	  #create variable report if binary variable is 1 or 0 for pcheck species
	  nsp=length(binTree$tip)
	  ledge=which(binTree$edge[,2]<=nsp)
	  status=binTree$edge.length[ledge]
	  names(status)=binTree$tip
	  
	  #add point for removed species that have a pseudogene. Add the max of the projection for that gene.
	  #iterate through ai
	  for ( i in 1:length(ai) ){
		#initialize gene-pseudogene projection vector "pproj"
		pproj=as.numeric(proj[ai[i], ,drop=F])
		maxpproj=max(pproj)

		#initialize binaryTreeUse Pseudogene vector "pbinedge"
		pbinedge=binTreeUse$edge.length
		
		#search pcheck for pseudogenes for that gene
		for (j in 1:length(pcheck)){
			#if pseudogene=0 then add as max in projection
			if ( length(pcheck)>0 ) {
				if ( !pseudo[ ainames[i],pcheck[j] ] ) {
					pproj=c( pproj , maxpproj )
					pbinedge=c( pbinedge , as.numeric(status[pcheck[j]]) )
				}
			}
		}

	 	#send to wilcox.test and cor()
	    corout[iiboth[ai[i]]]=cor(pproj,pbinedge)[1]
	    if ( sum(pbinedge == 1)>0 ) {
		    pout[iiboth[ai[i]]]=wilcox.test(pproj[pbinedge==1],pproj[pbinedge==0])$p.value
		} else {
		    pout[iiboth[ai[i]]]="NA"
		}	    
	  }

      done=done+length(ai)
      #show(length(ai))
        message(paste("Done with",done, "out of", todo))
      #  message(paste(sum(is.na(corout))," left"), appendLF = T)
      if(done>=maxDo){
        message("DONE")
        return(list(r=corout, p=pout))
      }
      #generate the projection
    }
  }
  
  return(list(r=corout, p=pout))
}

#modified to handle binary trees
correlateTreesProj=function(treeIn1, treeIn2, treeObj, residfun=residLN, plot=F, cutoff=-1, usePaths=T, tree1Bin=F, useIndex=F, species.list=NULL){
  if(is.character(treeIn1)){
    tree1=treeObj$trees[[treeIn1]]
  }
  else{
    tree1=treeIn1
  }
  if(is.character(treeIn2)){
    tree2=treeObj$trees[[treeIn2]]
  }
  else{    
    tree2=treeIn2
  }
  
  both=intersect(tree1$tip.label, tree2$tip.label)
  if(!is.null(species.list)){
    both=intersect(both, species.list)
  }
  if(is.character(treeIn1)){
    bothIndex=which(colSums(treeObj$report[c(treeIn1, treeIn2),])==2)
  }
  
  torm=setdiff(treeObj$masterTree$tip.label, both)
  tree1=pruneTree(tree1, both)
  tree1=unroot(tree1)
  if(tree1Bin){ #fix any edges that were created through pruning
    tree1$edge.length[tree1$edge.length>1]=1
  }
  tree2=pruneTree(tree2, both)
  tree2=unroot(tree2)
  allreport=treeObj$report[,both]
  ss=rowSums(allreport)
  iiboth=which(ss==length(both))
  if (! usePaths){
    allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
    for ( i in 1:length(iiboth)){            
      tmptree=rescaleTree(drop.tip(treeObj$trees[[iiboth[i]]], torm)) 
      allbranch[i, ]=tmptree$edge.length
    }
    
  }
  else{
    if(! useIndex){
      ee=edgeIndexRelativeMaster(tree1, treeObj$masterTree)
      ii= match(namePaths(ee,T), colnames(treeObj$paths)) 
      allbranch=treeObj$paths[iiboth,ii]
      allbranch=scaleMat_c(allbranch)
      nv=projection(t(allbranch), method="AVE", returnNV = T)
      mastertree=treeObj$master
      mastertree$edge.length=nv
      res=correlateTrees(tree1, tree2, mastertree, residfun=residfun, plot=plot, cutoff=cutoff, Tree1Bin=tree1Bin)
      
      res$nv=nv
      res$allbranch=allbranch
    }
    else{
      allbranch=getBranch(treeObj, bothIndex)
      show(rownames(allbranch)[1])
      allbranch=scaleMat_c(allbranch)
      nv=projection(t(allbranch), method="AVE",returnNV = T)
      rr=resid(allbranch, model.matrix(~0+nv))
      rownames(rr)=rownames(allbranch)
      show(dim(rr))
      show(rownames(allbranch)[1])
      plot(rr[name1,], rr[name2,])
      res=list()
    }
  }
  
  return(res)
}


# Older version. Doesn't work on binary trees.
#correlateTreesProj=function(name1, name2, treeList, residfun=residLN, plot=F, cutoff=-1, usePaths=F, species.list=NULL){
#  tree1=treeList$trees[[name1]]
#  tree2=treeList$trees[[name2]]
#  both=intersect(tree1$tip.label, tree2$tip.label)
#  if (! is.null(species.list)) {
#  	both=intersect(both,species.list)
#  }
#  torm=setdiff(treeList$masterTree$tip.label, both)
#  tree1=pruneTree(tree1, both)
#  tree2=pruneTree(tree2, both)
#  allreport=treeList$report[,both]
#  ss=rowSums(allreport)
#  iiboth=which(ss==length(both))
#  if (! usePaths){
#    allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
#    for ( i in 1:length(iiboth)){        
#      # tmptree=rescaleTree(pruneTree(treeList$trees[[iiboth[i]]], both))      
#      tmptree=rescaleTree(drop.tip(treeList$trees[[iiboth[i]]], torm)) 
#      #   show(c(length(tmptree$edge.length), ncol(allbranch)))
#      allbranch[i, ]=tmptree$edge.length
#    }
#    show(dim(allbranch))
#  }
#  else{
#    ee=edgeIndexRelativeMaster(tree1, treeList$masterTree)
#    ii= match(namePaths(ee,T), colnames(treeList$paths)) 
#    allbranch=treeList$paths[iiboth,ii]
#    allbranch=scaleMat_c(allbranch)
#  }
#  nv=projection(t(allbranch), method="AVE", returnNV = T)
#  #nv=colMeans(allbranch)
#  #  show(c(is.rooted(tree1), is.rooted(tree2), is.rooted(mastertree)))
#  mastertree=treeList$master
#  mastertree$edge.length=nv  
#  # par(mfrow=c(1,3))
#  res=correlateTrees(tree1, tree2, mastertree, residfun=residfun, plot=plot, cutoff=cutoff)
#  res$nv=nv 
#  res$allbranch=allbranch
#  return(res)
#}
#


correlateTreesAll1=function(name1, name2, treeList, residfun=residLN, plot=F, cutoff=-1, usePaths=F){
  tree1=treeList$trees[[name1]]
  tree2=treeList$trees[[name2]]
  both=intersect(tree1$tip.label, tree2$tip.label)
  torm=setdiff(treeList$mastertree$tip.lables, both)
  tree1=pruneTree(tree1, both)
  tree2=pruneTree(tree2, both)
  allreport=treeList$report[,both]
  ss=rowSums(allreport)
  
  iiboth=which(ss==length(both))
  if (! usePaths){
    allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
    for ( i in 1:length(iiboth)){        
      # tmptree=rescaleTree(pruneTree(treeList$trees[[iiboth[i]]], both))      
      tmptree=rescaleTree(drop.tip(treeList$trees[[iiboth[i]]], torm)) 
      #  show(c(length(tmptree$edge.lenght), ncol(allbranch)))
      allbranch[i, ]=tmptree$edge.length
    }
  }
  else{
    allbranch=treeList$paths[iiboth,getEdgeIndex(tree1, treeList$masterTree)]
    #allbranch=scaleMat_c(allbranch)
  }
  #nv=projection(t(allbranch), method="AVE", returnNV = T)
  nv=colMeans(allbranch)
  nn=names(iiboth)
  rownames(allbranch)=names(iiboth)
  allbranch=resid(allbranch, model.matrix(~1+nv))
  cc=cor(t(allbranch))
  show(dim(cc))
  show(length(nn))
  ii=which(trees$inter[nn,nn]==length(both))
  return(list(cc,ii))
}



edgeVars=function(tree, tip.vals){
  tip.vals=tip.vals[tree$tip.label]
  res=fastAnc(tree, x=tip.vals)
  
  res=c(tip.vals, res)
  names(res)[1:length(tip.vals)]=as.character(1:length(tip.vals))
  evals=matrix(nrow=nrow(tree$edge), ncol=2)

  evals[,1]=res[tree$edge[,1]]
  evals[,2]=res[tree$edge[,2]]
  newtree=tree
  newtree$edge.length=evals[,2]-evals[,1]
  newtree
}


resid=function(dat, lab){
  if (is.null(dim(lab))){
    mod=model.matrix(~1+lab);
  }
  else{
    mod=lab
  }
  n=dim(dat)[2]
  Id=diag(n)
  resid=dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% 
                   t(mod))
}

plotWtext=function(x,y, labels, ...){plot(x,y,...); text(x,y, labels=labels, ...)}

scaleMat=function(mat){t(apply(mat,1,scaleDist_c))}
scaleMat_c=cmpfun(scaleMat)


projection <- function(protein, rna=rvector, method=c("RNA","AVE","PCA"), returnNV=F)
{
  # *************************************************************************
  # FILE: projection.R
  # AUTHOR: Tetsuya Sato <sato@kuicr.kyoto-u.ac.jp> 
  # COPYRIGHT (C) 2005 
  # CREATE DATE: 3/28/2005
  # MEHOD: Projection operator method
  # VERSION: 0.01
  # UPDATE DATE: 3/28/2003
  # DESCRIPTION: Given a data matrix consisting of phylogenetic vectors,
  #              compute a data matrix by removing the phylogenetic
  #              relationship among organisms applying a projection operator.
  # INPUT(arguments):
  #   1st argument: data matrix consisting of phylogenetic vectors
  #   2nd argument: phylogenetic vector based on rRNA 16S
  #   3rd argument: type of the projection
  #       (if input RNA, use RNA vector)
  #       (if input AVE, use the mean of the data matrix)
  #       (if input PCA, use the PC 1 score in the PCA analysis)
  # OUTPUT:
  #   a data matrix as a result of removing the phylogenetic relationship
  #   among organisms from the original data matrix
  # *************************************************************************
  
  ###projection(ribosomal RNA)###
  if (match.arg(method) == "RNA") {
    rvector <- rna[,1]
    rvector.sc <- rvector/sqrt(sum(rvector^2))
    normv=rvector.sc
    result <- as.matrix(protein) - rvector.sc %*% t(rvector.sc) %*% as.matrix(protein)
  }
  
  ###projection(average vector)###
  if (match.arg(method) == "AVE") {
    n <- ncol(protein)
    m <- nrow(protein)
    protein.sc <- matrix(0,m,n)	
    for (i in 1:n) {
      unit <- protein[,i] / sqrt(sum(protein[,i]^2,na.rm=TRUE))
      protein.sc[,i] <- unit
    }
    mean.clean = function(x){mm=mean(x,na.rm=TRUE); mm }
  #  average <- apply(protein.sc, 1, mean.clean)
  average=rowMeans(protein.sc, na.rm=T) 
  avector.sc <- average/sqrt(sum(average^2))
    #    result <- average/sqrt(sum(average^2))
    normv=avector.sc
  Id=diag(nrow(protein))
  #result=(Id-avector.sc %*% t(avector.sc))%*%protein
    result <- as.matrix(protein) - avector.sc %*% t(avector.sc) %*% as.matrix(protein)
  }
  
  ###projection(PCA)###
  if (match.arg(method) == "PCA") {
    pvector <- as.vector(prcomp(protein)$x[,1])
    print(pvector)
    pvector.sc <- pvector/sqrt(sum(pvector^2))
    normv=pvector.sc
    result <- as.matrix(protein) - pvector.sc %*% t(pvector.sc) %*% as.matrix(protein)
  }
  if(returnNV){
    normv
  }
  else{
    result
  }
}

naresid=function(data, lab){
  if (is.null(dim(lab))){
    mod=model.matrix(~1+lab);
  }
  else{
    mod=lab
  }
  resid=matrix(nrow=nrow(data), ncol=ncol(data))
  for ( i in 1:nrow(data)){
    ii=which(!is.na(data[i,]))
    iiinv=which(is.na(data[i,]))
    if(length(ii) > 10) {
     dat=data[i,ii,drop=F]
     modtmp=mod[ii,]
     n=dim(dat)[2]
     Id=diag(n)
     resid[i,ii]=dat %*% (Id - modtmp %*% solve(t(modtmp) %*% modtmp) %*% 
                              t(modtmp))
    }
    resid[i, iiinv]=0;
  
  }
  rownames(resid)=rownames(data)
  colnames(resid)=colnames(data)
  resid
}

dupRM=function(x, max.rep=10){
  xa=round(x*10000)
  tt=table(xa)
  nn=names(which(tt>max.rep))
  for( n in nn){
    ii=which(as.factor(xa)==n)
    x[ii]=NA
  }
  x
}

projectRLM=function(data, nv){
  resid=matrix(nrow=nrow(data), ncol=ncol(data))
  rownames(resid)=rownames(data)
  colnames(resid)=colnames(data)
  for ( i in 1:nrow(data)){
    print(i)
    ii=which(!is.na(data[i,]))
    iir=which(is.na(data[i,]))
    y=data[i,ii]
    nvtmp=as.vector(nv[ii])
    lmres=rlm(y~1+nvtmp, psi=psi.huber)
    resid[i,ii]=lmres$residuals
    resid[i, iir]=0
  }
  resid
}

edgeVarsDiff=function(mastertree,tip.vals){
  
  cm=intersect(mastertree$tip.label, names(tip.vals))
  mastertree=pruneTree(mastertree, cm)
  
  #sets edge length to the difference between two nodes, the evolutionary change (i.e. the change between species A and its ancestral species)
  tip.vals=tip.vals[mastertree$tip.label]
  res=fastAnc(mastertree, x=tip.vals)
  
  res=c(tip.vals, res)
  names(res)[1:length(tip.vals)]=as.character(1:length(tip.vals))
  evals=matrix(nrow=nrow(mastertree$edge), ncol=2)
  
  evals[,1]=res[mastertree$edge[,1]]
  evals[,2]=res[mastertree$edge[,2]]
  newtree=mastertree
  newtree$edge.length=evals[,2]-evals[,1]
  newtree
}

plotContinuousCharXY=function(gene, treeListObj, tip.vals, tip.vals.ref=NULL,  col=NULL, residfun=residLO, useDiff=T, xlab="Char", cutoff=0.000004*3){
  #get the tree projection
  tip.vals=tip.vals[!is.na(tip.vals)]
  
  stopifnot(gene %in% names(treeListObj$trees))
  tree=treeListObj$trees[[gene]]
  stopifnot(!is.null(names(tip.vals)))
  both=intersect(tree$tip.label, names(tip.vals))
  
  stopifnot(length(both)>10)
  
  
  torm=setdiff(treeListObj$masterTree$tip.label, both)
  tree=pruneTree(tree, both)
  #  tip.vals=tip.vals[both]
  allreport=treeListObj$report[,both]
  ss=rowSums(allreport)
  iiboth=which(ss==length(both))
  
  
  ee=edgeIndexRelativeMaster(tree, treeListObj$masterTree)
  ii= match(namePaths(ee,T), colnames(treeListObj$paths)) 
  
  allbranch=treeListObj$paths[iiboth,ii]
  nv=projection(t(allbranch), method="AVE", returnNV = T)
  
 # iibad=which(allbranch<cutoff)
 # allbranch=scaleMat_c(allbranch)
  #allbranch[iibad]=NA
    

  
 
  
  
  tree=rescaleTree(tree)
  iigood=which(tree$edge.length>=cutoff)
  show(length(iigood))
  #tree$edge.length[iibad]=NA
  proj=double(length(tree$edge.length))
  proj[]=NA
  proj[iigood]=residfun(tree$edge.length[iigood], nv[iigood])
show(proj[1:10])
  cm=intersect(names(tip.vals), ucsctreeUse$tip)
  master.tree=pruneTree(ucsctreeUse, cm)
  treeChar=edgeVarsDiff(master.tree, tip.vals)
  treeChar=unroot(pruneTree(treeChar,both))
  
  
  nn=nameEdges(tree)
  
  
  nn[nn!=""]=speciesNames[nn[nn!=""], ]
  #par(mfrow=c(2,2), mai=rep(0.7,4))
  #plotWtext(sqrt(nv), sqrt(tree$edge.length), xlab="char", ylab="Gene branch length", labels = nn)
  if(!is.null(tip.vals.ref)){
    treeCharRef=edgeVarsDiff(master.tree, tip.vals.ref)
    show(treeCharRef)
    treeCharRef=unroot(pruneTree(treeCharRef, both))
    show(treeChar)
    proj=resid(rbind(proj), model.matrix(~1+treeCharRef$edge.length))[1,]
    treeChar$edge.length=resid(rbind(treeChar$edge.length), model.matrix(~1+treeCharRef$edge.length))[1,]
  }
  plotWtext(treeChar$edge.length, proj, xlab=xlab, ylab="relative gene branch length", labels = nn)
  stat=cor.test(treeChar$edge.length, proj, method="s")
  
  branch=tree$edge.length
  
  mtext(gene,side = 3, line=2, font=2)
  mtext(paste0("r=", round(stat$estimate,2), ";  p-value=", format.pval(stat$p.value)), side = 3, line=0.5, cex=.7)
  
  
}

correlateTreesChar2=function(treeListObj,  char, char.ref=NULL,  maxDo=NULL, cutoff=0.000004*3, start=NULL){
  maxT=treeListObj$numTrees
  if (is.null(maxDo)){
    
    maxDo=maxT
  }
  
  if(!is.null(char.ref)){
    cm=intersect(cm, names(char.ref))
    char.ref=char.ref[cm]
    
  }
  cm=intersect(names(char), ucsctreeUse$tip)
  master.tree=pruneTree(ucsctreeUse, cm)
  if(!is.null(char.ref)){
    refTree=edgeVarsDiff(master.tree, char.ref)
  }
  charTree=edgeVarsDiff(master.tree, char)
  corout=matrix(nrow=maxT, ncol=4)
  courout=as.data.frame(corout)
  rownames(corout)=names(treeListObj$trees)
  charReport=cbind(as.vector(as.numeric(cm %in% colnames(treeListObj$report[, cm]))))
  
  maxn=treeListObj$report[,cm]%*%(charReport)
  
  done=0
  todo=length(maxn>=10)
  corout[maxn<10]=0
  istart=1;
  if (!is.null(start) && ! is.numeric(start)){
    istart=match(start, names(treeListObj$trees))
  }
  else if (is.numeric(start)){
    istart=start
  }
  for (i in istart:maxT){
    message(paste("i=", i))
    if (is.na(corout[i,1])){
      tree1=treeListObj$trees[[i]]
      
      both=intersect(tree1$tip.label, cm)
      tree1=unroot(pruneTree(tree1,both))
      charTreeUse=unroot(pruneTree(charTree, both))
      if(!is.null(char.ref)){
        refTreeUse=unroot(pruneTree(refTree,both))
      }
      allreport=treeListObj$report[,both] 
      ss=rowSums(allreport)
      iiboth=which(ss==length(both))
      
      if(length(both)<10){
        next
      }
      
      
      torm=setdiff(treeListObj$masterTree$tip.label, both)
      
      #        allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
      #       for ( k in 1:length(iiboth)){        
      #         tmptree=rescaleTree(unroot(drop.tip(treeListObj$trees[[iiboth[k]]], torm)))         
      #         allbranch[k, ]=tmptree$edge.length           
      #       }
      
      ee=edgeIndexRelativeMaster(tree1, treeListObj$masterTree)
      ii= match(namePaths(ee,T), colnames(treeListObj$paths)) 
      
      allbranch=treeListObj$paths[iiboth,ii]
      nv=t(projection(t(allbranch), method="AVE", returnNV = T))
      iibad=which(allbranch<cutoff)
      allbranch=scaleMat_c(allbranch)
      allbranch[iibad]=NA
      
      if(cutoff>0){
        proj=naresid(allbranch, as.vector(nv))
      }
      else{
        proj=resid(allbranch, cbind(as.vector(nv)))
      }
      
      nb=length(both)
      ai=which(maxn[iiboth, 1]==nb)
      message("ai")
      show(length(ai))
      ee=charTreeUse$edge.length
      if(!is.null(char.ref)){
        ee=resid(rbind(ee), mod<-model.matrix(~1+refTreeUse$edge.length))
        proj=naresid(proj, refTreeUse$edge.length)
        
      }
      ee=as.vector(ee)
      
      for (j in ai){
        if(names(trees$trees)[iiboth[j]]=="APOE"){
          show(sum(!is.na(proj[j,])))
        }
        corout[iiboth[j],4]=paste(corout[iiboth[j],4], as.character((i)))
        if (sum(!is.na(proj[j,]))>0.5*ncol(proj)){
          if(names(trees$trees)[iiboth[j]]=="APOE"){
            show("Doing APOE")
          }
        # show(names(trees$trees)[iiboth[j]])
     ii2=which(!is.na(proj[j,]))

          cres=cor.test(proj[j,ii2], ee[ii2], method="s")
          
          corout[iiboth[j],1:3]=c(cres$estimate, cres$statistic, cres$p.value)
        
          }
      }
      
      done=done+length(ai)
      #  message(paste("Done with",done, "out of", todo))
      #  message(paste(sum(is.na(corout))," left"), appendLF = T)
      if(done>=maxDo){
        message("DONE")
        return(corout)
      }
      #generate the projection
    }
  }
  
  corout
}


# Modify drop.tip from "ape". 

# generalize the "ucsctreeUse" to a "globaltree" which will be a new argument in function.
# R Studio Server

correlateTreesChar3=function(treeListObj, globaltree, char, char.ref=NULL,  maxDo=NULL, cutoff=0.000004*3, start=NULL, sqrt=F){
  maxT=treeListObj$numTrees
  if (is.null(maxDo)){
    maxDo=maxT
  }
  
  if(!is.null(char.ref)){
    cm=intersect(cm, names(char.ref))
    char.ref=char.ref[cm]
  }
  cm=intersect(names(char), globaltree$tip)
  master.tree=pruneTree(globaltree, cm)
  if(!is.null(char.ref)){
    refTree=edgeVarsDiff(master.tree, char.ref)
  }
  charTree=edgeVarsDiff(master.tree, char)
  corout=matrix(nrow=maxT, ncol=4)
  courout=as.data.frame(corout)
  rownames(corout)=names(treeListObj$trees)
  charReport=cbind(as.vector(as.numeric(cm %in% colnames(treeListObj$report[, cm]))))
  
  maxn=treeListObj$report[,cm]%*%(charReport)
  
  done=0
  todo=length(maxn>=10)
  corout[maxn<10]=0
  istart=1;
  if (!is.null(start) && ! is.numeric(start)){
    istart=match(start, names(treeListObj$trees))
  }
  else if (is.numeric(start)){
    istart=start
  }
  for (i in istart:maxT){
    message(paste("i=", i))
    if (is.na(corout[i,4])){
      tree1=treeListObj$trees[[i]]
      both=intersect(tree1$tip.label, cm)
      tree1=unroot(pruneTree(tree1,both))
      charTreeUse=unroot(pruneTree(charTree, both))
      if(!is.null(char.ref)){
        refTreeUse=unroot(pruneTree(refTree,both))
      }
      allreport=treeListObj$report[,both] 
      ss=rowSums(allreport)
      iiboth=which(ss==length(both))
      
      if(length(both)<10){
        next
      }

      torm=setdiff(treeListObj$masterTree$tip.label, both)
      
      #        allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
      #       for ( k in 1:length(iiboth)){        
      #         tmptree=rescaleTree(unroot(drop.tip(treeListObj$trees[[iiboth[k]]], torm)))         
      #         allbranch[k, ]=tmptree$edge.length           
      #       }
      
      ee=edgeIndexRelativeMaster(tree1, treeListObj$masterTree)
      ii= match(namePaths(ee,T), colnames(treeListObj$paths)) 
      
      allbranch=treeListObj$paths[iiboth,ii]
      nv=t(projection(t(allbranch), method="AVE", returnNV = T))
    
      iibad=which(allbranch<cutoff)
      allbranch=scaleMat_c(allbranch)
      if(sqrt){
        nv=sqrt(nv)
        allbranch=sqrt(allbranch)
      }
      allbranch[iibad]=NA
      
      nb=length(both)
      ai=which(maxn[iiboth, 1]==nb)
      if(cutoff>0){
        proj=naresid(allbranch[ai, ,drop=F], as.vector(nv))	### <------ Here is the problem, in the "resid" and "naresid" functions above.
      }
      else{
        proj=resid(allbranch[ai, ,drop=F], cbind(as.vector(nv)))
      }

      message("ai")
      show(length(ai))
      ee=charTreeUse$edge.length
      if(!is.null(char.ref)){
        ee=resid(rbind(ee), mod<-model.matrix(~1+refTreeUse$edge.length))
        proj=naresid(proj, refTreeUse$edge.length)
        
      }
      ee=as.vector(ee)
      
      for (j in 1:length(ai)){
     
        jj=ai[j]
        corout[iiboth[jj],4]=paste(corout[iiboth[jj],4], as.character((i)))
        if (sum(!is.na(proj[j,]))>0.5*ncol(proj)){
          # show(names(trees$trees)[iiboth[j]])
          ii2=which(!is.na(proj[j,]))
          
          cres=cor.test(proj[j,ii2], ee[ii2], method="s")
       
          corout[iiboth[jj],1:3]=c(cres$estimate, cres$statistic, cres$p.value)
          
        }
      }
      
      done=done+length(ai)
      message(paste("done", done))
      #  message(paste("Done with",done, "out of", todo))
      #  message(paste(sum(is.na(corout))," left"), appendLF = T)
      if(done>=maxDo){
        message("DONE")
        return(corout)
      }
    }
  }
  
  corout
}


#Nathan
#assuming the binTree is already in canonical form
correlateTreesBinaryDiff=function(treeObj,  binTree, usePaths=T, maxDo=NULL, species.list=NULL, useSQ=F){
  maxT=treeObj$numTrees
  if (is.null(maxDo)){
    maxDo=maxT
  }
  corout=matrix(nrow=maxT, ncol=1)
  pout=matrix(nrow=maxT, ncol=1)
  diffout=matrix(nrow=maxT, ncol=4)
  colnames(diffout)=c("pos","neg","dif","z")
  rownames(corout)=rownames(pout)=rownames(diffout)=names(treeObj$trees)
  show(binTree$tip.label )
  binReport=as.vector(as.numeric(binTree$tip.label %in% colnames(treeObj$report)))
  show((binReport))
  names(binReport)=colnames(treeObj$report)
  maxn=treeObj$report[, species.list]%*%(binReport[species.list])
  
  done=0
  todo=length(maxn>=10)
  corout[maxn<10]=0
  
  for (i in 1:maxT){
    
    if (is.na(corout[i,1])){
      tree1=treeObj$trees[[i]]
      if(! is.null(species.list)){
        tree1=pruneTree(tree1, species.list)
      } 
      both=tree1$tip.label
      bothIndex=match(both, colnames(treeObj$report))
      allreport=treeObj$report[,bothIndex] 
      ss=rowSums(allreport)
      iiboth=which(ss==length(both))
      #  torm=setdiff(treeObj$masterTree$tip.label, both)
      binTreeUse=unroot(pruneTree(binTree,tree1$tip.label))  ###
###      binTreeUse=pruneTree(binTree,tree1$tip.label)		###
      allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
      if(length(both)<10){
        next
      }
      
      if(! usePaths){  
        torm=setdiff(treeObj$masterTree$tip.label, both)
        allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
        for ( k in 1:length(iiboth)){        
          tmptree=rescaleTree(unroot(drop.tip(treeObj$trees[[iiboth[k]]], torm)))         
          allbranch[k, ]=tmptree$edge.length           
        }
      }
      else{
        ii= match(namePaths(edgeIndexRelativeMaster(tree1, treeObj$masterTree),T), colnames(treeObj$paths)) 
        ii2=match(namePaths(edgeIndexRelativeMaster(binTreeUse, treeObj$masterTree),T), colnames(treeObj$paths)) 
        stopifnot(all(ii=ii2))
        allbranch=treeObj$paths[iiboth,ii]  
        allbranch=scaleMat_c(allbranch)
      }
      nb=length(both)
      if(!useSQ){
        proj=t(projection(t(allbranch), method="AVE", returnNV = F))
      }
      else{
        proj=projectionSQ(allbranch)
      }
      ai=which(maxn[iiboth, 1]==nb)      
      
      tmp=simpleAUCmat(binTreeUse$edge.length, (proj[ai, ,drop=F]))
      corout[iiboth[ai]]=tmp$auc
      pout[iiboth[ai]]=tmp$pp

      diff = simpleDiff(binTreeUse$edge.length, (proj[ai, ,drop=F]));
	  diffout[iiboth[ai],]=diff
	  
      done=done+length(ai)
      show(done)
      #  message(paste("Done with",done, "out of", todo))
      #  message(paste(sum(is.na(corout))," left"), appendLF = T)
      if(done>=maxDo){
        message("DONE")
        return(list(r=corout, p=pout, d=diffout))
      }
    }
  }
  return(list(r=corout, p=pout, d=diffout))
}


simpleDiff=function(lab, value){
  value=rbind(value)
  posn=sum(lab>0)
  negn=sum(lab<=0)
  if(posn<2||negn<2){
    meanpos=rep(NA, nrow(value))
    meanneg=rep(NA, nrow(value))
    diff=rep(NA, nrow(value))
    z=rep(NA, nrow(value))    
  }
  else{
    meanpos=apply(value[,lab>0, drop=F],1,mean)
    meanneg=apply(value[,lab<=0, drop=F],1,mean)
    varneg=apply(value[,lab<=0, drop=F],1,var)
    diff=meanpos-meanneg
    z=diff/sqrt(varneg)
  }
  mat=matrix(c(meanpos,meanneg,diff,z),ncol=4)
  colnames(mat)=c("pos","neg","dif","z")
  return(mat)
}

library(ggplot2)

norm01=function(x){
	x2=x-min(x)
	x2/max(x2)
}

plotBinary=function(corOutput, title="Species Plot", hjust=1.2){
  show(simpleAUCmat(lab=corOutput$l1, value=corOutput$e2))
  usenames=corOutput$names
  
  tmpn=speciesNames[usenames,]
  usenames[!is.na(tmpn)]=tmpn[!is.na(tmpn)]
  
  #usenames[usenames==""]="internal"
  theme_set(theme_bw( ))
  df=list()
  theme_set(theme_bw())
  df=list()
  iis=order(corOutput$e2)
  df$x=corOutput$e2[iis]
  df$y=usenames[iis]
  df$col=as.factor(corOutput$e1[iis])
  if(hjust>0){
    ll=c(min(df$x)-0.2, max(df$x)+0.05)
  }
  else{
    ll=c(min(df$x)-0.05, max(df$x)+0.2)
  }
  df=as.data.frame(df)
  ggplot(df, aes(x = x, y=y, col=col, label=y)) + scale_size_manual(values=c(3,7))+ geom_point(aes(size=col))+
    scale_color_manual(values = c("deepskyblue3", "brown1"))+
    scale_x_continuous(limits=ll)+
    geom_text(hjust=hjust, size=5)+
    ylab("Branches")+
    xlab("relative rate")+
    ggtitle(title)+
    geom_vline(xintercept=0, linetype="dotted")+
    theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.position="none")
}

# Unadorned version of plotBinary
plotRates=function(corOutput, title="Species Plot", hjust=1.2){
  usenames=corOutput$names
  tmpn=speciesNames[usenames,]
  usenames[!is.na(tmpn)]=tmpn[!is.na(tmpn)]
  
  theme_set(theme_bw( ))

  df=list()
  df$x=corOutput$e2
  df$y=usenames
  df=as.data.frame(df)

  ggplot(df, aes(x = x, y=y,  label=y)) +  geom_point() +
    geom_text(hjust=hjust, size=5)+
    ylab("Branches")+
    xlab("relative rate")+
    ggtitle(title)+
    geom_vline(xintercept=0, linetype="dotted")+
    theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.position="none")
}


plotRatesTree=function(corOutput, title="Species Plot", hjust=1.2){
  usenames=corOutput$names
  tmpn=speciesNames[usenames,]
  usenames[!is.na(tmpn)]=tmpn[!is.na(tmpn)]
  
  theme_set(theme_bw( ))

  df=list()
  df$x=corOutput$e2
  df$y=usenames
  df=as.data.frame(df)

  plot.phylo(tmp$tree2, edge.color=col[cut(vals, nlevels)], edge.width=3.5)

}


#modified to handle binary trees
correlateTreesProj=function(treeIn1, treeIn2, treeObj, residfun=residLN, plot=F, cutoff=-1, usePaths=T, tree1Bin=F, useIndex=F, species.list=NULL){
  if(is.character(treeIn1)){
    tree1=treeObj$trees[[treeIn1]]
  }
  else{
    tree1=treeIn1
  }
  if(is.character(treeIn2)){
    tree2=treeObj$trees[[treeIn2]]
  }
  else{    
    tree2=treeIn2
  }
  
  both=intersect(tree1$tip.label, tree2$tip.label)
  if(!is.null(species.list)){
    both=intersect(both, species.list)
  }
  if(is.character(treeIn1)){
    bothIndex=which(colSums(treeObj$report[c(treeIn1, treeIn2),])==2)
  }
  
  torm=setdiff(treeObj$masterTree$tip.label, both)
  tree1=pruneTree(tree1, both)
  tree1=unroot(tree1)
  if(tree1Bin){ #fix any edges that were created through pruning
    tree1$edge.length[tree1$edge.length>1]=1
  }
  tree2=pruneTree(tree2, both)
  tree2=unroot(tree2)
  allreport=treeObj$report[,both]
  ss=rowSums(allreport)
  iiboth=which(ss==length(both))
  if (! usePaths){
    allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
    for ( i in 1:length(iiboth)){            
      tmptree=rescaleTree(drop.tip(treeObj$trees[[iiboth[i]]], torm)) 
      allbranch[i, ]=tmptree$edge.length
    }
    
  }
  else{
    if(! useIndex){
      ee=edgeIndexRelativeMaster(tree1, treeObj$masterTree)
      ii= match(namePaths(ee,T), colnames(treeObj$paths)) 
      allbranch=treeObj$paths[iiboth,ii]
      allbranch=scaleMat_c(allbranch)
      nv=projection(t(allbranch), method="AVE", returnNV = T)
      mastertree=treeObj$master
      mastertree$edge.length=nv  
      res=correlateTrees(tree1, tree2, mastertree, residfun=residfun, plot=plot, cutoff=cutoff, Tree1Bin=tree1Bin)
      
      res$nv=nv 
      res$allbranch=allbranch
    }
    else{
      allbranch=getBranch(treeObj, bothIndex)
      show(rownames(allbranch)[1])
      allbranch=scaleMat_c(allbranch)
      nv=projection(t(allbranch), method="AVE",returnNV = T)
      rr=resid(allbranch, model.matrix(~0+nv))
      rownames(rr)=rownames(allbranch)
      show(dim(rr))
      show(rownames(allbranch)[1])
      plot(rr[name1,], rr[name2,])
      res=list()
    }
  }
  
  return(res)
}
