loadData <- function(fname) {
  nmc <- 100
  
  dhatMaxariList <- dhatBlock1List <- dhatBlock2List <- dhatBlock3List <- dhatBlock4List <- NULL
  dhatUsvtList <- dhatSvt1List <- dhatSvt2List <- dhatSvt3List <- NULL
  dhatSum1List <- dhatSum2List <- dhatSum3List <- dhatSum4List <- NULL
  
  ariTrueDList <- ariMaxariList <- ariBlock1List <- ariBlock2List <- ariBlock3List <- ariBlock4List <- NULL
  ariUsvtList <- ariSvt1List <- ariSvt2List <- ariSvt3List <- NULL
  ariSum1List <- ariSum2List <- ariSum3List <- ariSum4List <- NULL
  ariBlock1sList <- ariBlock2sList <- ariBlock3sList <- ariBlock4sList <- NULL
  
  khatTrueDList <- khatMaxariList <- khatBlock1List <- khatBlock2List <- khatBlock3List <- khatBlock4List <- NULL
  khatUsvtList <- khatSvt1List <- khatSvt2List <- khatSvt3List <- NULL
  khatSum1List <- khatSum2List <- khatSum3List <- khatSum4List <- NULL
  
  for (mc in 0 : (nmc - 1)) {
#    filename = paste(name, mc, ".RData", sep="")
#    load(filename)
    file <- drive_download(fname[mc+1,], overwrite = TRUE, verbose = FALSE)
    load(file$local_path)
    
    dhatMaxariList[[mc+1]] <- dhatMaxari
    dhatBlock1List[[mc+1]] <- dhatBlock1
    dhatBlock2List[[mc+1]] <- dhatBlock2
    dhatBlock3List[[mc+1]] <- dhatBlock3
    dhatBlock4List[[mc+1]] <- dhatBlock4
    dhatSvt1List[[mc+1]] <- dhatSvt1
    dhatSvt2List[[mc+1]] <- dhatSvt2
    dhatSvt3List[[mc+1]] <- dhatSvt3
    #dhatUsvtList[[mc+1]] <- dhatUsvt
    dhatSum1List[[mc+1]] <- dhatSum1
    dhatSum2List[[mc+1]] <- dhatSum2
    dhatSum3List[[mc+1]] <- dhatSum3
    dhatSum4List[[mc+1]] <- dhatSum4
    
    #ariTrueDList[[mc+1]] <- ariTrueD
    ariMaxariList[[mc+1]] <- ariMaxari
    ariBlock1List[[mc+1]] <- ariBlock1
    ariBlock2List[[mc+1]] <- ariBlock2
    ariBlock3List[[mc+1]] <- ariBlock3
    ariBlock4List[[mc+1]] <- ariBlock4
    #ariUsvtList[[mc+1]] <- ariUsvt
    ariSvt1List[[mc+1]] <- ariSvt1
    ariSvt2List[[mc+1]] <- ariSvt2
    ariSvt3List[[mc+1]] <- ariSvt3
    ariSum1List[[mc+1]] <- ariSum1
    ariSum2List[[mc+1]] <- ariSum2
    ariSum3List[[mc+1]] <- ariSum3
    ariSum4List[[mc+1]] <- ariSum4
    ariBlock1sList[[mc+1]] <- ariBlock1s
    ariBlock2sList[[mc+1]] <- ariBlock2s
    ariBlock3sList[[mc+1]] <- ariBlock3s
    ariBlock4sList[[mc+1]] <- ariBlock4s
    
    #khatTrueDList[[mc+1]] <- khatTrueD
    khatMaxariList[[mc+1]] <- khatMaxari
    khatBlock1List[[mc+1]] <- khatBlock1
    khatBlock2List[[mc+1]] <- khatBlock2
    khatBlock3List[[mc+1]] <- khatBlock3
    khatBlock4List[[mc+1]] <- khatBlock4
    #khatUsvtList[[mc+1]] <- khatUsvt
    khatSvt1List[[mc+1]] <- khatSvt1
    khatSvt2List[[mc+1]] <- khatSvt2
    khatSvt3List[[mc+1]] <- khatSvt3
    khatSum1List[[mc+1]] <- khatSum1
    khatSum2List[[mc+1]] <- khatSum2
    khatSum3List[[mc+1]] <- khatSum3
    khatSum4List[[mc+1]] <- khatSum4
    
    system(paste("rm ", file$local_path))
  }
  
  return(list(dhatMaxariList = dhatMaxariList,
              dhatBlock1List = dhatBlock1List,
              dhatBlock2List = dhatBlock2List,
              dhatBlock3List = dhatBlock3List,
              dhatBlock4List = dhatBlock4List,
              dhatSvt1List = dhatSvt1List,
              dhatSvt2List = dhatSvt2List,
              dhatSvt3List = dhatSvt3List,
              #dhatUsvtList = dhatUsvtList,
              dhatSum1List = dhatSum1List,
              dhatSum2List = dhatSum2List,
              dhatSum3List = dhatSum3List,
              dhatSum4List = dhatSum4List,
              
              #ariTrueDList = ariTrueDList,
              ariMaxariList = ariMaxariList,
              ariBlock1List = ariBlock1List,
              ariBlock2List = ariBlock2List,
              ariBlock3List = ariBlock3List,
              ariBlock4List = ariBlock4List,
              ariSvt1List = ariSvt1List,
              ariSvt2List = ariSvt2List,
              ariSvt3List = ariSvt3List,
              #ariUsvtList = ariUsvtList,
              ariSum1List = ariSum1List,
              ariSum2List = ariSum2List,
              ariSum3List = ariSum3List,
              ariSum4List = ariSum4List,
              ariBlock1sList = ariBlock1sList,
              ariBlock2sList = ariBlock2sList,
              ariBlock3sList = ariBlock3sList,
              ariBlock4sList = ariBlock4sList,
              
              khatTrueDList = khatTrueDList,
              khatMaxariList = khatMaxariList,
              khatBlock1List = khatBlock1List,
              khatBlock2List = khatBlock2List,
              khatBlock3List = khatBlock3List,
              khatBlock4List = khatBlock4List,
              #khatUsvtList = khatUsvtList,
              khatSvt1List = khatSvt1List,
              khatSvt2List = khatSvt2List,
              khatSvt3List = khatSvt3List,
              khatSum1List = khatSum1List,
              khatSum2List = khatSum2List,
              khatSum3List = khatSum3List,
              khatSum4List = khatSum4List))
}