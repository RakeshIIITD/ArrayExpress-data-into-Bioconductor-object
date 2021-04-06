# ---------------------------- Driver script --------------------------------

#--------------------------Libraries-----------------------------------------
require(affy)
require(Biobase)
require(limma)
require(oligo)
require(CCl4)
require(stringr)
require("XML")
require(RCurl)

#---------------------------------------------------------
#function to determine type of accession name and also specify the working directory
detectType =  function(accessionName,filePath = getwd()) 
{
  FileName = paste(accessionName,".sdrf", ".txt", sep = "")
  checkFileExist = paste(getwd(),FileName,sep = "/")
  #print(checkFileExist) : "/home/kd/Documents/Bdmh_Proj_files/E-ATMX-18.sdrf.txt"
  if(file.exists(checkFileExist)){
    #if(file.exits(paste(accessionName,".sdrf",".txt",sep = "")))
    #print("Entered")
    # sdrfFile = tryCatch({read.csv(paste(accessionName,"",".txt",sep = "") , header = TRUE , sep = "\t"  , row.names = 1 )} , warning = function(e) { return (NA)})
    
    sdrfFile = read.csv(paste(accessionName,".sdrf",".txt",sep = "") , header = TRUE , sep = "\t"  , row.names = 1 )
    
    ## check whether the file experiment is :
    ## (1) Affymetrix
    ## (2) One Color 
    ## (3) Two Color
    
    if(any(unique(sdrfFile$Label) %in% "biotin"))
    {
      return ("affy")
    }
    
    else if(length(unique(sdrfFile$Label))==2)
    {
      return ("twoColor")
    }
    else
    {
      return ("oneColor")
    }
    
  }
}

#--------------------------------------------------------------------------------
#detectType("E-ATMX-18")
#detectType("/home/rakesh/Desktop/E-ATMX-18/E-ATMX-18")

#--------------------------Extract raw data---------------------------------------

## function to create expressionSet object of one color experimen

## function to create AffyBatch object

extractDataFiles = function(fileList , formatType = "affy")
{
  extractedFiles = c()
  fileType = NULL
  zippedCount = 0 
  fileListLength = length(fileList)
  # 5 for E-ATMX-18
  for(i in seq( 1 , fileListLength ))
  {
    ## if file contains raw zipped files
    print(fileList[i])
    #if(grepl("raw[[:punct:]]zip",fileList[i])) # regex expression to match only raw file
    
    if(grepl("zip",fileList[i]))
    {
      zippedCount = zippedCount + 1
      #print(zippedCount)
      ## append file names in a single character vector from various raw zipped files
      
      files =  unzip( fileList[i] , exdir = paste(getwd(),"/extracted",sep = ""))
      #print(files)
      if(formatType %in% "affy" )
      {
        currentType = detectAffyFormat(files)
      }
      else
      {
        ## if format is not affy then detect file extension
        
        format = unlist(strsplit(files[1], split  = "\\."))
        currentType = format[2]
      }
      
      if(is.null(fileType)||fileType==currentType)
      {
        # if file type is same
        fileType = currentType
      }
      else
      {
        # else throw an error
        warning("File type in  all zip archives doesn't appear to be same. File formats in raw zipped file should be same ")
      }
      
      extractedFiles = append(extractedFiles , files)
      
    }
  }
  return (list(extractedFiles = extractedFiles , type = fileType))
}




##--------------- A function to detect the type of files-------------------------

detectAffyFormat = function(extractedFiles)
{
  affyType = "other"
  
  if(length( grep(".gpr", extractedFiles , ignore.case = TRUE ) )==length(extractedFiles))
  {
    affyType = "gpr" 
  }
  
  if(length( grep(".CEL", extractedFiles , ignore.case = TRUE ) )==length(extractedFiles))
  {
    affyType = "cel" 
  }
  
  
  return (affyType)
  ##-------------------giving other categories ------------------------------------
  ## ------------------------ TO-DO : Can add more categories----------------
}  

##Create affybatch object
createAffy = function(fileList, accessionName,columnName)
{
  ## extract from zipped format
  
  ## create affy either from CEL or gpr or from various file formats
  
  extractedFilesAndType = extractDataFiles(fileList)
  
  rawDataFiles = extractedFilesAndType$extractedFiles
  
  
  AffyType = extractedFilesAndType$type
  
  # abatch object
  abatch = NULL
  
  if(AffyType %in% "cel")
  {
    ##  TO-DO : logic to convert to .cel files
    
    
    FileName = paste(accessionName,".sdrf", ".txt", sep = "")
    print(FileName)
    checkFileExist = paste(getwd(),FileName,sep = "/")
    print(checkFileExist) # "/home/kd/Documents/Bdmh_Proj_files/E-AFMX-1.sdrf.txt"
    
   
    
    if(file.exists(FileName)){
      
      print(checkFileExist)
      sdrf = read.AnnotatedDataFrame(checkFileExist, row.names = NULL , fill = TRUE 
                                     , blank.lines.skip = TRUE , varMetadata.char = "$" , quote = "\"")
      
      # sdrf = read.AnnotatedDataFrame('/home/rakesh/Desktop/Data/E-MEXP-1422/E-MEXP-1422.sdrf.txt'
      #                                ,row.names = NULL , fill = TRUE 
      #                                , blank.lines.skip = TRUE , varMetadata.char = "$" , quote = "\"")
      # 
      #print(file.path())
      #setwd(paste(getwd(),"/extracted",sep = ""))
      #print(getwd())
      
      f  = list.files(paste(getwd(),"/extracted",sep = ""),full.names = TRUE)
      #print(f)
      
      #print("1")
      print("Constructing affybatch....")
      abatch = read.affybatch(f,phenoData = sdrf)
      print("Object successfully created!")
      #abatch = read.affybatch(f, phenoData = sdrf , description = exData)
    }
    else {
      #Show warning
      warning("File does not exists")
    }
  }
  
  if(AffyType %in% "gpr")
  {
    
    #gpr = read.maimages( files = "/home/rakesh/Desktop/AffyMetrix/extracted/2010-01-15_a11-993a_0532.gpr",source = "generic" ,columns = list(G = "GenePix:F532 Median",R = "GenePix:F532 Median - B532" ,Rb = "GenePix:B532 CV" , Gb = "GenePix:B532" ))
    gpr = read.maimages( files = f ,source = "generic" ,columns = columnName)
  }
  
  
  if(AffyType %in% "other")
  {
    #Other
  }
  ## first extract zipped files in form of folders
  return (abatch)
}

createOneColor = function()
{
  
}

createTwoColor = function(fileList,accessionName,columnName)
{
  # Need to prepare sdrf for n-channel set
  #phenoDataNChannelSet = read.AnnotatedDataFrame("E-ATMX-18.sdrf.txt", row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\"")
  # 
  FileName = paste(accessionName,".sdrf", ".txt", sep = "")
  #print(FileName)
  checkFileExist = paste(getwd(),FileName,sep = "/")
  #print(checkFileExist) # "/home/kd/Documents/Bdmh_Proj_files/E-ATMX-18.sdrf.txt"
  
  if(file.exists(FileName)){
    #print(list.files())
    
    print(checkFileExist)
    #Define sdrf file of two color experiment here
    sdrf = read.AnnotatedDataFrame(checkFileExist, row.names = NULL , fill = TRUE 
                                   , blank.lines.skip = TRUE , varMetadata.char = "$" , quote = "\"")
    
    print("sdrf read")
    extractedFilesAndType = extractDataFiles(fileList,"")
    rawDataFiles = extractedFilesAndType$extractedFiles
    
    # columns = list(G = "GenePix:F532 Median",R = "GenePix:F532 Median - B532" ,Rb = "GenePix:B532 CV" ,Gb = "GenePix:B532" )
    
    #assayDataFile = read.maimages( files = rawDataFiles, source = "generic" ,columns = NULL)
    #assayDataFile = read.maimages( files = rawDataFiles, source = "generic" ,columns = list(G = "GenePix:F532 Median",R = "GenePix:F532 Median - B532" ,Rb = "GenePix:F532 Mean - B532" ,Gb = "GenePix:B532 Mean" ))
    assayDataFile = read.maimages( files = rawDataFiles, source = "generic" ,columns = columnName)
    print(class(assayDataFile))
    assay = NULL
    if(class(assayDataFile) == "RGList" | class(assayDataFile) == "EListRaw"){
      #construct nchannelset
      if(class(assayDataFile) == "RGList"){
        if("Rb" %in% names(assayDataFile))
          assay = with(assayDataFile, assayDataNew(R = R, G = G, Rb = Rb, Gb = Gb)) #will not work if datacolumns where user specified
        else
          assay = with(assayDataFile, assayDataNew(G = G, R = R))
      }
    }
    cat("assayclass",class(assay))
    # # k = files with extracted data  //Done
    #featureData = new("AnnotatedDataFrame", data = assay$targets)
    
    mfeatureData = new("AnnotatedDataFrame", data = assayDataFile$targets )
    #mfeatureData = new("AnnotatedDataFrame", data = assayDataFile$targets)
    print(class(mfeatureData))
    
    #nset = new("NChannelSet",  )
    nset = new("NChannelSet", assayData = assay)
    
    
    nset@phenoData = sdrf
    nset@featureData = mfeatureData
    print("Saving NchannelSet Object......")
    save(nset, file = file.path(getwd(), "Data.r"))
  }
}

#---------------------------------------------------------------------------

# 2 parameters are passed to function, accessionName specified by user and path, if not specified then default be used
getFilesFromArrayExpress = function(accessionName, workingDirectory = getwd())
{
  #Reference : https://www.ebi.ac.uk/arrayexpress/help/programmatic_access.html
  baseUrl = "https://www.ebi.ac.uk/arrayexpress/xml/v3/files" ## Baseurl for files download
  xmlUrl = paste( baseUrl, accessionName, sep = "/") #String concatenation where accessionName is appended to baseurl
  #print(xmlUrl) #"https://www.ebi.ac.uk/arrayexpress/xml/v3/files/E-ATMX-18"
  #if(url.exists(xmlURL) == FALSE) stop("Url does not exists, Wrong accession name.....")
  xml = XML::xmlParse(suppressWarnings(readLines(xmlUrl)))
  #Used to read lines of xml from website.
  #Suppress warnings are used not to display warnings at user level.
  #print(xml)
  idfUrl = XML::xpathSApply(xml, "//*[kind='idf']/url", xmlValue) #Search for idf file url in xml obtained, xmlValue used to return R object instead of XpressNodeSet
  idfName = XML::xpathSApply(xml, "//*[kind='idf']/name", xmlValue) #Search for file name in xml
  #print(idfName) #"E-ATMX-18.idf.txt"
  #print(class(idfName)) #"character"
  #print(class(idfUrl)) #"character"
  #print(idfUrl) #"https://www.ebi.ac.uk/arrayexpress/files/E-ATMX-18/E-ATMX-18.idf.txt"
  DestinationIdfFile = paste(workingDirectory,idfName,sep = "/")
  #print(DestinationIdfFile) #"/home/kd/E-ATMX-18.idf.txt"
  #print(class(DestinationIdfFile)) #"character"
  #Below code checks for file existence in specified diretory. If available then prints a message else dpwnload it from internet
  if(!file.exists(DestinationIdfFile)){
    try(download.file(idfUrl,destfile = paste(workingDirectory,idfName,sep = "/"), mode = "wb")) #Used to download idf file from website
  }
  else warning(paste(idfName,"already exists in directory",sep = " "))
  
  #Used to download eSet file from website for the accession
  ESetRUrl = XML::xpathSApply(xml, "//*[kind='r-object']/url", xmlValue) #Search for esetR file url in xml obtained, xmlValue used to return R object instead of XpressNodeSet
  ESetRName = XML::xpathSApply(xml, "//*[kind='r-object']/name", xmlValue) #Search for file name in xml
  DestinationESetFile = paste(workingDirectory,ESetRName,sep = "/")
  if(!file.exists(DestinationESetFile)){
    try(download.file(ESetRUrl,destfile = DestinationESetFile, mode = "wb")) #Used to download idf file from website
  }
  else warning(paste(ESetRName,"already exists in directory",sep = " "))
  
  #Used to download idf file from website
  RawUrl = XML::xpathSApply(xml, "//*[kind='raw']/url", xmlValue) #Search for esetR file url in xml obtained, xmlValue used to return R object instead of XpressNodeSet
  RawName = XML::xpathSApply(xml, "//*[kind='raw']/name", xmlValue) #Search for file name in xml
  DestinationRawFile = paste(workingDirectory,RawName,sep = "/")
  if(!file.exists(DestinationRawFile)){
    try(download.file(RawUrl,destfile = DestinationRawFile, mode = "wb")) #Used to download idf file from website
  }
  else warning(paste(RawName,"already exists in directory",sep = " "))
  
  #Download sdrf file from website
  SdrfUrl = XML::xpathSApply(xml, "//*[kind='sdrf']/url", xmlValue) #Search for esetR file url in xml obtained, xmlValue used to return R object instead of XpressNodeSet
  SdrfName = XML::xpathSApply(xml, "//*[kind='sdrf']/name", xmlValue) #Search for file name in xml
  DestinationSdrfFile = paste(workingDirectory,SdrfName,sep = "/")
  if(!file.exists(DestinationSdrfFile)){
    try(download.file(SdrfUrl,destfile = DestinationSdrfFile, mode = "wb")) #Used to download idf file from website
  }
  else warning(paste(SdrfName,"already exists in directory",sep = " "))
  
  #Download adf file from website
  AdfUrl = XML::xpathSApply(xml, "//*[kind='adf' and extension='txt']/url", xmlValue) #Search for esetR file url in xml obtained, xmlValue used to return R object instead of XpressNodeSet
  AdfName = XML::xpathSApply(xml, "//*[kind='adf' and extension='txt']/name", xmlValue) #Search for file name in xml
  DestinationAdfFile = paste(workingDirectory,AdfName,sep = "/")
  if(!file.exists(DestinationAdfFile)){
    try(download.file(AdfUrl,destfile = DestinationAdfFile, mode = "wb")) #Used to download idf file from website
  }
  else warning(paste(AdfName,"already exists in directory",sep = " "))
  #Warnings are printed if file already exists in the directory
  # A vector of all destination file names are returned by creating a character vector.
  #Since last line of function is returned by default, so defined in this way.
  #c(DestinationAdfFile,DestinationESetFile,DestinationRawFile,DestinationIdfFile,DestinationSdrfFile)
  c(AdfName,ESetRName,RawName,idfName,SdrfName)
}
