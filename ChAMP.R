#Carregando o pacote
library("ChAMP")

#Subindo a sample sheet e os chips
testDir <- "C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\Camila"

#Realizando os filtros das pobres
myLoad <- champ.load(testDir,arraytype="EPIC")


myLoad$pd

champ.QC <- function(beta = myLoad$beta,
                     pheno=myLoad$pd$Sample_Group,
                     mdsPlot=TRUE,
                     densityPlot=TRUE,
                     dendrogram=TRUE,
                     PDFplot=TRUE,
                     Rplot=TRUE,
                     Feature.sel="None",
                     resultsDir="./CHAMP_QCimages/")
  
  champ.QC()

#Normalizando os dados
myNorm <- champ.norm(beta=myLoad$beta, arraytype="EPIC")


# Corrigir os efeitos de batch usando ComBat (corrigindo Slide e Array)
pheno = myLoad$pd

# Corrigir os efeitos de batch para o Slide
myCombat <- ComBat(dat = myNorm, batch = pheno$Slide, mod = NULL)

# Corrigir os efeitos de batch para o Array
myCombat <- ComBat(dat = myCombat, batch = pheno$Array, mod = NULL)

# Verificar a correção dos dados
champ.QC(beta = myCombat, pheno = pheno$Sample_Group, dendrogram = FALSE)


# Rodar o SVD após a correção de batch
champ.SVD(beta = as.data.frame(myCombat), pd = pheno)


champ.QC(myCombat) 


#Ajustando pela idade pq a idade não é nossa variavel de interesse, idade pode ser um vies por mudar a metilação
age <- champ.runCombat(beta=myNorm,
                       pd=myLoad$pd,
                       variablename="Age",
                       batchname=c("Slide"),
                       logitTrans=TRUE)

#Ajustando pelo tipo celular
myRefBase <- champ.refbase(beta = age, arraytype = "EPIC")
head (myRefBase)

myNorm <- myRefBase

#Importando os dados 
write.csv(myRefBase, file = "C:\\Users\\guihi\\OneDrive\\Área de Trabalho\\Resultados Camila\\DadosNoramalizados_Idade_Fracao.csv")

#Calculando as DMPs

DMP<- champ.DMP(beta = myRefBase$CorrectedBeta,
                pheno = myLoad$pd$Sample_Group,
                compare.group = c("I","II"),
                adjPVal = 0.05,
                adjust.method = NULL,
                arraytype = "EPIC")

#Calculando os DMRs

myDMR <- champ.DMR(beta = myRefBase$CorrectedBeta,
                   pheno = myLoad$pd$Sample_Group,
                   compare.group = c("I","II")
                   ,method="Bumphunter")


DMR.GUI(DMR=myDMR, beta = myRefBase$CorrectedBeta,
        pheno = myLoad$pd$Sample_Group,
        compare.group = c("I","II"))
