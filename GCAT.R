#
# Load required libraries
#
library('DSI')
library('DSOpal')
library('dsBaseClient')
library('dsOmicsClient')

#
# set the connection with Opal
#
builder <- DSI::newDSLoginBuilder()
builder$append(server = "gcat", url = "https://datashield.isglobal.org/repo",
               user =  "invited", password = "123456",
               profile = "rock-inma")
logindata <- builder$build()
conns <- DSI::datashield.login(logins = logindata)

#
# Load the resources and tables
#

DSI::datashield.assign.resource(conns, "resource", "GCAT.GCATCoreSpain_V2_small")
DSI::datashield.assign.expr(conns = conns, symbol = "geno",
                            expr = as.symbol("as.resource.object(resource)"))
DSI::datashield.assign.table(conns, "pheno", "GCAT.GCATCoreSpain_V2_phenotypes_treated")



# Put genotype and phenotypes in a single object (Genomic Data Storage)

ds.GenotypeData(x='geno', covars = 'pheno', columnId = "EGA_ID", sexId = "SEX",
                male_encoding = "MALE", female_encoding = "FEMALE",
                newobj.name = 'gds.object')

ds.genoDimensions('gds.object')

ds.getChromosomeNames('gds.object')
ds.varLabels('gds.object')

# Some descriptive analyses

ds.genoDimensions('geno')

pca <- ds.PCASNPS('geno')
plotPCASNPS(pca, group = 'gds.object$CMBD_ORIGIN')
ds.scatterPlot('pheno$PC1', 'pheno$PC2', datasources = conns)



# Test HWE

hwe <- ds.exactHWE('gds.object', controls_column = 'V_4019', controls = TRUE)

# Allele frequencies


freq <- ds.alleleFrequency('gds.object')


# dsOmicsClient::ds.getSNPSbyGen('gds.object', "BRCA2", name = 'gds.brca')

# 2721 Pure hyperglyceridemia
# 7330 Age-related osteoporosis without current pathological fracture
# 4019 Essential hypertension

ds.table('pheno$V_4019')

# GWAS
ans1 <- ds.fastGWAS('gds.object', 'V_4019 ~  PC1 + PC2 + PC3 + PC4 + SEX')
ans2 <- ds.metaGWAS('gds.object', 'V_4019 ~  PC1 + PC2 + PC3 + PC4 + SEX')
manhattan(ans2$gcat)
qqplot(ans2$gcat$p.value)
LocusZoom(ans2$gcat)

#
# Polygenic risk score: https://www.pgscatalog.org/
#


ds.colnames('pheno')

# HDL cholesterol
ds.PRS(resources = 'geno', pgs_id = "PGS000660", 
       table = 'pheno', table_id_column = "EGA_ID")

# Essential hypertension
ds.PRS(resources = 'geno', pgs_id = "PGS000958", 
       table = 'pheno', table_id_column = "EGA_ID")



ds.colnames('pheno')

ds.glm('V_4019 ~ SEX + prs_PGS000660', data='pheno', family="binomial")

ds.glm('V_4019 ~ SEX + prs_PGS000958', data='pheno', family="binomial")

