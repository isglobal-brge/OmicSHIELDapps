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
builder$append(server = "gcat", url = "https://datashield.isglobal.org",
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


# Some descriptive analyses

ds.genoDimensions('geno')

pca <- ds.PCASNPS('geno')
plotPCASNPS(pca)
ds.scatterPlot('pheno$PC1', 'pheno$PC2', datasources = conns)

# Put genotype and phenotypes in a single object (Genomic Data Storage)

ds.GenotypeData(x='geno', covars = 'pheno', columnId = "EGA_ID", sexId = "SEX",
                male_encoding = "MALE", female_encoding = "FEMALE",
                newobj.name = 'gds.object')

ds.genoDimensions('gds.object')

ds.getChromosomeNames('gds.object')
ds.varLabels('gds.object')

# Test HWE

hwe <- ds.exactHWE('gds.object', controls_column = 'V_4019', controls = TRUE)

# Allele frequencies


freq <- ds.alleleFrequency('gds.object')


# dsOmicsClient::ds.getSNPSbyGen('gds.object', "BRCA2", name = 'gds.brca')

# 2721 Pure hyperglyceridemia
# 7330 Age-related osteoporosis without current pathological fracture

ds.table('pheno$V_4019')

# GWAS

ans <- ds.metaGWAS('gds.object', 'V_4019 ~  PC1 + PC2 + PC3 + PC4 + SEX')
manhattan(ans$gcat)
qqplot(ans$gcat$p.value)
LocusZoom(ans$gcat)

#
# Polygenic risk score: https://www.pgscatalog.org/
#

# HDL cholesterol
ds.PRS(resources = 'resource', pgs_id = "PGS000660", 
       table = 'pheno', table_id_column = "EGA_ID")

ds.colnames('pheno')

ds.glm('V_4019 ~ SEX', data='pheno', family="gaussian")

datashield.logout(conns)


