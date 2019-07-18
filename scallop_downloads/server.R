





library("shiny")

#gene positions
load("~/srv/scallop-web-browser/2014-07-16 gene locations.rdata")
protein_pos_file<-"~/srv/scallop-web-browser/2017-04-07_protein_pos_data.rdata"
load(protein_pos_file)
p<-data.frame(
  row.names=data[,"hgnc_symbol"],
  pheno_id=data[,"No_in_GWAS_files"],
  CHR = data[,"trait_chr"],
  BP = data[,"trait_pos"],
  stringsAsFactors=F)
p<-p[order(rownames(p)),]
phenotypes_vector<-p[,"pheno_id"]
names(phenotypes_vector) <- rownames(p)






shinyServer(function(input, output) {
  
  
  
  output$text1 <- renderText({
    stop(safeError("The bulk data download function is unfortunately not available before publication, neither here nor at the zenodo.org location (DOI:10.5281/zenodo.2615265). This embargo is motivated because the review may change the final format of the sumstats and we do not wish for different versions to be used or implemented elsewhere. However, the full set of data is already browsable in the other modules implemented in this web-site, they are just not bulk-downloadable."))
  })
  
  
})



