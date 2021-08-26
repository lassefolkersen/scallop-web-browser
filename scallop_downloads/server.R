





library("shiny")

#gene positions
load("~/srv/scallop-web-browser/location_files/2014-07-16 gene locations.rdata")
protein_pos_file<-"~/srv/scallop-web-browser/location_files/2017-04-07_protein_pos_data.rdata"
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
    stop(safeError("The bulk data download function is unfortunately not available. However the bulk data has already been released at zenodo.org  (https://doi.org/10.5281/zenodo.2615265) and will soon be available at GWAS catalog. You can download them from there. Further, the full set of data is already browsable in the other modules implemented in this web-site, just not bulk-downloadable."))
  })
  
  
})



