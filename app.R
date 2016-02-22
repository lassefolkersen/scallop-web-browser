


#Generating random naming scheme until we figure out the actual names
cl<-c(247197891, 242713278, 199439629, 191246650, 180727832,
			170735623, 158630410, 146252219, 140191642, 135347681,
			134361903, 132289533, 114110907, 106354309, 100334282,
			88771793,  78646005,  76106388,  63802660,  62429769,
			46935585, 49396972, 154908521,  57767721)
p<-data.frame(
	row.names=paste("Phenotype",as.character(1:83)), 
	pheno_id=as.character(1:83), 
	CHR = sample(c(1,1,1,1,2,2,2,2,3,3,3,4,4,5,5,5,6,6,7,7,8:22),83,replace=T),
	stringsAsFactors=F)
for(i in 1:nrow(p)){
	p[i,"BP"]<-runif(min=100,max=cl[p[i,"CHR"]],1)
}


phenotypes_vector<-p[,"pheno_id"]
names(phenotypes_vector) <- rownames(p)

load("/home/ubuntu/misc/2014-07-16 gene locations.rdata")

server <- function(input, output) {
	
	output$mainPlot <- renderPlot({ 
		# Take a dependency on input$goButton
		
		if(input$goButton == 0){
			return("")
		}else{
			
			##################################
			#input-variables, log and register	
			##################################
			email <- isolate(input$email)
			type <- isolate(input$type)
			gene <- isolate(input$gene)
			distance <- isolate(input$distance)
			p_value_cutoff <- isolate(input$p_value_cutoff)
			top_label_count<-isolate(input$top_label_count)
			phenotype <- isolate(input$phenotype)

			if(!tolower(email) %in% c("lassefolkersen@gmail.com","daniel.ziemek@pfizer.com","anders.malarstig@pfizer.com") ){
				m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"plot",email)
				m<-paste(m,collapse="\t")
				write(m,file="/home/ubuntu/logs/illegal_access_log.txt",append=TRUE)
				Sys.sleep(2)
				stop("In the test-phase non-privileged users are not allowed")
			}else{
				
				m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"plot",email,gene,phenotype)
				m<-paste(m,collapse="\t")
				write(m,file="/home/ubuntu/logs/log.txt",append=TRUE)
			}
			
			
			
			##################################
			#DNA-aspect --- starting at a 
			#genetic location
			##################################
			if(type == "dna"){
				in_dir<-"~/data/2016-02-19_splits/"
				window<-1000000
				if(!gene%in%rownames(geneLocations)){stop(paste(gene,"not found"))}
				chr<-sub("^chr","",geneLocations[gene,"chr_name"])
				start<-geneLocations[gene,"start"] - distance
				end<-geneLocations[gene,"end"] + distance
				p1<-floor(start/window)
				p2<-floor(end/window)
				filename1<-paste(in_dir,"split_pheno_",phenotype,"_chr",chr,"_",p1,"_",p1+1,".txt.gz",sep="")
				if(!file.exists(filename1))stop(paste("Could not find file",filename1))
				dh<-read.table(filename1,sep="\t",header=T,stringsAsFactors=F)
				
				#in case we are at a window breakpoint
				if(p1!=p2){
					dh1<-dh
					filename2<-paste(in_dir,"split_pheno_",phenotype,"_chr",chr,"_",p2,"_",p2+1,".txt.gz",sep="")
					dh2<-read.table(filename2,sep="\t",header=T,stringsAsFactors=F)
					dh<-rbind(dh1,dh2)
				}
				
				GWAS<-dh[dh[,"BP"]>start & dh[,"BP"]<end,]
				if(nrow(GWAS)==0){stop(paste("No SNPs found around gene",gene))}
				
				GWAS[,"-log10(P)"] <- -log10(GWAS[,"P" ])
				rownames(GWAS)<-GWAS[,"SNP"]
				ylim<-c(0,max(c(5,GWAS[,"-log10(P)"])))
				xlim<-c(start,end)
				plot(NULL,
						 xlim=xlim,
						 ylim=ylim,
						 xlab=paste("chr",chr),
						 ylab="-log10(P)",
						 main=paste("SNPs around",gene)
				)
				
				points(
					x=GWAS[,"BP"],
					y=GWAS[,"-log10(P)"],
					pch=19,
					col="dodgerblue"
				)
				
				#highlight top-3 SNPs
				tooClose<-vector()
				tooCloseDist<-4000
				count<-min(c(nrow(GWAS),top_label_count))
				for(i in 1:count){
					snps<-rownames(GWAS)[order(GWAS[,"-log10(P)"],decreasing=T)]
					if(sum(!snps%in%tooClose)==0)break
					snp<-snps[!snps%in%tooClose][1]
					text(x=GWAS[snp,"BP"],y=GWAS[snp,"-log10(P)"],label=snp,adj=0,cex=0.9)
					tooCloseHere<-GWAS[snp,"BP"] - tooCloseDist < GWAS[,"BP"] & GWAS[snp,"BP"] + tooCloseDist > GWAS[,"BP"]
					tooClose<-c(tooClose,rownames(GWAS)[tooCloseHere])
					
				}
				lines(x=c(geneLocations[gene,"start"],geneLocations[gene,"end"]),y=c(0,0),lwd=4,col="black")
				
				
				##################################
				#Protein-aspect --- starting with 
				#a specific plasma protein
				##################################
				
			}else if(type == "protein"){
				in_dir <- "~/data/2016-02-20_significant_bits/"
				file_name<-paste(in_dir,"trimmed_pheno_",phenotype,".txt.gz",sep="")
				data<-read.table(file_name, header=T,sep="\t",stringsAsFactors=FALSE)
				data<-data[data[,"P"] < 10^-p_value_cutoff,]
				data[,"neglogp"] <- -log10(data[,"P"])
				
				
				data[,"trait_chr"]<-p[p[,"pheno_id"]%in%phenotype,"CHR"]
				data[,"trait_pos"]<-p[p[,"pheno_id"]%in%phenotype,"BP"]
				
				
				#calculate the "universal-BP" --- i.e. the BP-sum on chr + all-other-chr before it (=were on the circle to plot)
				#from here https://support.bioconductor.org/p/14766/
				cl<-c(247197891, 242713278, 199439629, 191246650, 180727832,
							170735623, 158630410, 146252219, 140191642, 135347681,
							134361903, 132289533, 114110907, 106354309, 100334282,
							88771793,  78646005,  76106388,  63802660,  62429769,
							46935585, 49396972, 154908521,  57767721)
				chrLengths<-data.frame(row.names=c(1:22,"X","Y"),sum=cl,number=1:length(cl),endsum=cumsum(cl))
				chrLengths[,"startsum"]<- chrLengths[,"endsum"] - chrLengths[,"sum"]
				data[,"snp_abs_pos"]<-chrLengths[data[,"CHR"],"startsum"] + data[,"BP"]
				data[,"trait_abs_pos"]<-chrLengths[data[,"trait_chr"],"startsum"] + data[,"trait_pos"]
				
				
				#This block transforms all the (chr-pos, -logP) positions from cartesian coordinates to polar coordinates (=this is the circularizing step)
				p_range<-1
				maxPos<-sum(cl)#max(data[,"snp_abs_pos"],na.rm=T)
				base<-2
				for(i in 1:nrow(data)){
					data[i,"r"] <- max(c(0,log10(data[i,"neglogp"]) / p_range))
					data[i,"theta"] <- (data[i,"snp_abs_pos"] / maxPos) * 2 * pi
					data[i,"x"] = (base+data[i,"r"]) * cos( data[i,"theta"] )
					data[i,"y"] = (base+data[i,"r"]) * sin( data[i,"theta"] )
					data[i,"x0"] = (base) * cos( data[i,"theta"] )
					data[i,"y0"] = (base) * sin( data[i,"theta"] )
				}
				data<-data[!is.na(data[,"x"]),]
				data<-data[order(data[,"snp_abs_pos"]),]
				data[,"colours"]<-rainbow(nrow(data), s = 1, v = 1, start = 0, end = 1, alpha = 1)
				
				
				#Open up the actual plotting mechanism with an empty canvas
				plot(NULL,xlim=c(-4,4),ylim=c(-4,4),ylab="",xlab="",xaxt="n",yaxt="n")
				
				#Drawing the lines for the manhattan plot (using polar coordinates from above)
				for(i in 1:nrow(data)){
					lines(
						x=c(data[i,"x"],data[i,"x0"]),
						y=c(data[i,"y"],data[i,"y0"]),
						col=data[i,"colours"])
				}
				
				
				#drawing the frame-work circles (we want a P 0, 10^-8 and a 10^-15 circe
				ch<-data.frame(offset=c(0,log10(15),log10(8)),lty=c(1,2,1),col=c("black","grey70","grey50"),stringsAsFactors=F)
				for(k in 1:nrow(ch)){
					res<-50
					b<-seq(0,2*pi,length.out=res)
					for(i in 1:(res-1)){
						x1 = (ch[k,"offset"] + base) * cos( b[i] )
						x2 = (ch[k,"offset"] + base) * cos( b[i+1] )
						y1 = (ch[k,"offset"] + base) * sin( b[i] )
						y2 = (ch[k,"offset"] + base) * sin( b[i+1] )
						lines(x=c(x1,x2),y=c(y1,y2),lty=ch[k,"lty"],col=ch[k,"col"])
					}
				}
				
				#drawing the internal lines
				for(j in 1:nrow(data)){
					orig<-data[j,"snp_abs_pos"]
					dest<-data[j,"trait_abs_pos"]
					colour<-rainbow(nrow(data), s = 1, v = 1, start = 0, end = 1, alpha = 1)[sum(data[,"snp_abs_pos"]<orig)+1]
					x1 = (base) * cos( 2*pi*(orig/maxPos) )
					y1 = (base) * sin( 2*pi*(orig/maxPos) )
					x2 = (base) * cos( 2*pi*(dest/maxPos) )
					y2 = (base) * sin( 2*pi*(dest/maxPos) )
					lines(x=c(x1,x2),y=c(y1,y2),col=colour,lwd=1)
					print(paste(x1,y1,x2,y2))
				}
				
				
				
				#drawing the chromosome sizes
				for(i in 1:nrow(chrLengths)){
					start<-chrLengths[i,"startsum"]
					# end<-chrLengths[i,"cumsum"]
					x1 = (base) * cos( 2*pi*(start/maxPos) )
					x2 = (base+0.1) * cos( 2*pi*(start/maxPos) )
					y1 = (base) * sin( 2*pi*(start/maxPos)  )
					y2 = (base+0.1) * sin( 2*pi*(start/maxPos)  )
					
					lines(x=c(x1,x2),y=c(y1,y2),lty=1,col="black")
					text(x1,y1,label=rownames(chrLengths)[i],cex=0.5,adj=-0.2,col="grey50")
				}
				

				#labelling top-X hits
				tooClose<-vector()
				tooCloseDist<-20000000
				count<-min(c(nrow(data),top_label_count))
				for(i in 1:count){
					snps<-data[order(data[,"neglogp"],decreasing=T),"SNP"]
					if(sum(!snps%in%tooClose)==0)break
					snp<-snps[!snps%in%tooClose][1]
					x<-data[data[,"SNP"]%in%snp,"x"]
					y<-data[data[,"SNP"]%in%snp,"y"]
					text(x=x,y=y,label=snp,adj=0,cex=0.8)
					pos<-data[data[,"SNP"]%in%snp,"snp_abs_pos"]
					tooCloseHere<-data[abs(data[,"snp_abs_pos"] - pos) <tooCloseDist,"SNP"]
					tooClose<-c(tooClose,tooCloseHere)
				}
			}else if(type == "download"){
				return(NULL)
				
			}else{ 
				stop("Aspect not implemented yet")
			}
			
		}
	})
}




ui <- fluidPage(
	titlePanel("Improve OLINK browser"),
	sidebarLayout(
		sidebarPanel(
			
			radioButtons("type", "Focus", 
									 c("Cis-effects"="dna",
									 	"Trans-effects" = "protein",
									 	"Data download"="download"), 
									 selected = "dna", inline = TRUE),
			
			textInput(inputId="email", label = "E-mail", value = ""),
			conditionalPanel(
				condition = "input.type == 'dna'",
				textInput(inputId="gene", label = "Gene", value = "")
			),
			# conditionalPanel(
				# condition = "input.type == 'protein'",
				selectInput("Protein", "Phenotype", choices = phenotypes_vector),
			# ),
			checkboxInput("advanced", "Advanced options", value = FALSE),
			conditionalPanel(
				condition = "input.type == 'dna'  & input.advanced",
				# selectInput("phenotype_dna", "Protein selection", choices = c("all",phenotypes_vector),selected="10"),
				sliderInput("distance","Distance from gene (bp)",min=100000,max=500000,value=200000)
			),
			conditionalPanel(
				condition = "input.type == 'protein'  & input.advanced",
				sliderInput("p_value_cutoff","P-value cutoff (log)",min=4,max=12,value=6,step=0.1)
			),
			conditionalPanel(
				condition = "input.advanced",
				sliderInput("top_label_count","#SNPs to label",min=3,max=30,value=3,step=1)
			),
			
			
			actionButton("goButton","Start"),
			width=4
		),
		mainPanel(
			plotOutput("mainPlot",width = "800px", height = "800px")
		)
	)
)

shinyApp(ui = ui, server = server)


