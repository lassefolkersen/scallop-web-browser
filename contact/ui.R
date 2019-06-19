source("../uifunctions.R")
initialize('con',TRUE)

shinyUI(bootstrapPage(
	head(),
	navigation(),
	titlePanel("About SCALLOP"),
	beginPage(),
	HTML("
<img src=../www/composite-banner-scallop.png alt=''>
<h1>SCALLOP genetics of the proteome</h1>
<p>The SCALLOP consortium is a collaborative framework for discovery and follow-up of genetic associations with proteins on the Olink Proteomics platform. To date, 19 PIs from 10 research institutions have joined the effort, which now comprises summary level data on SNP to protein level associations from almost 30,000 patients or controls. SCALLOP welcomes new members.</p>
	      <p>For more information please contact <u><a href='mailto:Anders.Malarstig@ki.se?subject=Scallop%20inquiry%20from%20scallopconsortium.com'>Anders Mälarstig</a></u></p>
	      <p>For the latest news about SCALLOP: <u><a href='https://www.olink.com/scallop/scallop-news/'>view this news page</a></u></p>
        <p>Or read more in our <u><a href=../www/Scallop_Flyer_v1.0.pdf>flyer</a></u>.</p>
        <br><br><br>
	      <h2>Introduction from Anders Mälarstig</h2>
        <iframe width='932' height='524' src='https://www.youtube.com/embed/mtUAz4uRODc' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe>
	      <br><br><br>
        <p><img src='../www/scallop_logo.png' alt='' width='200' height='121'></p>
        <h2>Current work</h2>
	      <p>Each SCALLOP member works on human study collections from the general population, clinical trials or patients with certain diseases such as coronary artery disease, rheumatoid arthritis, bipolar disease, heart failure, dementias or metabolic syndrome.</p>
	      <p>The aim of the SCALLOP consortium is to identify novel molecular connections and protein biomarkers that are causal in diseases.</p>
	      <p>This work starts with identification of so called protein quantitative trait loci, pQTLs, which are robust connections between a gene variant and the levels of a protein.</p>
	      <p>There are two types of pQTLs:</p>
	      <ul>
	      <li>cis-pQTLs are variants that are proximal to the gene encoding the protein under study whereas trans-pQTLs are distal regulation of proteins via an often unknown path.</li>
	      <li>Trans-pQTLs can provide unique insights of molecular connections in human biology.</li>
	      </ul>
	      <p>Cis-pQTLs are strong instruments for determining if a protein biomarker for disease is causing disease or elevated or suppressed as a consequence of it. The SCALLOP consortium is currently underway with mapping novel pQTLs for several 100s of proteins in unprecedented sample sizes, something which will yield much deeper insights into the trans-regulation of plasma proteins than what has been possible to date.</p>
	      <br><br><br>
        <h2>Identify causal protein biomarkers</h2>
	      <p><img src='../www/scallop-infographic.png'></p>
	      <br><br><br>
        <h2>Operations</h2>
	      <p>To be a member of the SCALLOP consortium you have to be the PI of a study collection with Olink proteomics and genome-wide genotyping data. We also expect members to sign up to the Consortium Agreement, which manages conduct and authorships. <u><a href='../www/SCALLOP_consortium_agreement_March_2019.pdf' target='_blank' rel='noopener'>Download the&nbsp; Consortium Agreement here</a></u>.</p>
	      <p>The leadership for subprojects within the SCALLOP consortium rotates and members can take new ideas and suggestions for additional subprojects to the monthly steering committee meetings.</p>
	      <p>SCALLOP uses a dedicated server for sharing of data. The server is set up under the Danish node of the TRYGGVE server structure. TRYGGVE allows sharing of sensitive data thanks to 2-step authorization procedures and high data security. Thanks to this structure SCALLOP is set up to move to individual-level data should the consortium wish to do so.</p>
	      <br><br><br>
        <h2>Repository browser</h2>
	      <table style='height: 559px;' width='1105'>
	      <tbody>
	      <tr bgcolor='#077183'>
	      <td width='208'><span style='color: #077183;'><strong><span style='color: #ffffff;'>Acronym</span></strong></span></td>
	        <td width='211'><span style='color: #ffffff;'><strong>Design</strong></span></td>
	        <td width='477'><span style='color: #ffffff;'><strong>Website</strong></span></td>
	        </tr>
	        <tr>
	        <td>ASAP</td>
	        <td>Aortic valve surgery</td>
	        <td><a href='http://www.ki.se/en/meds/team-per-eriksson'>www.ki.se/en/meds/team-per-eriksson</a></td>
	        </tr>
	        <tr>
	        <td>BioFinder</td>
	        <td>Dementia</td>
	        <td><a href='http://biofinder.se/'>http://biofinder.se</a></td>
	        </tr>
	        <tr>
	        <td>COMBINE</td>
	        <td>Rheumatoid arthritis</td>
	        <td><a href='http://www.combinesweden.se/'>www.combinesweden.se/</a></td>
	        </tr>
	        <tr>
	        <td>EpiHealth</td>
	        <td>Prospective observational</td>
	        <td><a href='http://www.epihealth.se/For-scientists/'>www.epihealth.se/For-scientists/</a></td>
	        </tr>
	        <tr>
	        <td>Estonian Biobank</td>
	        <td>Population based</td>
	        <td><a href='http://www.geenivaramu.ee/en/access-biobank'>www.geenivaramu.ee/en/access-biobank&nbsp;</a></td>
	        </tr>
	        <tr>
	        <td>HELIC MANOLIS</td>
	        <td>Population isolate</td>
	        <td><a href='http://www.helic.org/sites.html'>www.helic.org/sites.html</a></td>
	        </tr>
	        <tr>
	        <td>IMPROVE</td>
	        <td>Prospective, metabolic syndrome</td>
	        <td><a href='http://www.ki.se/en/meds/team-hamsten'>www.ki.se/en/meds/team-hamsten</a></td>
	        </tr>
	        <tr>
	        <td>INTERVAL</td>
	        <td>Blood donors</td>
	        <td><a href='http://www.intervalstudy.org.uk/'>www.intervalstudy.org.uk/</a></td>
	        </tr>
	        <tr>
	        <td>Kadoorie biobank</td>
	        <td>Pancreatic cancer</td>
	        <td><a href='http://www.ckbiobank.org/'>www.ckbiobank.org</a></td>
	        </tr>
	        <tr>
	        <td>LifeLines Deep</td>
	        <td>Population based</td>
	        <td><a href='http://www.lifelines.nl/researcher/biobank-lifelines/additional-studies/lifelines-deep'>www.lifelines.nl/researcher/biobank-lifelines/additional-studies/lifelines-deep</a></td>
	        </tr>
	        <tr>
	        <td>MPP-RES</td>
	        <td>Heart failure</td>
	        <td><a href='http://www.ludc.med.lu.se/malmoe-prevention-project-mpp/'>www.ludc.med.lu.se/malmoe-prevention-project-mpp/</a></td>
	        </tr>
	        <tr>
	        <td>NSPHS</td>
	        <td>Population isolate</td>
	        <td><a href='http://www.ncbi.nlm.nih.gov/pubmed/20568910'>www.ncbi.nlm.nih.gov/pubmed/20568910</a></td>
	        </tr>
	        <tr>
	        <td>ORCADES</td>
	        <td>Population isolate</td>
	        <td><a href='https://www.ed.ac.uk/viking/about-us/orcades'>www.ed.ac.uk/viking/about-us/orcades</a></td>
	        </tr>
	        <tr>
	        <td>Pfizer trials (RA, UC, psoriasis)</td>
	        <td>Clinical trials, RA</td>
	        <td><a href='http://www.pfizer.com/'>www.pfizer.com</a></td>
	        </tr>
	        <tr>
	        <td>PIVUS</td>
	        <td>Prospective observational</td>
	        <td><a href='http://www.medsci.uu.se/pivus/'>www.medsci.uu.se/pivus/</a></td>
	        </tr>
	        <tr>
	        <td>STABILITY</td>
	        <td>Acute coronary syndrome</td>
	        <td><a href='http://www.nejm.org/doi/full/10.1056/NEJMoa1315878#t=article%C2%A0'>www.nejm.org/doi/full/10.1056/NEJMoa1315878#t=article&nbsp;</a></td>
	      </tr>
	        <tr>
	        <td>STANLEY</td>
	        <td>Bipolar/depression</td>
	        <td><a href='http://www.ki.se/meb/stanleyswebic-studien'>www.ki.se/meb/stanleyswebic-studien</a></td>
	        </tr>
	        <tr>
	        <td>ULSAM</td>
	        <td>Prospective observational</td>
	        <td><a href='http://www.pubcare.uu.se/ulsam/'>www.pubcare.uu.se/ulsam/</a></td>
	        </tr>
	        <tr>
	        <td>VIS</td>
	        <td>Population isolate</td>
	        <td><a href='http://www.ed.ac.uk/mrc-human-genetics-unit/research/qtl-group/'>www.ed.ac.uk/mrc-human-genetics-unit/research/qtl-group/</a></td>
	        </tr>
	        </tbody>
	        </table>
	        <hr>
	        <h2>pQTL publications and data from SCALLOP members</h2>
	        <p>Folkersen L, Fauman E, Sabater-Lleal M, Strawbridge R, Frånberg M, Sennblad B, Baldassarre D, Veglia F, Humphries S, Rauramaa R, de Faire U, Smit A, Giral P, Kurl S, Mannarino E, Enroth S, Johansson A, Bosdotter Enroth S, Gustafsson S, Lind L, Lindgren C, Morris A, Giedraitis V, Silveira A, Franco-Cereceda A, Tremoli E, IMPROVE study group , Gyllensten U, Ingelsson E, Brunak S, Eriksson P, Ziemek D, Hamsten A and Mälarstig A.<strong>&nbsp;Mapping of 79 loci for 83 plasma protein biomarkers in cardiovascular disease.</strong>&nbsp;(2017) PLOS Genetics 13(4), doi.org/10.1371/journal.pgen.1006706<br>
	        
	        
	        <p>Enroth S, Johansson A, Bosdotter Enroth S and Ulf Gyllensten U.&nbsp;<strong>Strong effects of genetic and lifestyle factors on biomarker variation and use of personalized cutoffs</strong>. Nature Commun. (2014) Aug 22;5:4684. doi: 10.1038/ncomms5684.<br>
	        
	        <hr>
	        <p>&nbsp;</p>
          

	        </div>
	        

	      				"),
	footer()
))












