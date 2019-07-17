source("../uifunctions.R")
initialize('con',TRUE)

shinyUI(bootstrapPage(
	head(),
	navigation(),
	titlePanel("Terms of use"),
	beginPage(),
	HTML("
Everything on this site is available for use under the <u></u><a href='https://creativecommons.org/licenses/by/2.0/'>CC BY 2.0 license</a></u>.
<br><br>
	     All pQTL summary stats data has also been deposited at zenodo.org for long-term preservation using this DOI <a href='https://doi.org/10.5281/zenodo.2615265'>10.5281/zenodo.2615265</a>. 
	     <br><br>
	     Additionally, the code running this web-site is released at <a href='https://github.com/lassefolkersen/scallop-web-browser'>github</a> (GNU General Public License v3.0)."),
	footer()
))












