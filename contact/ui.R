source("../uifunctions.R")
initialize('con',TRUE)

shinyUI(bootstrapPage(
	head(),
	navigation(),
	titlePanel("Contact"),
	beginPage(),
	HTML(	"
<b>Analysis and design:</b><br><a href='http://www.cbs.dtu.dk/staff/show-staff.php?id=1202'>Lasse Folkersen</a> (<u><a href='http://www.google.com/recaptcha/mailhide/d?k=01pdzWyCfeU-_1PRAPdKlJfg==&amp;c=3eyQPG-VqkHu6ECGRBSHdRraKCXOUsVtLpuyWWt-dpY=' onclick='window.open('http://www.google.com/recaptcha/mailhide/d?k\07501pdzWyCfeU-_1PRAPdKlJfg\75\75\46c\0753eyQPG-VqkHu6ECGRBSHdRraKCXOUsVtLpuyWWt-dpY\075', '', 'toolbar=0,scrollbars=0,location=0,statusbar=0,menubar=0,resizable=0,width=500,height=300'); return false;' title='Reveal this e-mail address'>email</a></u>).<br><br>
				
<b>Project lead:</b><br>Anders Malarstig.<br><br>

<b>Analysis group:</b><br>
                Maria Sabater-Lleal<br>
                Rona Strawbridge<br>
                Eric Fauman<br>
                Daniel Ziemek<br>
<br><br>
<b>Web-page styling:</b><br>Stefan Delport.<br><br>

	      				"),
	footer()
))












