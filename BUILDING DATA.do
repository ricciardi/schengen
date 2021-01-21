adopath + "U:\My Documents\PAPERS PUBLICATION\ADOPATH"
sysdir set PLUS "U:\My Documents\PAPERS PUBLICATION\ADOPATH"
global datain "S:\databases\1_Projets\Schengen_AAlbanese\LFS\LFS microdata full set 25sept2020\Quarterly Files_95_2019\Quarterly Files"
global dataout "R:\CBW and Schengen\2. Datasets\quarterly"
global gdp "R:\CBW and Schengen\COVARIATES EUROSTAT\selected\xls"

 
 
 *BUILD GDP CONTROL VARIABLE
 cd "$gdp"
clear all
foreach var in   nama_10r_2gdp nama_10_pc {
	clear
import excel using "`var'.xls", firstrow  
drop if _n==1
cap: rename A i
drop if i=="Special value:"
drop if i==""
des _all, varlist
local lista=r(varlist)
foreach var2 in `lista' {
replace `var2'="" if `var2'==":" 
}
drop if i==""

destring _all, replace
reshape long j_, i(i) j(year) 
rename i country
local h=cond("`var'"=="nama_10_pc",  "GDPcapitaC"  , "GDPcapitaR" )
rename j_  `h'
compress
sort country year
save "R:\CBW and Schengen\COVARIATES EUROSTAT\selected/`var'", replace
}
foreach var in    nama_10r_2gdp nama_10_pc {
sort country year
	merge  country year using "R:\CBW and Schengen\COVARIATES EUROSTAT\selected/`var'"
	drop _merge
}
rename country REGION
gen country=substr(REGION, 1, 2 )
*country GDP to each regions
bys country year: egen GDPc_COUNTRY=mean(GDPcapitaC)
drop GDPcapitaC 
compress
sort REGION year
replace REGION=REGION+"0"  if REGION=="AT1" | REGION=="AT2" | REGION=="AT3"
replace REGION=REGION+"0"  if country=="DE"
replace REGION=REGION+"00"  if REGION=="DK"
replace REGION=REGION+"00"  if REGION=="NL"
replace REGION=REGION+"0"  if REGION=="HU1"
replace REGION="IE01" if REGION=="IE04" 
replace REGION= "FR24" if  REGION=="FRB0"
replace REGION= "FR26" if  REGION=="FRC1"
replace REGION= "FR43" if  REGION=="FRC2"
replace REGION= "FR25" if  REGION=="FRD1"
replace REGION= "FR23" if  REGION=="FRD2"
replace REGION= "FR30" if  REGION=="FRE1"
replace REGION= "FR22" if  REGION=="FRE2"
replace REGION= "FR42" if  REGION=="FRF1"
replace REGION= "FR21" if  REGION=="FRF2"
replace REGION= "FR41" if  REGION=="FRF3"
replace REGION= "FR51" if  REGION=="FRG0"
replace REGION= "FR52" if  REGION=="FRH0"
replace REGION= "FR61" if  REGION=="FRI1"
replace REGION= "FR63" if  REGION=="FRI2"
replace REGION= "FR53" if  REGION=="FRI3"
replace REGION= "FR81" if  REGION=="FRJ1"
replace REGION= "FR62" if  REGION=="FRJ2"
replace REGION= "FR72" if  REGION=="FRK1"
replace REGION= "FR71" if  REGION=="FRK2"
replace REGION= "FR82" if  REGION=="FRL0"
replace REGION= "FR83" if  REGION=="FRM0"
replace REGION= "PL11" if REGION=="PL71" 
replace REGION= "PL33" if REGION=="PL72" 
replace REGION= "PL31" if REGION=="PL81" 
replace REGION= "PL32" if REGION=="PL82" 
replace REGION= "PL34" if REGION=="PL84" 
replace REGION= "PL12" if REGION=="PL9" 

sort country REGION year
*use country information when regional is missing i.e. France (before 2015) and Switzerland
replace GDPcapitaR=GDPc_COUNTRY  if GDPcapitaR==.
drop GDPc_COUNTRY

**gdp ratio to adjacent country (highest GDP of foreign regions)
foreach var in  GDPcapitaR   { 
gen R_`var'=.
forvalues y=2005(1)2018 {
*******************
sum `var'  if (REGION=="ES11"  | REGION=="ES41" ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="PT11")  & year==`y'
sum `var'  if ( REGION=="ES43"  | REGION=="ES41" ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="PT16" )& year==`y'
sum `var'  if ( REGION=="ES61"  | REGION=="ES43" ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="PT18") & year==`y'
sum `var'  if ( REGION=="ES61")  & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="PT15" ) & year==`y'

sum `var'  if ( REGION=="PT11" ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="ES11" ) & year==`y'
sum `var'  if ( REGION=="PT11"  |  REGION=="PT16" ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="ES41" )& year==`y'
sum `var'  if ( REGION=="PT18"  |  REGION=="PT16"  )& year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="ES43" ) & year==`y'
sum `var'  if ( REGION=="PT18"  |  REGION=="PT15" ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="ES61" ) & year==`y'
sum `var'  if (   REGION=="FR61"  )& year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="ES21" | REGION=="ES22"  )& year==`y'
sum `var'  if (   REGION=="FR61" |  REGION=="FR62" ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="ES24" ) & year==`y'
sum `var'  if (   REGION=="FR81" |  REGION=="FR62"   ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="ES51"  ) & year==`y'

sum `var'  if ( REGION=="ES21"  | REGION=="ES22"  | REGION=="ES24"    ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="FR61"  ) & year==`y'
sum `var'  if ( REGION=="ES22"  | REGION=="ES24"    | REGION=="ES51"   ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="FR62"  ) & year==`y'
sum `var'  if ( REGION=="ES51"    ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="FR81"  ) & year==`y'
sum `var'  if ( REGION=="ITC3"   | REGION=="ITC1"   ) & year==`y' 
replace R_`var' =`var' /r(max) if ( REGION=="FR82"  ) & year==`y'
sum `var'  if (  country=="CH"    ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="FR71" | REGION=="FR43"  ) & year==`y'
sum `var'  if ( country=="CH"     ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="FR42"  ) & year==`y'
sum `var'  if (  REGION=="LU00"    | REGION=="DEC0"   | REGION=="BE34"  ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="FR41"  ) & year==`y'
sum `var'  if ( REGION=="BE35"    | REGION=="BE34"  | REGION=="BE32"  ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="FR21"  ) & year==`y'
sum `var'  if ( REGION=="BE25"    | REGION=="BE32"  | REGION=="BE35"  ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="FR30"  ) & year==`y'

sum `var'  if ( REGION=="FR30"    | REGION=="NL34"    ) & year==`y' 
replace R_`var' =`var' /r(max) if ( REGION=="BE25"  ) & year==`y'

sum `var'  if ( REGION=="NL34"     ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="BE23"  ) & year==`y'

sum `var'  if ( REGION=="NL41"     ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="BE21"  ) & year==`y'

sum `var'  if ( REGION=="NL41" | REGION=="NL42"      ) & year==`y'
replace R_`var' =`var' /r(max) if ( REGION=="BE22"  ) & year==`y'

sum `var'  if ( REGION=="NL42"   | REGION=="DEB0"   | REGION=="DEA0"   | REGION=="LU00"  ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="BE33"    ) & year==`y'
sum `var'  if ( REGION=="FR41"   | REGION=="FR21"   | REGION=="LU00"  ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="BE34"    ) & year==`y'
sum `var'  if ( REGION=="FR30"   | REGION=="FR21"    ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="BE32"    ) & year==`y'
sum `var'  if ( REGION=="FR21"    ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="BE35"  ) & year==`y'

sum `var'   if ( REGION=="BE23" | REGION=="BE21" | REGION=="BE22"     | REGION=="BE33"     | REGION=="DEA0"    | REGION=="DE90"     ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="NL00"  ) & year==`y'

sum `var'   if ( REGION=="DEF0" | REGION=="SE22"  ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="DK00"  ) & year==`y'

sum `var'   if ( REGION=="BE34" | REGION=="BE33"  | REGION=="DEB0" | REGION=="DEC0" | REGION=="FR41"   ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="LU00"  ) & year==`y'

sum `var'   if ( REGION=="IE01" |  REGION=="IE04"   ) & year==`y' 
replace R_`var' =`var' /r(max) if (  REGION=="UKN0"  ) & year==`y'

sum `var'   if ( REGION=="FR82"   ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="ITC3"  ) & year==`y'
sum `var'   if ( country=="CH"   ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="ITC1" |  REGION=="ITC2"  |  REGION=="ITC4" |   REGION=="ITH1"   ) & year==`y'
sum `var'   if ( REGION=="AT30" | REGION=="AT20"   ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="ITH3"  ) & year==`y'
sum `var'   if ( REGION=="SI04" | REGION=="AT20"   ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="ITH4"  ) & year==`y'

sum `var'   if ( REGION=="DE80" | REGION=="DE40"   ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="PL42"  ) & year==`y'

sum `var'   if ( REGION=="DED0" | REGION=="DE40"   ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="PL43"  ) & year==`y'
sum `var'   if ( REGION=="DED0"   ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="PL51"  ) & year==`y'

sum `var'   if ( REGION=="DED0"  |  REGION=="DE20"   ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="CZ04"  ) & year==`y'
sum `var'   if ( REGION=="AT30"  |  REGION=="DE20"  |  REGION=="AT10"   ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="CZ03" ) & year==`y'
sum `var'   if (  REGION=="AT10"  ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="CZ06" |  REGION=="SK01"  |  REGION=="SK02" |  REGION=="HU22" ) & year==`y'

sum `var'   if (   REGION=="AT20"  ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="SI03" ) & year==`y'
sum `var'   if (   REGION=="AT20"  |    REGION=="ITH4" ) & year==`y'
replace R_`var' =`var' /r(max) if (  REGION=="SI04" ) & year==`y'

*germany
sum `var'  if (  REGION=="FR42"  | country=="CH"   ) & year==`y'
replace R_`var' =`var' /r(max) if (   REGION=="DE10" ) & year==`y'
sum `var'  if (  REGION=="FR42"  | country=="LU"   | REGION=="BE33"   ) & year==`y'
replace R_`var' =`var' /r(max) if (   REGION=="DEB0" ) & year==`y'
sum `var'  if (  REGION=="FR42"  | country=="LU"   | REGION=="FR41"   ) & year==`y'
replace R_`var' =`var' /r(max) if (   REGION=="DEC0" ) & year==`y'
sum `var'  if (   REGION=="NL42" |   REGION=="NL22"  |  REGION=="NL21" | REGION=="BE33"   ) & year==`y'
replace R_`var' =`var' /r(max) if (   REGION=="DEA0" ) & year==`y'
sum `var'  if (  REGION=="NL11" |   REGION=="NL13" |  REGION=="NL21"   ) & year==`y'
replace R_`var' =`var' /r(max) if (   REGION=="DE90" ) & year==`y'

sum `var'  if (   REGION=="DK03"    ) & year==`y'
replace R_`var' =`var' /r(max) if (   REGION=="DEF0"  ) & year==`y'

}
}
keep if year>=2005
compress
sort REGION year

keep if country=="BE" |   country=="CZ" |   country=="DE" |   country=="DK" |   country=="ES" |   country=="FR" |   country=="HU" |   country=="IT" |   country=="LU" |   country=="NL" |   country=="PL" |   country=="SI" |   country=="SK" |   country=="UK"
sort year REGION
save  "R:\CBW and Schengen\COVARIATES EUROSTAT\selected\GDP_final.dta" , replace



*MOVE TO EUROSTAT DATABASE
clear all
foreach H in   BE CZ DE DK  ES  FR HU  IT  LU  NL  PL   SI SK UK {
forvalues Y=2005(1)2019 {
forvalues Q=1(1)4 {
capture {
import delimited "$datain/`H'`Y'Q`Q'.csv", clear
unab allvars: _all
local masterlist "age coeff    country   countryw        region     ilostat    "


local keeplist: list allvars & masterlist
keep `keeplist'
gen tq=yq(`Y', `Q')
tostring *, replace 
compress
save "$dataout/`H'`Y'Q`Q'.dta", replace
}
}
}
}

clear all
foreach H in   BE CZ DE DK  ES  FR HU  IT  LU  NL  PL   SI SK UK {
forvalues Y=2005(1)2019 {
forvalues Q=1(1)4 {
di "`H'`Y'Q`Q'"
capture {
append using "$dataout/`H'`Y'Q`Q'.dta"
}
	keep tq age coeff    country   countryw        region     ilostat 
}
}
}


destring tq age coeff                    ilostat  , replace
compress
format tq %tq
gen year=year(dofq(tq))
cap: tostring region, replace


*cleaning
keep if age>=22 & age<62
*drop overseas:
*france
drop if (region=="A1" | region=="A2" | region=="A3" | region=="A4" | region=="A5" | region=="Y1" | region=="Y2"  | region=="Y3"   | region=="Y4" ) & country=="FR"
*spain
drop if region=="70"   & country=="ES"
*ceita y melila 
drop if region=="63"   & country=="ES"
drop if region=="64"   & country=="ES"


*USE SAME DEFINITION OF REGIONS OVER TIME
*belgium
replace region = "A" if region=="10" & year<=2017 & country=="BE"
replace region = "B" if region=="24" & year<=2017 & country=="BE"
replace region = "C" if region=="31" & year<=2017 & country=="BE"
replace region = "D" if region=="21" & year<=2017 & country=="BE"
replace region = "E" if region=="22" & year<=2017 & country=="BE"
replace region = "F" if region=="23" & year<=2017 & country=="BE"
replace region = "G" if region=="25" & year<=2017 & country=="BE"
replace region = "H" if region=="32" & year<=2017 & country=="BE"
replace region = "I" if region=="33" & year<=2017 & country=="BE"
replace region = "L" if region=="34" & year<=2017 & country=="BE"
replace region = "M" if region=="35" & year<=2017 & country=="BE"
replace region="10" if  region=="A"   & country=="BE"
replace region = "24" if region=="B"  & country=="BE"
replace region = "31" if region=="C"  & country=="BE"
replace region = "21" if region=="D"  & country=="BE"
replace region = "22" if region=="E"  & country=="BE"
replace region = "23" if region=="F"  & country=="BE"
replace region = "25" if region=="G" & country=="BE"
replace region = "32" if region=="H"  & country=="BE"
replace region = "33" if region=="I"  & country=="BE"
replace region = "34" if region=="L"  & country=="BE"
replace region = "35" if region=="M"  & country=="BE"
replace region="10" if  region=="A"   & country=="BE"

*DK (NUTS AVAILABLE ONLY FROM 2007)
replace region="00"  if country=="DK"

*SI after 2001
replace region="3" if region=="1" & country=="SI" 
replace region="4" if region=="2" & country=="SI" 

*UK -  NORTHERN IRELAND, OTHERS NO BORDER
replace region="N0" if region=="B0" & country=="UK"



replace region="0"+region if country=="CZ"
replace region="0"+region if country=="SK"
replace region="0"+region if country=="LU"

compress
gen REGION=country+region
replace REGION="HU10" if REGION=="HU11" | REGION=="HU12" 
replace REGION="PL11" if REGION=="PL72"
replace REGION="PL33" if REGION=="PL81"
replace REGION="PL31" if REGION=="PL82"
replace REGION="PL32" if REGION=="PL84"
replace REGION="PL34" if REGION=="PL91"
replace REGION="PL12" if REGION=="PL92"
replace REGION= "FR24" if  REGION=="FRB0"
replace REGION= "FR26" if  REGION=="FRC1"
replace REGION= "FR43" if  REGION=="FRC2"
replace REGION= "FR25" if  REGION=="FRD1"
replace REGION= "FR23" if  REGION=="FRD2"
replace REGION= "FR30" if  REGION=="FRE1"
replace REGION= "FR22" if  REGION=="FRE2"
replace REGION= "FR42" if  REGION=="FRF1"
replace REGION= "FR21" if  REGION=="FRF2"
replace REGION= "FR41" if  REGION=="FRF3"
replace REGION= "FR51" if  REGION=="FRG0"
replace REGION= "FR52" if  REGION=="FRH0"
replace REGION= "FR61" if  REGION=="FRI1"
replace REGION= "FR63" if  REGION=="FRI2"
replace REGION= "FR53" if  REGION=="FRI3"
replace REGION= "FR81" if  REGION=="FRJ1"
replace REGION= "FR62" if  REGION=="FRJ2"
replace REGION= "FR72" if  REGION=="FRK1"
replace REGION= "FR71" if  REGION=="FRK2"
replace REGION= "FR82" if  REGION=="FRL0"
replace REGION= "FR83" if  REGION=="FRM0"
replace REGION= "PL11" if REGION=="PL71" 
replace REGION= "PL33" if REGION=="PL72" 
replace REGION= "PL31" if REGION=="PL81" 
replace REGION= "PL32" if REGION=="PL82" 
replace REGION= "PL34" if REGION=="PL84" 
replace REGION= "PL12" if REGION=="PL9" 
replace REGION="DE00" if country=="DE" & region=="0"
replace region="DE00" if country=="DE" & region=="0"
foreach var in NL   {
replace REGION=REGION+"0" if country=="`var'"
replace region=region+"0" if country=="`var'"
}
replace REGION="SI03" if REGION=="SI3"
replace region="03" if region=="3" & country=="SI"
replace REGION="SI04" if REGION=="SI4"  
replace region="04" if region=="4" &  country=="SI"
drop if REGION=="PL0"
drop if REGION=="HU0"


* SELECT ONLY BORDER REGION
*spain
cap: drop border
gen border=1 if REGION=="ES51" | REGION=="ES24"  | REGION=="ES22"  | REGION=="ES21"      | REGION=="ES11"  | REGION=="ES41"  | REGION=="ES43"  | REGION=="ES61"

*belgium
replace border=1 if REGION=="BE35" |  REGION=="BE34" | REGION=="BE33"  | REGION=="BE22"  | REGION=="BE21" | REGION=="BE23" | REGION=="BE25"  | REGION=="BE32"

*CZ 
replace border=1 if   REGION=="CZ03" | REGION=="CZ04"   | REGION=="CZ06" 

*dk
replace border=1 if country=="DK"

*GERMANY - NUTS 1 only from  2002. 
replace border=1       if REGION=="DE90" | REGION=="DEA0" | REGION=="DEB0" | REGION=="DEB0" | REGION=="DEC0" | REGION=="DEF0" | REGION=="DE10"

	   
*FRANCE
replace border=1 if REGION=="FR30" |  REGION=="FR21"   |  REGION=="FR41"  |  REGION=="FR42"   |  REGION=="FR43"   |  REGION=="FR71"  |  REGION=="FR82"    |  REGION=="FR81"   |  REGION=="FR62"  |  REGION=="FR61"  

*ITALY
replace border=1 if    REGION=="ITC3"  |  REGION=="ITH3" |  REGION=="ITC1" |  REGION=="ITC2" |  REGION=="ITC4" |  REGION=="ITH1"
  
*luxembourg
replace border=1 if country=="LU" 

*hungary 
replace border=1 if REGION=="HU22" 

*poland 
replace border=1 if 			  REGION=="PL42" |  REGION=="PL43"   |  REGION=="PL51"     

*NETHERLANDS -> NO VARIABLES ON REGIONS 
replace border=1 if country=="NL" 

*SLOVENIA
replace border=1 if REGION=="SI03" |  REGION=="SI04"   

*SLOVAKIA
replace border=1 if REGION=="SK01" |  REGION=="SK02"   

*uk - ONLY BORDER WITH IRELAND 
replace border=1 if REGION=="UKN0" 

replace border=0 if border==.

*drop internal regions
drop if border==0




*******OUTCOME
gen CBWbord=.
replace CBWbord=countryw=="PT"  if REGION=="ES61" | REGION=="ES43" |  REGION=="ES41" |  REGION=="ES11"
replace CBWbord=countryw=="FR"  if REGION=="ES21" | REGION=="ES22" |  REGION=="ES24" |  REGION=="ES51"
replace CBWbord=countryw=="ES"  if REGION=="FR61" | REGION=="FR62" |  REGION=="FR81" 
replace CBWbord=countryw=="BE"  if REGION=="FR30" | REGION=="FR21"
replace CBWbord=countryw=="BE" | countryw=="LU"   | countryw=="DE"  if REGION=="FR41" 
replace CBWbord=countryw=="DE" | countryw=="CH"  if REGION=="FR42" 
replace CBWbord=countryw=="CH"  if REGION=="FR43" 
replace CBWbord=countryw=="CH" | countryw=="IT"  if REGION=="FR71" 
replace CBWbord=countryw=="IT"  if REGION=="FR82" 
replace CBWbord=countryw=="FR"  if REGION=="BE32"  | REGION=="BE35" 
replace CBWbord=countryw=="FR"  | countryw=="LU"  if  REGION=="BE34" 
replace CBWbord=countryw=="NL"  | countryw=="LU"  | countryw=="DE" if  REGION=="BE33" 
replace CBWbord=countryw=="NL"   if  REGION=="BE22"  |  REGION=="BE21"  |  REGION=="BE23" 
replace CBWbord=countryw=="NL" | countryw=="FR"   if   REGION=="BE25" 
replace CBWbord=countryw=="BE"  | countryw=="DE"   if  country=="NL"
replace CBWbord=countryw=="DE"  | countryw=="SE"   if  country=="DK"
replace CBWbord=countryw=="IE" if  country=="UK"
replace CBWbord=countryw=="LT" if  REGION=="PL34"
replace CBWbord=countryw=="SK" if  REGION=="PL32" |  REGION=="PL21"
replace CBWbord=countryw=="SK" | countryw=="CZ" if  REGION=="PL22"
replace CBWbord=countryw=="CZ" if  REGION=="PL52"
replace CBWbord=countryw=="CZ" | countryw=="DE"  if  REGION=="PL51"
replace CBWbord=countryw=="DE"  if  REGION=="PL43" | REGION=="PL42"
replace CBWbord=countryw=="PL" | countryw=="DE"  if  REGION=="CZ05" 
replace CBWbord= countryw=="DE"  if  REGION=="CZ04" 
replace CBWbord= countryw=="DE"  | countryw=="AT"  if  REGION=="CZ03" 
replace CBWbord= countryw=="SK"  | countryw=="AT"  if  REGION=="CZ06" 
replace CBWbord= countryw=="SK"  | countryw=="PL"  if  REGION=="CZ07" |  REGION=="CZ08" 
replace CBWbord=countryw=="AT"   if  REGION=="SK01" 
replace CBWbord= countryw=="AT" | countryw=="CZ"  | countryw=="HU"   if  REGION=="SK02" 
replace CBWbord= countryw=="HU" | countryw=="PL" | countryw=="CZ"   if  REGION=="SK03" 
replace CBWbord= countryw=="HU" | countryw=="PL"   if   REGION=="SK04" 
replace CBWbord=countryw=="SI" | countryw=="AT" | countryw=="SK"     if  REGION=="HU22" 
replace CBWbord=countryw=="SK"  if  REGION=="HU31"  | REGION=="HU10"  | REGION=="HU21" 
* slovania does not have information on where exactly work - use anywhere in EU
replace CBWbord= countryw=="005-EU28"   if  REGION=="SI03" 
replace CBWbord= countryw=="005-EU28"  if   REGION=="SI04" 
replace CBWbord=countryw=="FR"   if   REGION=="ITC3"
replace CBWbord=countryw=="FR"  | countryw=="CH"   if   REGION=="ITC1" | REGION=="ITC2"
replace CBWbord=countryw=="CH"   if   REGION=="ITC4"
*bolzano
replace CBWbord=countryw=="AT"  | countryw=="CH"  if   REGION=="IT31" |  REGION=="ITH1"
replace CBWbord=countryw=="AT"   if   REGION=="ITH3"
replace CBWbord=countryw=="AT"  | countryw=="SI"  if   REGION=="ITH4"
*germany
replace CBWbord=countryw=="CH" | countryw=="FR"  if   REGION=="DE10"
replace CBWbord=countryw=="LU" | countryw=="FR"  if   REGION=="DEC0"
replace CBWbord=countryw=="LU" | countryw=="FR" | countryw=="BE"  if   REGION=="DEB0"
replace CBWbord= countryw=="NL" | countryw=="BE"  if   REGION=="DEA0"
replace CBWbord= countryw=="NL"  if   REGION=="DE90"
replace CBWbord= countryw=="DK"  if   REGION=="DEF0"
replace CBWbord= countryw=="PL"  if   REGION=="DE80" | REGION=="DE40"
replace CBWbord= countryw=="PL" | countryw=="CZ"   if   REGION=="DED0"
replace CBWbord= countryw=="AT" | countryw=="CZ"   if   REGION=="DE20"
replace CBWbord= countryw=="FR"   if   REGION=="DE10"
replace CBWbord= countryw=="FR" | countryw=="DE" | countryw=="BE"   if   country=="LU"

*clean IF MISSING INFO
replace CBWbord=. if countryw=="" | countryw=="." | countryw=="NO ANSWER" 
*conditional on employment
gen CBWbordEMPL=CBWbord if ilostat==1




*SCHENGEN AND FOM TREATMENTS
*clusters
gen tr_SCHENGEN=192 if REGION=="CZ03" | REGION=="CZ04" | REGION=="CZ06"  | REGION=="HU22"  | REGION=="PL42"  | REGION=="PL43"  | REGION=="PL51"  | REGION=="SI03"  | REGION=="SI04" |  REGION=="SK01" | REGION=="SK02"
replace tr_SCHENGEN=196 if REGION=="DE10" | REGION=="FR42" | REGION=="FR43"  | REGION=="FR71"  | REGION=="ITC1"  | REGION=="ITC2"  | REGION=="ITC4"  | REGION=="ITH1"
gen tr_FoM=204 if REGION=="CZ03" | REGION=="CZ04" | REGION=="CZ06"  | REGION=="HU22"  | REGION=="PL42"  | REGION=="PL43"  | REGION=="PL51"  | REGION=="SI03"  | REGION=="SI04" |  REGION=="SK01" | REGION=="SK02"
replace tr_FoM=189 if REGION=="DE10" | REGION=="FR42" | REGION=="FR43"  | REGION=="FR71"  | REGION=="ITC1"  | REGION=="ITC2"  | REGION=="ITC4"  | REGION=="ITH1"
format tr_SCHENGEN %tq
format tr_FoM %tq

gen treated_CBW=0 
replace treated_CBW=1 if tr_SCHENGEN==192
replace treated_CBW=2  if tr_SCHENGEN==196

label define tr222   0 "Controls (always-treated)" 1 "2008: Schengen; 2011: FoM"  2 "2007: FoM; 2009: Schengen"  , replace
label values  treated_CBW tr222

*treatment dummies
gen SCHENGEN_CBW=tq>=tr_SCHENGEN
gen FoM_CBW=tq>=tr_FoM
replace SCHENGEN_CBW=1 if treated_CBW==0
replace FoM_CBW=1 if treated_CBW==0

drop if tq==.
compress
gen quarter=quarter(dofq(tq))
compress


keep country REGION tq year    CBWbord CBWbordEMPL      countryw    age  treated_CBW SCHENGEN_CBW FoM_CBW    tr_SCHENGEN tr_FoM coeff
save "$dataout/LFS_redQ04_19_f.dta", replace



*from individual data to regional data

use "$dataout/LFS_redQ04_19_f.dta", replace
replace coef=coef*1000
sort REGION year
keep if tq>=tq(2005q1)
drop   country 
*calculate population size
gen CBWbordcoef=CBWbord*coef
bys REGION tq: egen N_CBWbord=total(CBWbordcoef)
replace N_CBWbord=round(N_CBWbord)
bys REGION tq: egen pop=total(coef)
replace pop=round(pop)
collapse pop N_CBWbord CBWbord CBWbordEMPL  treated_CBW SCHENGEN_CBW FoM_CBW    tr_SCHENGEN   year  [pw=coef],  by(REGION tq)
cap: drop _merge
sort year REGION
merge  year REGION using "R:\CBW and Schengen\COVARIATES EUROSTAT\selected\GDP_final.dta" 
*fix missing country in 2019(from gdp dataset)
replace country=substr(REGION, 1,2 ) if year==2019 & country==""
drop _merge
drop if tq==.
sort  REGION tq year 
compress


save "R:\CBW and Schengen\COVARIATES EUROSTAT\selected\FINAL_21.dta" , replace



