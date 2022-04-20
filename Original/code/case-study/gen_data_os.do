
clear all


/*##############################################################################
##############################################################################

	Inclusion/Exclusion Criteria

##############################################################################
##############################################################################*/


import delimited "~/path/f143_av3_os_pub.dat", clear 


*f134days: Days since randomization or enrollment
keep id elst2yr_3 plst2yr_3 hlst2yr_3 tepiwk 

*study population of postmenopausal women who in the questionnaire had reported no use of any hormone therapy during the prior 2-year period 
keep if hlst2yr_3==0 &  plst2yr_3==0 & elst2yr_3==0

save "~/path/stata_data/f143_os_l.dta", replace



/*##############################################################################
##############################################################################

	OUTCOME

##############################################################################
##############################################################################*/


import delimited "~/path/f134_os_pub.dat", clear 


*f134days: Days since randomization or enrollment
keep id f134days

save "~/path/stata_data/f134_os_pub_l.dta", replace



import delimited "~/path/fust_os_pub.dat", clear


/*

From fust_ct.pdf:

Usage Notes:
To determine the last/latest follow-up status recorded for a paticipant, 
use one of the following methods: 
	(1) select the participant record with a missing value for FUENDDY or 
	(2) select the participant record with the greatest FUSTARTDY value.

*/

keep id fustartdy fustatus 
sort id fustartdy 
by id: egen mmax = max(fustartdy)
keep if fustartdy == mmax


merge 1:1 id using "~/path/stata_data/f134_os_pub_l.dta"

drop if _merge!= 3
drop mmax
drop _merge

save "~/path/stata_data/fust_os_l.dta", replace





import delimited "~/path/outc_adj_os_pub.dat", clear

keep id chd chddy chdsrc colon colonsrc colorectal colorectalsrc pancreas colondy colorectaldy
merge 1:1 id using "~/path/stata_data/fust_os_l.dta"
drop if _merge!= 3


gen TIME_CHD = chddy 
replace TIME_CHD = fustartdy if chddy == .
replace TIME_CHD = f134days if fustatus == 1 & chd == 0

gen TIME_Colon = colondy
replace TIME_Colon = fustartdy if colondy == .
replace TIME_Colon = f134days if fustatus == 1 & colon == 0

gen TIME_ColonRectal = colorectaldy
replace TIME_ColonRectal = fustartdy if colorectaldy == .
replace TIME_ColonRectal = f134days if fustatus == 1 & colorectal == 0

drop chddy fustartdy _merge f134days 




save "~/path/stata_data/outc_adj_os_fust_os.dta", replace




import delimited "~/path/bmdhip_os_pub.dat", clear


keep if hipdays < 30 & hipdays > -30

keep id hipbmdcor hipdays 


keep id hipbmdcor

save "~/path/stata_data/bmdhip_os_l.dta", replace


import delimited "~/path/bmdspine_os_pub.dat", clear


keep id spnbmdcor spndays 

keep if spndays < 30 & spndays > -30

keep id spnbmdcor

save "~/path/stata_data/bmdspine_os_l.dta", replace


import delimited "~/path/bmdwbody_os_pub.dat", clear


keep id wbbmdcor wbdays 

keep if wbdays < 30 & wbdays > -30

keep id wbbmdcor

save "~/path/stata_data/bmdwbody_os_l.dta", replace

/*##############################################################################
##############################################################################

	CONFOUNDERS

##############################################################################
##############################################################################*/



/*******************************************************************************

	f45_os_pub

*******************************************************************************/

import delimited "~/path/f45_os_pub.dat", clear

keep if f45vtyp == 1

keep id f45multi f45mvmin

label variable f45multi "Multivitamin without Minerals"
label variable f45mvmin "Multivitamin with Minerals"

save "~/path/stata_data/f45_os_l.dta", replace




/*******************************************************************************

	f43_os_pub

*******************************************************************************/



import delimited "~/path/f43_os_pub.dat", clear

keep id totp totpcat totptime toth tothcat tothtime

label variable totp "Estrogen + Progestin Ever at baseline"
label define totpl 0 "No" 1 "Yes"
label values totp totpl

label variable toth "HRT use ever"
label define tothl 0 "No" 1 "Yes"
label values toth tothl


*recode totpcat 3/4 = 3
label variable totpcat "Estrogen + progesterone duration by category (Years)"
label define totpcatl 0 "None" 1 "<5" 2 "5- <10" 3 "10- <15" 4 "15+"
label values totpcat totpcatl

*recode tothcat 3/4 = 3
label variable tothcat "Estrogen and/or Estrogen+progesterone duration by category (Years)"
label define tothcatl 0 "None" 1 "<5" 2 "5- <10" 3 "10- <15" 4 "15+"
label values tothcat tothcatl


label variable totptime "Lifetime estrogen + progest duration"
replace totptime = 0 if totptime == .

label variable tothtime "Lifetime estrogen and/or Estrogen+progesterone duration"
replace tothtime = 0 if tothtime == .

save "~/path/stata_data/f43_os_l.dta", replace


/*******************************************************************************

	dem_os_pub

*******************************************************************************/



import delimited "~/path/dem_os_pub.dat", clear

keep id age ethnic educ income


label variable age "Age at baseline"

recode ethnic (2=5) (1=4) (4=3) (5=1) (3=2) (8=6)

label variable ethnic "Ethnicity"
label define ethnicl 1 "White" 2 "Black" 3 "Hispanic" 4 "American Indian" 5 "Asian/Pacific Islander" 6 "Unknown" 
label values ethnic ethnicl

label variable educ "Education"
label define educl 1 "Didn't go to school" 2 "Grade school (1-4 years)" 3 "Grade school (5-8 years)" 4 "Some high school (9-11 years)" 5 "High school diploma or GED" ///
6 "Vocational or training school" 7 "Some college or Associate Degree" 8 "College graduate or Baccalaureate Degree" 9 "Some post-graduate or professional" 10 "Master's Degree" ///
11 "Doctoral Degree (Ph.D,M.D.,J.D.,etc.)"
label values educ educl

label variable income "Family income"
label define incomel 1 "Less than $10,000" 2 "$10,000 to $19,999" 3 "$20,000 to $34,999" 4 "$35,000 to $49,999" 5 "$50,000 to $74,999" ///
6 "$75,000 to $99,999" 7 "$100,000 to $149,999" 8 "$150,000 or more" 9 "Don't know" ///
11 "Doctoral Degree (Ph.D,M.D.,J.D.,etc.)"
label values income incomel


save "~/path/stata_data/dem_os_l.dta", replace


/*******************************************************************************

	f80_os_pub

*******************************************************************************/



import delimited "~/path/f80_os_pub.dat", clear

*keep if F80VTYP == 1

sort id f80days
by id: gen nn = _n
by id: egen mmax = max(nn)

keep if nn == mmax

keep id bmix syst dias

save "~/path/stata_data/f80_os_l.dta", replace


import delimited "~/path/f80_os_pub.dat", clear

keep if f80vtyp == 1

rename bmix bmix_bl
rename syst syst_bl
rename dias dias_bl

keep id bmix_bl syst_bl dias_bl

save "~/path/stata_data/f80_os_l_bl.dta", replace


/*******************************************************************************

	f34_os_pub

*******************************************************************************/



import delimited "~/path/f34_os_pub.dat", clear


keep id smokevr cigsday alcswk alcohol

label variable smokevr "Smoke or smoked, cigarettes/day"
label define smokevrl 0 "No" 1 "Yes"   
label values smokevr smokevrl

replace cigsday = 0 if cigsday == .

label variable cigsday "Smoke or smoked, cigarettes/day"
label define cigsdayl 0 "0" 1 "<1" 2 "1-4" 3 "5-14" 4 "15-24" 5 "25-34" 6 "35-44" 7 "45+"  
label values cigsday cigsdayl


label variable alcohol "Alcohol intake"
label define alcoholl 1 "Non drinker" 2 "Past drinker" 3 "<1 drink per month" 4 "<1 drink per week" 5 "1 to <7 drinks per week" 6 "7+ drinks per week"  
label values alcohol alcoholl

label variable alcswk "Alcohol servings per week"



save "~/path/stata_data/f34_os_l.dta", replace


/*******************************************************************************

	f31_os_pub

*******************************************************************************/



import delimited "~/path/f31_os_pub.dat", clear

keep id parity agefbir booph meno

recode parity (-1=0) (0=1) (1=2) (2=3) (3=4) (4=5) (5=6)

label variable parity "Number of Term Pregnancies"
label define parityl 0 "Never pregnant" 1 "Never had term pregnancy" 2 "1" 3 "2" 4 "3" 5 "4" 6 "5+"
label values parity parityl


label variable agefbir "Age at First Birth"
label define agefbirl 0 "Never had term pregnancy" 1 "<20" 2 "20-29" 3 "30+" 
label values agefbir agefbirl

label variable booph "Bilateral Oophorectomy"
label define boophl 0 "No" 1 "Yes" 
label values booph boophl

label variable meno "Age at menopause"


save "~/path/stata_data/f31_os_l.dta", replace




/*******************************************************************************

	f2_os_pub

*******************************************************************************/



import delimited "~/path/f2_os_pub.dat", clear

keep id diabtrt mi stroke dvt colon_f2 endo_f2 skin_f2 melan_f2 othca10y brca_f2 diab

label variable diabtrt "Diabetes treated (pills or shots)"
label define diabtrtl 0 "No" 1 "Yes" 
label values diabtrt diabtrtl

label variable mi "MI ever"
label define mil 0 "No" 1 "Yes" 
label values mi mil

label variable stroke "Stroke ever"
label define strokel 0 "No" 1 "Yes" 
label values stroke strokel

label variable dvt "DVT ever"
label define dvtl 0 "No" 1 "Yes" 
label values dvt dvtl




save "~/path/stata_data/f2_os_l.dta", replace









/*******************************************************************************

	f30_os_pub

*******************************************************************************/



import delimited "~/path/f30_os_pub.dat", clear

keep id hyptpill hicholrp angina cabg hip55 numfalls cvd atrialfb cardrest aortican osteopor

label variable hyptpill "Pills for hypertension ever"
label define hyptpilll 0 "No" 1 "Yes" 
label values hyptpill hyptpilll

label variable hicholrp "High cholesterol requiring pills ever"
label define hicholrpl 0 "No" 1 "Yes" 
label values hicholrp hicholrpl

label variable angina "Angina ever"
label define anginal 0 "No" 1 "Yes" 
label values angina anginal

label variable cabg "Coronary bypass surgery ever"
label define cabgl 0 "No" 1 "Yes" 
label values cabg cabgl

label variable hip55 "Hip fracture age 55 or older"
label define hip55l 0 "No" 1 "Yes" 
label values hip55 hip55l

label variable numfalls "Times fell down last 12 months"
label define numfallsl 0 "None" 1 "1" 2 "2" 3 "3 or more time"
label values numfalls numfallsl


label variable atrialfb "Atrial fibrillation ever"
label define atrialfbl 0 "No" 1 "Yes" 
label values atrialfb atrialfbl

label variable cvd "Cardiovascular disease ever"
label define cvdl 0 "No" 1 "Yes" 
label values cvd cvdl


label variable cardrest "Cardiac arrest ever"
label define cardrestl 0 "No" 1 "Yes" 
label values cardrest cardrestl

label variable aortican "Aortic aneurysm ever"
label define aorticanl 0 "No" 1 "Yes" 
label values aortican aorticanl

label variable osteopor "Osteoporosis ever"
label define osteoporl 0 "No" 1 "Yes" 
label values osteopor osteoporl


save "~/path/stata_data/f30_os_l.dta", replace



/*******************************************************************************

	f60_os_pub

*******************************************************************************/



import delimited "~/path/f60_item_os_pub.dat", clear

keep id f60vtyp redmeat grndmeat steak bbqsand chili liver stew gravy spagmeat lnchham lnchbol fruits vegtabls

keep if  f60vtyp==1

label variable redmeat "Computed food group: Red meat, med serv/day"
label variable grndmeat "Item #44: Ground meat (hamburg/meatloaf/picadillo), med serv/day"
label variable steak "Item #45: Beef/pork/lamb m dish (steak/roast/ham), med serv/day"
label variable bbqsand "Item #46: Beef/pork/lamb as sandwich (steak/bbq), med serv/day"
label variable chili "Item #48: Chili with meat and beans, med serv/day"
label variable liver "Item #49: Liver (chicken liver & oth organ meat), med serv/day"
label variable stew "Item #47: Stew/potpie/casserole w/meat or chicken, med serv/day"
label variable gravy "Item #52: Gravy made w/meat dripping & wht sauce, med serv/day"
label variable spagmeat "Item #59: Spaghetti with meat sauce, med serv/day"
label variable lnchham "Item #71: Lunch meat (ham/turkey/lean meat), med serv/day"
label variable lnchbol "Item #72: Other lunch meat (bologna/salami/spam), med serv/day"

label variable fruits "Fruits, med serv/day"
label variable lnchbol "Vegetables, med serv/day"



save "~/path/stata_data/f60_item_os_pub.dta", replace


import delimited "~/path/f60_nutr_os_pub.dat", clear

keep if  f60vtyp==1

keep id f60enrgy f60enrgyj

label variable f60enrgy "Dietary Energy (kcal)"
label variable f60enrgyj "Dietary Energy (joules)"

save "~/path/stata_data/f60_nutr_os_pub.dta", replace



/*##############################################################################
##############################################################################

	MERGING

##############################################################################
##############################################################################*/




use "~/path/stata_data/outc_adj_os_fust_os.dta", clear


merge 1:1 id using "~/path/stata_data/f134_os_pub_l.dta"
drop if _merge!= 3
drop _merge


merge 1:1 id using "~/path/stata_data/f143_os_l.dta"
drop if _merge!= 3
drop _merge


merge 1:1 id using "~/path/stata_data/f43_os_l.dta"
drop if _merge!= 3
drop _merge

merge 1:1 id using "~/path/stata_data/f45_os_l.dta"
drop if _merge!= 3
drop _merge


merge 1:1 id using "~/path/stata_data/dem_os_l.dta"
drop if _merge!= 3
drop _merge

merge 1:1 id using "~/path/stata_data/f80_os_l.dta"
drop if _merge!= 3
drop _merge

merge 1:1 id using "~/path/stata_data/f80_os_l_bl.dta"
drop if _merge!= 3
drop _merge

merge 1:1 id using "~/path/stata_data/f31_os_l.dta"
drop if _merge!= 3
drop _merge

merge 1:1 id using "~/path/stata_data/f2_os_l.dta"
drop if _merge!= 3
drop _merge

merge 1:1 id using "~/path/stata_data/f30_os_l.dta"
drop if _merge!= 3
drop _merge

merge 1:1 id using "~/path/stata_data/f34_os_l.dta"
drop if _merge!= 3
drop _merge

merge 1:1 id using "~/path/stata_data/f60_item_os_pub.dta"
drop if _merge!= 3
drop _merge

merge 1:1 id using "~/path/stata_data/f60_nutr_os_pub.dta"
drop if _merge!= 3
drop _merge





gen time_since_menopause = age-meno

keep if time_since_menopause >=0


save "~/path/stata_data/final_os_l.dta", replace


use "~/path/stata_data/final_os_l.dta", clear


gen failure = 0
replace failure = 1 if chdsrc==1

gen failure_c = 0
replace failure_c = 1 if colonsrc==1

gen failure_cr = 0
replace failure_cr = 1 if colorectalsrc==1

gen HT_intervention = totp


saveold "~/path/stata_data/final_os_12.dta", replace



