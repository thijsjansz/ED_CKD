***** STATA CODE - I am using stata to select variables from the main UKBB dataset as this works faster than in R. 
***** I will then add primary care data and HES data and save this into a file I can further analyse using R.
***** before running this code, the code ED_CKD_mysql.R will have to be run.

****** Hospital episode statistics *********

foreach disease in  "ED""CKD" "DM" "HTN" "IHD" "CVD" "CHOL"{ 
clear 
set maxvar 64000

import delimited "/slade/projects/UKBB/phenotype_data/master/main_data/HES_data/hesin_diag150920.txt"
	

if `"`disease'"' == "ED" {
	local dis="ed"
	local icd10="F522 F529" /* F522 = failure of genital response, F529 unspecified sexual dysfunction not caused by organic disorder or disease; both are frequently used to code erectile dysfunction */
}
	
		if "`disease'" == "CKD" {
		local dis="ckd"
		local icd10 = "I120 I130 I131 I132 I139 N180 N188 N189 N19X Y841 Z491 Z492 Z992 T861 Z940 N165A Z490" /* the codes starting with I1 refer to hypertensive disease with renal failure; the codes starting with N1 refer to chronic renal failure. Y841 = kidney dialysis, Z491 = extracorporeal dialysis, Z492 = other dialysis, Z992 = dependence on renal dialysis, T861 = kidney transplant failure and rejection, Z940 = kidney transplant status, N165A = renal tubulointerstitial disorders in transplant rejection, Z490 = preparatory care for dialysis */
	}
	
		if "`disease'" == "DM" {
		local dis="diabetes"
		local icd10 = "O240 O241 O242 O243 E100 E101 E102 E103 E104 E105 E106 E107 E108 E109 E110 E111 E112 E113 E114 E115 E116 E117 E118 E119 E120 E121 E122 E123 E124 E125 E126 E127 E128 E129 E130 E131 E132 E133 E134 E135 E136 E137 E138 E139 E140 E141 E142 E143 E144 E145 E146 E147 E148 E149 G632A N083A" /* these ICD10 codes refer to diabetes - I have left out gestational diabetes */
	}
	
			if "`disease'" == "HTN" {
		local dis="hypertension"
		local icd10 = "I10X I150 I151 I152 I158 I159" /* I10X is primary hypertension, I15* is secondary hypertension */
	}
	
			if "`disease'" == "IHD" {
		local dis="ihd"
		local icd10 = "I200 I201 I208 I209 I210 I211 I212 I213 I214 I219 I220 I221 I228 I229 I230 I231 I232 I233 I234 I235 I236 I238 I241 I248 I249 I252 I255 I256 I258 I259" /* I20* angina, I21* acute MI, I22* further MI, I23* complications of MI, I24/I25* miscellaneous referring to ischaemic heart disease */
	}
	
			if "`disease'" == "CVD" {
		local dis="cvd"
		local icd10 = "I630 I631 I632 I633 I634 I635 I636 I638 I639 I64X G450 G451 G452 G453  G458 G459 G460A G461A G462A G463A G464A G465A G466A G467A G468A" /*  */
	}
	
			if "`disease'" == "CHOL" {
		local dis="chol"
		local icd10 = "E780 E782 E784 E785" /* pure hypercholesterolaemia, mixed, other, and unspecified hyperlipidaemia */
	}
	
	
	generate hes_flag_`dis' = .
	
	foreach i of local icd10 {
		display "Processing `i'"
		replace hes_flag_`dis' = 1 if regexm(diag_icd10,"^`i'.*")==1 /* Make flag to identify all records with code in user's list */
	}

	keep if hes_flag_`dis' == 1 

	save "/slade/home/tj358/pilot/`dis'_hes_intermediate_only.dta", replace 


***** Adding file contain date and selecting the earliest of dates *************************

clear 
import delimited "/slade/projects/UKBB/phenotype_data/master/main_data/HES_data/hesin150920.txt"


merge 1:m eid ins_index using "/slade/home/tj358/pilot/`dis'_hes_intermediate_only.dta"
keep if _merge == 3
keep if hes_flag_`dis' ==1

generate addate = date(admidate, "DMY") /* Generate admission date variable */
format addate %td /* Reformat the admission date to date format */
generate dischdate = date(disdate, "DMY") /* Generate discharge date variable */
format dischdate %td /* Reformat the discharge date to date format */
generate epistdate = date(epistart, "DMY") /* Generate epi start variable */
format epistdate %td /* Reformat the epi start to date format */
generate epienddate = date(epiend, "DMY") /* Generate epi end variable */
format epienddate %td /* Reformat the epi end to date format */

generate date_1st_`dis' = min(addate,dischdate,epistdate,epienddate) /* Find earliest of dates */
format date_1st_`dis' %td /* Reformat the epi end to date format */

sort eid date_1st_`dis' /*Order records by admission date*/
bysort eid: generate order=_n /* Generate a variable with the order of the admission date */
keep if order==1 /* Keep the first admission date */

*** Renaming variables for later merges with other files
rename diag_icd10 diag_icd10_`dis'
rename eid n_eid
rename dsource dsource_`dis'
rename mainspef mainspef_`dis'
rename mainspef_uni mainspef_uni_`dis'
rename carersi carersi_`dis'
rename date_1st_`dis' date_1st_hes_`dis'

keep n_eid date_1st_hes_`dis' diag_icd10_`dis' /* Keep relevant variables */

save "/slade/home/tj358/pilot/`dis'_hes_only.dta", replace
		
}

clear

******* GP records *********

foreach disease in  "ED""CKD" "DM" "HTN" "IHD" "CVD" "CHOL"{ 

set maxvar 64000

if "`disease'" == "ED" {
		local dis="ed"
	}
if "`disease'" == "CKD" {
		local dis="ckd"
	}
if "`disease'" == "DM" {
		local dis="diabetes"
	}
if "`disease'" == "HTN" {
		local dis="hypertension"
	}
if "`disease'" == "IHD" {
		local dis="ihd"
	}
if "`disease'" == "CVD" {
		local dis="cvd"
	}
if "`disease'" == "CHOL" {
		local dis="chol"
	}
import delimited using "/slade/home/tj358/pilot/`dis'_gpcases.tsv"

*/ date into readable format: */
gen event_date_`dis' = date(event_dt, "YMD")
format event_date_`dis' %td

sort  n_eid event_date_`dis'  /*Order records by admission date*/
bysort n_eid: generate order_gp = _n /* Generate a variable with the order of the admission date */
keep if order_gp==1 /* Keep earliest admission for each person */

drop value1 value2 value3 event_dt order_gp /* these values are empty and not needed */

rename read_2 read_2_`dis'_gp
rename read_3 read_3_`dis'_gp

save "/slade/home/tj358/pilot/`dis'_gp_only.dta", replace
clear

}

****** combining above with self-reported and first-occurrence data *******
clear

set maxvar 64000 

use n_eid ts_53_0_0 n_52_0_0 n_31_0_0 n_34_0_0 n_189_0_0 n_1647_* n_22011_0_0 n_22012_0_0 n_21001_0_0 n_20002_* n_20009_* n_21000_0_0 n_22001_0_0 n_20116_* s_40001_* s_40002_* ts_40000_* s_41204* n_22004_0_0 n_22005_0_0 n_22006_0_0 n_22010_0_0 using "/slade/local/UKBB/phenotype_data/master/main_data/raw_data_2019.dta"

* add cystatin C levels
merge 1:1 n_eid using "/slade/projects/UKBB/phenotype_data/master/main_data/biomarkers_2019.dta", keepusing(n_30720_0_0 ts_30721_0_0 n_30700_0_0 ts_30701_0_0)


***** add data on date of first occurrence
merge 1:1 n_eid using "/slade/local/UKBB/phenotype_data/master/main_data/ukb_first_occurrences.dta", keepusing(n_eid n_132031_0_0 ts_132030_0_0 n_132033_0_0 ts_132032_0_0 n_132035_0_0 ts_132034_0_0 n_130707_0_0 ts_130706_0_0 n_130709_0_0 ts_130708_0_0 n_130711_0_0 ts_130710_0_0 n_130713_0_0 ts_130712_0_0 n_130715_0_0 ts_130714_0_0 n_131287_0_0 ts_131286_0_0 n_131295_0_0 ts_131294_0_0 n_131297_0_0 ts_131296_0_0 n_131299_0_0 ts_131298_0_0 n_131301_0_0 ts_131300_0_0 n_131303_0_0 ts_131302_0_0 n_131305_0_0 ts_131304_0_0 n_131307_0_0 ts_131306_0_0 n_131367_0_0 ts_131366_0_0 n_131369_0_0 ts_131368_0_0 n_131057_0_0 ts_131056_0_0 n_131059_0_0 ts_131058_0_0) generate(first_occ_merge)
* There are no first occurrence variables for erectile dysfunction; the first 6 first occurrence variables (13203*) refer to renal failure, the others to diabetes.


drop if first_occ_merge ==2
drop first_occ_merge

***** add phenotypes from GP and HES records
foreach disease in  "ED""CKD" "DM" "HTN" "IHD" "CVD" "CHOL"{ 

if "`disease'" == "ED" {
		local dis="ed"
	}
if "`disease'" == "CKD" {
		local dis="ckd"
	}
if "`disease'" == "DM" {
		local dis="diabetes"
	}
if "`disease'" == "HTN" {
		local dis="hypertension"
	}
if "`disease'" == "IHD" {
		local dis="ihd"
	}
if "`disease'" == "CVD" {
		local dis="cvd"
	}
if "`disease'" == "CHOL" {
		local dis="chol"
	}
merge 1:1 n_eid using "/slade/home/tj358/pilot/`dis'_gp_only.dta", generate(`dis'_gp_merge)
drop if `dis'_gp_merge ==2
drop `dis'_gp_merge

merge 1:1 n_eid using /slade/home/tj358/pilot/`dis'_hes_only.dta, generate (`dis'_hes_merge)
drop if `dis'_hes_merge ==2 
drop `dis'_hes_merge
} 

*** Renaming variables

rename n_31_0_0 sex
rename n_22001_0_0 genetic_sex
rename ts_53_0_0 enrol_date
rename n_34_0_0 YOB
rename n_52_0_0 MOB
rename n_21001_0_0 BMI
rename n_22010_0_0 exclusion
rename n_22006_0_0 white_british
rename n_20116_0_0 smoking_status
rename n_21000_0_0 ethnic_background
rename n_189_0_0 TDI

rename n_20116_1_0 smoking_status_1
rename n_20116_2_0 smoking_status_2

rename n_1647_0_0 country_of_birth

** Dropping participants 
*drop if exclusion ==1 *we do not have to drop these participants per se if we are not doing genetic analyses yet
drop if sex ==. & genetic_sex==.

* Sorting by genetic relatedeness pairings and dropping one of each pair if kinship coefficient > 0.0844
sort n_22011_0_0, stable
by n_22011_0_0: generate order_kin= _n
*drop if order_kin >= 2 & n_22012_0_0 > 0.0844 & n_22011_0_0!=. * same here - only consider dropping if we are going to do genetic analyses
drop order_kin 


* removing participants who have withdrawn from the UKBB 
run "/slade/local/UKBB/phenotype_data/scripts/do_files/withdrawn_participants.do"
drop if withdrawn ==1

** Calculating dummy DOB which will be used in date analysis

gen DOB = mdy(MOB, 15, YOB)
format DOB %td
drop MOB YOB

generate enrol_age = datediff(DOB, enrol_date, "m")


* save file to be analysed in R instead:
drop country_of_birth n_1647_1_0 n_1647_2_0 s_40002_0_0 s_40002_0_1 s_40002_0_2 s_40002_0_3 s_40002_0_4 s_40002_0_5 s_40002_0_6 s_40002_0_7 s_40002_0_8 s_40002_0_9 s_40002_0_10 s_40002_0_11 s_40002_0_12 s_40002_0_13 s_40002_1_0 s_40002_1_1 s_40002_1_2 s_40002_1_3 s_40002_1_4 s_40002_1_5 s_40002_1_6 s_41204_0_0 s_41204_0_1 s_41204_0_2 s_41204_0_3 s_41204_0_4 s_41204_0_5 s_41204_0_6 s_41204_0_7 s_41204_0_8 s_41204_0_9 s_41204_0_10 s_41204_0_11 s_41204_0_12 s_41204_0_13 s_41204_0_14 s_41204_0_15 s_41204_0_16 s_41204_0_17 s_41204_0_18 s_41204_0_19 s_41204_0_20 s_41204_0_21 s_41204_0_22 s_41204_0_23 s_41204_0_24 s_41204_0_25 s_41204_0_26 s_41204_0_27 s_41204_0_28 s_41204_0_29 s_41204_0_30 s_41204_0_31 s_41204_0_32 s_41204_0_33 s_41204_0_34 s_41204_0_35 s_41204_0_36 s_41204_0_37 s_41204_0_38 s_41204_0_39 s_41204_0_40 s_41204_0_41 s_41204_0_42 s_41204_0_43 s_41204_0_44 s_41204_0_45 s_41204_0_46 s_41204_0_47 s_41204_0_48 s_41204_0_49 s_41204_0_50 s_41204_0_51 s_41204_0_52 s_41204_0_53 s_41204_0_54 s_41204_0_55 s_41204_0_56 s_41204_0_57 s_41204_0_58 s_41204_0_59 s_41204_0_60 s_41204_0_61 s_41204_0_62 s_41204_0_63 s_41204_0_64 s_41204_0_65 s_41204_0_66 s_41204_0_67 s_41204_0_68 s_41204_0_69 s_41204_0_70 s_41204_0_71 s_41204_0_72 s_41204_0_73 s_41204_0_74 s_41204_0_75 s_41204_0_76 s_41204_0_77 s_41204_0_78 s_41204_0_79 s_41204_0_80 s_41204_0_81 s_41204_0_82 s_41204_0_83 s_41204_0_84 s_41204_0_85 s_41204_0_86 s_41204_0_87 s_41204_0_88 s_41204_0_89 s_41204_0_90 s_41204_0_91 s_41204_0_92 s_41204_0_93 s_41204_0_94 s_41204_0_95 s_41204_0_96 s_41204_0_97 s_41204_0_98 s_41204_0_99 s_41204_0_100 s_41204_0_101 s_41204_0_102 s_41204_0_103 s_41204_0_104 s_41204_0_105 s_41204_0_106 s_41204_0_107 s_41204_0_108 s_41204_0_109 s_41204_0_110 s_41204_0_111 s_41204_0_112 s_41204_0_113 s_41204_0_114 s_41204_0_115 s_41204_0_116 s_41204_0_117 s_41204_0_118 s_41204_0_119 s_41204_0_120 s_41204_0_121 s_41204_0_122 s_41204_0_123 s_41204_0_124 s_41204_0_125 s_41204_0_126 s_41204_0_127 s_41204_0_128 s_41204_0_129 s_41204_0_130 s_41204_0_131 s_41204_0_132 s_41204_0_133 s_41204_0_134 s_41204_0_135 s_41204_0_136 s_41204_0_137 s_41204_0_138 s_41204_0_139 s_41204_0_140 s_41204_0_141 s_41204_0_142 s_41204_0_143 s_41204_0_144 s_41204_0_145 s_41204_0_146 s_41204_0_147 s_41204_0_148 s_41204_0_149 s_41204_0_150 s_41204_0_151 s_41204_0_152 s_41204_0_153 s_41204_0_154 s_41204_0_155 s_41204_0_156 s_41204_0_157 s_41204_0_158 s_41204_0_159 s_41204_0_160 s_41204_0_161 s_41204_0_162 s_41204_0_163 s_41204_0_164 s_41204_0_165 s_41204_0_166 s_41204_0_167 s_41204_0_168 s_41204_0_169 s_41204_0_170 s_41204_0_171 s_41204_0_172 s_41204_0_173 s_41204_0_174 s_41204_0_175 s_41204_0_176 s_41204_0_177 s_41204_0_178 s_41204_0_179 s_41204_0_180 s_41204_0_181 s_41204_0_182 s_41204_0_183 genetic_sex n_22004_0_0 n_22005_0_0 exclusion n_22011_0_0 n_22012_0_0 

save "/slade/home/tj358/pilot/R_ed_ckd_data.dta", replace
