infix using "C:\HAITISPA_FIN\ARCHIVE\HTSC6BFLSP.dct"

label variable facil    "Facility identification"
label variable provno   "Provider serial number"
label variable cltype   "Client type"
label variable clno     "Client number"
label variable cfacil   "Facility number"
label variable cprovcat "Provider category"
label variable cprovsex "Provider sex"
label variable cprovno  "Provider serial number"
label variable serv     "Service where observation occured (only for STI)"
label variable cday     "Client day"
label variable cmonth   "Client month"
label variable cyear    "Client year"
label variable ointv    "Observer code"
label variable c100     "Provider agreement"
label variable c101     "Client agreement"
label variable c102h    "Hour begin"
label variable c102m    "Minute begin"
label variable c103     "First observation for this provider of this type of client"
label variable clientwt "client weight"
label variable bothoe   "Both observation and exit completed"
label variable oc104    "Sex of child"
label variable oc105    "Client history"
label variable oc106    "General danger signs"
label variable oc107    "Symptomatic HIV suspected"
label variable oc108    "Physical exam"
label variable oc109    "Other assessments"
label variable oc110    "Counseling of caretaker"
label variable oc111    "Additional counseling"
label variable oc112    "Referrals and admissions"
label variable oc113    "Outcome of consultation"
label variable oc201    "Dehydration dx"
label variable oc202    "Respiratory dx"
label variable oc203    "Digestive system dx"
label variable oc204    "Malaria dx"
label variable oc205    "Fever/ measles dx"
label variable oc206    "Ear dx"
label variable oc207    "Throat dx"
label variable oc208    "Other dx"
label variable oc209    "Treatment prescribed"
label variable oc210    "General treatment"
label variable oc211    "Respiratory tx"
label variable oc212    "Malaria tx"
label variable oc213    "Dehydration tx"
label variable oc214    "Other tx/ advice"
label variable oc215    "First visit for child to facility"
label variable oc216    "Child immunized during visit"
label variable oc217h   "Ending  hour"
label variable oc217m   "Ending minute"
label variable xcday    "Client Day"
label variable xcmonth  "Client Month"
label variable xcyear   "Client Year"
label variable xcintv   "Interviewer code"
label variable xccode   "Client code"
label variable xc100    "Permission to intervIew"
label variable xc101h   "Time begins hour"
label variable xc101m   "Time begins minutes"
label variable xc102    "Nom de l'enfant malade"
label variable xc103d   "Day of birth"
label variable xc103m   "Month of birth"
label variable xc103y   "Year of birth"
label variable xc104    "Age in months"
label variable xc105    "Child had fever with this illness in past 2 days"
label variable xc106    "Child had convulsions with this illness"
label variable xc107    "Child has cough or difficult breathing with this illness"
label variable xc108    "Child can drink, eat and/or breastfeed"
label variable xc109    "Child vomits everything when eats or breastfeeds during this illness"
label variable xc110    "Child has had watery and frequent stools with this illness in past 2 days"
label variable xc111    "Child been excessively sleepy during illness"
label variable xc112    "Other reasons for visit"
label variable xc113    "Most recent visit for same illness"
label variable xc114    "Days ago illness began"
label variable xc115    "Provider gave diagnosis"
label variable xc116    "Action if no improvement"
label variable xc117    "Signs/symptoms for immediate return"
label variable xc118    "Provider mentioned reasons for follow-up visit"
label variable xc119    "Provider gave/prescribed medications for home"
label variable xc120    "Description of medications/prescriptions"
label variable xc121    "Provider gave instructions for giving medicines"
label variable xc122    "Confident knows how much medicine to give each day and how many days to give"
label variable xc123    "Child given 1st dose of any medicines at facility"
label variable xc124    "Child given an injection for sickness today"
label variable xc125    "Child weighed"
label variable xc126    "Child's growth mentioned during visit"
label variable xc127    "Provider asked about foods/amts normally given to child"
label variable xc128    "Provider directions regarding solids"
label variable xc129    "Provider directions regarding fluids"
label variable xc130    "Child given vaccination during visit"
label variable xc131    "Has vaccination card"
label variable xc132    "Vaccination card seen"
label variable xc134    "Provider instructed to go to other provider/lab/facility for further care"
label variable xc135_1  "Referral information"
label variable xc135_2  "Referral information"
label variable xc135_3  "Referral information"
label variable xc135_4  "Referral information"
label variable xc136    "Visited other provider or traditional healer before this visit"
label variable x201     "Time waited for provider"
label variable x202_01  "Problems encountered today"
label variable x202_02  "Problems encountered today"
label variable x202_03  "Problems encountered today"
label variable x202_04  "Problems encountered today"
label variable x202_05  "Problems encountered today"
label variable x202_06  "Problems encountered today"
label variable x202_07  "Problems encountered today"
label variable x202_08  "Problems encountered today"
label variable x202_09  "Problems encountered today"
label variable x202_10  "Problems encountered today"
label variable x202_11  "Problems encountered today"
label variable x203     "Part of prepayment plan"
label variable x204     "Fees for services for visit"
label variable x205     "Total fees"
label variable x206     "Nearest facillity to home"
label variable x207     "Main reason didn't go to nearest facility"
label variable x208     "Opinion of services given today"
label variable x209     "Would recommend to friend/family member"
label variable x301     "Relationship to child"
label variable x302     "Respondent's age"
label variable x303     "Ever attended school"
label variable x304     "Highest level attended"
label variable x305     "Literate"
label variable x306h    "End INTERVIew Hour"
label variable x306m    "End INTERVIew Minutes"

#delimit ;
label define CPROVCAT
     1 "Medical doctor (generalist)"
     2 "Surgeon (general)"
     3 "Medical specialist"
     4 "Lab technician"
     5 "Radiology technician"
     6 "Medical imagery technician"
     7 "Nurse"
     8 "Nurse/midwife"
     9 "Auxiliary nurse"
    12 "Dental lab technician"
    13 "Dental hygienist"
    14 "Auxiliary dentist"
    15 "Other community health workers"
    95 "Non-technical qualification"
    96 "Other"
;
label define CPROVSEX
     1 "Male"
     2 "Female"
;
label define SERV
     1 "ANC"
     2 "FP"
     3 "Sick child"
     4 "STI"
;
label define CMONTH
     1 "January"
     2 "February"
     3 "March"
     4 "April"
     5 "May"
     6 "June"
     7 "July"
     8 "August"
     9 "September"
    10 "October"
    11 "November"
    12 "December"
;
label define OINTV
     1 "Paul Brea"
     2 "Claudia Marchena"
     3 "Keith Purvis"
     4 "Albert Themme"
    10 "Team 1 superviseur 10"
    11 "Team 1 interviewer 11"
    12 "Team 1 interviewer 12"
    13 "Team 1 interviewer 13"
    14 "Team 1 interviewer 14"
    15 "Team 2 superviseur 15"
    16 "Team 2 interviewer 16"
    17 "Team 2 interviewer 17"
    18 "Team 2 interviewer 18"
    19 "Team 2 interviewer 19"
    20 "Team 3 superviseur 20"
    21 "Team 3 interviewer 21"
    22 "Team 3 interviewer 22"
    23 "Team 3 interviewer 23"
    24 "Team 3 interviewer 24"
    25 "Team 4 supervisor 25"
    26 "Team 4 interviewer 26"
    27 "Team 4 interviewer 27"
    28 "Team 4 interviewer 28"
    29 "Team 4 interviewer 29"
    30 "Team 5 supervisor 30"
    31 "Team 5 interviewer 31"
    32 "Team 5 interviewer 32"
    33 "Team 5 interviewer 33"
    34 "Team 5 interviewer 34"
    35 "Team 6 supervisor 35"
    36 "Team 6 interviewer 36"
    37 "Team 6 interviewer 37"
    38 "Team 6 interviewer 38"
    39 "Team 6 interviewer 39"
    40 "Team 7 supervisor 40"
    41 "Team 7 interviewer 41"
    42 "Team 7 interviewer 42"
    43 "Team 7 interviewer 43"
    44 "Team 7 interviewer 44"
    45 "Team 8 supervisor 45"
    46 "Team 8 interviewer 46"
    47 "Team 8 interviewer 47"
    48 "Team 8 interviewer 48"
    49 "Team 8 interviewer 49"
    50 "Team 9 supervisor 50"
    51 "Team 9 interviewer 51"
    52 "Team 9 interviewer 52"
    53 "Team 9 interviewer 53"
    54 "Team 9 interviewer 54"
    55 "Team 10 supervisor 55"
    56 "Team 10 interviewer 56"
    57 "Team 10 interviewer 57"
    58 "Team 10 interviewer 58"
    59 "Team 10 interviewer 59"
    60 "Team 11 supervisor 60"
    61 "Team 11 interviewer 61"
    62 "Team 11 interviewer 62"
    63 "Team 11 interviewer 63"
    64 "Team 11 interviewer 64"
    65 "Team 12 supervisor 65"
    66 "Team 12 interviewer 66"
    67 "Team 12 interviewer 67"
    68 "Team 12 interviewer 68"
    69 "Team 12 interviewer 69"
    70 "Team 13 supervisor 70"
    71 "Team 13 interviewer 71"
    72 "Team 13 interviewer 72"
    73 "Team 13 interviewer 73"
    74 "Team 13 interviewer 74"
    75 "Team 14 supervisor 75"
    76 "Team 14 interviewer 76"
    77 "Team 14 interviewer 77"
    78 "Team 14 interviewer 78"
    79 "Team 14 interviewer 79"
    80 "Team 15 supervisor 80"
    81 "Team 13 interviewer 81"
    82 "Team 15 interviewer 82"
    83 "Team 15 interviewer 83"
    84 "Team 15 interviewer 84"
    85 "Team 16 supervisor 85"
    86 "Team 16 interviewer 86"
    87 "Team 16 interviewer 87"
    88 "Team 16 interviewer 88"
    89 "Team 16 interviewer 89"
;
label define C100
     1 "Yes"
     2 "No"
;
label define C101
     1 "Yes"
     2 "No"
;
label define C103
     1 "Yes, first observation of this type of client for this provider"
     2 "No, not first of this type for this provider"
;
label define BOTHOE
     0 "No piece completed"
     1 "Both observation and exit completed"
     2 "Only observation completed"
     3 "Only exit completed"
;
label define OC104
     1 "Male"
     2 "Female"
;
label define OC113
     1 "Treated and sent home"
     2 "Child referred to provider same facility"
     3 "Child admitted to same facilty"
     4 "Child sent to lab"
     5 "Child referred to other facilty"
;
label define OC201
     1 "Severe dehydration"
     2 "Moderate dehydration"
     3 "Mild dehydration"
     4 "None of the above"
;
label define OC204
     1 "Malaria (clinical diagnosis)"
     2 "Malaria (blood smear)"
     3 "Malaria (rapid test)"
     4 "None of the above"
;
label define OC205
     1 "Fever of unknown origin"
     2 "Measles with no complications"
     3 "Measles with complicatins (mouth/eye or severe)"
     4 "None of the above"
;
label define OC207
     1 "Sore throat"
     2 "Other throat diagnosis (specify)"
     3 "None of the above"
;
label define OC208
     1 "Yes, other diagnosis (specify)"
     2 "No other diagnosis"
;
label define OC209
     1 "Yes"
     2 "No"
;
label define OC215
     1 "First visit"
     2 "Follow-up"
     8 "DK"
;
label define OC216
     1 "Yes, vaccinated child"
     2 "Yes, referred"
     3 "Not due for or completed vaccination"
     4 "Vaccine not available"
     5 "Child too sick"
     6 "Not day for vaccination"
     7 "Did not check for vaccination"
;
label define XCMONTH
     1 "January"
     2 "February"
     3 "March"
     4 "April"
     5 "May"
     6 "June"
     7 "July"
     8 "August"
     9 "September"
    10 "October"
    11 "November"
    12 "December"
;
label define XCINTV
     1 "Paul Brea"
     2 "Claudia Marchena"
     3 "Keith Purvis"
     4 "Albert Themme"
    10 "Team 1 superviseur 10"
    11 "Team 1 interviewer 11"
    12 "Team 1 interviewer 12"
    13 "Team 1 interviewer 13"
    14 "Team 1 interviewer 14"
    15 "Team 2 superviseur 15"
    16 "Team 2 interviewer 16"
    17 "Team 2 interviewer 17"
    18 "Team 2 interviewer 18"
    19 "Team 2 interviewer 19"
    20 "Team 3 superviseur 20"
    21 "Team 3 interviewer 21"
    22 "Team 3 interviewer 22"
    23 "Team 3 interviewer 23"
    24 "Team 3 interviewer 24"
    25 "Team 4 supervisor 25"
    26 "Team 4 interviewer 26"
    27 "Team 4 interviewer 27"
    28 "Team 4 interviewer 28"
    29 "Team 4 interviewer 29"
    30 "Team 5 supervisor 30"
    31 "Team 5 interviewer 31"
    32 "Team 5 interviewer 32"
    33 "Team 5 interviewer 33"
    34 "Team 5 interviewer 34"
    35 "Team 6 supervisor 35"
    36 "Team 6 interviewer 36"
    37 "Team 6 interviewer 37"
    38 "Team 6 interviewer 38"
    39 "Team 6 interviewer 39"
    40 "Team 7 supervisor 40"
    41 "Team 7 interviewer 41"
    42 "Team 7 interviewer 42"
    43 "Team 7 interviewer 43"
    44 "Team 7 interviewer 44"
    45 "Team 8 supervisor 45"
    46 "Team 8 interviewer 46"
    47 "Team 8 interviewer 47"
    48 "Team 8 interviewer 48"
    49 "Team 8 interviewer 49"
    50 "Team 9 supervisor 50"
    51 "Team 9 interviewer 51"
    52 "Team 9 interviewer 52"
    53 "Team 9 interviewer 53"
    54 "Team 9 interviewer 54"
    55 "Team 10 supervisor 55"
    56 "Team 10 interviewer 56"
    57 "Team 10 interviewer 57"
    58 "Team 10 interviewer 58"
    59 "Team 10 interviewer 59"
    60 "Team 11 supervisor 60"
    61 "Team 11 interviewer 61"
    62 "Team 11 interviewer 62"
    63 "Team 11 interviewer 63"
    64 "Team 11 interviewer 64"
    65 "Team 12 supervisor 65"
    66 "Team 12 interviewer 66"
    67 "Team 12 interviewer 67"
    68 "Team 12 interviewer 68"
    69 "Team 12 interviewer 69"
    70 "Team 13 supervisor 70"
    71 "Team 13 interviewer 71"
    72 "Team 13 interviewer 72"
    73 "Team 13 interviewer 73"
    74 "Team 13 interviewer 74"
    75 "Team 14 supervisor 75"
    76 "Team 14 interviewer 76"
    77 "Team 14 interviewer 77"
    78 "Team 14 interviewer 78"
    79 "Team 14 interviewer 79"
    80 "Team 15 supervisor 80"
    81 "Team 13 interviewer 81"
    82 "Team 15 interviewer 82"
    83 "Team 15 interviewer 83"
    84 "Team 15 interviewer 84"
    85 "Team 16 supervisor 85"
    86 "Team 16 interviewer 86"
    87 "Team 16 interviewer 87"
    88 "Team 16 interviewer 88"
    89 "Team 16 interviewer 89"
;
label define XC100
     1 "Yes, client agrees"
     2 "No, client refuses"
;
label define XC103D
    98 "NSP"
;
label define XC103M
     1 "January"
     2 "February"
     3 "March"
     4 "April"
     5 "May"
     6 "June"
     7 "July"
     8 "August"
     9 "September"
    10 "October"
    11 "November"
    12 "December"
    98 "NSP"
;
label define XC103Y
  9998 "Don't know"
;
label define XC104
    98 "Don't know"
;
label define XC105
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC106
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC107
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC108
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC109
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC110
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC111
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC113
     1 "Within past week"
     2 "Winth past 2-4 weeks"
     3 "More than 4 weeks ago"
     4 "No, child has not been to facility for this illness"
     8 "Don't know"
;
label define XC114
    98 "Don't know"
;
label define XC115
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC116
     1 "Return to facility"
     2 "Go to other facility"
     3 "Go to other health worker/pharmacy"
     4 "Go to traditional healer"
     5 "Nothing, just wait"
     8 "Don't know"
;
label define XC119
     1 "Yes, gave meds"
     2 "Yes, gave prescriptions"
     3 "Gave meds and prescriptions"
     4 "No"
;
label define XC120
     1 "Has all meds"
     2 "Has some meds, some unfilled prescriptions"
     3 "No medications seen, has only prescriptions"
;
label define XC121
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC122
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC123
     1 "Yes"
     2 "No"
;
label define XC124
     1 "Yes, received injection"
     2 "Yes, received prescription for injection"
     3 "No"
     8 "Don't know"
;
label define XC125
     1 "Yes"
     2 "No"
;
label define XC126
     1 "Yes"
     2 "No"
;
label define XC127
     1 "Yes"
     2 "No"
     8 "Can't remember"
;
label define XC128
     1 "Give less than usual"
     2 "Give same as usual"
     3 "Give more than usual"
     4 "Give nothing/don't feed"
     6 "Didn't discuss"
     8 "Don't know/not certain"
;
label define XC129
     1 "Give less than usual"
     2 "Give same as usual"
     3 "Give more than usual"
     4 "Give nothing/don't feed"
     6 "Didn't discuss"
     8 "Don't know/not certain"
;
label define XC130
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XC131
     1 "Yes"
     2 "No"
;
label define XC132
     1 "Yes"
     2 "No"
;
label define XC134
     1 "Yes"
     2 "No"
;
label define XC135
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define X201
     0 "Saw provider immediately"
   998 "Don't know"
;
label define X202
     1 "Major"
     2 "Minor"
     3 "No problem"
     8 "Don't know"
;
label define X203
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define X204
     1 "Yes"
     2 "No"
;
label define X205
999998 "Don't know"
;
label define X206
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define X207
     1 "Inconvenient operating hours"
     2 "Bad reputation"
     3 "Don't like personnel"
     4 "No medicine"
     5 "Prefers to remain anonymous"
     6 "It is more expensive"
     7 "Was referred"
    96 "Other"
    98 "Don't know"
;
label define X208
     1 "Very satisfied"
     2 "More or less satisfied"
     3 "Not satisfied"
;
label define X209
     1 "Yes"
     2 "No"
     8 "DK"
;
label define X301
     1 "Mother"
     2 "Father"
     3 "Sibling"
     4 "Aunt/Uncle"
     5 "Grand mom/ grand dad"
     6 "Other (specify)"
;
label define X302
    98 "Don't know"
;
label define X303
     1 "Yes"
     2 "No"
;
label define X304
     1 "Primary"
     2 "Post-primary/ vocational"
     3 "Secondary/ A-Level"
     4 "College (middle level)/University"
;
label define X305
     1 "Yes, read and write"
     2 "Yes, read only"
     3 "No"
;

#delimit cr
label values cprovcat CPROVCAT
label values cprovsex CPROVSEX
label values serv     SERV
label values cmonth   CMONTH
label values ointv    OINTV
label values c100     C100
label values c101     C101
label values c103     C103
label values bothoe   BOTHOE
label values oc104    OC104
label values oc113    OC113
label values oc201    OC201
label values oc204    OC204
label values oc205    OC205
label values oc207    OC207
label values oc208    OC208
label values oc209    OC209
label values oc215    OC215
label values oc216    OC216
label values xcmonth  XCMONTH
label values xcintv   XCINTV
label values xc100    XC100
label values xc103d   XC103D
label values xc103m   XC103M
label values xc103y   XC103Y
label values xc104    XC104
label values xc105    XC105
label values xc106    XC106
label values xc107    XC107
label values xc108    XC108
label values xc109    XC109
label values xc110    XC110
label values xc111    XC111
label values xc113    XC113
label values xc114    XC114
label values xc115    XC115
label values xc116    XC116
label values xc119    XC119
label values xc120    XC120
label values xc121    XC121
label values xc122    XC122
label values xc123    XC123
label values xc124    XC124
label values xc125    XC125
label values xc126    XC126
label values xc127    XC127
label values xc128    XC128
label values xc129    XC129
label values xc130    XC130
label values xc131    XC131
label values xc132    XC132
label values xc134    XC134
label values xc135_1  XC135
label values xc135_2  XC135
label values xc135_3  XC135
label values xc135_4  XC135
label values x201     X201
label values x202_01  X202
label values x202_02  X202
label values x202_03  X202
label values x202_04  X202
label values x202_05  X202
label values x202_06  X202
label values x202_07  X202
label values x202_08  X202
label values x202_09  X202
label values x202_10  X202
label values x202_11  X202
label values x203     X203
label values x204     X204
label values x205     X205
label values x206     X206
label values x207     X207
label values x208     X208
label values x209     X209
label values x301     X301
label values x302     X302
label values x303     X303
label values x304     X304
label values x305     X305
