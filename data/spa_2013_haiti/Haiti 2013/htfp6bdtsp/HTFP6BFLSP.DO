infix using "C:\HAITISPA_FIN\ARCHIVE\HTFP6BFLSP.dct"

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
label variable of104    "Sex of client"
label variable of105    "Client history (female clients only)"
label variable of106    "Client history (all clients)"
label variable of107    "Physical exam"
label variable of108    "Partners & STIs"
label variable of109    "Questions & concerns"
label variable of110    "Privacy/ confidentiality"
label variable of111a   "Methods prescribed"
label variable of111b   "Methods provided in facility"
label variable of113    "Questions for pills/injections"
label variable of115    "Questions for condoms"
label variable of117    "Questions for IUDs"
label variable of119    "Questions for Implants"
label variable of121    "Questions for PA or SDM"
label variable of123    "Questions for LAM"
label variable of125    "Questions for vasectomy"
label variable of127    "Questions for female sterilization"
label variable of129    "Questions for emergency contraception"
label variable of130    "Additional provider actions"
label variable of131a   "Client has had previous contact with provider in this facility"
label variable of131b   "Client has been pregnant"
label variable of201    "Clinical procedure"
label variable of202    "Same provider as consult"
label variable of203    "Permission from new provider"
label variable of204    "Provider type"
label variable of205    "Sex of provider"
label variable of206    "Pelvic exam conducted"
label variable of207    "Actions before procedure"
label variable of208    "Actions during procedure"
label variable of209    "Actions after procedure"
label variable of210    "IUD insertion or removal conducted"
label variable of211    "Actions before procedure"
label variable of212    "Actions during procedure"
label variable of213    "Actions after procedure"
label variable of214    "Client-provider interaction"
label variable of215    "Injectable contraceptive given"
label variable of216    "Actions before procedure"
label variable of217    "Actions during procedure"
label variable of219    "Needle/syringe provision NOW 218 on quest."
label variable of220    "Implant inserted or removed"
label variable of221    "Actions before procedure"
label variable of222    "Actions during procedure"
label variable of223    "Actions after procedure"
label variable of224    "Client-provider interaction"
label variable of225    "Needle/ syringe provision"
label variable of301    "FP status of client at beginning of consult"
label variable of302    "Client's principal reason for visit"
label variable of303    "Outcome of visit for FP user"
label variable of304    "Outcome of visit for non-user"
label variable of305    "Client left with method"
label variable of306    "Provider wrote on card after consult"
label variable of307h   "Hour consult ended"
label variable of307m   "Minute consult ended"
label variable xfpday   "Client Day"
label variable xfpmonth "Client Month"
label variable xfpyear  "Client Year"
label variable xfpintv  "Interviewer code"
label variable xfpccode "Client code"
label variable xf100    "Permission to continue"
label variable xf101h   "Time begins hour"
label variable xf101m   "Time begins minutes"
label variable xf102    "Sex of client"
label variable xf103    "Using FP when visited on day of survey"
label variable xf104    "Used FP in last 6 months"
label variable xf105    "Method using"
label variable xf106    "Provider asked about problems"
label variable xf107    "Client has had problems"
label variable xf107a   "Which problems have you had with the method?"
label variable xf108    "Client mentioned problem to provider"
label variable xf109    "Provider suggested action to resolve problem"
label variable xf110    "Outcome of visit"
label variable xf111    "Thought about switching before visit"
label variable xf112    "Thought about FP method to use before visit"
label variable xf113    "Method considered"
label variable xf114    "Provider discussed method in XF113"
label variable xf115    "FP methods provider discussed"
label variable xf116a   "FP method prescribed"
label variable xf116b   "FP method received"
label variable xf117_1  "Provider actions during consult"
label variable xf117_2  "Provider actions during consult"
label variable xf117_3  "Provider actions during consult"
label variable xf117_4  "Provider actions during consult"
label variable xf118a   "Question for pill"
label variable xf118b   "Question for male condom"
label variable xf118c   "Question for female condom"
label variable xf118d   "Question for IUCD"
label variable xf118e   "Question for Progestin injectable (Depo, 2-3 monthly)"
label variable xf118f   "Question for monthly injectable"
label variable xf118g   "Question for Implant"
label variable xf118h   "Question for rhythm"
label variable xf118i   "Question for vasectomy (male sterilization)"
label variable xf118j   "Question for tubal ligation (female sterilization)"
label variable xf118k   "Question for LAM"
label variable xf119    "Method protects against STIs/HIV"
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
label define OF104
     1 "Male"
     2 "Female"
;
label define OF131A
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define OF131B
     1 "Yes"
     2 "No"
     3 "Male client"
     8 "Don't know"
;
label define OF202
     1 "Yes"
     2 "No"
;
label define OF203
     1 "Yes"
     2 "No"
;
label define OF204
     1 "Generalist medical doctor"
     2 "Specialist medical doctor"
     3 "Nursing professional"
     4 "Auxiliary nurse"
     5 "Midwife professional"
     6 "Biologist"
     7 "Medical technologist"
     8 "Lab technician"
     9 "Other community health worker"
    96 "Other"
;
label define OF205
     1 "Male"
     2 "Female"
;
label define OF206
     1 "Yes"
     2 "No"
;
label define OF215
     1 "Yes"
     2 "No"
;
label define OF219
     1 "Provided by facility"
     2 "Provided by client"
     8 "Don't know"
;
label define OF225
     1 "Provided by facility"
     2 "Provided by client"
     8 "Don't know"
;
label define OF301
     1 "Current user"
     2 "Nonuser, used in past"
     3 "Nonuser, no past use"
     8 "Not determined"
;
label define OF302
     1 "Resupply or routine follow-up"
     2 "Discuss problem with method"
     3 "Desire to change method (no problem)"
     4 "Desire to discontinue method (no problem)"
     5 "Discuss other problem"
;
label define OF303
     1 "Continued with current method"
     2 "Switched method"
     3 "Planned method switch, not received today, continued use of current method"
     4 "Planned method switch, not received today, discontinued current method"
     5 "Decided to stop using family planning"
;
label define OF304
     1 "Accepted to start method"
     2 "Did not decide on method"
;
label define OF305
     1 "Yes, left with method"
     2 "No, method not in stock"
     3 "No, requires appointment"
     4 "No, delay receiving due to health problem"
     5 "No, pregnancy status uncertain"
     6 "Other"
;
label define OF306
     1 "Yes"
     2 "No"
     3 "No individual card used"
     8 "Don't know"
;
label define OF307H
    97 "inconsistent"
;
label define OF307M
    97 "inconsistent"
;
label define XFPMONTH
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
label define XFPINTV
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
label define XF100
     1 "Yes, client agrees"
     2 "No, client refuses"
;
label define XF102
     1 "Male"
     2 "Female"
;
label define XF103
     1 "Yes"
     2 "No"
;
label define XF104
     1 "Yes"
     2 "No"
;
label define XF106
     1 "Yes, asked"
     2 "No, did not asked"
;
label define XF107
     1 "Yes"
     2 "No"
;
label define XF108
     1 "Yes"
     2 "No"
;
label define XF109
     1 "Yes"
     2 "No"
;
label define XF110
     1 "Continue or restart same method"
     2 "Switch methods"
     3 "Stop using method (due to problems)"
     4 "Stop using method (elective-no problems)"
;
label define XF111
     1 "Yes"
     2 "No"
;
label define XF112
     1 "Yes"
     2 "No"
;
label define XF114
     1 "Yes"
     2 "No"
;
label define XF117
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XF118A
     1 "Once a day"
     2 "Other"
     8 "Don't know"
;
label define XF118B
     1 "Once"
     2 "Other"
     8 "Don't know"
;
label define XF118C
     1 "Any oil or lubricant"
     2 "Other"
     8 "Don't know"
;
label define XF118D
     1 "Check string"
     2 "Other"
     8 "Don't know"
;
label define XF118E
     1 "2-3 months"
     2 "Other"
     8 "Don't know"
;
label define XF118F
     1 "One month"
     2 "Other"
     8 "Don't know"
;
label define XF118G
     1 "3-7 years"
     2 "Other"
     8 "Don't know"
;
label define XF118I
     1 "Yes, definitely can make woman pregnant"
     2 "Yes, only slight risk of making a woman pregnant"
     3 "No, cannot make a woman pregnant"
     8 "Don't know"
;
label define XF118J
     1 "Yes, definitely can get pregnant"
     2 "Yes, only slight risk of getting pregnant"
     3 "No, cannot become pregnant"
     8 "Don't know"
;
label define XF118K
     1 "Yes"
     2 "No"
     8 "Don't know"
;
label define XF119
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
label values of104    OF104
label values of131a   OF131A
label values of131b   OF131B
label values of202    OF202
label values of203    OF203
label values of204    OF204
label values of205    OF205
label values of206    OF206
label values of215    OF215
label values of219    OF219
label values of225    OF225
label values of301    OF301
label values of302    OF302
label values of303    OF303
label values of304    OF304
label values of305    OF305
label values of306    OF306
label values of307h   OF307H
label values of307m   OF307M
label values xfpmonth XFPMONTH
label values xfpintv  XFPINTV
label values xf100    XF100
label values xf102    XF102
label values xf103    XF103
label values xf104    XF104
label values xf106    XF106
label values xf107    XF107
label values xf108    XF108
label values xf109    XF109
label values xf110    XF110
label values xf111    XF111
label values xf112    XF112
label values xf114    XF114
label values xf117_1  XF117
label values xf117_2  XF117
label values xf117_3  XF117
label values xf117_4  XF117
label values xf118a   XF118A
label values xf118b   XF118B
label values xf118c   XF118C
label values xf118d   XF118D
label values xf118e   XF118E
label values xf118f   XF118F
label values xf118g   XF118G
label values xf118i   XF118I
label values xf118j   XF118J
label values xf118k   XF118K
label values xf119    XF119
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
