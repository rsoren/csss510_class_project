infix using "C:\HAITISPA_FIN\ARCHIVE\HTAN6BFLSP.dct"

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
label variable oa104    "Client history"
label variable oa105    "Aspects of prior pregnancies"
label variable oa106a   "Danger signs of current pregnancy asked by provider or client mentioned"
label variable oa106b   "Danger signs of current pregnancy counselled by provider"
label variable oa107    "Physical examination"
label variable oa108_1  "Routine tests"
label variable oa108_2  "Routine tests"
label variable oa108_3  "Routine tests"
label variable oa108_4  "Routine tests"
label variable oa109    "HIV counseling and testing"
label variable oa110    "Healthy pregnancy"
label variable oa111    "Iron prophylaxis"
label variable oa112    "Tetanus Toxoid injection"
label variable oa113    "Deworming"
label variable oa114    "Malaria"
label variable oa115    "Preparation for delivery"
label variable oa116    "Newborn/ postpartum recommendations"
label variable oa117    "Encouraged questions"
label variable oa118    "Used visual aids"
label variable oa119    "Provider looked at client card"
label variable oa120    "Provider wrote on card"
label variable oa121    "Outcome of consult"
label variable oa122    "Weeks pregnant"
label variable oa123    "Visit number"
label variable oa124    "Number of pregnancies for client"
label variable oa125h   "Hour ended observation"
label variable oa125m   "Minute ended observation"
label variable xanday   "Client Day"
label variable xanmonth "Client Month"
label variable xanyear  "Client Year"
label variable xanintv  "Interviewer code"
label variable xanccode "Client code"
label variable xa100    "Permission to continue"
label variable xa101h   "Time begins hour"
label variable xa101m   "Time begins minutes"
label variable xa102    "Antenatal care card"
label variable xa103    "TT recorded on card"
label variable xa104    "Weeks pregnant recorded on ANC card"
label variable xa105    "Card indicated client received IPT - NOT USED"
label variable xa106    "First pregnancy"
label variable xa107    "Visit number to this facility for this pregnancy"
label variable xa108    "Ever given iron tablets"
label variable xa109    "Provider ever explained how to take iron tablets"
label variable xa110    "Provider ever discussed side effects of iron tablets"
label variable xa111    "Client's knowledge of side effects of iron tablets"
label variable xa112    "Given or prescribed anti-malarial"
label variable xa113    "Swallowed pills or have them"
label variable xa114    "Provider ever advise use of mosquito net treated with insecticide"
label variable xa115    "Provider ever offered ITN free of charge"
label variable xa116    "Provider offered to sell ITN or recommended place to purchase"
label variable xa117    "Provider ever discussed diet/nutrition during pregnancy"
label variable xa118    "Provider talked about danger signs in pregnancy"
label variable xa119    "Client knowledge of danger signs"
label variable xa120    "Provider advice about danger signs"
label variable xa121    "Provider ever discussed preparatin for delivery"
label variable xa122    "Client knowledge of preparatory needs"
label variable xa123    "Money available for delivery"
label variable xa124    "Provider ever discussed delivery plans"
label variable xa125    "Where client plans to deliver"
label variable xa126    "Client knowledge of danger signs during/after delivery"
label variable xa127    "Provider ever advised about importance of exclusive breastfeeding"
label variable xa128    "Number of months provider advised for exclusive breastfeeding"
label variable xa129    "Provider talked about FP"
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
label define OA117
     1 "Yes, asked questions"
     2 "No, did not ask questions"
;
label define OA118
     1 "Yes, used visual aids"
     2 "No, did not use visual aids"
;
label define OA119
     1 "Yes, looked at card"
     2 "No, did not look at card"
     3 "Card not used"
;
label define OA120
     1 "Yes"
     2 "No"
     3 "Card not used"
     8 "Don't know"
;
label define OA121
     1 "Client goes home"
     2 "Client referred (to lab or other provider) at same facility"
     3 "Client admitted to same facility"
     4 "Client referred to other faciltiy"
;
label define OA123
     1 "First visit"
     2 "Second visit"
     3 "Third visit"
     4 "Fourth visit"
     5 "Fifth or more visit"
     8 "Don't know"
;
label define OA124
     1 "First pregnancy"
     2 "Not first pregnancy"
     8 "Don't know"
;
label define XANINTV
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
label define XA100
     1 "Yes, agrees"
     2 "No, client refuses"
;
label define XA102
     1 "Yes"
     2 "No, card kept at facility"
     3 "No card/book used"
;
label define XA103
     1 "Yes 1 time"
     2 "Yes 2 times"
     3 "Yes, 3+ times"
     4 "No record"
;
label define XA104
    95 "Not available"
;
label define XA105
     1 "Yes, 1 dose"
     2 "Yes, 2 doses"
     3 "Yes, 3 doses"
     4 "Yes, 4 doses"
     5 "No"
;
label define XA106
     1 "First pregnancy"
     2 "Not first pregancy"
;
label define XA107
     1 "First visit"
     2 "Second visit"
     3 "Third visit"
     4 "Fourth visit"
     5 "More than 4 visits"
;
label define XA108
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA109
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA110
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA112
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA113
     1 "Yes"
     2 "No"
;
label define XA114
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA115
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA116
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA117
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA118
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA121
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Ne sait pas"
;
label define XA123
     1 "Yes, enough"
     2 "Yes, but not enough"
     3 "No"
;
label define XA124
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA125
     1 "At this facility"
     2 "At other health facility"
     3 "At own home"
     4 "At TBAs home"
     6 "Other location (specify)"
     8 "No/ Don't know/ undecided"
;
label define XA127
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
     8 "Don't know"
;
label define XA128
     1 "Between 4-6 months"
     2 "6  months"
     6 "Other"
     8 "Don't know"
;
label define XA129
     1 "Yes, this visit only"
     2 "Yes, this and previous visit"
     3 "Yes, previous visit only"
     4 "No"
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
label values oa117    OA117
label values oa118    OA118
label values oa119    OA119
label values oa120    OA120
label values oa121    OA121
label values oa123    OA123
label values oa124    OA124
label values xanintv  XANINTV
label values xa100    XA100
label values xa102    XA102
label values xa103    XA103
label values xa104    XA104
label values xa105    XA105
label values xa106    XA106
label values xa107    XA107
label values xa108    XA108
label values xa109    XA109
label values xa110    XA110
label values xa112    XA112
label values xa113    XA113
label values xa114    XA114
label values xa115    XA115
label values xa116    XA116
label values xa117    XA117
label values xa118    XA118
label values xa121    XA121
label values xa123    XA123
label values xa124    XA124
label values xa125    XA125
label values xa127    XA127
label values xa128    XA128
label values xa129    XA129
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
