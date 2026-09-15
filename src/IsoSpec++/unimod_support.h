// Generated from data/unimod.csv and its Unimod OBO source by scripts/build_unimod_table.py.
// Index = original Unimod ID; true = supported by the embedded IsoSpec table.
#pragma once

#include <cstddef>
#include <cstdint>

namespace IsoSpec {

inline constexpr bool unimod_supported[] = {
    false,  // UNIMOD:0: ontology root, not a modification
    true,   // UNIMOD:1
    true,   // UNIMOD:2
    true,   // UNIMOD:3
    true,   // UNIMOD:4
    true,   // UNIMOD:5
    true,   // UNIMOD:6
    true,   // UNIMOD:7
    true,   // UNIMOD:8
    false,  // UNIMOD:9: isotope-labeled token 2H(8); pinned-isotope conversion not supported
    true,   // UNIMOD:10
    true,   // UNIMOD:11
    false,  // UNIMOD:12: isotope-labeled token 2H(8); pinned-isotope conversion not supported
    true,   // UNIMOD:13
    false,  // UNIMOD:14: not present in source Unimod snapshot
    false,  // UNIMOD:15: not present in source Unimod snapshot
    false,  // UNIMOD:16: not present in source Unimod snapshot
    true,   // UNIMOD:17
    false,  // UNIMOD:18: not present in source Unimod snapshot
    false,  // UNIMOD:19: not present in source Unimod snapshot
    true,   // UNIMOD:20
    true,   // UNIMOD:21
    false,  // UNIMOD:22: not present in source Unimod snapshot
    true,   // UNIMOD:23
    true,   // UNIMOD:24
    true,   // UNIMOD:25
    true,   // UNIMOD:26
    true,   // UNIMOD:27
    true,   // UNIMOD:28
    true,   // UNIMOD:29
    true,   // UNIMOD:30
    true,   // UNIMOD:31
    false,  // UNIMOD:32: not present in source Unimod snapshot
    false,  // UNIMOD:33: not present in source Unimod snapshot
    true,   // UNIMOD:34
    true,   // UNIMOD:35
    true,   // UNIMOD:36
    true,   // UNIMOD:37
    false,  // UNIMOD:38: not present in source Unimod snapshot
    true,   // UNIMOD:39
    true,   // UNIMOD:40
    false,  // UNIMOD:41: non-element token Hex; group-to-formula expansion not supported
    true,   // UNIMOD:42
    false,  // UNIMOD:43: non-element token HexNAc; group-to-formula expansion not supported
    true,   // UNIMOD:44
    true,   // UNIMOD:45
    true,   // UNIMOD:46
    true,   // UNIMOD:47
    true,   // UNIMOD:48
    true,   // UNIMOD:49
    true,   // UNIMOD:50
    true,   // UNIMOD:51
    true,   // UNIMOD:52
    true,   // UNIMOD:53
    false,  // UNIMOD:54: non-element token HexA; group-to-formula expansion not supported
    true,   // UNIMOD:55
    false,  // UNIMOD:56: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    false,  // UNIMOD:57: not present in source Unimod snapshot
    true,   // UNIMOD:58
    false,  // UNIMOD:59: isotope-labeled token 13C(3); pinned-isotope conversion not supported
    true,   // UNIMOD:60
    false,  // UNIMOD:61: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    false,  // UNIMOD:62: isotope-labeled token 2H(6); pinned-isotope conversion not supported
    false,  // UNIMOD:63: isotope-labeled token 2H(9); pinned-isotope conversion not supported
    true,   // UNIMOD:64
    false,  // UNIMOD:65: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    false,  // UNIMOD:66: isotope-labeled token 13C(4); pinned-isotope conversion not supported
    false,  // UNIMOD:67: not present in source Unimod snapshot
    false,  // UNIMOD:68: not present in source Unimod snapshot
    false,  // UNIMOD:69: not present in source Unimod snapshot
    false,  // UNIMOD:70: not present in source Unimod snapshot
    false,  // UNIMOD:71: not present in source Unimod snapshot
    false,  // UNIMOD:72: not present in source Unimod snapshot
    false,  // UNIMOD:73: not present in source Unimod snapshot
    false,  // UNIMOD:74: not present in source Unimod snapshot
    false,  // UNIMOD:75: not present in source Unimod snapshot
    false,  // UNIMOD:76: not present in source Unimod snapshot
    false,  // UNIMOD:77: not present in source Unimod snapshot
    false,  // UNIMOD:78: not present in source Unimod snapshot
    false,  // UNIMOD:79: not present in source Unimod snapshot
    false,  // UNIMOD:80: not present in source Unimod snapshot
    false,  // UNIMOD:81: not present in source Unimod snapshot
    false,  // UNIMOD:82: not present in source Unimod snapshot
    false,  // UNIMOD:83: not present in source Unimod snapshot
    false,  // UNIMOD:84: not present in source Unimod snapshot
    false,  // UNIMOD:85: not present in source Unimod snapshot
    false,  // UNIMOD:86: not present in source Unimod snapshot
    false,  // UNIMOD:87: not present in source Unimod snapshot
    false,  // UNIMOD:88: not present in source Unimod snapshot
    true,   // UNIMOD:89
    true,   // UNIMOD:90
    false,  // UNIMOD:91: isotope-labeled token 2H(10); pinned-isotope conversion not supported
    true,   // UNIMOD:92
    true,   // UNIMOD:93
    true,   // UNIMOD:94
    false,  // UNIMOD:95: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    false,  // UNIMOD:96: not present in source Unimod snapshot
    false,  // UNIMOD:97: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    false,  // UNIMOD:98: not present in source Unimod snapshot
    false,  // UNIMOD:99: not present in source Unimod snapshot
    false,  // UNIMOD:100: not present in source Unimod snapshot
    false,  // UNIMOD:101: not present in source Unimod snapshot
    false,  // UNIMOD:102: not present in source Unimod snapshot
    false,  // UNIMOD:103: not present in source Unimod snapshot
    false,  // UNIMOD:104: not present in source Unimod snapshot
    true,   // UNIMOD:105
    false,  // UNIMOD:106: isotope-labeled token 13C(9); pinned-isotope conversion not supported
    true,   // UNIMOD:107
    true,   // UNIMOD:108
    false,  // UNIMOD:109: not present in source Unimod snapshot
    false,  // UNIMOD:110: not present in source Unimod snapshot
    false,  // UNIMOD:111: not present in source Unimod snapshot
    true,   // UNIMOD:112
    true,   // UNIMOD:113
    true,   // UNIMOD:114
    true,   // UNIMOD:115
    true,   // UNIMOD:116
    true,   // UNIMOD:117
    true,   // UNIMOD:118
    true,   // UNIMOD:119
    false,  // UNIMOD:120: not present in source Unimod snapshot
    true,   // UNIMOD:121
    true,   // UNIMOD:122
    true,   // UNIMOD:123
    false,  // UNIMOD:124: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:125: not present in source Unimod snapshot
    true,   // UNIMOD:126
    true,   // UNIMOD:127
    true,   // UNIMOD:128
    true,   // UNIMOD:129
    true,   // UNIMOD:130
    true,   // UNIMOD:131
    false,  // UNIMOD:132: not present in source Unimod snapshot
    false,  // UNIMOD:133: not present in source Unimod snapshot
    true,   // UNIMOD:134
    true,   // UNIMOD:135
    true,   // UNIMOD:136
    false,  // UNIMOD:137: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:138: not present in source Unimod snapshot
    true,   // UNIMOD:139
    true,   // UNIMOD:140
    true,   // UNIMOD:141
    false,  // UNIMOD:142: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:143: non-element token HexNAc(2); group-to-formula expansion not supported
    false,  // UNIMOD:144: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:145: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:146: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:147: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:148: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:149: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:150: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:151: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:152: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:153: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:154: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:155: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:156: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:157: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:158: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:159: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:160: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:161: non-element token Hex(3); group-to-formula expansion not supported
    true,   // UNIMOD:162
    false,  // UNIMOD:163: not present in source Unimod snapshot
    false,  // UNIMOD:164: not present in source Unimod snapshot
    false,  // UNIMOD:165: not present in source Unimod snapshot
    false,  // UNIMOD:166: not present in source Unimod snapshot
    false,  // UNIMOD:167: not present in source Unimod snapshot
    false,  // UNIMOD:168: not present in source Unimod snapshot
    false,  // UNIMOD:169: not present in source Unimod snapshot
    false,  // UNIMOD:170: isotope-labeled token 18O; pinned-isotope conversion not supported
    false,  // UNIMOD:171: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    true,   // UNIMOD:172
    false,  // UNIMOD:173: not present in source Unimod snapshot
    false,  // UNIMOD:174: not present in source Unimod snapshot
    false,  // UNIMOD:175: not present in source Unimod snapshot
    true,   // UNIMOD:176
    false,  // UNIMOD:177: not present in source Unimod snapshot
    true,   // UNIMOD:178
    false,  // UNIMOD:179: not present in source Unimod snapshot
    false,  // UNIMOD:180: not present in source Unimod snapshot
    false,  // UNIMOD:181: not present in source Unimod snapshot
    false,  // UNIMOD:182: not present in source Unimod snapshot
    false,  // UNIMOD:183: not present in source Unimod snapshot
    false,  // UNIMOD:184: isotope-labeled token 13C(9); pinned-isotope conversion not supported
    false,  // UNIMOD:185: isotope-labeled token 13C(9); pinned-isotope conversion not supported
    true,   // UNIMOD:186
    true,   // UNIMOD:187
    false,  // UNIMOD:188: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:189: not present in source Unimod snapshot
    false,  // UNIMOD:190: not present in source Unimod snapshot
    false,  // UNIMOD:191: not present in source Unimod snapshot
    false,  // UNIMOD:192: not present in source Unimod snapshot
    false,  // UNIMOD:193: isotope-labeled token 18O(2); pinned-isotope conversion not supported
    true,   // UNIMOD:194
    true,   // UNIMOD:195
    false,  // UNIMOD:196: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    true,   // UNIMOD:197
    false,  // UNIMOD:198: isotope-labeled token 2H(5); pinned-isotope conversion not supported
    false,  // UNIMOD:199: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    true,   // UNIMOD:200
    false,  // UNIMOD:201: not present in source Unimod snapshot
    false,  // UNIMOD:202: not present in source Unimod snapshot
    false,  // UNIMOD:203: not present in source Unimod snapshot
    false,  // UNIMOD:204: not present in source Unimod snapshot
    true,   // UNIMOD:205
    true,   // UNIMOD:206
    true,   // UNIMOD:207
    true,   // UNIMOD:208
    true,   // UNIMOD:209
    false,  // UNIMOD:210: not present in source Unimod snapshot
    true,   // UNIMOD:211
    false,  // UNIMOD:212: isotope-labeled token 2H(5); pinned-isotope conversion not supported
    false,  // UNIMOD:213: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:214: isotope-labeled token 13C(3); pinned-isotope conversion not supported
    false,  // UNIMOD:215: not present in source Unimod snapshot
    false,  // UNIMOD:216: not present in source Unimod snapshot
    false,  // UNIMOD:217: not present in source Unimod snapshot
    false,  // UNIMOD:218: not present in source Unimod snapshot
    false,  // UNIMOD:219: not present in source Unimod snapshot
    false,  // UNIMOD:220: not present in source Unimod snapshot
    false,  // UNIMOD:221: not present in source Unimod snapshot
    false,  // UNIMOD:222: not present in source Unimod snapshot
    false,  // UNIMOD:223: not present in source Unimod snapshot
    false,  // UNIMOD:224: not present in source Unimod snapshot
    false,  // UNIMOD:225: not present in source Unimod snapshot
    false,  // UNIMOD:226: not present in source Unimod snapshot
    false,  // UNIMOD:227: not present in source Unimod snapshot
    false,  // UNIMOD:228: not present in source Unimod snapshot
    false,  // UNIMOD:229: not present in source Unimod snapshot
    false,  // UNIMOD:230: not present in source Unimod snapshot
    false,  // UNIMOD:231: not present in source Unimod snapshot
    false,  // UNIMOD:232: not present in source Unimod snapshot
    false,  // UNIMOD:233: not present in source Unimod snapshot
    false,  // UNIMOD:234: not present in source Unimod snapshot
    false,  // UNIMOD:235: not present in source Unimod snapshot
    false,  // UNIMOD:236: not present in source Unimod snapshot
    false,  // UNIMOD:237: not present in source Unimod snapshot
    false,  // UNIMOD:238: not present in source Unimod snapshot
    false,  // UNIMOD:239: not present in source Unimod snapshot
    false,  // UNIMOD:240: not present in source Unimod snapshot
    false,  // UNIMOD:241: not present in source Unimod snapshot
    false,  // UNIMOD:242: not present in source Unimod snapshot
    true,   // UNIMOD:243
    false,  // UNIMOD:244: not present in source Unimod snapshot
    false,  // UNIMOD:245: not present in source Unimod snapshot
    false,  // UNIMOD:246: not present in source Unimod snapshot
    false,  // UNIMOD:247: not present in source Unimod snapshot
    false,  // UNIMOD:248: not present in source Unimod snapshot
    false,  // UNIMOD:249: not present in source Unimod snapshot
    false,  // UNIMOD:250: not present in source Unimod snapshot
    false,  // UNIMOD:251: not present in source Unimod snapshot
    false,  // UNIMOD:252: not present in source Unimod snapshot
    true,   // UNIMOD:253
    true,   // UNIMOD:254
    true,   // UNIMOD:255
    true,   // UNIMOD:256
    false,  // UNIMOD:257: not present in source Unimod snapshot
    false,  // UNIMOD:258: isotope-labeled token 18O; pinned-isotope conversion not supported
    false,  // UNIMOD:259: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    true,   // UNIMOD:260
    true,   // UNIMOD:261
    false,  // UNIMOD:262: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    false,  // UNIMOD:263: not present in source Unimod snapshot
    true,   // UNIMOD:264
    false,  // UNIMOD:265: not present in source Unimod snapshot
    false,  // UNIMOD:266: not present in source Unimod snapshot
    false,  // UNIMOD:267: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:268: isotope-labeled token 13C(5); pinned-isotope conversion not supported
    false,  // UNIMOD:269: isotope-labeled token 13C(9); pinned-isotope conversion not supported
    true,   // UNIMOD:270
    true,   // UNIMOD:271
    true,   // UNIMOD:272
    false,  // UNIMOD:273: not present in source Unimod snapshot
    false,  // UNIMOD:274: not present in source Unimod snapshot
    true,   // UNIMOD:275
    true,   // UNIMOD:276
    false,  // UNIMOD:277: not present in source Unimod snapshot
    true,   // UNIMOD:278
    false,  // UNIMOD:279: not present in source Unimod snapshot
    true,   // UNIMOD:280
    true,   // UNIMOD:281
    false,  // UNIMOD:282: not present in source Unimod snapshot
    false,  // UNIMOD:283: not present in source Unimod snapshot
    false,  // UNIMOD:284: isotope-labeled token 2H(2); pinned-isotope conversion not supported
    true,   // UNIMOD:285
    false,  // UNIMOD:286: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:287: not present in source Unimod snapshot
    true,   // UNIMOD:288
    true,   // UNIMOD:289
    true,   // UNIMOD:290
    true,   // UNIMOD:291
    true,   // UNIMOD:292
    true,   // UNIMOD:293
    true,   // UNIMOD:294
    false,  // UNIMOD:295: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:296: not present in source Unimod snapshot
    false,  // UNIMOD:297: not present in source Unimod snapshot
    false,  // UNIMOD:298: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    true,   // UNIMOD:299
    false,  // UNIMOD:300: not present in source Unimod snapshot
    true,   // UNIMOD:301
    true,   // UNIMOD:302
    true,   // UNIMOD:303
    false,  // UNIMOD:304: not present in source Unimod snapshot
    false,  // UNIMOD:305: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:306: not present in source Unimod snapshot
    false,  // UNIMOD:307: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:308: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:309: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:310: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:311: non-element token Hex(5); group-to-formula expansion not supported
    true,   // UNIMOD:312
    true,   // UNIMOD:313
    true,   // UNIMOD:314
    false,  // UNIMOD:315: not present in source Unimod snapshot
    true,   // UNIMOD:316
    false,  // UNIMOD:317: not present in source Unimod snapshot
    true,   // UNIMOD:318
    true,   // UNIMOD:319
    true,   // UNIMOD:320
    false,  // UNIMOD:321: not present in source Unimod snapshot
    false,  // UNIMOD:322: not present in source Unimod snapshot
    true,   // UNIMOD:323
    true,   // UNIMOD:324
    true,   // UNIMOD:325
    false,  // UNIMOD:326: not present in source Unimod snapshot
    true,   // UNIMOD:327
    false,  // UNIMOD:328: not present in source Unimod snapshot
    false,  // UNIMOD:329: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    false,  // UNIMOD:330: isotope-labeled token 2H(6); pinned-isotope conversion not supported
    false,  // UNIMOD:331: not present in source Unimod snapshot
    true,   // UNIMOD:332
    true,   // UNIMOD:333
    false,  // UNIMOD:334: not present in source Unimod snapshot
    true,   // UNIMOD:335
    false,  // UNIMOD:336: not present in source Unimod snapshot
    true,   // UNIMOD:337
    false,  // UNIMOD:338: not present in source Unimod snapshot
    false,  // UNIMOD:339: not present in source Unimod snapshot
    true,   // UNIMOD:340
    false,  // UNIMOD:341: not present in source Unimod snapshot
    true,   // UNIMOD:342
    true,   // UNIMOD:343
    true,   // UNIMOD:344
    true,   // UNIMOD:345
    false,  // UNIMOD:346: not present in source Unimod snapshot
    false,  // UNIMOD:347: not present in source Unimod snapshot
    true,   // UNIMOD:348
    true,   // UNIMOD:349
    true,   // UNIMOD:350
    true,   // UNIMOD:351
    true,   // UNIMOD:352
    true,   // UNIMOD:353
    true,   // UNIMOD:354
    false,  // UNIMOD:355: not present in source Unimod snapshot
    false,  // UNIMOD:356: not present in source Unimod snapshot
    true,   // UNIMOD:357
    false,  // UNIMOD:358: not present in source Unimod snapshot
    true,   // UNIMOD:359
    true,   // UNIMOD:360
    true,   // UNIMOD:361
    true,   // UNIMOD:362
    true,   // UNIMOD:363
    false,  // UNIMOD:364: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    true,   // UNIMOD:365
    false,  // UNIMOD:366: isotope-labeled token 18O; pinned-isotope conversion not supported
    false,  // UNIMOD:367: not present in source Unimod snapshot
    true,   // UNIMOD:368
    true,   // UNIMOD:369
    false,  // UNIMOD:370: not present in source Unimod snapshot
    true,   // UNIMOD:371
    true,   // UNIMOD:372
    false,  // UNIMOD:373: not present in source Unimod snapshot
    true,   // UNIMOD:374
    true,   // UNIMOD:375
    true,   // UNIMOD:376
    true,   // UNIMOD:377
    true,   // UNIMOD:378
    true,   // UNIMOD:379
    true,   // UNIMOD:380
    true,   // UNIMOD:381
    true,   // UNIMOD:382
    false,  // UNIMOD:383: not present in source Unimod snapshot
    false,  // UNIMOD:384: not present in source Unimod snapshot
    true,   // UNIMOD:385
    false,  // UNIMOD:386: not present in source Unimod snapshot
    true,   // UNIMOD:387
    true,   // UNIMOD:388
    true,   // UNIMOD:389
    true,   // UNIMOD:390
    true,   // UNIMOD:391
    true,   // UNIMOD:392
    false,  // UNIMOD:393: non-element token Hex(2); group-to-formula expansion not supported
    true,   // UNIMOD:394
    true,   // UNIMOD:395
    true,   // UNIMOD:396
    true,   // UNIMOD:397
    true,   // UNIMOD:398
    false,  // UNIMOD:399: not present in source Unimod snapshot
    true,   // UNIMOD:400
    true,   // UNIMOD:401
    true,   // UNIMOD:402
    true,   // UNIMOD:403
    false,  // UNIMOD:404: not present in source Unimod snapshot
    true,   // UNIMOD:405
    false,  // UNIMOD:406: not present in source Unimod snapshot
    true,   // UNIMOD:407
    false,  // UNIMOD:408: non-element token Hex; group-to-formula expansion not supported
    true,   // UNIMOD:409
    true,   // UNIMOD:410
    true,   // UNIMOD:411
    false,  // UNIMOD:412: isotope-labeled token 2H(5); pinned-isotope conversion not supported
    true,   // UNIMOD:413
    true,   // UNIMOD:414
    true,   // UNIMOD:415
    true,   // UNIMOD:416
    true,   // UNIMOD:417
    false,  // UNIMOD:418: not present in source Unimod snapshot
    true,   // UNIMOD:419
    true,   // UNIMOD:420
    true,   // UNIMOD:421
    true,   // UNIMOD:422
    true,   // UNIMOD:423
    true,   // UNIMOD:424
    true,   // UNIMOD:425
    true,   // UNIMOD:426
    false,  // UNIMOD:427: not present in source Unimod snapshot
    false,  // UNIMOD:428: non-element token HexNAc; group-to-formula expansion not supported
    false,  // UNIMOD:429: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:430: not present in source Unimod snapshot
    true,   // UNIMOD:431
    true,   // UNIMOD:432
    true,   // UNIMOD:433
    true,   // UNIMOD:434
    true,   // UNIMOD:435
    true,   // UNIMOD:436
    true,   // UNIMOD:437
    true,   // UNIMOD:438
    true,   // UNIMOD:439
    true,   // UNIMOD:440
    false,  // UNIMOD:441: not present in source Unimod snapshot
    true,   // UNIMOD:442
    true,   // UNIMOD:443
    true,   // UNIMOD:444
    true,   // UNIMOD:445
    false,  // UNIMOD:446: not present in source Unimod snapshot
    true,   // UNIMOD:447
    true,   // UNIMOD:448
    true,   // UNIMOD:449
    true,   // UNIMOD:450
    true,   // UNIMOD:451
    true,   // UNIMOD:452
    true,   // UNIMOD:453
    false,  // UNIMOD:454: non-element token HexN; group-to-formula expansion not supported
    true,   // UNIMOD:455
    false,  // UNIMOD:456: not present in source Unimod snapshot
    true,   // UNIMOD:457
    false,  // UNIMOD:458: not present in source Unimod snapshot
    false,  // UNIMOD:459: not present in source Unimod snapshot
    false,  // UNIMOD:460: not present in source Unimod snapshot
    false,  // UNIMOD:461: not present in source Unimod snapshot
    false,  // UNIMOD:462: not present in source Unimod snapshot
    false,  // UNIMOD:463: not present in source Unimod snapshot
    false,  // UNIMOD:464: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:465: not present in source Unimod snapshot
    false,  // UNIMOD:466: not present in source Unimod snapshot
    false,  // UNIMOD:467: not present in source Unimod snapshot
    false,  // UNIMOD:468: not present in source Unimod snapshot
    false,  // UNIMOD:469: not present in source Unimod snapshot
    false,  // UNIMOD:470: not present in source Unimod snapshot
    false,  // UNIMOD:471: not present in source Unimod snapshot
    true,   // UNIMOD:472
    false,  // UNIMOD:473: not present in source Unimod snapshot
    false,  // UNIMOD:474: not present in source Unimod snapshot
    false,  // UNIMOD:475: not present in source Unimod snapshot
    true,   // UNIMOD:476
    false,  // UNIMOD:477: isotope-labeled token 2H(9); pinned-isotope conversion not supported
    true,   // UNIMOD:478
    false,  // UNIMOD:479: not present in source Unimod snapshot
    false,  // UNIMOD:480: not present in source Unimod snapshot
    false,  // UNIMOD:481: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    false,  // UNIMOD:482: not present in source Unimod snapshot
    false,  // UNIMOD:483: not present in source Unimod snapshot
    false,  // UNIMOD:484: not present in source Unimod snapshot
    false,  // UNIMOD:485: not present in source Unimod snapshot
    false,  // UNIMOD:486: not present in source Unimod snapshot
    false,  // UNIMOD:487: not present in source Unimod snapshot
    true,   // UNIMOD:488
    false,  // UNIMOD:489: not present in source Unimod snapshot
    false,  // UNIMOD:490: non-element token Hep; group-to-formula expansion not supported
    false,  // UNIMOD:491: not present in source Unimod snapshot
    false,  // UNIMOD:492: not present in source Unimod snapshot
    true,   // UNIMOD:493
    true,   // UNIMOD:494
    true,   // UNIMOD:495
    false,  // UNIMOD:496: not present in source Unimod snapshot
    false,  // UNIMOD:497: not present in source Unimod snapshot
    true,   // UNIMOD:498
    false,  // UNIMOD:499: isotope-labeled token 13C(2); pinned-isotope conversion not supported
    true,   // UNIMOD:500
    true,   // UNIMOD:501
    false,  // UNIMOD:502: not present in source Unimod snapshot
    true,   // UNIMOD:503
    true,   // UNIMOD:504
    true,   // UNIMOD:505
    true,   // UNIMOD:506
    false,  // UNIMOD:507: not present in source Unimod snapshot
    false,  // UNIMOD:508: not present in source Unimod snapshot
    false,  // UNIMOD:509: not present in source Unimod snapshot
    false,  // UNIMOD:510: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    false,  // UNIMOD:511: not present in source Unimod snapshot
    false,  // UNIMOD:512: non-element token Hex(2); group-to-formula expansion not supported
    true,   // UNIMOD:513
    true,   // UNIMOD:514
    true,   // UNIMOD:515
    false,  // UNIMOD:516: not present in source Unimod snapshot
    false,  // UNIMOD:517: not present in source Unimod snapshot
    true,   // UNIMOD:518
    true,   // UNIMOD:519
    true,   // UNIMOD:520
    false,  // UNIMOD:521: not present in source Unimod snapshot
    true,   // UNIMOD:522
    true,   // UNIMOD:523
    false,  // UNIMOD:524: not present in source Unimod snapshot
    false,  // UNIMOD:525: isotope-labeled token 13C; pinned-isotope conversion not supported
    true,   // UNIMOD:526
    false,  // UNIMOD:527: not present in source Unimod snapshot
    true,   // UNIMOD:528
    true,   // UNIMOD:529
    true,   // UNIMOD:530
    true,   // UNIMOD:531
    false,  // UNIMOD:532: isotope-labeled token 13C(2); pinned-isotope conversion not supported
    false,  // UNIMOD:533: isotope-labeled token 13C; pinned-isotope conversion not supported
    true,   // UNIMOD:534
    true,   // UNIMOD:535
    false,  // UNIMOD:536: isotope-labeled token 13C; pinned-isotope conversion not supported
    false,  // UNIMOD:537: isotope-labeled token 13C; pinned-isotope conversion not supported
    true,   // UNIMOD:538
    true,   // UNIMOD:539
    true,   // UNIMOD:540
    true,   // UNIMOD:541
    true,   // UNIMOD:542
    true,   // UNIMOD:543
    true,   // UNIMOD:544
    true,   // UNIMOD:545
    true,   // UNIMOD:546
    true,   // UNIMOD:547
    true,   // UNIMOD:548
    true,   // UNIMOD:549
    true,   // UNIMOD:550
    true,   // UNIMOD:551
    true,   // UNIMOD:552
    true,   // UNIMOD:553
    true,   // UNIMOD:554
    true,   // UNIMOD:555
    true,   // UNIMOD:556
    true,   // UNIMOD:557
    true,   // UNIMOD:558
    true,   // UNIMOD:559
    true,   // UNIMOD:560
    true,   // UNIMOD:561
    true,   // UNIMOD:562
    true,   // UNIMOD:563
    true,   // UNIMOD:564
    true,   // UNIMOD:565
    true,   // UNIMOD:566
    true,   // UNIMOD:567
    true,   // UNIMOD:568
    true,   // UNIMOD:569
    true,   // UNIMOD:570
    true,   // UNIMOD:571
    true,   // UNIMOD:572
    true,   // UNIMOD:573
    true,   // UNIMOD:574
    true,   // UNIMOD:575
    true,   // UNIMOD:576
    true,   // UNIMOD:577
    true,   // UNIMOD:578
    false,  // UNIMOD:579: not present in source Unimod snapshot
    true,   // UNIMOD:580
    true,   // UNIMOD:581
    true,   // UNIMOD:582
    false,  // UNIMOD:583: not present in source Unimod snapshot
    true,   // UNIMOD:584
    true,   // UNIMOD:585
    false,  // UNIMOD:586: not present in source Unimod snapshot
    false,  // UNIMOD:587: not present in source Unimod snapshot
    true,   // UNIMOD:588
    true,   // UNIMOD:589
    true,   // UNIMOD:590
    false,  // UNIMOD:591: not present in source Unimod snapshot
    false,  // UNIMOD:592: not present in source Unimod snapshot
    false,  // UNIMOD:593: not present in source Unimod snapshot
    true,   // UNIMOD:594
    true,   // UNIMOD:595
    true,   // UNIMOD:596
    true,   // UNIMOD:597
    true,   // UNIMOD:598
    true,   // UNIMOD:599
    true,   // UNIMOD:600
    true,   // UNIMOD:601
    true,   // UNIMOD:602
    true,   // UNIMOD:603
    true,   // UNIMOD:604
    true,   // UNIMOD:605
    true,   // UNIMOD:606
    true,   // UNIMOD:607
    true,   // UNIMOD:608
    true,   // UNIMOD:609
    true,   // UNIMOD:610
    true,   // UNIMOD:611
    false,  // UNIMOD:612: not present in source Unimod snapshot
    true,   // UNIMOD:613
    true,   // UNIMOD:614
    true,   // UNIMOD:615
    true,   // UNIMOD:616
    true,   // UNIMOD:617
    true,   // UNIMOD:618
    true,   // UNIMOD:619
    true,   // UNIMOD:620
    true,   // UNIMOD:621
    true,   // UNIMOD:622
    true,   // UNIMOD:623
    true,   // UNIMOD:624
    true,   // UNIMOD:625
    true,   // UNIMOD:626
    true,   // UNIMOD:627
    true,   // UNIMOD:628
    true,   // UNIMOD:629
    true,   // UNIMOD:630
    true,   // UNIMOD:631
    true,   // UNIMOD:632
    true,   // UNIMOD:633
    true,   // UNIMOD:634
    true,   // UNIMOD:635
    true,   // UNIMOD:636
    true,   // UNIMOD:637
    true,   // UNIMOD:638
    true,   // UNIMOD:639
    true,   // UNIMOD:640
    true,   // UNIMOD:641
    true,   // UNIMOD:642
    true,   // UNIMOD:643
    true,   // UNIMOD:644
    true,   // UNIMOD:645
    true,   // UNIMOD:646
    true,   // UNIMOD:647
    true,   // UNIMOD:648
    true,   // UNIMOD:649
    true,   // UNIMOD:650
    true,   // UNIMOD:651
    true,   // UNIMOD:652
    true,   // UNIMOD:653
    true,   // UNIMOD:654
    true,   // UNIMOD:655
    true,   // UNIMOD:656
    true,   // UNIMOD:657
    true,   // UNIMOD:658
    true,   // UNIMOD:659
    true,   // UNIMOD:660
    true,   // UNIMOD:661
    true,   // UNIMOD:662
    true,   // UNIMOD:663
    true,   // UNIMOD:664
    true,   // UNIMOD:665
    true,   // UNIMOD:666
    true,   // UNIMOD:667
    true,   // UNIMOD:668
    true,   // UNIMOD:669
    true,   // UNIMOD:670
    true,   // UNIMOD:671
    true,   // UNIMOD:672
    true,   // UNIMOD:673
    true,   // UNIMOD:674
    true,   // UNIMOD:675
    true,   // UNIMOD:676
    true,   // UNIMOD:677
    true,   // UNIMOD:678
    true,   // UNIMOD:679
    true,   // UNIMOD:680
    true,   // UNIMOD:681
    true,   // UNIMOD:682
    true,   // UNIMOD:683
    true,   // UNIMOD:684
    true,   // UNIMOD:685
    true,   // UNIMOD:686
    false,  // UNIMOD:687: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    false,  // UNIMOD:688: not present in source Unimod snapshot
    false,  // UNIMOD:689: not present in source Unimod snapshot
    false,  // UNIMOD:690: not present in source Unimod snapshot
    false,  // UNIMOD:691: not present in source Unimod snapshot
    false,  // UNIMOD:692: not present in source Unimod snapshot
    false,  // UNIMOD:693: not present in source Unimod snapshot
    false,  // UNIMOD:694: not present in source Unimod snapshot
    false,  // UNIMOD:695: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:696: isotope-labeled token 2H(9); pinned-isotope conversion not supported
    true,   // UNIMOD:697
    false,  // UNIMOD:698: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    false,  // UNIMOD:699: not present in source Unimod snapshot
    false,  // UNIMOD:700: not present in source Unimod snapshot
    false,  // UNIMOD:701: not present in source Unimod snapshot
    false,  // UNIMOD:702: not present in source Unimod snapshot
    false,  // UNIMOD:703: not present in source Unimod snapshot
    false,  // UNIMOD:704: not present in source Unimod snapshot
    false,  // UNIMOD:705: not present in source Unimod snapshot
    false,  // UNIMOD:706: not present in source Unimod snapshot
    false,  // UNIMOD:707: not present in source Unimod snapshot
    false,  // UNIMOD:708: not present in source Unimod snapshot
    false,  // UNIMOD:709: not present in source Unimod snapshot
    false,  // UNIMOD:710: not present in source Unimod snapshot
    false,  // UNIMOD:711: not present in source Unimod snapshot
    false,  // UNIMOD:712: not present in source Unimod snapshot
    false,  // UNIMOD:713: not present in source Unimod snapshot
    false,  // UNIMOD:714: not present in source Unimod snapshot
    false,  // UNIMOD:715: not present in source Unimod snapshot
    false,  // UNIMOD:716: not present in source Unimod snapshot
    false,  // UNIMOD:717: not present in source Unimod snapshot
    false,  // UNIMOD:718: not present in source Unimod snapshot
    false,  // UNIMOD:719: not present in source Unimod snapshot
    true,   // UNIMOD:720
    true,   // UNIMOD:721
    false,  // UNIMOD:722: not present in source Unimod snapshot
    true,   // UNIMOD:723
    true,   // UNIMOD:724
    true,   // UNIMOD:725
    true,   // UNIMOD:726
    true,   // UNIMOD:727
    true,   // UNIMOD:728
    true,   // UNIMOD:729
    false,  // UNIMOD:730: isotope-labeled token 13C(7); pinned-isotope conversion not supported
    false,  // UNIMOD:731: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:732: not present in source Unimod snapshot
    false,  // UNIMOD:733: not present in source Unimod snapshot
    true,   // UNIMOD:734
    true,   // UNIMOD:735
    true,   // UNIMOD:736
    false,  // UNIMOD:737: isotope-labeled token 13C(4); pinned-isotope conversion not supported
    false,  // UNIMOD:738: isotope-labeled token 13C; pinned-isotope conversion not supported
    true,   // UNIMOD:739
    false,  // UNIMOD:740: isotope-labeled token 13C(12); pinned-isotope conversion not supported
    false,  // UNIMOD:741: isotope-labeled token 13C(12); pinned-isotope conversion not supported
    false,  // UNIMOD:742: not present in source Unimod snapshot
    true,   // UNIMOD:743
    true,   // UNIMOD:744
    false,  // UNIMOD:745: not present in source Unimod snapshot
    true,   // UNIMOD:746
    true,   // UNIMOD:747
    true,   // UNIMOD:748
    false,  // UNIMOD:749: not present in source Unimod snapshot
    true,   // UNIMOD:750
    true,   // UNIMOD:751
    false,  // UNIMOD:752: not present in source Unimod snapshot
    false,  // UNIMOD:753: not present in source Unimod snapshot
    false,  // UNIMOD:754: not present in source Unimod snapshot
    false,  // UNIMOD:755: not present in source Unimod snapshot
    false,  // UNIMOD:756: not present in source Unimod snapshot
    false,  // UNIMOD:757: not present in source Unimod snapshot
    false,  // UNIMOD:758: not present in source Unimod snapshot
    false,  // UNIMOD:759: not present in source Unimod snapshot
    false,  // UNIMOD:760: not present in source Unimod snapshot
    false,  // UNIMOD:761: not present in source Unimod snapshot
    true,   // UNIMOD:762
    false,  // UNIMOD:763: isotope-labeled token 2H(6); pinned-isotope conversion not supported
    false,  // UNIMOD:764: isotope-labeled token 2H(6); pinned-isotope conversion not supported
    true,   // UNIMOD:765
    true,   // UNIMOD:766
    true,   // UNIMOD:767
    false,  // UNIMOD:768: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    false,  // UNIMOD:769: not present in source Unimod snapshot
    false,  // UNIMOD:770: not present in source Unimod snapshot
    true,   // UNIMOD:771
    false,  // UNIMOD:772: isotope-labeled token 13C(5); pinned-isotope conversion not supported
    true,   // UNIMOD:773
    true,   // UNIMOD:774
    false,  // UNIMOD:775: isotope-labeled token 13C(2); pinned-isotope conversion not supported
    false,  // UNIMOD:776: isotope-labeled token 2H(5); pinned-isotope conversion not supported
    false,  // UNIMOD:777: not present in source Unimod snapshot
    false,  // UNIMOD:778: not present in source Unimod snapshot
    false,  // UNIMOD:779: not present in source Unimod snapshot
    false,  // UNIMOD:780: not present in source Unimod snapshot
    false,  // UNIMOD:781: not present in source Unimod snapshot
    false,  // UNIMOD:782: not present in source Unimod snapshot
    false,  // UNIMOD:783: not present in source Unimod snapshot
    false,  // UNIMOD:784: not present in source Unimod snapshot
    false,  // UNIMOD:785: not present in source Unimod snapshot
    false,  // UNIMOD:786: not present in source Unimod snapshot
    false,  // UNIMOD:787: not present in source Unimod snapshot
    false,  // UNIMOD:788: not present in source Unimod snapshot
    false,  // UNIMOD:789: not present in source Unimod snapshot
    false,  // UNIMOD:790: not present in source Unimod snapshot
    false,  // UNIMOD:791: not present in source Unimod snapshot
    false,  // UNIMOD:792: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    false,  // UNIMOD:793: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:794: not present in source Unimod snapshot
    false,  // UNIMOD:795: not present in source Unimod snapshot
    false,  // UNIMOD:796: not present in source Unimod snapshot
    false,  // UNIMOD:797: not present in source Unimod snapshot
    false,  // UNIMOD:798: not present in source Unimod snapshot
    false,  // UNIMOD:799: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    true,   // UNIMOD:800
    true,   // UNIMOD:801
    false,  // UNIMOD:802: not present in source Unimod snapshot
    false,  // UNIMOD:803: not present in source Unimod snapshot
    false,  // UNIMOD:804: not present in source Unimod snapshot
    false,  // UNIMOD:805: not present in source Unimod snapshot
    false,  // UNIMOD:806: not present in source Unimod snapshot
    false,  // UNIMOD:807: not present in source Unimod snapshot
    false,  // UNIMOD:808: not present in source Unimod snapshot
    false,  // UNIMOD:809: not present in source Unimod snapshot
    false,  // UNIMOD:810: not present in source Unimod snapshot
    true,   // UNIMOD:811
    false,  // UNIMOD:812: not present in source Unimod snapshot
    false,  // UNIMOD:813: not present in source Unimod snapshot
    false,  // UNIMOD:814: not present in source Unimod snapshot
    false,  // UNIMOD:815: not present in source Unimod snapshot
    false,  // UNIMOD:816: not present in source Unimod snapshot
    false,  // UNIMOD:817: not present in source Unimod snapshot
    false,  // UNIMOD:818: not present in source Unimod snapshot
    false,  // UNIMOD:819: not present in source Unimod snapshot
    false,  // UNIMOD:820: not present in source Unimod snapshot
    true,   // UNIMOD:821
    true,   // UNIMOD:822
    false,  // UNIMOD:823: not present in source Unimod snapshot
    true,   // UNIMOD:824
    true,   // UNIMOD:825
    false,  // UNIMOD:826: not present in source Unimod snapshot
    true,   // UNIMOD:827
    false,  // UNIMOD:828: not present in source Unimod snapshot
    false,  // UNIMOD:829: not present in source Unimod snapshot
    true,   // UNIMOD:830
    false,  // UNIMOD:831: not present in source Unimod snapshot
    false,  // UNIMOD:832: not present in source Unimod snapshot
    false,  // UNIMOD:833: not present in source Unimod snapshot
    false,  // UNIMOD:834: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    false,  // UNIMOD:835: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:836: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    true,   // UNIMOD:837
    false,  // UNIMOD:838: not present in source Unimod snapshot
    false,  // UNIMOD:839: not present in source Unimod snapshot
    false,  // UNIMOD:840: not present in source Unimod snapshot
    false,  // UNIMOD:841: not present in source Unimod snapshot
    false,  // UNIMOD:842: not present in source Unimod snapshot
    false,  // UNIMOD:843: not present in source Unimod snapshot
    false,  // UNIMOD:844: not present in source Unimod snapshot
    false,  // UNIMOD:845: not present in source Unimod snapshot
    true,   // UNIMOD:846
    false,  // UNIMOD:847: not present in source Unimod snapshot
    true,   // UNIMOD:848
    true,   // UNIMOD:849
    false,  // UNIMOD:850: not present in source Unimod snapshot
    true,   // UNIMOD:851
    false,  // UNIMOD:852: not present in source Unimod snapshot
    false,  // UNIMOD:853: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    false,  // UNIMOD:854: not present in source Unimod snapshot
    false,  // UNIMOD:855: not present in source Unimod snapshot
    false,  // UNIMOD:856: not present in source Unimod snapshot
    false,  // UNIMOD:857: not present in source Unimod snapshot
    false,  // UNIMOD:858: not present in source Unimod snapshot
    true,   // UNIMOD:859
    true,   // UNIMOD:860
    true,   // UNIMOD:861
    false,  // UNIMOD:862: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    false,  // UNIMOD:863: not present in source Unimod snapshot
    false,  // UNIMOD:864: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:865: not present in source Unimod snapshot
    false,  // UNIMOD:866: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    false,  // UNIMOD:867: not present in source Unimod snapshot
    false,  // UNIMOD:868: not present in source Unimod snapshot
    false,  // UNIMOD:869: not present in source Unimod snapshot
    false,  // UNIMOD:870: not present in source Unimod snapshot
    false,  // UNIMOD:871: not present in source Unimod snapshot
    false,  // UNIMOD:872: not present in source Unimod snapshot
    false,  // UNIMOD:873: not present in source Unimod snapshot
    false,  // UNIMOD:874: not present in source Unimod snapshot
    false,  // UNIMOD:875: not present in source Unimod snapshot
    true,   // UNIMOD:876
    true,   // UNIMOD:877
    false,  // UNIMOD:878: not present in source Unimod snapshot
    false,  // UNIMOD:879: not present in source Unimod snapshot
    false,  // UNIMOD:880: not present in source Unimod snapshot
    false,  // UNIMOD:881: not present in source Unimod snapshot
    false,  // UNIMOD:882: not present in source Unimod snapshot
    false,  // UNIMOD:883: not present in source Unimod snapshot
    true,   // UNIMOD:884
    false,  // UNIMOD:885: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    true,   // UNIMOD:886
    true,   // UNIMOD:887
    true,   // UNIMOD:888
    false,  // UNIMOD:889: isotope-labeled token 13C(3); pinned-isotope conversion not supported
    true,   // UNIMOD:890
    true,   // UNIMOD:891
    false,  // UNIMOD:892: not present in source Unimod snapshot
    true,   // UNIMOD:893
    true,   // UNIMOD:894
    true,   // UNIMOD:895
    true,   // UNIMOD:896
    false,  // UNIMOD:897: isotope-labeled token 15N(4); pinned-isotope conversion not supported
    true,   // UNIMOD:898
    true,   // UNIMOD:899
    false,  // UNIMOD:900: not present in source Unimod snapshot
    true,   // UNIMOD:901
    true,   // UNIMOD:902
    true,   // UNIMOD:903
    true,   // UNIMOD:904
    true,   // UNIMOD:905
    true,   // UNIMOD:906
    false,  // UNIMOD:907: non-element token Hex; group-to-formula expansion not supported
    true,   // UNIMOD:908
    false,  // UNIMOD:909: not present in source Unimod snapshot
    false,  // UNIMOD:910: non-element token dHex; group-to-formula expansion not supported
    true,   // UNIMOD:911
    true,   // UNIMOD:912
    false,  // UNIMOD:913: not present in source Unimod snapshot
    true,   // UNIMOD:914
    false,  // UNIMOD:915: not present in source Unimod snapshot
    false,  // UNIMOD:916: not present in source Unimod snapshot
    false,  // UNIMOD:917: not present in source Unimod snapshot
    false,  // UNIMOD:918: not present in source Unimod snapshot
    false,  // UNIMOD:919: not present in source Unimod snapshot
    false,  // UNIMOD:920: not present in source Unimod snapshot
    false,  // UNIMOD:921: not present in source Unimod snapshot
    false,  // UNIMOD:922: not present in source Unimod snapshot
    false,  // UNIMOD:923: isotope-labeled token 13C(4); pinned-isotope conversion not supported
    false,  // UNIMOD:924: not present in source Unimod snapshot
    false,  // UNIMOD:925: not present in source Unimod snapshot
    true,   // UNIMOD:926
    false,  // UNIMOD:927: not present in source Unimod snapshot
    true,   // UNIMOD:928
    false,  // UNIMOD:929: not present in source Unimod snapshot
    false,  // UNIMOD:930: not present in source Unimod snapshot
    true,   // UNIMOD:931
    true,   // UNIMOD:932
    true,   // UNIMOD:933
    true,   // UNIMOD:934
    true,   // UNIMOD:935
    true,   // UNIMOD:936
    true,   // UNIMOD:937
    true,   // UNIMOD:938
    true,   // UNIMOD:939
    true,   // UNIMOD:940
    true,   // UNIMOD:941
    true,   // UNIMOD:942
    true,   // UNIMOD:943
    false,  // UNIMOD:944: isotope-labeled token 2H(9); pinned-isotope conversion not supported
    false,  // UNIMOD:945: not present in source Unimod snapshot
    true,   // UNIMOD:946
    true,   // UNIMOD:947
    true,   // UNIMOD:948
    true,   // UNIMOD:949
    true,   // UNIMOD:950
    true,   // UNIMOD:951
    true,   // UNIMOD:952
    true,   // UNIMOD:953
    true,   // UNIMOD:954
    true,   // UNIMOD:955
    true,   // UNIMOD:956
    true,   // UNIMOD:957
    true,   // UNIMOD:958
    true,   // UNIMOD:959
    true,   // UNIMOD:960
    true,   // UNIMOD:961
    false,  // UNIMOD:962: not present in source Unimod snapshot
    false,  // UNIMOD:963: not present in source Unimod snapshot
    false,  // UNIMOD:964: not present in source Unimod snapshot
    false,  // UNIMOD:965: not present in source Unimod snapshot
    false,  // UNIMOD:966: not present in source Unimod snapshot
    true,   // UNIMOD:967
    false,  // UNIMOD:968: not present in source Unimod snapshot
    false,  // UNIMOD:969: not present in source Unimod snapshot
    false,  // UNIMOD:970: not present in source Unimod snapshot
    true,   // UNIMOD:971
    true,   // UNIMOD:972
    true,   // UNIMOD:973
    false,  // UNIMOD:974: not present in source Unimod snapshot
    false,  // UNIMOD:975: not present in source Unimod snapshot
    false,  // UNIMOD:976: not present in source Unimod snapshot
    true,   // UNIMOD:977
    true,   // UNIMOD:978
    true,   // UNIMOD:979
    false,  // UNIMOD:980: not present in source Unimod snapshot
    true,   // UNIMOD:981
    false,  // UNIMOD:982: not present in source Unimod snapshot
    false,  // UNIMOD:983: not present in source Unimod snapshot
    true,   // UNIMOD:984
    false,  // UNIMOD:985: isotope-labeled token 13C(4); pinned-isotope conversion not supported
    false,  // UNIMOD:986: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:987: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:988: not present in source Unimod snapshot
    true,   // UNIMOD:989
    false,  // UNIMOD:990: not present in source Unimod snapshot
    true,   // UNIMOD:991
    false,  // UNIMOD:992: not present in source Unimod snapshot
    true,   // UNIMOD:993
    false,  // UNIMOD:994: isotope-labeled token 15N; pinned-isotope conversion not supported
    false,  // UNIMOD:995: isotope-labeled token 15N(2); pinned-isotope conversion not supported
    false,  // UNIMOD:996: isotope-labeled token 15N(3); pinned-isotope conversion not supported
    true,   // UNIMOD:997
    false,  // UNIMOD:998: not present in source Unimod snapshot
    false,  // UNIMOD:999: not present in source Unimod snapshot
    true,   // UNIMOD:1000
    true,   // UNIMOD:1001
    true,   // UNIMOD:1002
    true,   // UNIMOD:1003
    false,  // UNIMOD:1004: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:1005: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:1006: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    false,  // UNIMOD:1007: isotope-labeled token 2H(6); pinned-isotope conversion not supported
    true,   // UNIMOD:1008
    true,   // UNIMOD:1009
    true,   // UNIMOD:1010
    false,  // UNIMOD:1011: not present in source Unimod snapshot
    true,   // UNIMOD:1012
    false,  // UNIMOD:1013: not present in source Unimod snapshot
    true,   // UNIMOD:1014
    true,   // UNIMOD:1015
    false,  // UNIMOD:1016: not present in source Unimod snapshot
    true,   // UNIMOD:1017
    true,   // UNIMOD:1018
    false,  // UNIMOD:1019: isotope-labeled token 2H(6); pinned-isotope conversion not supported
    true,   // UNIMOD:1020
    true,   // UNIMOD:1021
    true,   // UNIMOD:1022
    true,   // UNIMOD:1023
    true,   // UNIMOD:1024
    false,  // UNIMOD:1025: not present in source Unimod snapshot
    false,  // UNIMOD:1026: not present in source Unimod snapshot
    true,   // UNIMOD:1027
    true,   // UNIMOD:1028
    false,  // UNIMOD:1029: not present in source Unimod snapshot
    false,  // UNIMOD:1030: not present in source Unimod snapshot
    true,   // UNIMOD:1031
    true,   // UNIMOD:1032
    true,   // UNIMOD:1033
    false,  // UNIMOD:1034: isotope-labeled token 2H(5); pinned-isotope conversion not supported
    true,   // UNIMOD:1035
    true,   // UNIMOD:1036
    true,   // UNIMOD:1037
    true,   // UNIMOD:1038
    true,   // UNIMOD:1039
    false,  // UNIMOD:1040: not present in source Unimod snapshot
    true,   // UNIMOD:1041
    true,   // UNIMOD:1042
    true,   // UNIMOD:1043
    true,   // UNIMOD:1044
    true,   // UNIMOD:1045
    true,   // UNIMOD:1046
    true,   // UNIMOD:1047
    true,   // UNIMOD:1048
    true,   // UNIMOD:1049
    true,   // UNIMOD:1050
    true,   // UNIMOD:1051
    true,   // UNIMOD:1052
    true,   // UNIMOD:1053
    true,   // UNIMOD:1054
    true,   // UNIMOD:1055
    true,   // UNIMOD:1056
    true,   // UNIMOD:1057
    true,   // UNIMOD:1058
    true,   // UNIMOD:1059
    true,   // UNIMOD:1060
    true,   // UNIMOD:1061
    true,   // UNIMOD:1062
    true,   // UNIMOD:1063
    true,   // UNIMOD:1064
    true,   // UNIMOD:1065
    true,   // UNIMOD:1066
    true,   // UNIMOD:1067
    true,   // UNIMOD:1068
    true,   // UNIMOD:1069
    true,   // UNIMOD:1070
    true,   // UNIMOD:1071
    true,   // UNIMOD:1072
    true,   // UNIMOD:1073
    true,   // UNIMOD:1074
    true,   // UNIMOD:1075
    true,   // UNIMOD:1076
    true,   // UNIMOD:1077
    true,   // UNIMOD:1078
    true,   // UNIMOD:1079
    true,   // UNIMOD:1080
    true,   // UNIMOD:1081
    true,   // UNIMOD:1082
    true,   // UNIMOD:1083
    true,   // UNIMOD:1084
    true,   // UNIMOD:1085
    true,   // UNIMOD:1086
    true,   // UNIMOD:1087
    true,   // UNIMOD:1088
    true,   // UNIMOD:1089
    true,   // UNIMOD:1090
    true,   // UNIMOD:1091
    true,   // UNIMOD:1092
    true,   // UNIMOD:1093
    true,   // UNIMOD:1094
    true,   // UNIMOD:1095
    true,   // UNIMOD:1096
    true,   // UNIMOD:1097
    true,   // UNIMOD:1098
    true,   // UNIMOD:1099
    true,   // UNIMOD:1100
    true,   // UNIMOD:1101
    true,   // UNIMOD:1102
    true,   // UNIMOD:1103
    true,   // UNIMOD:1104
    true,   // UNIMOD:1105
    true,   // UNIMOD:1106
    true,   // UNIMOD:1107
    true,   // UNIMOD:1108
    true,   // UNIMOD:1109
    true,   // UNIMOD:1110
    true,   // UNIMOD:1111
    true,   // UNIMOD:1112
    true,   // UNIMOD:1113
    true,   // UNIMOD:1114
    true,   // UNIMOD:1115
    true,   // UNIMOD:1116
    true,   // UNIMOD:1117
    false,  // UNIMOD:1118: not present in source Unimod snapshot
    true,   // UNIMOD:1119
    true,   // UNIMOD:1120
    true,   // UNIMOD:1121
    true,   // UNIMOD:1122
    true,   // UNIMOD:1123
    true,   // UNIMOD:1124
    true,   // UNIMOD:1125
    true,   // UNIMOD:1126
    true,   // UNIMOD:1127
    true,   // UNIMOD:1128
    true,   // UNIMOD:1129
    true,   // UNIMOD:1130
    true,   // UNIMOD:1131
    true,   // UNIMOD:1132
    true,   // UNIMOD:1133
    true,   // UNIMOD:1134
    true,   // UNIMOD:1135
    true,   // UNIMOD:1136
    true,   // UNIMOD:1137
    true,   // UNIMOD:1138
    true,   // UNIMOD:1139
    true,   // UNIMOD:1140
    true,   // UNIMOD:1141
    true,   // UNIMOD:1142
    true,   // UNIMOD:1143
    true,   // UNIMOD:1144
    true,   // UNIMOD:1145
    true,   // UNIMOD:1146
    true,   // UNIMOD:1147
    true,   // UNIMOD:1148
    true,   // UNIMOD:1149
    true,   // UNIMOD:1150
    true,   // UNIMOD:1151
    true,   // UNIMOD:1152
    true,   // UNIMOD:1153
    true,   // UNIMOD:1154
    true,   // UNIMOD:1155
    true,   // UNIMOD:1156
    true,   // UNIMOD:1157
    true,   // UNIMOD:1158
    true,   // UNIMOD:1159
    true,   // UNIMOD:1160
    true,   // UNIMOD:1161
    true,   // UNIMOD:1162
    true,   // UNIMOD:1163
    true,   // UNIMOD:1164
    true,   // UNIMOD:1165
    true,   // UNIMOD:1166
    true,   // UNIMOD:1167
    true,   // UNIMOD:1168
    true,   // UNIMOD:1169
    true,   // UNIMOD:1170
    true,   // UNIMOD:1171
    true,   // UNIMOD:1172
    true,   // UNIMOD:1173
    true,   // UNIMOD:1174
    true,   // UNIMOD:1175
    true,   // UNIMOD:1176
    true,   // UNIMOD:1177
    true,   // UNIMOD:1178
    true,   // UNIMOD:1179
    true,   // UNIMOD:1180
    true,   // UNIMOD:1181
    true,   // UNIMOD:1182
    true,   // UNIMOD:1183
    true,   // UNIMOD:1184
    true,   // UNIMOD:1185
    true,   // UNIMOD:1186
    true,   // UNIMOD:1187
    true,   // UNIMOD:1188
    true,   // UNIMOD:1189
    true,   // UNIMOD:1190
    true,   // UNIMOD:1191
    true,   // UNIMOD:1192
    true,   // UNIMOD:1193
    true,   // UNIMOD:1194
    true,   // UNIMOD:1195
    true,   // UNIMOD:1196
    true,   // UNIMOD:1197
    true,   // UNIMOD:1198
    true,   // UNIMOD:1199
    true,   // UNIMOD:1200
    true,   // UNIMOD:1201
    true,   // UNIMOD:1202
    true,   // UNIMOD:1203
    true,   // UNIMOD:1204
    true,   // UNIMOD:1205
    true,   // UNIMOD:1206
    true,   // UNIMOD:1207
    true,   // UNIMOD:1208
    true,   // UNIMOD:1209
    true,   // UNIMOD:1210
    true,   // UNIMOD:1211
    true,   // UNIMOD:1212
    true,   // UNIMOD:1213
    true,   // UNIMOD:1214
    true,   // UNIMOD:1215
    true,   // UNIMOD:1216
    true,   // UNIMOD:1217
    true,   // UNIMOD:1218
    true,   // UNIMOD:1219
    true,   // UNIMOD:1220
    true,   // UNIMOD:1221
    true,   // UNIMOD:1222
    true,   // UNIMOD:1223
    true,   // UNIMOD:1224
    true,   // UNIMOD:1225
    true,   // UNIMOD:1226
    true,   // UNIMOD:1227
    true,   // UNIMOD:1228
    true,   // UNIMOD:1229
    true,   // UNIMOD:1230
    true,   // UNIMOD:1231
    true,   // UNIMOD:1232
    true,   // UNIMOD:1233
    true,   // UNIMOD:1234
    true,   // UNIMOD:1235
    true,   // UNIMOD:1236
    true,   // UNIMOD:1237
    true,   // UNIMOD:1238
    true,   // UNIMOD:1239
    true,   // UNIMOD:1240
    true,   // UNIMOD:1241
    true,   // UNIMOD:1242
    true,   // UNIMOD:1243
    true,   // UNIMOD:1244
    true,   // UNIMOD:1245
    true,   // UNIMOD:1246
    true,   // UNIMOD:1247
    true,   // UNIMOD:1248
    true,   // UNIMOD:1249
    true,   // UNIMOD:1250
    true,   // UNIMOD:1251
    true,   // UNIMOD:1252
    true,   // UNIMOD:1253
    true,   // UNIMOD:1254
    true,   // UNIMOD:1255
    true,   // UNIMOD:1256
    true,   // UNIMOD:1257
    true,   // UNIMOD:1258
    false,  // UNIMOD:1259: not present in source Unimod snapshot
    false,  // UNIMOD:1260: not present in source Unimod snapshot
    true,   // UNIMOD:1261
    true,   // UNIMOD:1262
    true,   // UNIMOD:1263
    true,   // UNIMOD:1264
    false,  // UNIMOD:1265: not present in source Unimod snapshot
    false,  // UNIMOD:1266: isotope-labeled token 13C(4); pinned-isotope conversion not supported
    false,  // UNIMOD:1267: isotope-labeled token 13C(4); pinned-isotope conversion not supported
    false,  // UNIMOD:1268: not present in source Unimod snapshot
    false,  // UNIMOD:1269: not present in source Unimod snapshot
    true,   // UNIMOD:1270
    true,   // UNIMOD:1271
    false,  // UNIMOD:1272: not present in source Unimod snapshot
    false,  // UNIMOD:1273: not present in source Unimod snapshot
    false,  // UNIMOD:1274: not present in source Unimod snapshot
    false,  // UNIMOD:1275: not present in source Unimod snapshot
    true,   // UNIMOD:1276
    true,   // UNIMOD:1277
    true,   // UNIMOD:1278
    true,   // UNIMOD:1279
    false,  // UNIMOD:1280: not present in source Unimod snapshot
    true,   // UNIMOD:1281
    true,   // UNIMOD:1282
    true,   // UNIMOD:1283
    false,  // UNIMOD:1284: not present in source Unimod snapshot
    false,  // UNIMOD:1285: not present in source Unimod snapshot
    false,  // UNIMOD:1286: non-element token Hex(2); group-to-formula expansion not supported
    true,   // UNIMOD:1287
    true,   // UNIMOD:1288
    true,   // UNIMOD:1289
    true,   // UNIMOD:1290
    false,  // UNIMOD:1291: isotope-labeled token 2H(6); pinned-isotope conversion not supported
    true,   // UNIMOD:1292
    true,   // UNIMOD:1293
    false,  // UNIMOD:1294: not present in source Unimod snapshot
    false,  // UNIMOD:1295: not present in source Unimod snapshot
    false,  // UNIMOD:1296: isotope-labeled token 13C(3); pinned-isotope conversion not supported
    false,  // UNIMOD:1297: isotope-labeled token 13C(3); pinned-isotope conversion not supported
    false,  // UNIMOD:1298: isotope-labeled token 13C(4); pinned-isotope conversion not supported
    false,  // UNIMOD:1299: isotope-labeled token 2H(10); pinned-isotope conversion not supported
    false,  // UNIMOD:1300: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    true,   // UNIMOD:1301
    false,  // UNIMOD:1302: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    false,  // UNIMOD:1303: non-element token NeuAc; group-to-formula expansion not supported
    false,  // UNIMOD:1304: non-element token NeuGc; group-to-formula expansion not supported
    true,   // UNIMOD:1305
    false,  // UNIMOD:1306: isotope-labeled token 2H(6); pinned-isotope conversion not supported
    false,  // UNIMOD:1307: not present in source Unimod snapshot
    false,  // UNIMOD:1308: not present in source Unimod snapshot
    false,  // UNIMOD:1309: not present in source Unimod snapshot
    true,   // UNIMOD:1310
    false,  // UNIMOD:1311: not present in source Unimod snapshot
    true,   // UNIMOD:1312
    true,   // UNIMOD:1313
    true,   // UNIMOD:1314
    true,   // UNIMOD:1315
    false,  // UNIMOD:1316: not present in source Unimod snapshot
    true,   // UNIMOD:1317
    false,  // UNIMOD:1318: not present in source Unimod snapshot
    false,  // UNIMOD:1319: not present in source Unimod snapshot
    true,   // UNIMOD:1320
    false,  // UNIMOD:1321: isotope-labeled token 13C; pinned-isotope conversion not supported
    false,  // UNIMOD:1322: isotope-labeled token 2H(2); pinned-isotope conversion not supported
    false,  // UNIMOD:1323: isotope-labeled token 2H(2); pinned-isotope conversion not supported
    false,  // UNIMOD:1324: isotope-labeled token 2H(4); pinned-isotope conversion not supported
    false,  // UNIMOD:1325: not present in source Unimod snapshot
    true,   // UNIMOD:1326
    true,   // UNIMOD:1327
    true,   // UNIMOD:1328
    false,  // UNIMOD:1329: not present in source Unimod snapshot
    true,   // UNIMOD:1330
    true,   // UNIMOD:1331
    false,  // UNIMOD:1332: not present in source Unimod snapshot
    false,  // UNIMOD:1333: not present in source Unimod snapshot
    false,  // UNIMOD:1334: not present in source Unimod snapshot
    false,  // UNIMOD:1335: not present in source Unimod snapshot
    false,  // UNIMOD:1336: not present in source Unimod snapshot
    false,  // UNIMOD:1337: not present in source Unimod snapshot
    false,  // UNIMOD:1338: not present in source Unimod snapshot
    false,  // UNIMOD:1339: not present in source Unimod snapshot
    true,   // UNIMOD:1340
    true,   // UNIMOD:1341
    false,  // UNIMOD:1342: isotope-labeled token 13C(4); pinned-isotope conversion not supported
    false,  // UNIMOD:1343: not present in source Unimod snapshot
    true,   // UNIMOD:1344
    true,   // UNIMOD:1345
    false,  // UNIMOD:1346: not present in source Unimod snapshot
    false,  // UNIMOD:1347: not present in source Unimod snapshot
    true,   // UNIMOD:1348
    true,   // UNIMOD:1349
    true,   // UNIMOD:1350
    false,  // UNIMOD:1351: not present in source Unimod snapshot
    false,  // UNIMOD:1352: not present in source Unimod snapshot
    false,  // UNIMOD:1353: not present in source Unimod snapshot
    false,  // UNIMOD:1354: not present in source Unimod snapshot
    true,   // UNIMOD:1355
    true,   // UNIMOD:1356
    false,  // UNIMOD:1357: not present in source Unimod snapshot
    false,  // UNIMOD:1358: isotope-labeled token 2H(5); pinned-isotope conversion not supported
    false,  // UNIMOD:1359: not present in source Unimod snapshot
    false,  // UNIMOD:1360: not present in source Unimod snapshot
    false,  // UNIMOD:1361: not present in source Unimod snapshot
    false,  // UNIMOD:1362: not present in source Unimod snapshot
    true,   // UNIMOD:1363
    true,   // UNIMOD:1364
    true,   // UNIMOD:1365
    false,  // UNIMOD:1366: not present in source Unimod snapshot
    false,  // UNIMOD:1367: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1368: isotope-labeled token 2H(6); pinned-isotope conversion not supported
    false,  // UNIMOD:1369: not present in source Unimod snapshot
    false,  // UNIMOD:1370: isotope-labeled token 2H(3); pinned-isotope conversion not supported
    false,  // UNIMOD:1371: isotope-labeled token 2H(9); pinned-isotope conversion not supported
    false,  // UNIMOD:1372: isotope-labeled token 13C(2); pinned-isotope conversion not supported
    false,  // UNIMOD:1373: not present in source Unimod snapshot
    false,  // UNIMOD:1374: not present in source Unimod snapshot
    false,  // UNIMOD:1375: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1376: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1377: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1378: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1379: non-element token dHex; group-to-formula expansion not supported
    true,   // UNIMOD:1380
    true,   // UNIMOD:1381
    true,   // UNIMOD:1382
    true,   // UNIMOD:1383
    true,   // UNIMOD:1384
    true,   // UNIMOD:1385
    false,  // UNIMOD:1386: not present in source Unimod snapshot
    true,   // UNIMOD:1387
    true,   // UNIMOD:1388
    true,   // UNIMOD:1389
    true,   // UNIMOD:1390
    true,   // UNIMOD:1391
    false,  // UNIMOD:1392: isotope-labeled token 13C(4); pinned-isotope conversion not supported
    false,  // UNIMOD:1393: isotope-labeled token 13C(3); pinned-isotope conversion not supported
    false,  // UNIMOD:1394: isotope-labeled token 2H(2); pinned-isotope conversion not supported
    false,  // UNIMOD:1395: isotope-labeled token 2H(2); pinned-isotope conversion not supported
    false,  // UNIMOD:1396: isotope-labeled token 2H(2); pinned-isotope conversion not supported
    true,   // UNIMOD:1397
    false,  // UNIMOD:1398: isotope-labeled token 13C(6); pinned-isotope conversion not supported
    true,   // UNIMOD:1399
    false,  // UNIMOD:1400: non-element token NeuAc; group-to-formula expansion not supported
    false,  // UNIMOD:1401: not present in source Unimod snapshot
    false,  // UNIMOD:1402: isotope-labeled token 2H(7); pinned-isotope conversion not supported
    false,  // UNIMOD:1403: isotope-labeled token 2H(6); pinned-isotope conversion not supported
    false,  // UNIMOD:1404: not present in source Unimod snapshot
    true,   // UNIMOD:1405
    true,   // UNIMOD:1406
    false,  // UNIMOD:1407: not present in source Unimod snapshot
    false,  // UNIMOD:1408: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1409: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1410: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1411: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1412: non-element token HexNAc; group-to-formula expansion not supported
    false,  // UNIMOD:1413: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1414: isotope-labeled token 2H(9); pinned-isotope conversion not supported
    false,  // UNIMOD:1415: not present in source Unimod snapshot
    false,  // UNIMOD:1416: not present in source Unimod snapshot
    false,  // UNIMOD:1417: not present in source Unimod snapshot
    false,  // UNIMOD:1418: not present in source Unimod snapshot
    false,  // UNIMOD:1419: isotope-labeled token 15N(-1); pinned-isotope conversion not supported
    true,   // UNIMOD:1420
    true,   // UNIMOD:1421
    false,  // UNIMOD:1422: not present in source Unimod snapshot
    true,   // UNIMOD:1423
    false,  // UNIMOD:1424: not present in source Unimod snapshot
    false,  // UNIMOD:1425: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:1426: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:1427: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1428: non-element token Pent(2); group-to-formula expansion not supported
    false,  // UNIMOD:1429: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1430: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1431: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1432: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1433: non-element token HexNAc(3); group-to-formula expansion not supported
    false,  // UNIMOD:1434: non-element token HexNAc; group-to-formula expansion not supported
    false,  // UNIMOD:1435: non-element token HexNAc; group-to-formula expansion not supported
    false,  // UNIMOD:1436: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1437: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1438: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1439: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1440: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1441: non-element token Pent(3); group-to-formula expansion not supported
    false,  // UNIMOD:1442: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:1443: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1444: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1445: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1446: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1447: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1448: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1449: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1450: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1451: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1452: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1453: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1454: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1455: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1456: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1457: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1458: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1459: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1460: non-element token Hex(7); group-to-formula expansion not supported
    false,  // UNIMOD:1461: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1462: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1463: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1464: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1465: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1466: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1467: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1468: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1469: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1470: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1471: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1472: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1473: non-element token Hex(8); group-to-formula expansion not supported
    false,  // UNIMOD:1474: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1475: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1476: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1477: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1478: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1479: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1480: non-element token Hex(7); group-to-formula expansion not supported
    false,  // UNIMOD:1481: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1482: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1483: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1484: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1485: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1486: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1487: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1488: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1489: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1490: non-element token Hex(7); group-to-formula expansion not supported
    false,  // UNIMOD:1491: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1492: non-element token Pent(3); group-to-formula expansion not supported
    false,  // UNIMOD:1493: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1494: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1495: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1496: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1497: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1498: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1499: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1500: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1501: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1502: non-element token Hex(7); group-to-formula expansion not supported
    false,  // UNIMOD:1503: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1504: non-element token Hex(8); group-to-formula expansion not supported
    false,  // UNIMOD:1505: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1506: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1507: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1508: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1509: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1510: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1511: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1512: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1513: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1514: non-element token Hex(7); group-to-formula expansion not supported
    false,  // UNIMOD:1515: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1516: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1517: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1518: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1519: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1520: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1521: non-element token Hex(7); group-to-formula expansion not supported
    false,  // UNIMOD:1522: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1523: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1524: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1525: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1526: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1527: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1528: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1529: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1530: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1531: non-element token Hex(9); group-to-formula expansion not supported
    false,  // UNIMOD:1532: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1533: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1534: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1535: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1536: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1537: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1538: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1539: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1540: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1541: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1542: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1543: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1544: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1545: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1546: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1547: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1548: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1549: non-element token Hex(7); group-to-formula expansion not supported
    false,  // UNIMOD:1550: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1551: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1552: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1553: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1554: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1555: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1556: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1557: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1558: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1559: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1560: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1561: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1562: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1563: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1564: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1565: non-element token HexNAc(3); group-to-formula expansion not supported
    false,  // UNIMOD:1566: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1567: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1568: non-element token HexNAc(2); group-to-formula expansion not supported
    false,  // UNIMOD:1569: not present in source Unimod snapshot
    false,  // UNIMOD:1570: non-element token HexNAc; group-to-formula expansion not supported
    false,  // UNIMOD:1571: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1572: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:1573: non-element token HexNAc(2); group-to-formula expansion not supported
    false,  // UNIMOD:1574: not present in source Unimod snapshot
    false,  // UNIMOD:1575: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1576: not present in source Unimod snapshot
    false,  // UNIMOD:1577: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1578: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1579: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1580: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1581: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1582: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1583: non-element token HexNAc(2); group-to-formula expansion not supported
    false,  // UNIMOD:1584: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1585: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1586: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1587: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1588: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1589: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1590: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1591: non-element token HexNAc(4); group-to-formula expansion not supported
    false,  // UNIMOD:1592: non-element token HexNAc(1); group-to-formula expansion not supported
    false,  // UNIMOD:1593: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1594: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1595: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1596: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1597: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1598: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1599: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1600: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1601: not present in source Unimod snapshot
    false,  // UNIMOD:1602: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1603: not present in source Unimod snapshot
    false,  // UNIMOD:1604: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1605: not present in source Unimod snapshot
    false,  // UNIMOD:1606: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1607: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1608: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1609: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1610: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1611: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1612: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1613: not present in source Unimod snapshot
    false,  // UNIMOD:1614: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1615: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1616: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1617: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1618: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1619: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1620: token Ac has no supported element symbol or group expansion
    false,  // UNIMOD:1621: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1622: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1623: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:1624: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1625: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1626: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1627: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1628: non-element token HexNAc(5); group-to-formula expansion not supported
    false,  // UNIMOD:1629: not present in source Unimod snapshot
    false,  // UNIMOD:1630: token Ac(2) has no supported element symbol or group expansion
    false,  // UNIMOD:1631: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1632: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1633: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1634: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1635: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1636: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1637: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1638: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1639: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1640: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1641: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1642: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1643: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1644: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1645: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1646: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1647: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1648: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1649: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1650: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1651: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1652: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1653: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1654: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1655: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1656: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1657: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1658: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1659: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1660: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1661: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1662: non-element token Hex; group-to-formula expansion not supported
    false,  // UNIMOD:1663: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1664: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1665: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1666: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1667: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1668: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1669: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1670: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1671: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1672: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1673: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1674: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1675: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1676: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1677: not present in source Unimod snapshot
    false,  // UNIMOD:1678: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1679: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1680: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1681: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1682: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1683: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1684: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1685: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1686: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1687: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1688: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1689: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1690: non-element token Hex(7); group-to-formula expansion not supported
    false,  // UNIMOD:1691: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1692: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1693: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1694: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1695: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1696: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1697: non-element token dHex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1698: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1699: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1700: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1701: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1702: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1703: non-element token dHex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1704: not present in source Unimod snapshot
    false,  // UNIMOD:1705: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1706: non-element token dHex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1707: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1708: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1709: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1710: not present in source Unimod snapshot
    false,  // UNIMOD:1711: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1712: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1713: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1714: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1715: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1716: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1717: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1718: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1719: non-element token dHex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1720: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1721: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1722: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1723: non-element token Hex(8); group-to-formula expansion not supported
    false,  // UNIMOD:1724: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1725: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1726: non-element token dHex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1727: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1728: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1729: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1730: non-element token dHex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1731: not present in source Unimod snapshot
    false,  // UNIMOD:1732: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1733: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1734: not present in source Unimod snapshot
    false,  // UNIMOD:1735: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1736: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1737: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1738: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1739: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1740: non-element token dHex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1741: not present in source Unimod snapshot
    false,  // UNIMOD:1742: non-element token Hex(9); group-to-formula expansion not supported
    false,  // UNIMOD:1743: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1744: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1745: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1746: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1747: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1748: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1749: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1750: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1751: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1752: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1753: non-element token Hex(10); group-to-formula expansion not supported
    false,  // UNIMOD:1754: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1755: non-element token Hex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1756: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1757: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1758: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1759: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1760: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1761: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1762: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1763: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1764: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1765: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1766: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1767: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1768: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1769: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1770: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1771: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1772: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1773: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1774: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1775: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1776: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1777: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1778: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1779: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1780: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1781: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1782: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1783: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1784: non-element token dHex(1); group-to-formula expansion not supported
    false,  // UNIMOD:1785: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1786: token Ac has no supported element symbol or group expansion
    false,  // UNIMOD:1787: isotope-labeled token 13C(2); pinned-isotope conversion not supported
    false,  // UNIMOD:1788: not present in source Unimod snapshot
    true,   // UNIMOD:1789
    false,  // UNIMOD:1790: not present in source Unimod snapshot
    false,  // UNIMOD:1791: not present in source Unimod snapshot
    false,  // UNIMOD:1792: not present in source Unimod snapshot
    false,  // UNIMOD:1793: not present in source Unimod snapshot
    false,  // UNIMOD:1794: not present in source Unimod snapshot
    false,  // UNIMOD:1795: not present in source Unimod snapshot
    false,  // UNIMOD:1796: not present in source Unimod snapshot
    false,  // UNIMOD:1797: not present in source Unimod snapshot
    false,  // UNIMOD:1798: not present in source Unimod snapshot
    true,   // UNIMOD:1799
    true,   // UNIMOD:1800
    true,   // UNIMOD:1801
    false,  // UNIMOD:1802: not present in source Unimod snapshot
    false,  // UNIMOD:1803: not present in source Unimod snapshot
    false,  // UNIMOD:1804: not present in source Unimod snapshot
    false,  // UNIMOD:1805: not present in source Unimod snapshot
    false,  // UNIMOD:1806: not present in source Unimod snapshot
    false,  // UNIMOD:1807: not present in source Unimod snapshot
    false,  // UNIMOD:1808: not present in source Unimod snapshot
    false,  // UNIMOD:1809: not present in source Unimod snapshot
    false,  // UNIMOD:1810: not present in source Unimod snapshot
    false,  // UNIMOD:1811: not present in source Unimod snapshot
    false,  // UNIMOD:1812: not present in source Unimod snapshot
    false,  // UNIMOD:1813: not present in source Unimod snapshot
    false,  // UNIMOD:1814: not present in source Unimod snapshot
    false,  // UNIMOD:1815: not present in source Unimod snapshot
    false,  // UNIMOD:1816: not present in source Unimod snapshot
    false,  // UNIMOD:1817: not present in source Unimod snapshot
    false,  // UNIMOD:1818: not present in source Unimod snapshot
    false,  // UNIMOD:1819: not present in source Unimod snapshot
    false,  // UNIMOD:1820: not present in source Unimod snapshot
    false,  // UNIMOD:1821: not present in source Unimod snapshot
    false,  // UNIMOD:1822: not present in source Unimod snapshot
    false,  // UNIMOD:1823: not present in source Unimod snapshot
    true,   // UNIMOD:1824
    true,   // UNIMOD:1825
    true,   // UNIMOD:1826
    false,  // UNIMOD:1827: isotope-labeled token 2H(2); pinned-isotope conversion not supported
    true,   // UNIMOD:1828
    true,   // UNIMOD:1829
    true,   // UNIMOD:1830
    true,   // UNIMOD:1831
    true,   // UNIMOD:1832
    true,   // UNIMOD:1833
    true,   // UNIMOD:1834
    true,   // UNIMOD:1835
    true,   // UNIMOD:1836
    true,   // UNIMOD:1837
    true,   // UNIMOD:1838
    true,   // UNIMOD:1839
    false,  // UNIMOD:1840: non-element token dHex; group-to-formula expansion not supported
    true,   // UNIMOD:1841
    false,  // UNIMOD:1842: not present in source Unimod snapshot
    true,   // UNIMOD:1843
    false,  // UNIMOD:1844: not present in source Unimod snapshot
    true,   // UNIMOD:1845
    true,   // UNIMOD:1846
    false,  // UNIMOD:1847: not present in source Unimod snapshot
    true,   // UNIMOD:1848
    true,   // UNIMOD:1849
    false,  // UNIMOD:1850: not present in source Unimod snapshot
    false,  // UNIMOD:1851: not present in source Unimod snapshot
    false,  // UNIMOD:1852: not present in source Unimod snapshot
    false,  // UNIMOD:1853: not present in source Unimod snapshot
    false,  // UNIMOD:1854: not present in source Unimod snapshot
    false,  // UNIMOD:1855: not present in source Unimod snapshot
    false,  // UNIMOD:1856: not present in source Unimod snapshot
    false,  // UNIMOD:1857: not present in source Unimod snapshot
    false,  // UNIMOD:1858: not present in source Unimod snapshot
    false,  // UNIMOD:1859: not present in source Unimod snapshot
    false,  // UNIMOD:1860: not present in source Unimod snapshot
    false,  // UNIMOD:1861: not present in source Unimod snapshot
    false,  // UNIMOD:1862: not present in source Unimod snapshot
    false,  // UNIMOD:1863: not present in source Unimod snapshot
    false,  // UNIMOD:1864: not present in source Unimod snapshot
    false,  // UNIMOD:1865: not present in source Unimod snapshot
    false,  // UNIMOD:1866: not present in source Unimod snapshot
    false,  // UNIMOD:1867: not present in source Unimod snapshot
    true,   // UNIMOD:1868
    false,  // UNIMOD:1869: not present in source Unimod snapshot
    true,   // UNIMOD:1870
    true,   // UNIMOD:1871
    true,   // UNIMOD:1872
    true,   // UNIMOD:1873
    false,  // UNIMOD:1874: not present in source Unimod snapshot
    true,   // UNIMOD:1875
    false,  // UNIMOD:1876: not present in source Unimod snapshot
    true,   // UNIMOD:1877
    true,   // UNIMOD:1878
    true,   // UNIMOD:1879
    true,   // UNIMOD:1880
    true,   // UNIMOD:1881
    true,   // UNIMOD:1882
    true,   // UNIMOD:1883
    false,  // UNIMOD:1884: not present in source Unimod snapshot
    true,   // UNIMOD:1885
    true,   // UNIMOD:1886
    true,   // UNIMOD:1887
    true,   // UNIMOD:1888
    true,   // UNIMOD:1889
    false,  // UNIMOD:1890: not present in source Unimod snapshot
    false,  // UNIMOD:1891: not present in source Unimod snapshot
    false,  // UNIMOD:1892: not present in source Unimod snapshot
    false,  // UNIMOD:1893: not present in source Unimod snapshot
    false,  // UNIMOD:1894: not present in source Unimod snapshot
    false,  // UNIMOD:1895: not present in source Unimod snapshot
    true,   // UNIMOD:1896
    true,   // UNIMOD:1897
    true,   // UNIMOD:1898
    true,   // UNIMOD:1899
    true,   // UNIMOD:1900
    true,   // UNIMOD:1901
    true,   // UNIMOD:1902
    true,   // UNIMOD:1903
    false,  // UNIMOD:1904: not present in source Unimod snapshot
    true,   // UNIMOD:1905
    true,   // UNIMOD:1906
    true,   // UNIMOD:1907
    true,   // UNIMOD:1908
    false,  // UNIMOD:1909: not present in source Unimod snapshot
    true,   // UNIMOD:1910
    true,   // UNIMOD:1911
    true,   // UNIMOD:1912
    true,   // UNIMOD:1913
    true,   // UNIMOD:1914
    true,   // UNIMOD:1915
    true,   // UNIMOD:1916
    true,   // UNIMOD:1917
    true,   // UNIMOD:1918
    false,  // UNIMOD:1919: not present in source Unimod snapshot
    true,   // UNIMOD:1920
    false,  // UNIMOD:1921: not present in source Unimod snapshot
    true,   // UNIMOD:1922
    true,   // UNIMOD:1923
    true,   // UNIMOD:1924
    true,   // UNIMOD:1925
    true,   // UNIMOD:1926
    true,   // UNIMOD:1927
    true,   // UNIMOD:1928
    true,   // UNIMOD:1929
    false,  // UNIMOD:1930: non-element token Pent(2); group-to-formula expansion not supported
    false,  // UNIMOD:1931: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:1932: non-element token Hex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1933: non-element token Pent(2); group-to-formula expansion not supported
    false,  // UNIMOD:1934: non-element token HexNAc(2); group-to-formula expansion not supported
    false,  // UNIMOD:1935: non-element token Pent(3); group-to-formula expansion not supported
    false,  // UNIMOD:1936: non-element token Pent(2); group-to-formula expansion not supported
    false,  // UNIMOD:1937: non-element token Pent(2); group-to-formula expansion not supported
    false,  // UNIMOD:1938: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1939: non-element token Pent; group-to-formula expansion not supported
    false,  // UNIMOD:1940: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1941: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1942: non-element token HexA(2); group-to-formula expansion not supported
    false,  // UNIMOD:1943: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1944: non-element token Hex(5); group-to-formula expansion not supported
    false,  // UNIMOD:1945: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1946: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1947: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:1948: non-element token Sulf; group-to-formula expansion not supported
    false,  // UNIMOD:1949: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1950: non-element token dHex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1951: non-element token Sulf; group-to-formula expansion not supported
    false,  // UNIMOD:1952: non-element token Sulf(2); group-to-formula expansion not supported
    false,  // UNIMOD:1953: non-element token Hex(9); group-to-formula expansion not supported
    false,  // UNIMOD:1954: non-element token Sulf; group-to-formula expansion not supported
    false,  // UNIMOD:1955: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1956: non-element token Sulf(2); group-to-formula expansion not supported
    false,  // UNIMOD:1957: non-element token Hex(9); group-to-formula expansion not supported
    false,  // UNIMOD:1958: non-element token Sulf(2); group-to-formula expansion not supported
    false,  // UNIMOD:1959: non-element token Hex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1960: non-element token dHex(4); group-to-formula expansion not supported
    false,  // UNIMOD:1961: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1962: non-element token Hex(10); group-to-formula expansion not supported
    false,  // UNIMOD:1963: non-element token dHex; group-to-formula expansion not supported
    false,  // UNIMOD:1964: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1965: non-element token dHex(2); group-to-formula expansion not supported
    false,  // UNIMOD:1966: non-element token Sulf; group-to-formula expansion not supported
    false,  // UNIMOD:1967: token Ac has no supported element symbol or group expansion
    false,  // UNIMOD:1968: non-element token Hex(3); group-to-formula expansion not supported
    false,  // UNIMOD:1969: token Ac(2) has no supported element symbol or group expansion
    true,   // UNIMOD:1970
    true,   // UNIMOD:1971
    true,   // UNIMOD:1972
    true,   // UNIMOD:1973
    true,   // UNIMOD:1974
    true,   // UNIMOD:1975
    true,   // UNIMOD:1976
    true,   // UNIMOD:1977
    true,   // UNIMOD:1978
    true,   // UNIMOD:1979
    false,  // UNIMOD:1980: not present in source Unimod snapshot
    false,  // UNIMOD:1981: not present in source Unimod snapshot
    false,  // UNIMOD:1982: not present in source Unimod snapshot
    false,  // UNIMOD:1983: not present in source Unimod snapshot
    false,  // UNIMOD:1984: not present in source Unimod snapshot
    false,  // UNIMOD:1985: not present in source Unimod snapshot
    true,   // UNIMOD:1986
    true,   // UNIMOD:1987
    false,  // UNIMOD:1988: not present in source Unimod snapshot
    true,   // UNIMOD:1989
    true,   // UNIMOD:1990
    true,   // UNIMOD:1991
    true,   // UNIMOD:1992
    false,  // UNIMOD:1993: isotope-labeled token 13C(9); pinned-isotope conversion not supported
    false,  // UNIMOD:1994: not present in source Unimod snapshot
    false,  // UNIMOD:1995: not present in source Unimod snapshot
    false,  // UNIMOD:1996: not present in source Unimod snapshot
    false,  // UNIMOD:1997: not present in source Unimod snapshot
    false,  // UNIMOD:1998: not present in source Unimod snapshot
    true,   // UNIMOD:1999
    false,  // UNIMOD:2000: not present in source Unimod snapshot
    true,   // UNIMOD:2001
    false,  // UNIMOD:2002: not present in source Unimod snapshot
    false,  // UNIMOD:2003: not present in source Unimod snapshot
    false,  // UNIMOD:2004: not present in source Unimod snapshot
    false,  // UNIMOD:2005: not present in source Unimod snapshot
    true,   // UNIMOD:2006
    true,   // UNIMOD:2007
    true,   // UNIMOD:2008
    false,  // UNIMOD:2009: not present in source Unimod snapshot
    false,  // UNIMOD:2010: not present in source Unimod snapshot
    false,  // UNIMOD:2011: not present in source Unimod snapshot
    false,  // UNIMOD:2012: not present in source Unimod snapshot
    false,  // UNIMOD:2013: not present in source Unimod snapshot
    true,   // UNIMOD:2014
    false,  // UNIMOD:2015: isotope-labeled token 13C(9); pinned-isotope conversion not supported
    false,  // UNIMOD:2016: isotope-labeled token 13C(7); pinned-isotope conversion not supported
    true,   // UNIMOD:2017
    false,  // UNIMOD:2018: not present in source Unimod snapshot
    false,  // UNIMOD:2019: not present in source Unimod snapshot
    false,  // UNIMOD:2020: not present in source Unimod snapshot
    false,  // UNIMOD:2021: not present in source Unimod snapshot
    true,   // UNIMOD:2022
    false,  // UNIMOD:2023: not present in source Unimod snapshot
    false,  // UNIMOD:2024: not present in source Unimod snapshot
    true,   // UNIMOD:2025
    false,  // UNIMOD:2026: not present in source Unimod snapshot
    true,   // UNIMOD:2027
    false,  // UNIMOD:2028: non-element token Hex(6); group-to-formula expansion not supported
    false,  // UNIMOD:2029: non-element token Hex(7); group-to-formula expansion not supported
    false,  // UNIMOD:2030: not present in source Unimod snapshot
    false,  // UNIMOD:2031: not present in source Unimod snapshot
    false,  // UNIMOD:2032: not present in source Unimod snapshot
    true,   // UNIMOD:2033
    true,   // UNIMOD:2034
    true,   // UNIMOD:2035
    true,   // UNIMOD:2036
    true,   // UNIMOD:2037
    false,  // UNIMOD:2038: not present in source Unimod snapshot
    true,   // UNIMOD:2039
    true,   // UNIMOD:2040
    true,   // UNIMOD:2041
    true,   // UNIMOD:2042
    false,  // UNIMOD:2043: not present in source Unimod snapshot
    true,   // UNIMOD:2044
    false,  // UNIMOD:2045: not present in source Unimod snapshot
    false,  // UNIMOD:2046: not present in source Unimod snapshot
    false,  // UNIMOD:2047: not present in source Unimod snapshot
    false,  // UNIMOD:2048: not present in source Unimod snapshot
    false,  // UNIMOD:2049: not present in source Unimod snapshot
    false,  // UNIMOD:2050: isotope-labeled token 13C(15); pinned-isotope conversion not supported
    false,  // UNIMOD:2051: not present in source Unimod snapshot
    true,   // UNIMOD:2052
    true,   // UNIMOD:2053
    true,   // UNIMOD:2054
    true,   // UNIMOD:2055
    false,  // UNIMOD:2056: not present in source Unimod snapshot
    true,   // UNIMOD:2057
    true,   // UNIMOD:2058
    true,   // UNIMOD:2059
    true,   // UNIMOD:2060
    true,   // UNIMOD:2061
    true,   // UNIMOD:2062
    false,  // UNIMOD:2063: not present in source Unimod snapshot
    false,  // UNIMOD:2064: not present in source Unimod snapshot
    false,  // UNIMOD:2065: not present in source Unimod snapshot
    false,  // UNIMOD:2066: not present in source Unimod snapshot
    true,   // UNIMOD:2067
    true,   // UNIMOD:2068
    true,   // UNIMOD:2069
    true,   // UNIMOD:2070
    false,  // UNIMOD:2071: not present in source Unimod snapshot
    true,   // UNIMOD:2072
    true,   // UNIMOD:2073
    true,   // UNIMOD:2074
    false,  // UNIMOD:2075: not present in source Unimod snapshot
    false,  // UNIMOD:2076: not present in source Unimod snapshot
    false,  // UNIMOD:2077: not present in source Unimod snapshot
    false,  // UNIMOD:2078: not present in source Unimod snapshot
    true,   // UNIMOD:2079
    true,   // UNIMOD:2080
    true,   // UNIMOD:2081
    true,   // UNIMOD:2082
    true,   // UNIMOD:2083
    true,   // UNIMOD:2084
    true,   // UNIMOD:2085
    true,   // UNIMOD:2086
    false,  // UNIMOD:2087: not present in source Unimod snapshot
    false,  // UNIMOD:2088: isotope-labeled token 13C(2); pinned-isotope conversion not supported
    false,  // UNIMOD:2089: not present in source Unimod snapshot
    false,  // UNIMOD:2090: not present in source Unimod snapshot
    false,  // UNIMOD:2091: not present in source Unimod snapshot
    false,  // UNIMOD:2092: not present in source Unimod snapshot
    false,  // UNIMOD:2093: not present in source Unimod snapshot
    false,  // UNIMOD:2094: not present in source Unimod snapshot
    false,  // UNIMOD:2095: not present in source Unimod snapshot
    false,  // UNIMOD:2096: not present in source Unimod snapshot
    false,  // UNIMOD:2097: not present in source Unimod snapshot
    false,  // UNIMOD:2098: not present in source Unimod snapshot
    false,  // UNIMOD:2099: not present in source Unimod snapshot
    false,  // UNIMOD:2100: not present in source Unimod snapshot
    false,  // UNIMOD:2101: not present in source Unimod snapshot
    false,  // UNIMOD:2102: not present in source Unimod snapshot
    false,  // UNIMOD:2103: not present in source Unimod snapshot
    false,  // UNIMOD:2104: not present in source Unimod snapshot
    false,  // UNIMOD:2105: not present in source Unimod snapshot
    true,   // UNIMOD:2106
    true,   // UNIMOD:2107
    true,   // UNIMOD:2108
    true,   // UNIMOD:2109
    true,   // UNIMOD:2110
    true,   // UNIMOD:2111
    true,   // UNIMOD:2112
    true,   // UNIMOD:2113
    true,   // UNIMOD:2114
    true,   // UNIMOD:2115
    true,   // UNIMOD:2116
    true,   // UNIMOD:2117
    true,   // UNIMOD:2118
    true,   // UNIMOD:2119
    true,   // UNIMOD:2120
    true,   // UNIMOD:2121
    false,  // UNIMOD:2122: isotope-labeled token 13C(10); pinned-isotope conversion not supported
    false,  // UNIMOD:2123: isotope-labeled token 13C(13); pinned-isotope conversion not supported
    false,  // UNIMOD:2124: not present in source Unimod snapshot
    false,  // UNIMOD:2125: not present in source Unimod snapshot
    true,   // UNIMOD:2126
    true,   // UNIMOD:2127
    true,   // UNIMOD:2128
    true,   // UNIMOD:2129
    true,   // UNIMOD:2130
    true,   // UNIMOD:2131
    true,   // UNIMOD:2132
    false,  // UNIMOD:2133: not present in source Unimod snapshot
    false,  // UNIMOD:2134: not present in source Unimod snapshot
    true,   // UNIMOD:2135
    true,   // UNIMOD:2136
    false,  // UNIMOD:2137: not present in source Unimod snapshot
    true,   // UNIMOD:2138
    true,   // UNIMOD:2139
    false,  // UNIMOD:2140: isotope-labeled token 2H(2); pinned-isotope conversion not supported
    true,   // UNIMOD:2141
    true,   // UNIMOD:2142
    true,   // UNIMOD:2143
    true,   // UNIMOD:2144
    false,  // UNIMOD:2145: not present in source Unimod snapshot
    false,  // UNIMOD:2146: isotope-labeled token 2H(8); pinned-isotope conversion not supported
    true,   // UNIMOD:2147
};

//! Tests support in the shipped table without loading or parsing it.
//! Unsupported and unknown IDs return false. Override CSVs do not affect this
//! ledger. A supported modification still needs a compatible base composition.
inline constexpr bool is_unimod_supported(std::uint64_t id) noexcept {
    return id < sizeof(unimod_supported) / sizeof(unimod_supported[0])
        && unimod_supported[static_cast<std::size_t>(id)];
}

}  // namespace IsoSpec
