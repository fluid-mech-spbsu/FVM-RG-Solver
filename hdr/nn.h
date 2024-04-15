
/// ------------------------------------------------------------------------------ ///
/// ------------------------ NN for energy calculation --------------------------- ///
/// ------------------------------ (p,T) -> rho*U/1000  -------------------------- ///

double T_MIN = 250;
double T_MAX = 1999;

double P_MIN = 5;
double P_MAX = 499;

double enlayers0weights[2][50] = {
    {9.7095096111297607e-01,  1.6720622777938843e+00,
        -6.8972796201705933e-01,  2.9848214983940125e-01,
         2.1246190071105957e+00, -3.3672651648521423e-01,
        -8.2278043031692505e-01,  3.6726310849189758e-01,
         4.3615421652793884e-01,  1.0418766736984253e+00,
         8.9897692203521729e-01,  4.6442168951034546e-01,
        -1.4565324783325195e+00,  3.5252955555915833e-01,
         8.0000913143157959e-01, -1.4310457706451416e+00,
         6.0025829076766968e-01,  3.4324774146080017e-01,
         3.4862998127937317e-01, -1.4499273300170898e+00,
         4.3922644853591919e-01,  3.3114135265350342e-01,
         3.8790634274482727e-01, -5.1621079444885254e-01,
        -8.8808405399322510e-01,  2.6918002963066101e-01,
         2.1852126717567444e-01,  4.9526995420455933e-01,
         9.1483926773071289e-01, -8.2394492626190186e-01,
         4.4216388463973999e-01,  1.1740294694900513e+00,
         1.3746597766876221e+00, -3.6556455492973328e-01,
         8.8404172658920288e-01, -5.9465909004211426e-01,
         3.0086815357208252e-01, -1.0434908866882324e+00,
         7.2619986534118652e-01, -2.4692287445068359e+00,
         8.0421787500381470e-01, -4.2879715561866760e-01,
         8.7324976921081543e-01, -2.8993734717369080e-01,
         5.7638508081436157e-01, -1.3071371316909790e+00,
        -1.0249700546264648e+00, -1.1512367725372314e+00,
        -1.2222436666488647e+00, -3.5645744204521179e-01 },
    {-3.7314903736114502e-01,  9.1708536148071289e+00,
         6.9856479763984680e-02,  8.2866638898849487e-02,
         1.5199012756347656e+01, -5.7550504803657532e-02,
        -2.8307456970214844e+00,  5.9811372309923172e-02,
         4.8092972487211227e-02, -3.9142253398895264e+00,
        -4.2727322578430176e+00,  4.0087964385747910e-02,
         2.2354290008544922e+01,  8.2503646612167358e-02,
        -1.4798313379287720e-01, -5.7451395988464355e+00,
         2.2916032001376152e-02,  8.3738580346107483e-02,
         7.0102773606777191e-02, -2.5503039360046387e+00,
        -2.4647932052612305e+00,  6.7403107881546021e-02,
         4.3153040111064911e-02, -2.8115211054682732e-02,
         6.6467437744140625e+00,  9.3327268958091736e-02,
         1.2385966628789902e-01, -3.4557660110294819e-03,
        -7.0055012702941895e+00,  3.8528337478637695e+00,
         3.0030494555830956e-02, -5.6207376718521118e-01,
         5.1636590957641602e+00, -3.9864189922809601e-02,
         3.5319766998291016e+00,  1.6502998769283295e-03,
         1.0263830423355103e-01,  9.0239953994750977e+00,
        -2.9896779060363770e+00, -1.9030084609985352e+01,
         2.9237918853759766e+00, -1.8387908115983009e-02,
         3.3974370956420898e+00, -8.2789547741413116e-02,
         1.2468426488339901e-02,  1.7438796758651733e+00,
        -4.4902467727661133e+00,  1.3721853256225586e+01,
         1.1621583938598633e+01, -6.9442160427570343e-02 } };

double enlayers0bias[1][50] = { 
        -0.45583322644233703613,  0.12476086616516113281,
        -0.23756675422191619873,  1.10412633419036865234,
         0.49291768670082092285, -1.13258635997772216797,
        -1.15310716629028320312,  0.75866520404815673828,
         0.44312319159507751465, -0.45851644873619079590,
        -0.61183345317840576172,  0.39597401022911071777,
         1.97273850440979003906,  0.65266621112823486328,
        -0.53593230247497558594,  0.08894436061382293701,
         1.03463256359100341797,  0.66691994667053222656,
         0.81949180364608764648,  0.41844636201858520508,
         0.19167479872703552246,  1.12334847450256347656,
         0.95032387971878051758, -0.21746501326560974121,
         0.53808820247650146484,  1.30103969573974609375,
         1.35725605487823486328,  0.69470465183258056641,
        -0.74671953916549682617,  0.53752815723419189453,
         0.92633962631225585938, -0.84959924221038818359,
        -0.08246970921754837036, -1.16538882255554199219,
         1.10661745071411132812, -0.00449737021699547768,
         0.77182358503341674805,  1.07115805149078369141,
        -0.25089383125305175781, -0.44669628143310546875,
         1.49848914146423339844, -0.89548516273498535156,
         1.07344985008239746094, -1.22289693355560302734,
         1.12342965602874755859,  0.64610612392425537109,
        -0.82020413875579833984,  1.51448237895965576172,
         1.46232032775878906250, -0.71622782945632934570
};

double enlayers1weights[50][1] = {
        1.08522474765777587891,  1.53735280036926269531,
         -1.04821848869323730469,  0.80898332595825195312,
          0.93198990821838378906, -0.95304083824157714844,
         -1.02355802059173583984,  0.52799016237258911133,
          0.73226249217987060547,  1.02936494350433349609,
          0.98963052034378051758,  0.86609619855880737305,
         -3.26262354850769042969,  1.20017242431640625000,
          1.45039856433868408203, -1.17457783222198486328,
          0.82496637105941772461,  1.16071116924285888672,
          0.90450662374496459961, -1.30912804603576660156,
          1.01487112045288085938,  0.59253454208374023438,
          0.74286371469497680664, -1.32134878635406494141,
         -2.05862569808959960938,  0.84414905309677124023,
          1.15009021759033203125,  0.65777164697647094727,
          2.20019364356994628906, -1.66092157363891601562,
          1.22171044349670410156,  0.98953890800476074219,
          0.83221799135208129883, -0.97536331415176391602,
          0.73297524452209472656, -1.07044494152069091797,
          0.99677985906600952148, -2.74497032165527343750,
          0.99707865715026855469, -0.96116721630096435547,
          0.79099279642105102539, -1.06953179836273193359,
          0.82346808910369873047, -0.76188915967941284180,
          0.77346211671829223633, -0.83839589357376098633,
         -0.69429451227188110352, -3.13439917564392089844,
         -2.63587832450866699219, -0.53120851516723632812
};

double enlayers1bias = 0.45196458697319030762;

/// ------------------------------------------------------------------------------ ///
/// ------------------------ NN for zeta calculation ----------------------------- ///
/// --------------------------- (p,T) -> -log10(zeta)  --------------------------- ///


double zT_MIN = 250; // i generated a little bit less data than for Evibr(p,T)
double zT_MAX = 1998.;

double zP_MIN = 5;
double zP_MAX = 499.;

double zetalayers0weights[2][50] = { {-9.37648863e-02, -5.99425256e-01, -3.07872057e-01,
        -1.26079440e-01, -2.52141523e+00,  1.60194731e+00,
        -2.57682763e-02,  1.94326624e-01,  6.76833963e+00,
         1.54433742e-01, -5.11013642e-02,  7.17732608e-01,
        -2.36252844e-01,  4.17091519e-01,  2.54197866e-01,
        -3.00195843e-01, -4.38170016e-01,  7.93419220e-03,
        -4.47483182e-01,  1.01837003e+00,  1.18852966e-01,
        -2.26717025e-01,  1.61883879e+00, -5.86007655e-01,
         4.64691162e-01,  1.40785053e-01,  9.21683073e-01,
        -1.80505514e-01,  8.54420960e-02,  1.51350963e+00,
         1.75410479e-01,  3.99253339e-01,  4.02582139e-02,
        -1.63309649e-01,  7.54685178e-02, -4.33143139e-01,
         2.46845797e-01,  8.50038975e-03,  1.75314248e-01,
        -3.54863793e-01,  4.36040372e-01,  2.98143148e-01,
         9.18047130e-02, -2.39618763e-01,  1.17350519e-01,
         1.81364131e+00,  9.26165506e-02, -7.22320557e-01,
        -9.58946869e-02,  8.07082653e-02 },
         {2.24077143e-03,  1.26313889e+00, -2.40907833e-01,
         3.34416270e-01, -2.26521587e+00,  1.69530046e+00,
        -1.45371750e-01, -5.32392412e-02, -2.44282886e-01,
         3.55031580e-01, -2.42964864e-01,  1.04507160e+00,
        -9.84826088e-02,  6.41863823e-01, -3.57391536e-01,
        -2.11804986e-01, -7.50287771e-01, -1.87618867e-01,
         1.94646925e-01, -7.21312940e-01,  2.82259703e-01,
        -1.11425065e-01,  1.77433610e+00, -8.48071039e-01,
         7.91788757e-01,  3.41486692e-01,  1.20403731e+00,
        -1.46590129e-01,  2.63724536e-01,  8.13482821e-01,
        -5.73252773e+00,  6.44861579e-01, -8.86700004e-02,
        -2.65189648e-01,  3.60388249e-01, -6.67734206e-01,
         6.87483326e-02,  3.79927270e-02, -3.39547873e-01,
        -5.21461785e-01,  5.79355657e-01,  4.30076152e-01,
         2.62117516e-02, -3.62450302e-01,  6.36757091e-02,
         1.87549698e+00,  1.53540969e-02, -1.09925210e+00,
         3.57846141e-01, -1.16268769e-01 } };

double zetalayers0bias[1][50] = {{
        -1.02756094932556152344, -0.29870748519897460938,
        -0.66525459289550781250, -0.06200504675507545471,
        -0.38191878795623779297,  0.70538127422332763672,
         0.67578822374343872070,  0.25331139564514160156,
         0.56663203239440917969, -0.85684460401535034180,
         0.67507296800613403320,  0.77759552001953125000,
        -0.40008288621902465820,  0.58673542737960815430,
         0.26268830895423889160, -0.12336005270481109619,
        -0.39630213379859924316,  0.64052534103393554688,
         0.26071909070014953613,  0.18840198218822479248,
        -0.62157851457595825195, -0.38586997985839843750,
         0.43116647005081176758, -0.04028697684407234192,
         0.40372729301452636719, -0.92187410593032836914,
         0.84807229042053222656, -0.62651491165161132812,
        -0.14570054411888122559, -0.53557294607162475586,
        -0.23154276609420776367,  0.78466922044754028320,
         0.34254661202430725098,  0.68648207187652587891,
        -0.06456223130226135254, -0.70691263675689697266,
        -0.13643346726894378662,  0.39762762188911437988,
         0.22882552444934844971, -0.69820362329483032227,
         0.43460237979888916016,  0.94540214538574218750,
         0.55929177999496459961, -0.53665041923522949219,
         0.97807353734970092773,  0.48854133486747741699,
         0.51745438575744628906, -0.38095211982727050781,
        -0.44661819934844970703,  0.93205767869949340820 }
};

double zetalayers1weights[50][1] = {
         -0.19960658252239227295,  0.56992202997207641602,
         -0.26582461595535278320, -0.09982536733150482178,
         -0.44911971688270568848,  0.25213906168937683105,
          0.41086167097091674805,  0.28653848171234130859,
         -1.60343170166015625000, -0.48586964607238769531,
          0.48740464448928833008,  0.20114022493362426758,
         -0.32826337218284606934,  0.30328190326690673828,
          0.03835103288292884827, -0.09140224009752273560,
         -0.23669573664665222168,  0.19331817328929901123,
          0.09430427849292755127, -0.39834225177764892578,
         -0.21526955068111419678, -0.39068588614463806152,
          0.30379083752632141113, -0.13795542716979980469,
          0.24918039143085479736, -0.40506091713905334473,
          0.31451532244682312012, -0.30634111166000366211,
         -0.05891290307044982910, -0.72682183980941772461,
          1.39453220367431640625,  0.34692019224166870117,
          0.10433116555213928223,  0.24707099795341491699,
         -0.02319119125604629517, -0.28775000572204589844,
         -0.00511640030890703201,  0.41391021013259887695,
          0.10358305275440216064, -0.33957934379577636719,
          0.30348896980285644531,  0.29789525270462036133,
          0.36810266971588134766, -0.09731773287057876587,
          0.19714540243148803711,  0.33938893675804138184,
          0.42877137660980224609, -0.29537957906723022461,
         -0.09849094599485397339,  0.33548015356063842773 };

double zetalayers1bias = 0.25275662541389465332;


/// ------------------------------------------------------------------------------ ///
/// ------------------------ NN for temp calculation ----------------------------- ///
/// ---------------------------- rho*U/1000 -> (p,T)  --------------------------- ///


double tlayers0weights[1][50] = { {
    -6.04180992e-03, -3.83959571e-03, -2.88097382e+00,
        -4.56410553e-03,  5.03548145e+00, -4.51865292e+00,
        -1.20262611e+00,  7.00696290e-01, -1.82791686e+00,
         2.81509042e-01, -3.08088398e+00,  1.84900433e-01,
         2.62791419e+00, -5.11569977e-01, -8.06707685e-05,
         9.19996892e-05, -4.79483277e-01, -4.80032742e-01,
         4.65027004e-01,  5.15201151e-01,  2.56633669e-01,
        -9.33535218e-01,  1.30222067e-02, -9.89991203e-02,
         1.63434923e-03,  1.16455078e-03,  2.96085440e-02,
         1.25965336e-02,  1.09422743e+00, -7.46648312e-01,
        -3.24858265e-04, -1.13013339e+00, -1.72288402e-03,
        -7.16837402e-03,  6.77977741e-01,  8.02758157e-01,
        -6.69447836e-05, -2.25350974e-04,  2.53848100e+00,
        -4.54610169e-01, -8.81023152e-05,  3.45065966e-02,
         2.40758729e+00,  6.26266585e-04, -1.30983487e-01,
         3.96048635e-01, -6.05776119e+00,  1.30912349e-01,
         7.22190261e-01, -7.33687341e-01
         } };

double tlayers0biases[1][50] = { 1.17618043441325426102e-03,  7.15832517016679048538e-04,
         8.93882140517234802246e-02,  6.40956102870404720306e-04,
        -6.17041707038879394531e-01,  9.69992995262145996094e-01,
         1.37560456991195678711e-01, -4.49563384056091308594e-01,
         6.02368488907814025879e-02, -4.94886264204978942871e-02,
         2.42396071553230285645e-01, -6.83195114135742187500e-01,
        -1.29935383796691894531e-01,  8.33953246474266052246e-02,
         1.49938696267781779170e-05, -1.48962963066878728569e-05,
         5.99043607711791992188e-01,  1.39259114861488342285e-01,
        -7.67664015293121337891e-02,  5.59995651245117187500e-01,
         5.03825664520263671875e-01,  2.91708499193191528320e-01,
        -2.44802562519907951355e-03,  5.92090606689453125000e-01,
        -1.40680553158745169640e-04, -2.78371298918500542641e-04,
        -5.49948215484619140625e-03, -7.67266228795051574707e-02,
        -2.54465192556381225586e-01,  1.18398867547512054443e-01,
         4.38459683209657669067e-05,  1.89647465944290161133e-01,
         4.80702219647355377674e-05,  1.39569968450814485550e-03,
        -7.00103044509887695312e-02, -1.28728747367858886719e-01,
         1.57758204295532777905e-05,  4.46273843408562242985e-05,
        -3.30651521682739257812e-01, -2.69706677645444869995e-02,
         1.41977552630123682320e-05,  1.18493676185607910156e-01,
        -8.18421661853790283203e-01,  1.10188469989225268364e-04,
         2.21550725400447845459e-02, -6.26316443085670471191e-02,
         7.54881262779235839844e-01,  4.41790580749511718750e-01,
         2.06697478890419006348e-01,  1.13581568002700805664e-01 };

double tlayers1weights[2][50] = { { -1.99904534383676946163e-04, -1.68538870639167726040e-04,
1.04049235582351684570e-01, -7.64937649364583194256e-05,
6.33394896984100341797e-01, 7.53042161464691162109e-01,
-2.27088063955307006836e-01, 6.45481348037719726562e-02,
-4.76022586226463317871e-02, 7.73298647254705429077e-03,
-1.59334912896156311035e-02, 5.25896064937114715576e-03,
-3.00182551145553588867e-01, -7.64501392841339111328e-02,
8.96464007382746785879e-06, 2.40436056628823280334e-05,
-3.02500128746032714844e-02, 2.28662285953760147095e-02,
2.73178424686193466187e-02, 3.64625044167041778564e-02,
-6.10213540494441986084e-02, -9.05532985925674438477e-02,
7.85947544500231742859e-04, 8.95126089453697204590e-02,
-2.49171935138292610645e-05, 1.03213696274906396866e-03,
7.18535564374178647995e-04, -2.65876892954111099243e-02,
1.86820328235626220703e-01, -4.94806505739688873291e-02,
4.75764863949734717607e-07, -5.81776276230812072754e-02,
6.97998402756638824940e-05, -2.90170079097151756287e-04,
1.26035079360008239746e-01, 1.40090033411979675293e-01,
1.10778364614816382527e-05, -1.03379279607906937599e-05,
-2.95714139938354492188e-01, -7.56188407540321350098e-02,
-2.34188628382980823517e-05, 1.05103895068168640137e-01,
9.64588746428489685059e-02, -4.06716062570922076702e-05,
-6.04457920417189598083e-03, 9.66391488909721374512e-02,
-6.26537442207336425781e-01, 1.59740269184112548828e-01,
9.27792116999626159668e-02, -4.04924042522907257080e-02},
{2.04571770154871046543e-05, 2.53334601438837125897e-05,
-5.78432679176330566406e-02, -4.47811216872651129961e-05,
3.48931849002838134766e-01, 1.15461480617523193359e+00,
1.90880000591278076172e-02, -1.19449444115161895752e-01,
1.13868586719036102295e-01, -4.02331911027431488037e-03,
2.64320462942123413086e-01, -1.61735698580741882324e-01,
-4.02472048997879028320e-01, 3.79129014909267425537e-02,
-1.55487978190649300814e-05, -1.83429147000424563885e-05,
-3.94304618239402770996e-02, -1.52574712410569190979e-02,
-1.39049747958779335022e-02, 1.50479897856712341309e-01,
1.42828568816184997559e-01, 5.52737303078174591064e-02,
-1.41859549330547451973e-04, 8.57863649725914001465e-02,
-2.00979070541507098824e-06, -8.45971808303147554398e-04,
-8.19233246147632598877e-05, -3.82493101060390472412e-02,
1.31935805082321166992e-01, 1.30620347335934638977e-02,
-9.79886772256577387452e-06, -3.60114723443984985352e-02,
5.33384154550731182098e-05, -9.47349617490544915199e-05,
1.82953953742980957031e-01, 4.25821319222450256348e-02,
-1.48598128362209536135e-05, -6.10436700299032963812e-06,
1.00489497184753417969e-01, -2.64202011749148368835e-03,
1.79956077772658318281e-05, -2.68145166337490081787e-02,
3.86875152587890625000e-01, -1.02006095403339713812e-04,
2.76625109836459159851e-03, -4.81164306402206420898e-02,
-5.22036433219909667969e-01, 1.36001691222190856934e-01,
1.27045422792434692383e-01, 3.27382073737680912018e-03} }; 

double tlayers1biases[2][2] = { { 0.18312291800975799561, 0.04472917690873146057} };