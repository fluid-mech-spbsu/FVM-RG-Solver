#pragma once
/*
/// ------------------------------------------------------------------------------ ///
/// ------------------------ NN for energy calculation --------------------------- ///
/// ------------------------------ (p,T) -> rho*U/1000  -------------------------- ///

double T_MIN = 250.;
double T_MAX = 1999.;

double P_MIN = 5.;
double P_MAX = 499.;

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


double zT_MIN = 250.; // i generated a little bit less data than for Evibr(p,T)
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

double tlayers0bias[1][50] = { 
         1.17618043441325426102e-03,  7.15832517016679048538e-04,
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

double tlayers1biases[1][2] = { 0.18312291800975799561, 0.04472917690873146057 };


/// ------------------------------------------------------------------------------ ///
/// ------------------------ NN for lambda calculation ----------------------------- ///
/// --------------------------- (p,T) -> -log10(lambda)  --------------------------- ///


double llayers0weights[2][50] = {
        { 1.89327136e-01,  1.52862206e-01,  4.36335891e-01,
        -6.34332776e-01,  1.85923954e-03,  4.79507565e-01,
         4.67219830e-01, -1.68965086e-01, -1.53981894e-01,
        -1.92868352e-01,  2.76274711e-01,  3.08225185e-01,
        -1.67778939e-01,  7.90728509e-01, -1.98719308e-01,
         1.03945740e-01,  7.70158246e-02, -4.86927092e-01,
         1.13151439e-01, -4.06081051e-01,  1.23465762e-01,
        -3.62581849e-01, -6.33604527e-02, -2.43579984e-01,
         1.13958925e-01,  3.28026801e-01, -8.70528594e-02,
        -1.84360772e-01,  3.92506570e-01,  2.30460972e-01,
        -1.33766476e-02, -4.68701236e-02,  3.98116605e-03,
        -1.05997038e+00,  5.61046600e-01, -2.47636676e-01,
        -3.74905050e-01,  1.08730830e-01, -2.73888439e-01,
        -1.69356406e-01,  1.50909498e-01, -1.45420760e-01,
         2.60670036e-01, -1.33838177e-01, -1.75040334e-01,
        -2.46181805e-02, -2.99846064e-02,  1.66471884e-01,
         2.00141802e-01, -6.27913028e-02},
        {-2.84964796e-02,  2.65689909e-01,  9.12979618e-02,
         2.62193501e-01,  4.08219546e-03, -1.30411297e-01,
         5.06971657e-01, -1.23100114e+00, -1.56175882e-01,
        -2.00696170e-01,  2.36448944e-01,  2.77360827e-01,
         1.18458439e-02,  1.28630614e+00, -2.17459336e-01,
         8.92909825e-01,  2.24040031e+00, -4.46039379e-01,
         1.52172633e-02,  4.06882346e-01, -2.37857282e-01,
        -3.52991879e-01, -9.23736468e-02,  1.64510465e+00,
        -6.73002720e-01,  9.78778377e-02, -1.20621867e-01,
         2.26525068e-01,  4.00546253e-01,  5.29016435e-01,
        -3.09177991e-02, -1.22125261e-02,  3.67185974e+00,
        -3.68458599e-01, -1.91756308e-01, -1.10173273e+00,
        -1.26717389e-01,  4.60287444e-02, -3.51372510e-01,
        -2.17800036e-01,  3.13572198e-01,  1.93379056e-02,
        -2.46192083e-01,  5.31159528e-02, -1.97064072e-01,
        -9.66664627e-02,  5.09020090e-02,  3.61452609e-01,
        -4.05622721e-02, -4.01200473e-01} };

double llayers0bias[1][50] = { {
        0.39302831888198852539,  0.76095491647720336914,
        -0.51882410049438476562,  0.08305519819259643555,
         0.57761657238006591797, -0.01894748210906982422,
         0.82360404729843139648,  0.26801350712776184082,
         0.74895876646041870117, -0.33436629176139831543,
        -0.76039564609527587891,  0.74755126237869262695,
        -0.25586193799972534180, -0.11179817467927932739,
        -0.73620998859405517578, -0.32138773798942565918,
         0.25895592570304870605, -0.38464513421058654785,
         0.41815757751464843750,  0.21434396505355834961,
         0.37207141518592834473, -0.16551820933818817139,
        -0.99586176872253417969,  0.17169044911861419678,
         0.18308568000793457031,  0.74456328153610229492,
        -0.17514449357986450195, -0.19894334673881530762,
         0.76498490571975708008,  0.26126915216445922852,
        -0.73152446746826171875, -0.78945791721343994141,
         0.47093388438224792480,  0.01511394605040550232,
        -0.33045816421508789062,  0.07286833226680755615,
        -0.15001051127910614014, -0.85984987020492553711,
        -0.62532532215118408203, -0.42258939146995544434,
         0.35638940334320068359,  0.47395417094230651855,
        -0.00655508134514093399,  0.74355906248092651367,
        -0.84927439689636230469,  0.36830711364746093750,
         0.71281415224075317383, -0.54927104711532592773,
         0.60226237773895263672,  0.65142685174942016602} };

double llayers2weights[1][50] = { {
        0.26339396834373474121,  0.26533851027488708496,
         -0.07538151741027832031,  0.03726078569889068604,
          0.26834005117416381836, -0.04223218932747840881,
          0.18197999894618988037,  0.27285045385360717773,
          0.49703106284141540527, -0.30349639058113098145,
         -0.25827738642692565918,  0.20693373680114746094,
         -0.15028336644172668457, -0.29910263419151306152,
         -0.20111373066902160645, -0.22561708092689514160,
         -0.49750542640686035156, -0.20614430308341979980,
          0.38642793893814086914, -0.04971105232834815979,
          0.28190740942955017090, -0.02088427916169166565,
         -0.15926386415958404541, -0.57544577121734619141,
          0.14729703962802886963,  0.27678227424621582031,
         -0.29839691519737243652, -0.01790418662130832672,
          0.25340676307678222656, -0.00136217952240258455,
         -0.29386615753173828125, -0.19074733555316925049,
         -1.28328192234039306641,  0.08512758463621139526,
         -0.01432672329246997833,  0.16793838143348693848,
         -0.00182083400432020426, -0.16354949772357940674,
         -0.14037890732288360596, -0.10419515520334243774,
          0.15291371941566467285,  0.52130317687988281250,
         -0.02314160391688346863,  0.11186494678258895874,
         -0.20050972700119018555,  0.29729220271110534668,
          0.09937122464179992676, -0.43684437870979309082,
          0.31647661328315734863,  0.21715416014194488525} };

double llayers2bias = 0.10583660751581192017;
*/

double zlT_MIN = 250.; // i generated a little bit less data than for Evibr(p,T)
double zlT_MAX = 1998.;

double zlP_MIN = 5;
double zlP_MAX = 499.;

double zllayers0weight[2][50] = {{        
-1.92654657e+00, -2.77099043e-01, -8.79129171e-02,
2.33522519e-01, 1.08874214e+00, 4.26189005e-02,
4.66783606e-02, -1.65077746e-01, 1.61720161e+01,
-3.62863392e-01, 8.82563069e-02, -1.35827041e+00,
5.33097267e-01, 1.34700462e-01, -2.94291586e-01,
-9.70418072e+00, -9.93295765e+00, -2.76893348e-01,
7.44624019e-01, 6.45287752e-01, 2.96161454e-02,
-3.73180079e+00, -7.15624452e-01, -1.17333829e+00,
-2.62197196e-01, -4.08124179e-01, 8.21605921e-01,
-1.41589195e-01, -8.15379918e-02, 2.24036314e-02,
5.15928209e-01, -4.94329840e-01, 7.27807581e-01,
4.26056767e+00, 5.12335002e-01, -3.26404542e-01,
-9.33327913e-01, -5.52368581e-01, 4.47492376e-02,
-4.15556692e-02, 1.99507475e-02, 3.67646724e-01,
5.28832197e-01, -1.12258089e+00, -1.24655724e-01,
-6.67710602e-01, 3.41551423e-01, 4.15353119e-01,
2.44453281e-01, -3.35405850e+00},
{-1.87120700e+00, 2.93855369e-02, -2.06160069e-01,
4.03824538e-01, 1.42376781e+00, 1.68917429e+00,
1.66204154e-01, -2.25755572e-03, 2.84141272e-01,
-4.17349428e-01, -4.40831408e-02, -3.75048876e-01,
6.22313499e-01, 1.95699975e-01, -3.58588099e-01,
-3.53798962e+00, -3.59485006e+00, -3.80317926e-01,
8.59285176e-01, 7.22268879e-01, -4.11672384e-01,
-4.02361989e-01, -7.56942511e-01, -1.58091044e+00,
-3.00684273e-01, -4.77788329e-01, 3.15385401e-01,
4.02834743e-01, -1.16720021e-01, -9.89575684e-02,
5.94067693e-01, -5.58759987e-01, 3.79225820e-01,
2.39446878e+00, 5.62086105e-01, -3.67004514e-01,
-1.17653131e+00, -5.66789448e-01, -5.84202480e+00,
-1.25382051e-01, -8.08334351e-02, 4.25798625e-01,
-1.29091293e-01, -1.76013482e+00, 4.99833226e-02,
-8.93086553e-01, -9.57194045e-02, 4.48255539e-01,
4.06332523e-01, -2.15468550e+00}};

double zllayers0bias[1][50] = { 
    -0.26304721832275390625, 0.51912444829940795898,
    -0.72981357574462890625, -0.29251429438591003418,
    0.26017591357231140137, 0.19593632221221923828,
    0.60714036226272583008, 0.94195455312728881836,
    0.48026159405708312988, -0.84667938947677612305,
    -0.89447790384292602539, 0.38649165630340576172,
    0.72627019882202148438, 0.21481493115425109863,
    -0.46387118101119995117, -0.21044062077999114990,
    -0.17914871871471405029, -0.66577774286270141602,
    0.84760004281997680664, 0.79382193088531494141,
    -0.26329776644706726074, 0.10819828510284423828,
    0.45865976810455322266, -0.25527858734130859375,
    -0.42044502496719360352, -0.15785746276378631592,
    -0.43555620312690734863, -0.10910969972610473633,
    0.36864736676216125488, -0.40394392609596252441,
    0.59041583538055419922, -0.73928987979888916016,
    -0.19730249047279357910, 0.56662011146545410156,
    0.97389727830886840820, -0.47055703401565551758,
    -0.92157924175262451172, -0.26250097155570983887,
    -0.07698321342468261719, -0.86952799558639526367,
    -1.04729819297790527344, 0.74410194158554077148,
    -0.21602410078048706055, -0.08917310833930969238,
    0.72368574142456054688, -0.57488399744033813477,
    -0.36609876155853271484, 0.73629838228225708008,
    0.04097023233771324158, -0.51260232925415039062 };

double zllayers2weight[2][50] = { {
    -2.72109091281890869141e-01, 3.19180935621261596680e-01,
    -1.61882728338241577148e-01, 5.83370178937911987305e-02,
    3.84644478559494018555e-01, 1.79975256323814392090e-02,
    1.75658613443374633789e-01, 3.77842664718627929688e-01,
    -1.72117745876312255859e+00, -3.17058116197586059570e-01,
    -1.97787642478942871094e-01, 4.78431999683380126953e-01,
    3.69668453931808471680e-01, 4.65353690087795257568e-02,
    -2.85812705755233764648e-01, -5.24899780750274658203e-01,
    -6.31338715553283691406e-01, -2.78480172157287597656e-01,
    3.56306016445159912109e-01, 2.26108297705650329590e-01,
    -1.83124281466007232666e-02, 8.86450648307800292969e-01,
    2.45271742343902587891e-01, -3.93756002187728881836e-01,
    -3.92337888479232788086e-01, -3.40860933065414428711e-01,
    -4.37015444040298461914e-01, 2.61658787727355957031e-01,
    2.74296365678310394287e-02, -7.30871036648750305176e-02,
    3.08211982250213623047e-01, -1.48325368762016296387e-01,
    -4.34642061591148376465e-02, 3.68761718273162841797e-01,
    2.70357400178909301758e-01, -3.40439170598983764648e-01,
    -2.34572634100914001465e-01, -2.22920298576354980469e-01,
    1.40972363948822021484e+00, -3.75186383724212646484e-01,
    -2.73018568754196166992e-01, 2.46187642216682434082e-01,
    -1.78010344505310058594e-01, -2.21055582165718078613e-01,
    3.96017879247665405273e-01, -1.13472409546375274658e-01,
    -1.20535068213939666748e-01, 2.48986110091209411621e-01,
    1.05917125940322875977e-01, -4.18973505496978759766e-01},
    {-1.01322503760457038879e-02, -3.24506200850009918213e-02,
    8.27905088663101196289e-02, -2.39178866147994995117e-01,
    -3.10708969831466674805e-01, -9.17156517505645751953e-01,
    2.25548557937145233154e-02, 9.96742025017738342285e-02,
    9.09607764333486557007e-03, -2.14000880718231201172e-01,
    -1.28265783190727233887e-01, -3.96351469680666923523e-03,
    1.71760946512222290039e-01, 2.72018108516931533813e-02,
    -9.74122360348701477051e-02, 6.87665713485330343246e-04,
    4.35131341218948364258e-02, 1.96375966072082519531e-01,
    8.52757245302200317383e-02, -2.89660897105932235718e-02,
    1.95958882570266723633e-01, 1.83068979531526565552e-02,
    1.15952089428901672363e-01, 4.18049283325672149658e-02,
    -1.50505751371383666992e-01, 5.47598302364349365234e-02,
    -2.06201970577239990234e-01, -2.49215915799140930176e-01,
    1.80862590670585632324e-01, -9.86364036798477172852e-02,
    1.55634611845016479492e-01, -4.75015910342335700989e-03,
    9.19781997799873352051e-02, 7.72011950612068176270e-02,
    -1.34350195527076721191e-01, -2.98010528087615966797e-01,
    -2.29861646890640258789e-01, -2.63750493526458740234e-01,
    1.97323396801948547363e-01, -1.46606057882308959961e-01,
    -1.80113181471824645996e-01, 1.87506750226020812988e-01,
    7.45202004909515380859e-02, 1.50642031803727149963e-02,
    1.91313892602920532227e-01, -1.73073917627334594727e-01,
    9.94940251111984252930e-02, 9.52403470873832702637e-02,
    -1.25124633312225341797e-01, 2.17552599497139453888e-03} };

double zllayers2bias[2][1] = {0.14365679025650024414, 0.21037971973419189453};
