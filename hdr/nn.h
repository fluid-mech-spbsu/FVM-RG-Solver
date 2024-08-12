#pragma once

// macroparameters range used to train NN:
double T_MIN = 200.01;
double T_MAX = 1200.0099999990905;

double T_MINc = 200.001;
double T_MAXc = 1200.0009999779447;

/// --------------------------------------------------------------------------- ///
/// ------------------------ NN for U calculation ----------------------------- ///
/// --------------------------- (T) -> (U)  ----------------------------------- ///

// layers weights and biases:
double U_layers0weight[50] = { 
-1.045205116271972656e+00,
- 4.474741458892822266e+00,
- 1.062483072280883789e+00,
1.031360030174255371e+00,
- 1.041225075721740723e+00,
- 1.010295510292053223e+00,
4.286369800567626953e+00,
1.019177436828613281e+00,
- 1.027483940124511719e+00,
1.010388851165771484e+00,
1.011942267417907715e+00,
1.012321352958679199e+00,
- 1.016682147979736328e+00,
- 1.037028193473815918e+00,
- 3.653058767318725586e+00,
1.011606693267822266e+00,
1.052469968795776367e+00,
- 1.057688593864440918e+00,
- 1.076019644737243652e+00,
1.009980440139770508e+00,
1.017671346664428711e+00,
1.008183240890502930e+00,
1.170516848564147949e+00,
- 1.014241099357604980e+00,
- 1.137836694717407227e+00,
- 5.325095653533935547e+00,
1.021373629570007324e+00,
- 1.010689139366149902e+00,
- 1.012288331985473633e+00,
4.015087127685546875e+00,
1.053568482398986816e+00,
1.023804903030395508e+00,
1.009511232376098633e+00,
- 1.037393569946289062e+00,
4.180373668670654297e+00,
1.084039211273193359e+00,
1.012554049491882324e+00,
1.019868850708007812e+00,
1.011994957923889160e+00,
1.024624347686767578e+00,
- 1.021714091300964355e+00,
- 1.041794419288635254e+00,
- 1.021737694740295410e+00,
- 1.013175845146179199e+00,
1.013158440589904785e+00,
- 5.113942623138427734e+00,
- 1.011857748031616211e+00,
1.029571890830993652e+00,
- 1.009374976158142090e+00,
1.014750123023986816e+00
};

double U_layers0bias[50] = {
-1.483242988586425781e+00,
3.679351329803466797e+00,
- 1.275290608406066895e+00,
1.735049605369567871e+00,
- 1.544365525245666504e+00,
- 4.130483150482177734e+00,
- 4.002218723297119141e+00,
2.138432979583740234e+00,
- 1.833636879920959473e+00,
3.816214561462402344e+00,
2.768021821975708008e+00,
2.700391530990600586e+00,
- 2.280236244201660156e+00,
- 1.617747902870178223e+00,
1.134980440139770508e+00,
3.424731254577636719e+00,
1.386350750923156738e+00,
- 1.326136112213134766e+00,
- 1.152704358100891113e+00,
3.710908651351928711e+00,
2.219212055206298828e+00,
3.724607706069946289e+00,
6.442546844482421875e-01,
- 2.472927570343017578e+00,
- 7.802239060401916504e-01,
5.743694782257080078e+00,
2.038328647613525391e+00,
- 3.050732612609863281e+00,
- 2.707691431045532227e+00,
- 2.093955516815185547e+00,
1.373166799545288086e+00,
1.947587847709655762e+00,
4.200872898101806641e+00,
- 1.610949158668518066e+00,
- 2.850971460342407227e+00,
1.090004682540893555e+00,
2.667131662368774414e+00,
2.104908227920532227e+00,
2.748806238174438477e+00,
1.919087052345275879e+00,
- 2.024462938308715820e+00,
- 1.535134434700012207e+00,
- 2.023525714874267578e+00,
- 2.585970401763916016e+00,
2.589997053146362305e+00,
5.512524127960205078e+00,
- 2.779947996139526367e+00,
1.777921795845031738e+00,
- 3.384214639663696289e+00,
2.426962852478027344e+00
} ;

double U_layers2weight[50] = { 
	-2.257354164123535156e+01,-2.450259208679199219e+01,
	-2.258590316772460938e+01,2.242447662353515625e+01,
	-2.243202209472656250e+01,-2.257923889160156250e+01,
	2.409734535217285156e+01,2.267354965209960938e+01,
	-2.273462104797363281e+01,2.259887886047363281e+01,
	2.280170631408691406e+01,2.278064727783203125e+01,
	-2.274802017211914062e+01,-2.251764488220214844e+01,
	-2.462505912780761719e+01,2.243091201782226562e+01,
	2.245374298095703125e+01,-2.293915176391601562e+01,
	-2.300812530517578125e+01,2.244137573242187500e+01,
	2.232524490356445312e+01,2.294830894470214844e+01,
	2.306818962097167969e+01,-2.236000442504882812e+01,
	-2.249052429199218750e+01,-2.464269256591796875e+01,
	2.234285163879394531e+01,-2.296242332458496094e+01,
	-2.247395706176757812e+01,2.537311935424804688e+01,
	2.255511474609375000e+01,2.294081878662109375e+01,
	2.274678039550781250e+01,-2.244305419921875000e+01,
	2.460834121704101562e+01,2.265431213378906250e+01,
	2.244418144226074219e+01,2.267352294921875000e+01,
	2.253345870971679688e+01,2.241946983337402344e+01,
	-2.240127944946289062e+01,-2.250557327270507812e+01,
	-2.249547576904296875e+01,-2.236046409606933594e+01,
	2.293047714233398438e+01,-2.417689704895019531e+01,
	-2.263146591186523438e+01,2.237266731262207031e+01,
	-2.292586898803710938e+01,2.291049957275390625e+01
};

double U_layers2bias = 2.236671829223632812e+01;


/// --------------------------------------------------------------------------- ///




double C_layers0weight[50] =
{ -1.146090745925903320e+00,
-1.308274459838867188e+01,
1.334185695648193359e+01,
1.759899711608886719e+01,
-1.794814705848693848e+00,
1.377716541290283203e+01,
-3.018612623214721680e+00,
1.561042904853820801e+00,
-1.146209478378295898e+00,
1.757652282714843750e+01,
1.147100329399108887e+00,
1.386612224578857422e+01,
-1.147094607353210449e+00,
-1.331074047088623047e+01,
-1.146094799041748047e+00,
-1.146203637123107910e+00,
-1.477080917358398438e+01,
1.763002395629882812e+01,
-1.145826816558837891e+00,
1.664743041992187500e+01,
2.501905441284179688e+00,
-1.520271492004394531e+01,
1.147057414054870605e+00,
-1.578381538391113281e+01,
-1.146091222763061523e+00,
-1.376730728149414062e+01,
-1.428048515319824219e+01,
-2.721378803253173828e+00,
-1.280115509033203125e+01,
2.501838684082031250e+00,
-1.336396408081054688e+01,
-1.146032333374023438e+00,
-1.431108951568603516e+01,
-1.317889690399169922e+01,
2.501876115798950195e+00,
-1.588150596618652344e+01,
-1.147091627120971680e+00,
-1.503754234313964844e+01,
-1.397862815856933594e+01,
1.146970868110656738e+00,
1.314689064025878906e+01,
-1.146965146064758301e+00,
-1.284422874450683594e+01,
1.146209239959716797e+00,
-1.147048354148864746e+00,
-1.304532337188720703e+01,
-1.147055745124816895e+00,
1.624380874633789062e+01,
1.146147370338439941e+00,
-1.146040558815002441e+00
};

double C_layers0bias[50] = {
-7.378516197204589844e+00,
4.895833969116210938e+00,
-5.855842590332031250e+00,
-1.786227798461914062e+01,
-7.919993400573730469e+00,
-2.103714704513549805e+00,
5.268118530511856079e-02,
7.810003757476806641e+00,
-7.378392696380615234e+00,
-1.604640007019042969e+01,
7.377569675445556641e+00,
-7.409589767456054688e+00,
-7.377572536468505859e+00,
4.347427845001220703e+00,
-7.378526210784912109e+00,
-7.378398418426513672e+00,
1.011956119537353516e+01,
-1.710246849060058594e+01,
-7.378777980804443359e+00,
-1.436084079742431641e+01,
8.287686347961425781e+00,
1.336921691894531250e+00,
7.377624988555908203e+00,
1.146970653533935547e+01,
-7.378515243530273438e+00,
5.471219062805175781e+00,
8.653259277343750000e+00,
-8.034616470336914062e+00,
2.791937589645385742e+00,
8.287709236145019531e+00,
2.282803297042846680e+00,
-7.378611087799072266e+00,
8.191262245178222656e+00,
6.104412078857421875e+00,
8.287691116333007812e+00,
1.224692916870117188e+01,
-7.377571105957031250e+00,
9.733752250671386719e+00,
6.998837471008300781e+00,
7.377734184265136719e+00,
-4.184040546417236328e+00,
-7.377709388732910156e+00,
3.417123556137084961e+00,
7.378398418426513672e+00,
-7.377615451812744141e+00,
3.222812175750732422e+00,
-7.377624988555908203e+00,
-1.325002670288085938e+01,
7.378444194793701172e+00,
-7.378608703613281250e+00
};

double C_layers2weight[50] = {
-5.573278045654296875e+01, -5.572343444824218750e+01, 5.623962402343750000e+01,
5.523085021972656250e+01, -5.588820266723632812e+01, 5.556607055664062500e+01,
-5.551027679443359375e+01, 5.574328231811523438e+01, -5.626499557495117188e+01,
5.636742401123046875e+01, 5.612387084960937500e+01, 5.583558273315429688e+01,
-5.608723068237304688e+01, -5.585348129272460938e+01, -5.589969635009765625e+01,
-5.571704864501953125e+01, -5.619874191284179688e+01, 5.447074890136718750e+01,
-5.591801834106445312e+01, 5.579175186157226562e+01, 5.571110916137695312e+01,
-5.540037155151367188e+01, 5.583506774902343750e+01, -5.609511566162109375e+01,
-5.582115936279296875e+01, -5.565476608276367188e+01, -5.604417419433593750e+01,
-5.601683044433593750e+01, -5.581887435913085938e+01, 5.586734008789062500e+01,
-5.591412353515625000e+01, -5.588951492309570312e+01, -5.566689300537109375e+01,
-5.583782577514648438e+01, 5.590239715576171875e+01, -5.566028976440429688e+01,
-5.581641769409179688e+01, -5.610017013549804688e+01, -5.597616195678710938e+01,
5.584181594848632812e+01, 5.606116867065429688e+01, -5.588084793090820312e+01,
-5.600899887084960938e+01, 5.614326477050781250e+01, -5.573280334472656250e+01,
-5.626875305175781250e+01, -5.598905563354492188e+01, 5.613317871093750000e+01,
5.578377914428710938e+01, -5.590976715087890625e+01 };

double C_layers2bias = 5.548314666748046875e+01;
