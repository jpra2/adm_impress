import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

x_zC1_ENO3 = np.array([0, 0.0286839510401, 0.0787598235742, 0.139457850888, 0.160702160448, 0.163583545402, 0.169570685451, 0.175628679227, 0.183204123687, 0.19382037399, 0.210506427024, 0.231738927629, 0.252971428234, 0.275721379522, 0.306076297657, 0.333390409948, 0.374361578385, 0.424437450919, 0.474513323453, 0.495757633013, 0.507909047431, 0.533705709039, 0.564054722696, 0.638409806156, 0.69910783347, 0.70820663309, 0.718816978915, 0.726386518897, 0.730921157514, 0.736967342336, 0.739984530269, 0.742989909248, 0.749024285115, 0.752029664094, 0.756546589279, 0.759551968258, 0.764062988965, 0.767062463466, 0.771573484173, 0.774578863152, 0.779101692814, 0.783606809043, 0.791164540071, 0.797204820415, 0.803239196283, 0.812308473515, 0.822918819341, 0.833523260689, 0.848685958562, 0.853232406134, 0.857778853705, 0.859272686478, 0.860766519252, 0.860707474478, 0.863600668387, 0.868011312979, 0.869516954707, 0.872540047118, 0.887708649469, 0.910470409712, 0.9468892261, 0.993930197269, 1])

zC1_ENO3 = np.array([0.902723735409, 0.902723735409, 0.902723735409, 0.902723735409, 0.902723735409, 0.801556420233, 0.747081712062, 0.739299610895, 0.731517509728, 0.727626459144, 0.72373540856, 0.715953307393, 0.708171206226, 0.700389105058, 0.704280155642, 0.704280155642, 0.704280155642, 0.704280155642, 0.704280155642, 0.704280155642, 0.712062256809, 0.712062256809, 0.712062256809, 0.712062256809, 0.712062256809, 0.708171206226, 0.700389105058, 0.688715953307, 0.677042801556, 0.661478599222, 0.649805447471, 0.630350194553, 0.607003891051, 0.587548638132, 0.56420233463, 0.544747081712, 0.517509727626, 0.494163424125, 0.466926070039, 0.447470817121, 0.428015564202, 0.396887159533, 0.377431906615, 0.357976653696, 0.334630350195, 0.311284046693, 0.303501945525, 0.291828793774, 0.284046692607, 0.280155642023, 0.27626459144, 0.260700389105, 0.24513618677, 0.206225680934, 0.112840466926, 0.0194552529183, 0.011673151751, 0.00389105058366, 0, 0, 0, 0, 0])

x_zC2_ENO3 = np.array([0, 0.0333333333333, 0.110606060606, 0.151515151515, 0.157575757576, 0.163636363636, 0.172727272727, 0.181818181818, 0.19696969697, 0.216666666667, 0.25, 0.278787878788, 0.309090909091, 0.351515151515, 0.480303030303, 0.568181818182, 0.657575757576, 0.687878787879, 0.69696969697, 0.706060606061, 0.715151515152, 0.719696969697, 0.722727272727, 0.727272727273, 0.730303030303, 0.733333333333, 0.737878787879, 0.742424242424, 0.748484848485, 0.75303030303, 0.757575757576, 0.762121212121, 0.766666666667, 0.772727272727, 0.777272727273, 0.780303030303, 0.787878787879, 0.792424242424, 0.79696969697, 0.8, 0.801515151515, 0.804545454545, 0.810606060606, 0.815151515152, 0.822727272727, 0.831818181818, 0.837878787879, 0.843939393939, 0.851515151515, 0.857575757576, 0.863636363636, 0.863636363636, 0.871212121212, 0.883333333333, 0.90303030303, 0.930303030303, 0.963636363636, 1])

zC2_ENO3 = np.array([0.201550387597, 0.201679586563, 0.201979093258, 0.202137655626, 0.198285177355, 0.186680761099, 0.182840028189, 0.182875264271, 0.182933991074, 0.183010335917, 0.183139534884, 0.183251115809, 0.183368569415, 0.183533004463, 0.184032182288, 0.184372797745, 0.184719285882, 0.184836739488, 0.18487197557, 0.192659149636, 0.204322292694, 0.215967817712, 0.231483439042, 0.243128964059, 0.254768616397, 0.270284237726, 0.309061545689, 0.35946676063, 0.402125910265, 0.464159032182, 0.518440216115, 0.572721400047, 0.627002583979, 0.6774136716, 0.72781888654, 0.770466290815, 0.809255344139, 0.838652807141, 0.853050270143, 0.873565891473, 0.886827578107, 0.902343199436, 0.913994597134, 0.92176415316, 0.933421423538, 0.937332628612, 0.933480150341, 0.92962767207, 0.921905097486, 0.883168898285, 0.662262156448, 0.530479210712, 0.515004698144, 0.511175710594, 0.507376086446, 0.507376086446, 0.507376086446, 0.507376086446])

x_zC3_ENO3 = np.array([0, 0.0242056052588, 0.0635397138043, 0.105899523007, 0.122540876623, 0.14825933221, 0.160362134839, 0.161720733762, 0.166176226299, 0.170696979046, 0.176736514867, 0.187320534422, 0.196391703647, 0.206969790455, 0.219066660338, 0.235702081207, 0.255357269987, 0.268966990199, 0.282576710411, 0.299218064026, 0.314352500059, 0.347629274544, 0.383937682432, 0.414194689005, 0.432366691189, 0.46866916633, 0.508003274876, 0.532208880135, 0.613902797883, 0.677442511688, 0.756110728779, 0.772746149648, 0.786349937113, 0.802973492489, 0.812038728968, 0.825654381926, 0.836238401481, 0.843796720378, 0.849830323453, 0.854357008947, 0.858830299722, 0.861470371865, 0.865771612995, 0.876314103325, 0.894468307269, 0.917161062199, 0.951956619758, 0.986752177318, 1])

zC3 = np.array([0, 0, 0, 0, 0, 0, 0, 0.102600916016, 0.15752034932, 0.169302783645, 0.177169605354, 0.18113267994, 0.18508982178, 0.192974441729, 0.196943449061, 0.200930254633, 0.20885047106, 0.212825411139, 0.216800351219, 0.216865611429, 0.213003393531, 0.217055459313, 0.217197845227, 0.217316500154, 0.205623057026, 0.209686988301, 0.209841239707, 0.209936163649, 0.210256531954, 0.210505707302, 0.210814210114, 0.214801015686, 0.222697501127, 0.234527397423, 0.242406084625, 0.242459479342, 0.246422553929, 0.250373763022, 0.262162130093, 0.270023019056, 0.313177816275, 0.568090130283, 0.724969742993, 0.756383635112, 0.756454828069, 0.756543819265, 0.756680272432, 0.756816725599, 0.756858254823])

x_S = np.array([0, 0.155824508321, 0.160563231044, 0.16525486087, 0.169811098618, 0.177398941587, 0.186493757248, 0.201645896737, 0.215267517086, 0.230413769963, 0.247072882144, 0.262207361797, 0.27584075537, 0.298527758319, 0.31515743744, 0.334830495005, 0.452827634112, 0.46038604402, 0.483078933581, 0.523932021404, 0.584440506955, 0.719084985019, 0.7357440972, 0.752414982605, 0.76152157149, 0.769103527847, 0.781247608564, 0.79794792703, 0.813123612967, 0.828293412292, 0.841944465702, 0.84953230867, 0.855760344249, 0.860822830636, 1])

x_S_ENO3 = np.array([1, 1, 0.867704280156, 0.766536964981, 0.75486381323, 0.739299610895, 0.727626459144, 0.712062256809, 0.708171206226, 0.696498054475, 0.684824902724, 0.68093385214, 0.669260700389, 0.673151750973, 0.68093385214, 0.677042801556, 0.68093385214, 0.684824902724, 0.684824902724, 0.68093385214, 0.684824902724, 0.684824902724, 0.673151750973, 0.653696498054, 0.634241245136, 0.622568093385, 0.5953307393, 0.556420233463, 0.525291828794, 0.498054474708, 0.474708171206, 0.459143968872, 0.342412451362, 0, 0])


plt.figure(1)
plt.plot(zC1_x, zC1, '-g')
plt.grid()
plt.legend(('ENO3'))
plt.ylabel('C1 global molar fraction ')
plt.title('Mallison_ex1_C1')
plt.xlabel('Distance')
plt.savefig('C1_ENO.png')

plt.figure(2)
plt.plot(zC2_x, zC2, '-g')
plt.grid()
plt.legend(('ENO3'))
plt.ylabel('C2 global molar fraction ')
plt.title('Mallison_ex1_C2')
plt.xlabel('Distance')
plt.savefig('C2_ENO.png')

plt.figure(3)
plt.plot(zC3_x, zC3, '-g')
plt.grid()
plt.legend(('ENO3'))
plt.ylabel('C3 global molar fraction ')
plt.title('Mallison_ex1_C3')
plt.xlabel('Distance')
plt.savefig('C3_ENO.png')

plt.figure(4)
plt.plot(S_x, S, '-g')
plt.grid()
plt.legend(('ENO3'))
plt.ylabel('Saturation ')
plt.title('Mallison_ex1_S')
plt.xlabel('Distance')
plt.savefig('S_ENO.png')
import pdb; pdb.set_trace()