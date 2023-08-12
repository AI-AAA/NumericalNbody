import numpy as np

def lagrange_points(G,masses, positions):
    Ms, Me = masses[0],masses[1]
    R = np.linalg.norm(positions[1] - positions[0])
    Rv=positions[1] - positions[0]

    # L1
    L1 = np.array([R * (1-(Me / (3 * Ms))**(1 / 3)),0])
    # L2
    L2 = np.array([R * (1+(Me / (3 * Ms))**(1 / 3)),0])

    # L3
    L3 = np.array([-R * (1-(7*Me )/ (12*Ms)),0])

    # L4 and L5
    r=np.linalg.norm(R)
    L4 = np.array([r/2,np.sqrt(3)/2*r])
    L5 = np.array([r/2,-np.sqrt(3)/2*r])
    # Velocities
    L_list=np.array([L1, L2, L3, L4, L5])
    velocity=[]
    for i in range(5):
        θ = np.arctan2(L_list[i][1],L_list[i][0])
        v0 = np.sqrt( np.abs((G*Ms/(np.linalg.norm(positions[0]-L_list[i])**2)-G*Me/(np.linalg.norm(positions[1]-L_list[i])**2))*np.linalg.norm(L_list[i])))
        v=[-v0*np.sin(θ),v0*np.cos(θ)]
        velocity.append(v)

    return L_list, np.array(velocity)
def loadlagrangeData(Ln1,Ln2,N=20):
    G = 6.67e-11
    Msun = 2e30
    Me=5.97219e24
    AU = 1.5e11
    position=np.array([[0,0],[AU,0]])
    
    θ = np.arctan2(position[1][1],position[1][0])
    v0 = np.sqrt(G*Msun/AU) 
    #v0,[-v0*np.sin(θ),v0*np.cos(θ)]

    planet_name=["Sun","Earth"]
    planet_name+=[str(i) for i in range(N*5)]
    mass=[Msun,Me]
    
    velocity=np.array([[0,0],[-v0*np.sin(θ),v0*np.cos(θ)]])
    for i in range(N):
        pos_list,vel_list = lagrange_points(G,mass,position)
        pos_list +=np.random.rand(5, 2)
        position = np.vstack((position, pos_list[Ln1:Ln2]))
        velocity = np.vstack((velocity, vel_list[Ln1:Ln2]))
            
    mass = np.array(mass+[Me/1e24 for i in range(len(position)-2)])       
    momenta=np.array([mass[i]*velocity[i] for i in range(len(velocity))])
    T=365*3600*24
    def loadNBodyData():return G,planet_name,mass,position,momenta,T
    return loadNBodyData
    
def loadInnerSolarData():
    #mass unit: solar mass
    #veloctiy unit: au/day
    G = 2.95912208286E-04
    planet_name=["Sun","Mercury","Venus","Earth","Mars"]
    mass=np.array([1,1.6605481518732714e-07,
                   2.448327885340709e-06,3.003364344983656e-06,
                   3.2271058586874533e-07
                   ])
    position=np.array([[-3.524143431880468E-03 , 2.534181028743003E-03 , 1.674467975283666E-03],
          [-9.391450761268468E-03 , 2.748013882823011E-01 , 1.477277269382629E-01],
          [-2.221385884278196E-01 , 6.217431385999292E-01 , 2.941256985324626E-01],
          [-9.787816485876855E-01 ,-1.984028709388082E-01 ,-8.542683592683444E-02],
          [-1.078438465580282E+00 , 1.128512796230021E+00 , 5.471392912518396E-01]
          ])
    velocity=np.array([[-7.386335895962630E-04 ,-8.340058001867163E-04 ,-2.923490967948519E-04],
          [ -3.451180438095160E-02 ,-1.635655215876665E-03 , 2.779884358317152E-03],
          [-2.007431882024974E-02 ,-6.993009687456287E-03 ,-1.840309313507213E-03],
          [2.753735677460387E-03 ,-1.629015646703982E-02 ,-6.992089963868113E-03],
          [-1.082074452728994E-02 ,-8.144131358535762E-03 ,-3.373324468051556E-03]
          ])
    position = np.array([ (position[i, :] - position[0, :]) for i in range(len(position))])
    velocity = np.array([ (velocity[i, :] - velocity[0, :]) for i in range(len(velocity))])
    momenta = np.array([mass[i]*velocity[i] for i in range(len(velocity))])
    T=365
    return G,planet_name,mass,position,momenta,T

def loadNBodyChoreographedData3(N=3):
    G,planet_name,mass=1,["A","B","C"],np.array(np.ones(N))

    position=np.array([[0.7548983603205632,-0.3514110697402242],
                   [-0.6817800937262974,-0.4780556224427804],
                   [-0.07311826659426574,0.8294666921830046]
                   ])
    momenta=np.array([[-0.35141106974012876,-0.7548983603206323],
                   [-0.4780556224427675,0.6817800937264157],
                   [0.8294666921828963,0.07311826659421666]
                   ])
    T = 8*3.14
   
    return G,planet_name,mass,position,momenta,T
def loadNBodyChoreographedData(case):
    G=1
    init = [[0.7548983603205632,-0.3514110697402242,-0.35141106974012876,-0.7548983603206323,
            -0.6817800937262974,-0.4780556224427804,-0.4780556224427675,0.6817800937264157,
            -0.07311826659426574,0.8294666921830046,0.8294666921828963,0.07311826659421666],
        
            [-0.699242298027974,-0.6944467629710684,0.6944467629708051,-0.6992422980282715,
            0.6944467629708471,-0.6992422980282327,0.6992422980279377,0.694446762971098,
            0.699242298027974,0.6944467629710683,-0.6944467629708043,0.6992422980282715,
            -0.6944467629708471,0.6992422980282328,-0.6992422980279385,-0.694446762971098],
            
            [-0.07880227334416882,0.5570371897354746,0.15998292728488323,1.1593418791674066,
            0.5940359608209828,0.383319210563721,-0.5557289806160467,-0.9029539156799118,
            -0.5152336874768139,-0.9403564002991956,0.39574605333116347,-0.2563879634874948],
            
            [-1.1844475691074707,0.5898810975762604,0.4087319255099535,0.39755580006528823,
            0.20705693570351705,-0.20720580654667664,1.3649438991518443,0.2819370010443289,
            1.1835106293686215,-0.5912376198019829,-0.4103624168353448,-0.39612196466777105,
            -0.20611999596466776,0.20856232877239902,-1.363313407826453,-0.2833708364418462]
            ]
    N = int(len(init[case])/4)
    mass=np.array(np.ones(N))
    planet_name=[f"{i+1}" for i in range(N)]
    position=np.array(init[case]).reshape(N,4)[:, :2]
    momenta=np.array(init[case]).reshape(N,4)[:, 2:]
    T = [8*np.pi,8*np.pi,8*np.pi,9.9908e-2]
    def loaddata(): return G,planet_name,mass,position,momenta,T[case]
    return loaddata

def loadThreeBodyequalmassData(case=0):
    G,planet_name,mass=1,["A","B","C"],np.array([1,1,1])
    #L1,2,3,68
    init=np.array([[0.3471168881,0.5327249454],
                   [0.3068934205,0.1255065670],
                   [0.6150407229,0.5226158545],
                   [0.4160674674,0.2971499303]
                   ])
    period = np.array([6.3259139829,6.2346748391,37.3205235945,74.7892110946
                       ])
    init=init[case]
    position=np.array([[-1,0,0],[1,0,0],[0,0,0]])
    momenta=np.array([[init[0],init[1],0],[init[0],init[1],0],[-2*init[0],-2*init[1],0]])
    
    def loadThreeBodyequalmass(): return G,planet_name,mass,position,momenta,period[case]
    return loadThreeBodyequalmass
    
def loadPeriodicThreeBodyData(case):
    G= 1
    planet_name=["A","B","C"]
    mass=np.array([1,1,1])
    inits=[[-0.0347,1.1856,0.2693,-1.0020,-0.2328,-0.5978,0.2495,-0.1076,0.2059,-0.9396,-0.4553,1.0471],
    [  -0.602885898116520,1.059162128863347-1,0.252709795391000,1.058254872224370-1,-0.355389016941814,1.038323764315145-1,0.122913546623784,0.747443868604908,-0.019325586404545,1.369241993562101,-0.103587960218793,-2.116685862168820],
    [-0.97000436,0.24308753,  0.97000436,-0.24308753, 0, 0,-0.5*0.93240737,-0.5*0.86473146,-0.5*0.93240737,-0.5*0.86473146,0.93240737, 0.86473146],
    [0.716248295712871,0.384288553041130,0.086172594591232,1.342795868576616,0.538777980807643,0.481049882655556,1.245268230895990,2.444311951776573,-0.675224323690062,-0.962879613630031,-0.570043907205925,-1.481432338146543], 
    [1,0, -0.5,np.sqrt(3)/2, -0.5,-np.sqrt(3)/2,0,1, -np.sqrt(3)/2,-0.5,  np.sqrt(3)/2,-0.5],
    [0.486657678894505,0.755041888583519,-0.681737994414464,0.293660233197210,-0.022596327468640,-0.612645601255358,-0.182709864466916,0.363013287999004,-0.579074922540872,-0.748157481446087,0.761784787007641,0.385144193447218],
    [0.536387073390469,0.054088605007709,-0.252099126491433,0.694527327749042,-0.275706601688421,-0.335933589317989,-0.569379585580752,1.255291102530929,0.079644615251500,-0.458625997341406,0.489734970329286,-0.796665105189482],
    [0.517216786720872,0.556100331579180,0.002573889407142,0.116484954113653,-0.202555349022110,-0.731794952123173,0.107632564012758,0.681725256843756,-0.534918980283418,-0.854885322576851,0.427286416269208,0.173160065733631],
    [0.419698802831451,1.190466261251569,0.076399621770974,0.296331688995343,0.100310663855700,-0.729358656126973,0.102294566002840,0.687248445942575,0.148950262064203,0.240179781042517,-0.251244828059670,-0.927428226977280],
    [0.906009977920936,0.347143444586515,-0.263245299491018,0.140120037699958,-0.252150695248079,-0.661320078798829,0.242474965162160,1.045019736387070,-0.360704684300259,-0.807167979921595,0.118229719138103,-0.237851756465475],
    [1.666163752077218-1,-1.081921852656887+1,0.974807336315507-1,-0.545551424117481+1,0.896986706257760-1,-1.765806200083609+1,0.841202975403070,0.029746212757039,0.142642469612081,-0.492315648524683,-0.983845445011510,0.462569435774018],
    [1.451145020734434,-0.209755165361865,-0.729818019566695,0.408242931368610,0.509179927131025,0.050853900894748,0.424769074671482,-0.201525344687377,0.074058478590899,0.054603427320703,-0.498827553247650,0.146921917372517],
    [0.654504904492769,1.135463234751087,-0.008734462769570,-1.481635584087405,-0.487782115417060,0.236845409927442,0.294339951803385,-0.605376046698418,0.038678643994549,0.434105189236202,-0.333018595797923,0.171270857462206],
    [0.708322208567308,0.473311928206857,0.167739458699109,-0.057913961029346,-0.506578687023757,-0.306825234531109,0.824266639919723,0.522197827478425,-0.077017015655090,-0.167288552679470,-0.747249624264592,-0.354909274799606],
    [0.865355422048074,0.629488893636092,0.085036793558805,-0.013305780703954,-0.090983494772311,-0.892179296799402,0.288687787606716,0.171289709266963,-0.220256752038501,0.090375753070602,-0.068431035568003,-0.261665462337436],
    [0.335476420318203,-0.243208301824394,0.010021708193205,0.363104062311693,0.030978712523174,0.423035485079015,1.047838171160758,0.817404215288346,-0.847200907807940,-0.235749148338353,-0.200636552532016,-0.581655492859626],
    [-0.015958157830872,0.829387672711051,-0.023912093175726,0.211010340157083,-0.035305026035053,-0.728544174709096,0.008894140366887,0.735934914230558,-0.012693330102595,-1.001493752025633,0.003743585058501,0.265560077637935]]
    init=inits[case]
    position = np.array(init[:int(len(init)/2)]).reshape(3,2)
    velocity = np.array(init[int(len(init)/2):]).reshape(3,2)
    position = np.array([[x, y, 0] for x, y in position])
    velocity = np.array([[x, y, 0] for x, y in velocity])
    momenta = np.array([mass[i]*velocity[i] for i in range(len(velocity))])
    T=[2.9623,2.246101255307486,5,8.094721472532424,12,
       4.012156415940929,5.026055004808709,5.179484537709100,
       5.095053913455357,6.868155929188273,6.868155929188273,
       6.868155929188273,6.868155929188273,6.868155929188273,
       6.868155929188273, 5.956867234319493,1.253088248568183]
    def loadThreeBodySpecial(): return G,planet_name,mass,position,momenta,T[case]
    return loadThreeBodySpecial
def random_coordinates(num_points, min_norm, max_norm):
    coordinates = []

    for _ in range(num_points):
        while True:
            x = np.random.uniform(-max_norm, max_norm)
            y = np.random.uniform(-max_norm, max_norm)
            z = np.random.uniform(-max_norm, max_norm)
            norm = np.linalg.norm([x, y, z])

            if min_norm <= norm <= max_norm:
                coordinates.append((x, y, z))
                break
                
    return np.array(coordinates)
def generate_random_coordinates(num_points, min_length, max_length):
    coordinates = []
    for _ in range(num_points):
        theta = np.random.uniform(0, np.pi)
        phi = np.random.uniform(0, 2 * np.pi)
        r = np.random.uniform(min_length, max_length)

        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        coordinates.append((x, y, z))

    return np.array(coordinates)
def loadNBodyData(N=100):
    G = 2.95912208286E-4
    #mass unit: solor mass
    #velocity unit: AU
    planet_name=[str(i) for i in range(N)]
    mass =     np.vstack((np.ones((1, 1)), np.random.uniform(low=1/1000, high=1/1000, size=(N-1, 1))))
    position = np.vstack((np.zeros((1, 3)), generate_random_coordinates(N-1, 0.2, 2)))
    velocity = [np.sqrt((G*1)/np.linalg.norm(position[i+1])) for i in range(len(position)-1)]
    velocity = np.vstack((np.zeros((1, 3)), np.array([generate_random_coordinates(1, vel, vel)for vel in velocity]).reshape(N-1,3)))           
    momenta = np.array([mass[i]*velocity[i] for i in range(len(velocity))])
    T=365
    def loadNBodyData():return G,planet_name,mass,position,momenta,T
    return loadNBodyData
def loadOuterSolarData():
    #mass unit: solor mass
    #velocity unit: au/day
    G= 2.95912208286 * (10**-4)
    planet_name=["Sun","Jupiter","Saturn","Uranus","Neptune","Pluto"]
    mass=np.array([1.00000597682,0.000954786104043,
                   0.000954786104043,0.0000437273164546,
                   0.0000517759138449,1/(1.3*10**8)])
    position=np.array([[0,0,0],
          [-3.5023653,-3.8169847,-1.5507963],
          [9.0755314,-3.0458353,-1.6483708],
          [8.3101420,-16.2901086,-7.2521278],
          [11.4707666,-25.7294829,-10.8169456],
          [-15.5387357,-25.2225594,-3.1902382]
          ])
    velocity=np.array([[0,0,0],
          [0.00565429,-0.00412490,-0.00190589],
          [0.00168318,0.00483525,0.00192462],
          [0.00354178,0.00137102,0.00055029],
          [0.00288930,0.00114527,0.00039677],
          [0.00276725,-0.00170702,-0.00136504]
          ])
    momenta = np.array([mass[i]*velocity[i] for i in range(len(velocity))])
    T=365
    return G,planet_name,mass,position,momenta,T


def loadThreeBodyData():
    G = 6.673997487650961e-20
    planet_name=["A","B","C"]
    mass=np.array([1,1,1])
    position=np.array([[-0.97000436,0.24208753,0],
                       [0,0,0],
                       [0.97000436,-0.24208753,0]
                    ])
    velocity=np.array([[0.4662036850,0.4323657300,0],
                       [-0.933240737,-0.86473146,0],
                       [0.4662036850,0.4323657300,0]
                    ])
    momenta = np.array([mass[i]*velocity[i] for i in range(len(velocity))])
    return G,planet_name,mass,position,momenta
