import numpy as np
import copy

def norm(vector):
    return np.sqrt(np.array(vector).dot(np.array(vector)))

def negH_q(G,mass,current_position):
    result=copy.deepcopy(current_position)
    for i in range(len(mass)):
        hq=0
        for j in range(len(mass)):
            if i ==j: continue
            q_i=current_position[i]
            q_j=current_position[j]
            r=q_j-q_i
            hq+=G*mass[i]*mass[j]*r/norm(r)**3
        result[i]=copy.deepcopy(hq)
    #print(result,type(result))
    return np.array(result)
   
def H_p(mass,current_momenta):
    return np.array([current_momenta[i]/mass[i] for i in range(len(current_momenta))])

def energy_H(G,masses,momenta,coords):
    N=len(masses)
    #G=2.95912208286 * (10**-4)
    def my_norm(p): return np.sum(p**2, axis = 0)
    
    coords_0 = coords[0, :]
    momenta_0 = momenta[0, :]#NEW
    value    = 0.5 / masses[0] * my_norm(momenta_0)**2
    for i in range(1, N):
        for j in range(0, 1):
            coords_i = coords[i, :] 
            coords_j = coords[j, :] 
            T = 0.5 / masses[i]* my_norm (coords_i)**2
            U = - G * masses[i] * masses[j] * my_norm(coords_i - coords_j)
            value = value + T + U
    return value
def Totalenergy(G,mass, momenta, position):
    N = len(mass)
    # Calculate the first term
    term1=0.5*np.sum([(1/mass[i])*np.dot(momenta[i],momenta[i])for i in range(N)])
    # Calculate the second term
    term2 = 0
    for i in range(1, N):
        for j in range(i):
            distance = np.linalg.norm(position[i] - position[j])
            term2 += (mass[i] * mass[j]) / distance
    term2 *= -G
    # Calculate the total Energy
    H_value = term1 + term2

    return H_value

def Numerical_Method(G,method,h,mass,position0,momenta0):
    
    position0 = np.asarray(position0, dtype=np.float64)
    momenta0 = np.asarray(momenta0, dtype=np.float64)
    if method=="Explicity Euler":
        position = position0+h*H_p(mass,momenta0) 
        momenta  = momenta0+h*negH_q(G,mass,position0) 
    elif method=="Symplectic Euler":
        position =position0+ h*H_p(mass,momenta0)
        momenta  =momenta0+ h*negH_q(G,mass,position) 
    elif method=="Stormer–Verlet":
        momenta_mid = momenta0 + h/2 * negH_q(G,mass,position0)
        position    = position0 + h * H_p(mass,momenta_mid)
        momenta     = momenta_mid + h/2 * negH_q(G,mass,position)
    elif method=="4th order Runge-Kutta":
        momenta1 = h*negH_q(G,mass,position0)
        position1 = h*H_p(mass,momenta0)
        
        momenta2 = h*negH_q(G,mass,position0 + position1/2)
        position2 = h*H_p(mass,momenta0  +momenta1/2)
            
        momenta3 = h*negH_q(G,mass,position0 +position2/2)
        position3 = h*H_p(mass,momenta0  +momenta2/2)
            
        momenta4 = h*negH_q(G,mass,position0 +position3)
        position4 = h*H_p(mass,momenta0  +momenta3)
        
        momenta  = momenta0+(momenta1+2*momenta2+2*momenta3+momenta4)/6
        position = position0+(position1+2*position2+2*position3+position4)/6
    elif method == "yoshida6_step":
        w1 = 1.3512071919596578
        w0 = -1.7024143839193157
        w = [w1, w0, w1]

        position = copy.deepcopy(position0)
        momenta = copy.deepcopy(momenta0)

        for wi in w:
            half_h = 0.5 * h * wi
            full_h = h * wi

            # First half-step for momenta
            momenta += half_h * negH_q(G, mass, position)

            # Full step for position
            position += full_h * H_p(mass, momenta)

            # Second half-step for momenta
            momenta += half_h * negH_q(G, mass, position)
    elif method == "8":
        def k_value(position, momenta, factor):
            return h * negH_q(G, mass, position + factor[0] * position)

        def l_value(position, momenta, factor):
            return h * H_p(mass, momenta + factor * momenta)
        c = [0, 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200, 549002324/971916982, 13/20, 1201146811/1299019798, 1, 1]
        a = [        [1/18],
        [1/48, 1/16],
        [1/32, 0, 3/32],
        [5/16, 0, -75/64, 75/64],
        [3/80, 0, 0, 3/16, 3/20],
        [29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 23124283/1800000000],
        [16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555],
        [39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287],
        [246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789],
        [-1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653],
        [185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083],
        [403863854/491063109, 0, 0, -5068492393/434740067, -411421997/543043805, 652783627/914296604, 11173962825/925320556, -13158990841/6184727034, 3936647629/1978049680, -160528059/685178525, 248638103/1413531060],
    ]
        b = [14005451/335480064, 0, 0, 0, -59238493/1068277825, 181606767/758867731, 561292985/797845732, -1041891430/1371343529, 760417239/1151165299, 118820643/751138087, -528747749/2220607170, 1/4]
        position = copy.deepcopy(position0)
        momenta = copy.deepcopy(momenta0)
        k = [k_value(position0, momenta0, [sum(a[i][j] * position[j] for j in range(i)), sum(a[i][j] * momenta[j] for j in range(i))]) for i in range(len(a))]
        l = [l_value(position0, momenta0, sum(a[i][j] * momenta[j] for j in range(i))) for i in range(len(a))]

        momenta = momenta0 + sum(b[i] * k[i] for i in range(len(b)))
        position = position0 + sum(b[i] * l[i] for i in range(len(b)))
        pass
    else:
        print("Method Error")
    totoalenergy = Totalenergy(G,mass, momenta, position)
    energy= energy_H(G,mass, momenta, position)
    return position,momenta,energy,totoalenergy

def SimulateOribt(datafunction,Ts,step,method):#data,number of perid,number total step
    all_positions=list([])
    allmomenta=list([])
    allenergy=list([])
    alltotalenergy=list([])
    G,planet_name,mass,position,momenta,T=datafunction()
    all_positions.append(position)
    allmomenta.append(momenta)
    allenergy.append(energy_H(G,mass, momenta, position))
    alltotalenergy.append(Totalenergy(G,mass, momenta, position))
    h=Ts*T/step
    print(h,step*Ts)
    for i in range(int(step*Ts)):
        position,momenta,energy,totoalenergy = Numerical_Method(G,method,h,mass,all_positions[-1],allmomenta[-1])
        all_positions.append(position)
        allmomenta.append(momenta)
        allenergy.append(energy)
        alltotalenergy.append(totoalenergy)
        #if i == 20:
            #print(all_positions)
    return planet_name,all_positions,allenergy,alltotalenergy,method,Ts*T
    #plot_energy(Ts,h,allenergy)