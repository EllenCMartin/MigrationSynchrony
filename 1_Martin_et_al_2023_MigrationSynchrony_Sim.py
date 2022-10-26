
# -*- coding: utf-8 -*-
"""
Created on May 25th, 2021. 
Updated October 26, 2022.

@author: Aline Lee / Ellen Martin
"""
## Set working directory:
import os
if os.environ['HOMEPATH'] == "\\Users\\ellencm":
  os.chdir('//home.ansatt.ntnu.no/ellencm/Desktop/Dissertation/PAPER_I_ MIGRATION/Code_Saving')

if os.environ['HOMEPATH'] == "\\Users\\lee":
  os.chdir('//home.ansatt.ntnu.no/lee/Documents/Aline_work/Students/Ellen/Migration_paper/May22')


from math import exp, pi, sqrt, log
import numpy as np
from scipy import integrate
from scipy.special import j0, logit, expit
from scipy.stats import norm
from itertools import combinations, chain
import pickle
rng = np.random.default_rng()



# Function for calculating distances from each point to each of the other (thinned) points
def PointDist(gridsize, thin):  

    # Total number of x-values (=y-values)
    xmax = round(gridsize/thin) 
    # Indices of full grid
    gridallY, gridallX = np.indices((gridsize, gridsize))
    gridallY = np.flip(gridallY, 0)
    #Indices of thinned grid
    gridthinY, gridthinX = np.indices((xmax, xmax))*thin
    gridthinY = np.flip(gridthinY, 0)

    # Split into x an y values
    gridallXX = np.concatenate(gridallX, axis=None)
    gridallYY = np.concatenate(gridallY, axis=None)
    gridthinXX = np.concatenate(gridthinX, axis=None)
    gridthinYY = np.concatenate(gridthinY, axis=None)

    # For focal point in grid, x-x and y-y for each other point in grid (distance on each axis)
    xminy = [np.stack((gridthinXX-gridallXX[i], gridthinYY-gridallYY[i]), axis=-1) for i in range(0, gridsize**2)]
    distlist =  np.asarray([[e for e in np.round(np.sqrt(np.einsum('ij,ij->i', xminy[k], xminy[k])))] for k in range(0, gridsize**2)], dtype="int64")
    # distlist[i] gives the distances from point i to all thinned points. All points (not just thinned ones) are included as i
    # Storing files for various grid sizes 
    with open('PDt%s_thin%s.pickle' %(gridsize, thin),'wb') as f:
        pickle.dump(distlist, f)      

    return("Distlist stored as PDt%s_thin%s" %(gridsize, thin))


# Function for calculating xi as a funtion of distance - this only needs to be done once
def xi_dist(gridsize, le, sigma):   
    def xi_calc(dist, le, sigma):
        def intxi(u, dist, le):
            return exp(-u**2*le**2/4)*j0(dist*u)*u 
        xi_int = integrate.quad(intxi, a=0, b=np.inf, args=(dist, le), limit=400)
        xi = le*sigma/sqrt(2*pi)*xi_int[0]
        return(xi)
    iterator = (xi_calc(k, le, sigma) for k in range(0, round(sqrt(gridsize**2+gridsize**2))))
    xi_list = np.fromiter(iterator, float)
    return(xi_list)


# Function for creating sddv object which is applied to the nonbreeding ground survival.
def survDensdep(abundances, sdd):
    sddv = -abundances*sdd
    return(sddv)


### Function for drawing environmental field given xi-values.
def env_f(gridsize, thin, distlist, xi_list):
    env_field = np.zeros((int(gridsize/thin)**2))
    rand = rng.normal(size=gridsize**2)
    iterable = ((xi_list[f] for f in distlist[i]) for i in range(gridsize**2))
    env_field = np.einsum('ij,i->j', np.fromiter(chain.from_iterable(iterable), float).reshape(gridsize**2, (int(gridsize/thin)**2)), rand)
    return(env_field)


### Function for creating dispersal field
def Dispchange(new_abundances, gridsize, thin, thinned_indices, distlist, disprate, lg, numbermigroutes, thinnedpoints, migfec):  
            
    disp_probs_t = norm.pdf(range(1, round(sqrt(gridsize**2+gridsize**2))+1), 0, lg)  
    dispnums = new_abundances*disprate
    nodisp = new_abundances - dispnums
    newAbundances = np.zeros((1, thinnedpoints), dtype="float64")
    migfecfield = np.zeros((1, thinnedpoints), dtype="float64")
    for i in range(thinnedpoints):
       dprobs = disp_probs_t[distlist[thinned_indices,i][np.arange(thinnedpoints)!=i]]/sum(disp_probs_t[distlist[thinned_indices,i][np.arange(thinnedpoints)!=i]])       
       newAbundances += np.insert(dispnums[i]*dprobs, i, nodisp[i])  
       migfecfield += np.insert(dispnums[i]*dprobs, i, nodisp[i])*migfec[i]
    migfecfield = np.where(newAbundances>0, migfecfield/newAbundances, 0)
    return newAbundances, migfecfield  


### Function for creating distances files
def distances(gridsize, thin, zonewidth, midsection):
    distlist = distlistload(150, 2)
    thinnedpoints = int((gridsize/thin)**2) 
    gridsize_th = int(gridsize/thin)  
    thinned_indices = (np.ravel_multi_index(np.indices((round(gridsize/thin), round(gridsize/thin)))*thin, (gridsize, gridsize))+gridsize).reshape(thinnedpoints)
    full_to_thinned = dict(map(reversed, enumerate(thinned_indices))) # Use to easily convert indices of full grid to thinned ones
    indices_random_th = [full_to_thinned[e] for e in thinned_indices] 
    y_indices, x_indices = np.unravel_index(indices_random_th, (gridsize, gridsize))
    coords = list(zip(x_indices,y_indices))
 
    midsection_th = int(midsection/thin) 
    side = (gridsize_th-midsection_th)/2 

    midpoints_row1_full = np.int64(np.arange(gridsize*side*thin+side*thin, gridsize*side*thin+side*thin+midsection))
    midpoints_full = np.int64(np.concatenate([midpoints_row1_full+gridsize*a for a in range(0, midsection)]))
    midpoints_th = midpoints_full[np.isin(midpoints_full, thinned_indices)] 
    numthinnedpoints = len(midpoints_th)
   
    y_indicesmid, x_indicesmid = np.unravel_index(midpoints_th, (gridsize, gridsize))
    coords = list(zip(x_indicesmid,y_indicesmid))

    midpoint_indices_thinned = [i for i, x in 
         enumerate(np.isin(thinned_indices, midpoints_th)) if x]
    distancesmid = [distlist[midpoints_th[a]][midpoint_indices_thinned[b]] for a,b in combinations(np.arange(numthinnedpoints),2)] 

    with open('Midpoint_Distances_Gridsize%s_Thin%s.csv' %(gridsize, thin), 'w')  as output:
        output.write(str(distancesmid))
        
    with open('Coords_Gridsize%s_Thin%s.csv' %(gridsize, thin), 'w')  as output:
         output.write(str(coords))       




# Parameters:
    # nit - number of iterations for the simulation
    # abundances - one-dimensional numpy.ndarray of abundances (initial size, (gridsize**2,)) All start at the same value. 
    # xi_list - output from xi_dist function. Gives the xi-values that will produce the appropriate environmental field associated with different distances.
    # distlist - distlist[i] gives the distances from point i to all thinned points. All points (not just thinned ones) are included as i.
    # gridsize - side length of spatial grid.
    # thin - reducing number of points in grid for saving purposes (gridsize divided by thin).
    # midsection - one side length of inner grid which we are saving.
    # le - spatial scale of environmental noise.
    # sigma - standard deviation of environmental noise.
    # disprate - dispersal rate (probability of dispersing at breeding ground).
    # lg - spatial scale of dispersal distribution.
    # breedsurv - survival probability at breeding ground. 
    # nonbreedsurv - survival probability at nonbreeding ground.
    # survenv - relative effect of environmental noise on survival. 
    # fec - fecundity parameter (mean number of offspring).
    # fecenv - relative effect of environmental noise on fecundity.
    # sdd - density dependence acting on survival (see survDensdep function).
    # randommigroutes - 0 = not random, 1 = random.
    # correlation_between - correlation between different migraiton routes.
    # correlation_within - correlation within the different migraion routes.
    # weight - weighting of the effect of migration and nonbreeding ground. Default = 1 on breeding ground.
    # numbermigroutes - different number of migration routes in simulation (1, 2, or 4).
    # zonewidth - width of transition zone between migration routes where individuals have a non-zero probability of migrating in a different migration route.
    # vsplit - location of vertical division of migration route assignment. Used when changing sizes of migration routes.  Number refers to location in midsection.
    # hsplit - location of horizontal division of migration route assignment. Used when changing sizes of migration routes. Number refers to location in midsection.
    # carryover - 0 = no carryover effect, 1 = carry over effect.
    # burn_in: number of iterations at which to start saving output.
    

def breed_nonbreed_update(nit, abundances, xi_list, distlist, gridsize, thin, midsection, le, sigma, disprate, lg, breedsurv, nonbreedsurv, survenv, fec, fecenv, sdd, randommigroutes, correlation_between, correlation_within, weight, numbermigroutes, zonewidth, vsplit, hsplit, carryover, burn_in=20):
    counter = 0      
    new_abundances = abundances
    thinnedpoints = int((gridsize/thin)**2)  
    gridsize_th = int(gridsize/thin)
    zonewidth_th = int(zonewidth/thin)
    midsection_th = int(midsection/thin)
    side = (gridsize_th-midsection_th)/2
    hsplit_th_full = int(hsplit/thin)+side
    vsplit_th_full = int(vsplit/thin)+side

    if (zonewidth % 2) != 0:
        raise ValueError('Please select an even zonewidth') 
    if (zonewidth % 4) != 0:
        raise ValueError('Please select a zonewidth divisible by 4') 
    if (gridsize-max(vsplit, (gridsize-vsplit)))<(zonewidth/2) or (gridsize-max(hsplit, (gridsize-hsplit)))<(zonewidth/2):
        raise ValueError('Zonewidth is too wide for size of migration area')   

            
    # If simulating random migration routes(randommigroutes == 1), set up the migration routes:
    if randommigroutes == 1:
        if numbermigroutes == 0:
            raise ValueError('Random migration route cannot be selected with no migration')

        elif numbermigroutes == 1:
            raise ValueError('Random migration route cannot be selected with one migration route')

        elif numbermigroutes == 2:
            migroute = np.random.choice(2, thinnedpoints, replace=True) 
         

        elif numbermigroutes == 4:
            migroute = np.random.choice(4, thinnedpoints, replace=True)
                   
    # If simulating proximity grid migration (randommigroutes == 0), set up the migration routes:
    else:    
        migrouteprobs = np.zeros((thinnedpoints, numbermigroutes))     
          
        if numbermigroutes == 0:
            migroute = np.zeros(thinnedpoints, dtype=int)  

        elif numbermigroutes == 1:  
            migroute = np.zeros(thinnedpoints, dtype=int)

        elif numbermigroutes == 2:  
            migroute = np.empty(thinnedpoints, dtype=int)
            migrouteprobs[np.arange(gridsize_th*hsplit_th_full-gridsize_th*zonewidth_th/2, dtype=int), 0] = 1
            migrouteprobs[np.arange(gridsize_th*hsplit_th_full-gridsize_th*zonewidth_th/2, gridsize_th*hsplit_th_full+gridsize_th*zonewidth_th/2, dtype=int), 0] = 0.5            
            migrouteprobs[np.arange(gridsize_th*hsplit_th_full-gridsize_th*zonewidth_th/2, gridsize_th*hsplit_th_full+gridsize_th*zonewidth_th/2, dtype=int), 1] = 0.5            
            migrouteprobs[np.arange(gridsize_th*hsplit_th_full+gridsize_th*zonewidth_th/2, thinnedpoints, dtype=int), 1] = 1
            
            for i in np.arange(thinnedpoints):
                migroute[i] = rng.choice(np.arange(2), p=migrouteprobs[i,])  
       
        elif numbermigroutes == 4:      
            migroute = np.empty(thinnedpoints, dtype=int)
            # Upper section (migration route 1 and 2)
            for i in np.arange(hsplit_th_full-zonewidth_th/2, dtype=int):
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, dtype=int)+gridsize_th*i,0] = 1
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, vsplit_th_full+zonewidth_th/2, dtype=int)+gridsize_th*i,0] = 0.5
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, vsplit_th_full+zonewidth_th/2, dtype=int)+gridsize_th*i,1] = 0.5
                migrouteprobs[np.arange(vsplit_th_full+zonewidth_th/2,gridsize_th, dtype=int)+gridsize_th*i,1] = 1
        
            # Mid section (between 1 and 3, and between 2 and 4, including middle square)
            for i in np.arange(hsplit_th_full-zonewidth_th/2, hsplit_th_full+zonewidth_th/2, dtype=int):
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, dtype=int)+gridsize_th*i,0] = 0.5
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, dtype=int)+gridsize_th*i,2] = 0.5
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, vsplit_th_full+zonewidth_th/2, dtype=int)+gridsize_th*i,0] = 0.25
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, vsplit_th_full+zonewidth_th/2, dtype=int)+gridsize_th*i,1] = 0.25
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, vsplit_th_full+zonewidth_th/2, dtype=int)+gridsize_th*i,2] = 0.25
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, vsplit_th_full+zonewidth_th/2, dtype=int)+gridsize_th*i,3] = 0.25
                migrouteprobs[np.arange(vsplit_th_full+zonewidth_th/2,gridsize_th, dtype=int)+gridsize_th*i,1] = 0.5        
                migrouteprobs[np.arange(vsplit_th_full+zonewidth_th/2,gridsize_th, dtype=int)+gridsize_th*i,3] = 0.5  
            
            # Bottom section (migration route 3 and 4)
            for i in np.arange(hsplit_th_full+zonewidth_th/2, gridsize_th, dtype=int):
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, dtype=int)+gridsize_th*i,2] = 1
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, vsplit_th_full+zonewidth_th/2, dtype=int)+gridsize_th*i,2] = 0.5
                migrouteprobs[np.arange(vsplit_th_full-zonewidth_th/2, vsplit_th_full+zonewidth_th/2, dtype=int)+gridsize_th*i,3] = 0.5
                migrouteprobs[np.arange(vsplit_th_full+zonewidth_th/2,gridsize_th, dtype=int)+gridsize_th*i,3] = 1  
        
            for i in np.arange(thinnedpoints):
                migroute[i] = rng.choice(np.arange(4), p=migrouteprobs[i,])   
    # migroute - one-dimensional numpy array indicating the migration route for each grid cell
    
    # Removal of edge effects and thinning of indices to save: 
    thinned_indices = (np.ravel_multi_index(np.indices((round(gridsize/thin), round(gridsize/thin)))*thin, (gridsize, gridsize))+gridsize).reshape(thinnedpoints)
    midpoints_row1_full = np.int64(np.arange(gridsize*side*thin+side*thin, gridsize*side*thin+side*thin+midsection))
    midpoints_full = np.int64(np.concatenate([midpoints_row1_full+gridsize*a for a in range(0, midsection)]))
    midpoints_th = midpoints_full[np.isin(midpoints_full, thinned_indices)]

    midpoint_indices_thinned = [i for i, x in
        enumerate(np.isin(thinned_indices, midpoints_th)) if x]
    
    for i in range(nit):
        counter = counter + 1
    
## Setting up nonbreeding environmental conditions and spatial fluctuations in environmental noise: 
        if numbermigroutes >= 1:
            covariance_between = correlation_between * (sigma**2) 
            covariance_within = correlation_within * (sigma**2)            
            cov = np.ndarray([thinnedpoints, thinnedpoints], dtype=float)
            for j in np.arange(numbermigroutes):
                cov[np.ix_(migroute==j,migroute==j)] = covariance_within
                cov[np.ix_(migroute==j,migroute!=j)] = covariance_between
            np.fill_diagonal(cov, sigma**2, wrap=False)
            mean = np.zeros(thinnedpoints, dtype="float64")
            envnoise_nonbreed = rng.multivariate_normal(mean, cov, 1)[0]
        
        else:
            envnoise_nonbreed = env_f(gridsize, thin, distlist, xi_list)
        migfec = envnoise_nonbreed.copy()  
        
          
 ## NIT1 ## Send populations to nonbreeding ground (migrants) or simulate winter conditions on breeding ground (nonmigrants):
        new_abundances = new_abundances*expit(logit(nonbreedsurv) + survDensdep(new_abundances, sdd) + weight*envnoise_nonbreed)       

    
## Populations return to breeding ground, apply dispersal: Take the populations composed of individuals who survived over winter and apply dispersal, breeding ground survival, fecundity, and carryover fecundity to them.

        ## NIT2 ## Dispersal on breeding ground when returning from migration OR before the breeding season (nonmigrants) + create grid of carry over effects of migration on fecundity.
        disp_abundance_migfec = Dispchange(new_abundances, gridsize, thin, thinned_indices, distlist, disprate, lg, numbermigroutes, thinnedpoints, migfec) 
        
        # Draw environmental field for summer season on breeding ground
        env_fieldsummer = env_f(gridsize, thin, distlist, xi_list)

        ## NIT3 ## Survival on breeding ground 
        surv_ab = disp_abundance_migfec[0] * expit(logit(breedsurv) + env_fieldsummer*survenv)
        
        ## NIT4 ## Fecundity on breeding ground
        offsp = (np.exp(log(fec) + env_fieldsummer*fecenv  + carryover*disp_abundance_migfec[1]))*surv_ab   
            
        new_abundances = surv_ab[0] + offsp[0] 
        new_abundances_midpoints = [new_abundances[e] for e in midpoint_indices_thinned] 
   
        if counter == burn_in:
            with open('Output_S%s_F%s_Corr%s_Migroutes%s_Nit%s_Carryover%s_HighCorWithin.csv' %(nonbreedsurv, fec, correlation_between, numbermigroutes, nit, carryover), 'w')  as output:
                output.write(str(new_abundances_midpoints) + '\n')
                
        
        if counter > burn_in:
            with open('Output_S%s_F%s_Corr%s_Migroutes%s_Nit%s_Carryover%s_HighCorrWithin.csv' %(nonbreedsurv, fec, correlation_between, numbermigroutes, nit, carryover), 'a')  as output:
                output.write(str(new_abundances_midpoints) + '\n')
          
         

###############################################################################################################
### RUNNING SIMULATIONS #######################################################################################
###############################################################################################################


## 1) Step 1: Run the distance calculations (just need to do this once for every combo of distances & thin)
PointDist(150, 2)  ##(gridsize, thin), creates a pickle point files

xi_list = xi_dist(150, 9, 0.10) ##(gridsize, le, sigma)

## 2) Step 2: Load saved distances
def distlistload(gridsize, thin):
    with open('PDt%s_thin%s.pickle' %(gridsize, thin), 'rb')  as f: 
        return(pickle.load(f)) 
    
distlist = distlistload(150, 2) ##distlist(gridsize, thin), load the saved pickle distance file
   
    
distances(150, 2, 8, 50) ##(gridsize, thin, zonewidth, midsection), creates distances and coordinates csv files

abundances = np.repeat(20, 75**2) ##(abundance, number of populations)

breed_nonbreed_update(nit=500, abundances=abundances, xi_list=xi_list, distlist=distlist, gridsize=150, thin=2, midsection=50, le=9, sigma=0.10, disprate=0.02, lg=5, breedsurv=1, nonbreedsurv=0.5, survenv=1, fec=2.0, fecenv=1, sdd=0.03, randommigroutes=0, correlation_between=0, correlation_within=0.75, weight=1, numbermigroutes=0, zonewidth=8, vsplit=25, hsplit=25, carryover=1)

