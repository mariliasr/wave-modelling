def spec2d(Hs, Tp, pdir, ms=2, units='deg', normalize=True):
    
    import numpy as np
    from oceanwaves.utils import trapz_and_repeat

    """ This function returns the 2d-spectrum from Hs, Tp, pdir and spread

    -----------------------------------------------------

    REFERENCE:
    - https://www.orcina.com/webhelp/OrcaFlex/Content/html/Waves,Wavespectra.htm
    - http://research.dnv.com/hci/ocean/bk/c/a28/s3.htm
    - https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/applications/waves/

    -----------------------------------------------------

    Parameters: 

    Hs = []; # Sig. higth 
    Tp = []; # Peak Period  
    pdir = [];  # Peak Direction  
    ms = []; # Cossene power. Default value is 2.
    units = [deg] or [rad]; # Direction unit """

    # Setup parameters:

    freq = np.array([0.0500, 0.0566, 0.0642, 0.0727, 0.0824, 0.0933, 0.1057, 0.1198,
                     0.1357, 0.1538, 0.1742, 0.1974, 0.2236, 0.2533, 0.2870, 0.3252,
                     0.3684, 0.4174, 0.4729, 0.5357, 0.6070, 0.6877, 0.7791, 0.8827,
                     1.0000]) # frequncy bin in swan

    directions =  np.arange(0., 360., 10) # Directional bin swan
    
    if Tp == 0:
        spec_2d = np.zeros((len(freq), len(directions)))
    else:
        fp = 1/Tp # peak frequency

        gam = 3.3 # default value according to SWAN manual

        g = 9.81 # gravity
        
        
        sigma = freq * 0
        sigma[freq<fp] = 0.07
        sigma[freq>=fp]= 0.09
        sigma = np.array(sigma)

        # Pierson-Moskowitz Spectrum

        alpha = 1 / (0.06533 * gam ** 0.8015 + 0.13467)/16; # Yamaguchi (1984), used in SWAN

        pm = alpha * Hs ** 2 * Tp **-4 * freq ** -5 * np.exp(-1.25 * (Tp * freq)**-4) 

        # apply JONSWAP shape

        jon = pm * gam ** np.exp(-0.5 * (Tp * freq - 1) ** 2. / (sigma ** 2.))

        jon[np.isnan(jon)] = 0

        # Optionally correct total energy of user-discretized spectrum to match Hm0, 
        # as above methods are only an approximation

        eps = np.finfo(float).eps

        if normalize is True:
            corr = Hs ** 2/(16*trapz_and_repeat(jon, freq))
            jon = jon*corr


        # Directional Spreading

        '''Generate wave spreading'''

        from math import gamma

        directions = np.asarray(directions, dtype=np.float)

        # convert units to radians
        if units.lower().startswith('deg'):
            directions = np.radians(directions)
            pdir = np.radians(pdir)
        elif units.lower().startswith('rad'):
            pass
        else:
            raise ValueError('Unknown units: %s')

        # compute directional spreading
        A1 = (2.**ms) * (gamma(ms / 2 + 1))**2. / (np.pi * gamma(ms + 1))
        cdir = A1 * np.maximum(0., np.cos(directions - pdir))
        #cdir = np.maximum(0., np.cos(directions - pdir))**s

        # convert to original units
        if units.lower().startswith('deg'):
            directions = np.degrees(directions)
            pdir = np.degrees(pdir)
            cdir = np.degrees(cdir)

        # normalize directional spreading
        if normalize:
            cdir /= trapz_and_repeat(cdir, directions - pdir, axis=-1)

        cdir = np.array([list(cdir)] * len(jon))

        jon_list = list(jon)
        
        jon = np.array([ele for ele in jon_list for i in range(len(directions))]
                      ).reshape(len(freq), len(directions))

        jon2 = jon * cdir
        jon2[np.isnan(jon2)] = 0
        spec_2d = jon2
    
    return spec_2d