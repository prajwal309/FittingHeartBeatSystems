 
    #Use the timing of the primary and secondary eclipse in order to fit for the transits.
    def CreateSyntheticLightCurve():

        import batman
        TStep = 1.0/(60*24)
        Time = np.arange(0,30, TStep)

        # Define the parameters of the star and planet
        params = batman.TransitParams()
        params.t0 = 0.0  # Transit center time
        params.per = 5.0  # Orbital period (days)
        params.rp = 0.1  # Planet radius (in units of the stellar radius)
        params.a = 10.0  # Semi-major axis (in units of stellar radius)
        params.inc = 89.0  # Orbital inclination (in degrees)
        params.ecc = 0.4  # Eccentricity
        params.w = 90.0  # Argument of periastron (in degrees)
        params.u = [0.1, 0.3]  # Limb darkening coefficients (linear and quadratic)
        params.limb_dark = "quadratic"  # Limb darkening model

        # Create a BATMAN model
        m = batman.TransitModel(params, Time)

        # Compute the eccentric planet light curve
        Flux = m.light_curve(params)

        # Compute the eccentric planet secondary transit light curve
        params.fp = 0.011
        params.t_secondary = 1.5

        m = batman.TransitModel(params, Time, transittype="secondary")
        FluxSecondary = m.light_curve(params) - 0.011


        # Plot both the primary and secondary transit light curves
        plt.figure(figsize=(10, 6))
        plt.plot(Time, Flux, 'b-', label='Primary Transit', lw=2)
        plt.plot(Time, FluxSecondary, 'r--', label='Secondary Transit', lw=2)
        plt.xlabel('Time (days)')
        plt.ylabel('Relative Flux')
        plt.title('Eccentric Planet Primary and Secondary Transit Light Curves')
        plt.legend()
        plt.grid(True)
        plt.show()

        
        
        Flux+= np.random.normal(len(Time))*1e-3

        return Time, Flux, params
