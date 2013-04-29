.. _parameters:

Parameters and Data Files
=========================

Parameters
----------

For all on/off integer flags, 0 is off and 1 is on.

``use_chemistry`` (int)
    Flag to activate the grackle machinery.  Default: 0.

``with_radiative_cooling`` (int)
    Flag to include radiative cooling and actually update the thermal energy during the chemistry solver.  If off, the chemistry species will still be updated.  The most common reason to set this to off is to iterate the chemistry network to an equilibrium state.  Default: 1.

``primordial_chemistry`` (int)
    Flag to control which primordial chemistry network is used.  Default: 0.

    - 0: no chemistry network.
    - 1: 6-species atomic H and He.  Active species: H, H\ :sup:`+`, He, He\ :sup:`+`, \ :sup:`++`, e\ :sup:`-`.
    - 2: 9-species network including atomic species above and species for molecular hydrogen formation.  This network includes formation from the H\ :sup:`-` and H\ :sub:`2`\ :sup:`+` channels, three-body formation (H+H+H and H+H+H\ :sub:`2`), H\ :sub:`2` rotational transitions, chemical heating, and collision-induced emission (optional).  Active species: above + H\ :sup:`-`, H\ :sub:`2`, H\ :sub:`2`\ :sup:`+`.
    - 3: 12-species network include all above plus HD rotation cooling.  Active species: above plus D, D\ :sup:`+`, HD.

``h2_on_dust`` (int)
    - Flag to enable H\ :sub:`2` formation on dust grains, dust cooling, and dust-gas heat transfer follow `Omukai (2000) <http://adsabs.harvard.edu/abs/2000ApJ...534..809O>`_.  This assumes that the dust to gas ratio scales with the metallicity.  **This is not extensively tested and should probably not be used for now.** Default: 0.

``metal_cooling`` (int)
    Flag to enable metal cooling using the Cloudy tables.  If enabled, the cooling table to be used must be specified with the ``grackle_data_file`` parameter.  Default: 0.

``cmb_temperature_floor`` (int)
    Flag to enable an effective CMB temperature floor.  This is implemented by subtracting the value of the cooling rate at T\ :sub:`CMB` from the total cooling rate.  Default: 1.

``UVbackground`` (int)
    Flag to enable a UV background.  If enabled, the cooling table to be used must be specified with the ``grackle_data_file`` parameter.  Default: 0.

``include_metal_heating`` (int)
    Flag to include heating of the metals by the UV background.  Default: 0.

``grackle_data_file`` (string)
    Path to the data file containing the metal cooling and UV background tables.  Default: "".

``Gamma`` (float)
    The ratio of specific heats for an ideal gas.  A direct calculation for the molecular component is used if ``primordial_chemistry`` > 1.  Default:  5/3.

``three_body_rate`` (int)
    Flag to control which three-body H\ :sub:`2` formation rate is used.  0: `Abel, Bryan & Norman (2002) <http://adsabs.harvard.edu/abs/2002Sci...295...93A>`_, 1: `Palla, Salpeter & Stahler (1983) <http://adsabs.harvard.edu/abs/1983ApJ...271..632P>`_, 2: `Cohen & Westberg (1983) <http://adsabs.harvard.edu/abs/1983JPCRD..12..531C>`_, 3: `Flower & Harris (2007) <http://adsabs.harvard.edu/abs/2007MNRAS.377..705F>`_, 4: `Glover (2008) <http://adsabs.harvard.edu/abs/2008AIPC..990...25G>`_.  These are discussed in `Turk et. al. (2011) <http://adsabs.harvard.edu/abs/2011ApJ...726...55T>`_.  Default: 0.

``cie_cooling`` (int)
    Flag to enable H\ :sub:`2` collision-induced emission cooling from `Ripamonti & Abel (2004) <http://adsabs.harvard.edu/abs/2004MNRAS.348.1019R>`_.  Default: 0.

``h2_optical_depth_approximation`` (int)
    Flag to enable H\ :sub:`2` cooling attenuation from `Ripamonti & Abel (2004) <http://adsabs.harvard.edu/abs/2004MNRAS.348.1019R>`_.  Default: 0.

``photoelectric_heating`` (int)
    Flag to enable a spatially uniform heating term approximating photo-electric heating from dust from `Tasker & Bryan (2008) <http://adsabs.harvard.edu/abs/2008ApJ...673..810T>`_.  Default: 0.

``photoelectric_heating_rate`` (float)
    If ``photoelectric_heating`` enabled, the heating rate in units of erg cm\ :sup:`-3` s\ :sup:`-1`.  Default: 8.5e-26.

``Compton_xray_heating`` (int)
   Flag to enable Compton heating from an X-ray background following `Madau & Efstathiou (1999) <http://adsabs.harvard.edu/abs/1999ApJ...517L...9M>`_.  Default: 0.

``LWbackground_intensity`` (float)
    Intensity of a constant Lyman-Werner H\ :sub:`2` photo-dissociating radiation field in units of 10\ :sup:`-21` erg s\ :sup:`-1` cm\ :sup:`-2` Hz\ :sup:`-1` sr\ :sup:`-1`.  Default: 0.

``LWbackground_sawtooth_suppression`` (int)
    Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption (giving a sawtooth pattern), taken from `Haiman & Abel, & Rees (2000) <http://adsabs.harvard.edu/abs/2000ApJ...534...11H>`_.  Default: 0.

Data Files
----------

These files contain the metal heating and cooling rates and the UV background photo-heating and photo-ionization rates.  For all three files, the number density range is -10 < log\ :sub:`10` (n\ :sub:`H` / cm\ :sup:`-3`) < 4 and the temperature range is 1 < log\ :sub:`10` (T / K) < 9.  Extrapolation is performed when outside of the data range.  All data files are located in the **input** directory in the source.

 - **CloudyData_noUVB.h5** - metal cooling rates for collisional ionization equilibrium.

 - **CloudyData_UVB=FG2011.h5** - metal heating and cooling rates and UV background rates from the work of `Faucher-Giguere et. al. (2009) <http://adsabs.harvard.edu/abs/2009ApJ...703.1416F>`_, updated in 2011.  The maxmimum redshift is 10.6.  Above that, collisional ionization equilibrium is assumed.

 - **CloudyData_UVB=HM2012.h5** - metal heating and cooling rates and UV background rates from the work of `Haardt & Madau (2012) <http://adsabs.harvard.edu/abs/2012ApJ...746..125H>`_.  The maximum redshift is 15.13.  Above that, collisional ionization equilibrium is assumed.