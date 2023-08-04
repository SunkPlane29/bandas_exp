A NICER View of PSR J0740+6620:
    Nested Samples for Millisecond Pulsar Parameter Estimation

Authors:
    Riley T.E., Watts A.L., Ray P.S., Bogdanov S., Guillot S., Morsink S.M.,
    Bilous A.V., Arzoumanian Z., C. Devarshi, Deneva J.S., Gendreau K.C.,
    Harding A.K., Ho W.C.G., Lattimer J.M., Loewenstein M., Ludlam R.M.,
    Markwardt C.B., Okajima T., Prescod-Weinstein C., Remillard R.A.,
    Wolff M.T., Fonseca E., Cromartie H.T., Kerr M., Pennucci T.T.,
    Parthasarathy A., Ransom S., Stairs I., Guillemot L., Cognard I.

Accompanying article:
    Riley et al. 2021

This repository is also associated with the articles:
    * Fonseca et al. 2021
    * Raaijmakers et al. 2021
    * Wolff et al. 2021


Deposit contents
================

Tree
----

Here is a directory tree of important files with short descriptions.

* NICER_x_XMM_J0740__XPSI_v0p7_STU__NSX_FIH__NANOGrav_x_CHIME.ipynb
    = Python analysis notebook for production analysis of the ST-U model
* NICER_x_XMM_J0740__XPSI_v0p7_STU__NSX_FIH__NANOGrav_x_CHIME__skymaps.ipynb
    = Python analysis notebook for production analysis ST-U model skymaps
* NICER_x_XMM_J0740__XPSI_v0p7_STU__NSX_FIH__NANOGrav_x_CHIME__exploratory.ipynb
    = Python analysis notebook for exploratory analysis of the ST-U model
* NICER_x_XMM_J0740__XPSI_v0p7_ST_PST__NSX_FIH__NANOGrav_x_CHIME.ipynb
    = Python analysis notebook for exploratory analysis of ST+PST model
* J0740_products
    > NICER
        - nixtiref20170601v002_rmf_full.txt
            = RMF
        - rmf_nicer_channel_energies.txt
            = channel nominal energy bounds
        - nixtiaveonaxis20170601v004_offaxis_d51_arf.txt
            = ARF (on-axis is column three)
        - PSRJ0740p6620.Aug172020.Ztest.Det34Out.GTIOptimalEvents.txt
            = J0740+6620 XTI event list
    > XMM
        - pn
            + j0740_pn_spec_src_evt_obs123_0.2-15keV.txt
                = J0740+6620 XMM EPIC pn event list
            + j0740_pn_arf.txt
                = ARF
            + pntffg_blanksky_spec.txt
                = blank-sky spectrum
            + pn_rmf_full.txt
                = RMF
            + pn_energy_bounds.txt
                = channel nominal energy bounds
        - MOS1
            + j0740_mos1_spec_src_evt_obs123_0.2-12keV.txt
                = J0740+6620 XMM EPIC mos1 event list
            + j0740_mos1_arf.txt
                = ARF
            + m1tffg_blanksky_spec.txt
                = blank-sky spectrum
            + mos1_rmf_full.txt
                = RMF
            + mos1_energy_bounds.txt
                = channel nominal energy bounds
        - MOS2
            + j0740_mos2_spec_src_evt_obs123_0.2-12keV.txt
                = J0740+6620 XMM EPIC MOS1 event list
            + j0740_mos2_arf.txt
                = ARF
            + m2tffg_blanksky_spec.txt
                = blank-sky spectrum
            + mos2_rmf_full.txt
                = RMF
            + mos2_energy_bounds.txt
                = channel nominal energy bounds
* model_data
    > nsx_H_v200804.out
        = NSX fully-ionized non-magnetic hydrogen (FI-H) atmosphere emergent
          intensity table
    > README_v200804.txt
        = README for atmosphere table
    > interstellar_phot_frac.txt
        = interstellar attenuation fraction table at a fiducial column density
* STU
    = event data <- X-PSI ST-U model
    > CustomPhotosphere.py
        = module for photon specific intensity emergent from hot region with
          fully-ionized hydrogen atmosphere
    > CustomPhotosphere_He.py
        = module for photon specific intensity emergent from hot region with
          fully-ionized helium atmosphere
    > CustomInterstellar.py
        = module for interstellar attenuation model
    > CustomInstrument.py
        = module for instrument models (NICER XTI and XMM-Newton EPIC cameras)
    > CustomSignal.py
        = module for model per-signal likelihood evaluation
    > CustomPrior.py
        = module for exploratory analysis prior with deprecated NANOGrav x CHIME
          radio-timing information
    > CustomPrior_diffuse.py
        = module for exploratory analysis prior with diffuse mass, distance, and
          inclination prior information
    > CustomPrior_Cromartie.py
        = module for exploratory analysis prior with deprecated Cromartie+2019
          radio-timing information
    > CustomPrior_GLS.py
        = module for production analysis prior with updated NANOGrav x CHIME
          radio-timing information
    > CustomPrior_GLS_compressed_scaling.py
        = module for production analysis prior with updated NANOGrav x CHIME
          radio-timing information and a compressed effective area scaling prior
    > NICERxXMM
        - FI_H
            = NICER & XMM event data <- ST-U with fully-ionized hydrogen
              atmosphere (FI-H)
            + <exploratory runs>
            + run10
                = headline posterior conditional on NICER & XMM event data
                ++ NICER_x_XMM_J0740__XPSI_STU_NSX_FIH__mass_radius.txt
                    = headline mass-radius posterior samples (equally-weighted,
                      using perplexity to calculate effective sample size and
                      thus the number of samples to draw from the weighted set)
                ++ NICER_x_XMM_J0740__XPSI_STU_NSX_FIH__weight_mass_radius.txt
                    = headline mass-radius posterior samples (weighted)
                ++ NICER_x_XMM_J0740__XPSI_STU_NSX_FIH__contour_radius_mass_68.txt
                    = headline mass-radius posterior 68% iso-density contour
                ++ NICER_x_XMM_J0740__XPSI_STU_NSX_FIH__contour_radius_mass_95.txt
                    = headline mass-radius posterior 95% iso-density contour
                ++ NICER_x_XMM__XPSI_posterior_PDF__columns__radius_km__probability_density.txt
                    = headline marginal posterior PDF of radius
                ++ NICER_x_XMM__EAPC__XPSI_posterior_PDF__columns__radius_km__probability_density.txt
                    = headline marginal posterior PDF of radius with effective
                      area scaling factor prior compression
                ++ main.py
                    = model script that combines model modules and calls sampler
                ++ config.ini
                    = configuration file for the main script
                ++ init
                    = start of nested sampling process
                    +++ job.sh
                        = job script for Slurm batch scheduler
                    +++ out
                        = standard output file (X-PSI and MultiNest output)
                    +++ err
                        = standard error file (usually empty)
                ++ resume1
                    = resume nested sampling process (if necessary)
                    +++ job.sh
                        = slightly different job script to handle samples
                ++ samples
                    +++ <sample files>
        - FI_He
            = NICER & XMM event data <- ST-U with fully-ionized helium
              atmosphere (FI-He)
            + <exploratory runs>
    > main_NxX.py
        = main script/module for NICER x XMM <- FI-H headline posterior
    > main_NxX_IS.py
        = main script/module for NICER x XMM <- FI-H posterior with compressed
          effective area scaling prior, called by importance sampler
    > NICER
        - FI_H
            = NICER event data <- ST-U with fully-ionized hydrogen atmosphere
            + <exploratory runs>
            + run12
                = headline posterior conditional on NICER event data
                ++ NICER_J0740__XPSI_STU_NSX_FIH__mass_radius.txt
                    = headline mass-radius posterior samples (equally-weighted,
                      using perplexity to calculate effective sample size and
                      thus the number of samples to draw from the weighted set)
                ++ NICER_J0740__XPSI_STU_NSX_FIH__weight_mass_radius.txt
                    = headline mass-radius posterior samples (weighted)
                ++ NICER_J0740__XPSI_STU_NSX_FIH__contour_radius_mass_68.txt
                    = headline mass-radius posterior 68% iso-density contour
                ++ NICER_J0740__XPSI_STU_NSX_FIH__contour_radius_mass_95.txt
                    = headline mass-radius posterior 95% iso-density contour
                ++ NICER__XPSI_posterior_PDF__columns__radius_km__probability_density.txt
                    = headline marginal posterior PDF of radius
                ++ <run files>
        - FI_He
            = NICER event data <- ST-U with fully-ionized helium atmosphere
            + <exploratory runs>
    > main_NICER.py
        = main script/module for NICER <- FI-H posterior with deprecated
          NANOGrav x CHIME radio-timing prior
    > main_NICER_IS.py
        = main script/module for NICER <- FI-H headline posterior (with updated
          NANOGrav x CHIME radio-timing prior), called by importance sampler
    > XMM
        - FI_H
            = XMM event data <- ST-U with fully-ionized hydrogen atmosphere
            + <exploratory runs>
            + run3
                = headline posterior conditional on XMM event data
                ++ XMM__XPSI_posterior_PDF__columns__radius_km__probability_density.txt
                    = headline marginal posterior PDF of radius
                ++ <run files>
        - FI_He
            = XMM event data <- ST-U with fully-ionized helium atmosphere
            + <exploratory runs>
    > main_XMM.py
        = main script/module for XMM <- FI-H headline posterior
    > <main scripts/modules for exploratory analysis>
* ST_PST
    = event data <- X-PSI ST+PST model; similar directory structure to STU
* importance_sampling
    > NICER_x_XMM_importance_sample.py
        = instrument effective area prior compression
    > NICER_importance_sample.py
        = update from deprecated NANOGrav x CHIME radio timing prior


Sample files
------------

In each model directory there are one or more run directories.
In each run directory there are a set of files generated by MultiNest v3.12.
The files contain information about one or more nested sampling processes.

The file names contain conventional identifiers appended by MultiNest.
Other information can be found in the filename too.
As an example, consider the ST-U model path:
    STU/NICERxXMM/FI_H/run10/nlive40000_eff0.1_noCONST_noMM_noIS_tol-1.txt

The separated elements here are:
    * the run ID "run10"
    * the number of live points "nlive"=40000
    * the efficiency setting "eff"=0.1
    * the state of constant efficiency variant "CONST"=False
    * the state of mode-separation variant "MM"=False
    * the state of the importance nested sampling augmentation "IS"=False
    * the termination tolerance condition "tol"=10^{-1}

Other files generated by MultiNest on path "STU/NICERxXMM/FI_H/run10/" have the
same prefix "nlive40000_eff0.1_noCONST_noMM_noIS_tol-1", but have been
automatically appended with identifiers and file extensions:
    * .txt = weighted nested sample file
    * ev.dat = accumulated dead points (parameter vectors) ordered in likelihood
    * resume.dat = information for run resume, written to disk at higher cadence
    * stats.dat = estimators (mode breakdown if "MM" variant activated)
    * summary.txt = largely similar to stats.dat but not as human readable
    * post_equal_weights.dat = equally weighted posterior samples (more MC noise)
    * live.points = current live points in the unit hypercube (sampling space)
    * phys_live.points = current live points in the parameter space
    * phys_live-birth.txt = birth iteration numbers of current live points
    * dead-birth.txt = birth iteration numbers of dead points
    * __importance_sampled = importance sampled posterior using X-PSI tool
    * files with "transformed" in the name were generated during post-processing

The "birth" files are generated by MultiNest >= v3.11, and are needed for
thread decomposition and thus Monte Carlo error analysis via bootstrapping, and
for run combination.

Note that importance nested sampling is not compatible for evidence estimation
because of the way we typically implement the joint prior density. We do not
activate it for any inferences in the journal article, but for headline
posteriors at high resolution, errors are thought to be small. Moreover, the
augmented sampling process writes to disk all the rejected points too, leading
to additional information albeit requiring much larger disk space of O(10) GBs.
These files are archived but some of the other files refer to additional
estimators based on importance nested sampling.

Many of these files contain vectors as a subset of a row, with other elements
of a row being associated with a given vector. The vectors are in some cases
parameter vectors in a physical space, whilst in other cases they are points
in a unit hypercube. The vectors are ordered, and a vector in the unit
hypercube maps to a vector in the physical space. The ordering is defined
by the likelihood and prior callbacks. To retrieve the ordering for
each model, you may consult the custom prior subclass for each model (see
the "Parameters of common interest" and "Model & job scripts" section below),
wherein the parameter vector is given in the class docstring.

Information about these sampling process settings can be located in the
accompanying Letter and in other literature, first and foremost:
    * README @ https://github.com/farhanferoz/MultiNest
    * Feroz+ (2009) for MultiNest
    * Feroz+ (2013) for importance nested sampling with MultiNest
    * Higson+ (2018a,b) for error analysis and run combination by
        thread decomposition, for diagnostic tests, and for open-source
        software "nestcheck" to facilitate this post-processing
    * Riley (2019), A NICER perspective of neutron star parameter estimation
Other recommended reading includes:
    * Skilling (2006) for seminal work on nested sampling
    * Handley+ (2015) for digestible nested sampling and prior implementation
        theory, and for the PolyChord nested sampling algorithm variant
    * Speagle (2019) for a near comprehensive review of nested sampling
        together with a pure Python implementation "dynesty" with extensive API
        documentation, tutorials, and theory
    * Handley (2019) for open-source software "anesthetic" was not used in
        this work but can be a useful tool
See the end of this file for references.

Typically in each model model directory there is an "output" directory that
contains verbose output generated by MultiNest during each run. This
information is closely monitored, along with a subset of the above files,
during a run to gain understanding of how the process is behaving.
(This is vague, but we must refer to the literature to learn more about the
nested sampling processes and the MultiNest implementation.)


Parameters of interest
----------------------

In each likelihood function we have separated out the headline samples in the
joint space of mass and radius (M-R samples) for easier usage. These samples
are weighted. For example, in the file
    STU/NICERxXMM/FI_H/NICER_J0740__XPSI_STU_NSX_FIH__weight_mass_radius.txt,
each row contains a sample, and the columns are as follows:
    * importance weight
    * -2 . log(likelihood)
    * mass in solar masses ("total" or "gravitational")
    * equatorial radius in km (given as Schwarzschild radial coordinate)
If you want equally-weighted samples you can generate them manually from these
weighted samples, or extract them from a file such as
    STU/NICERxXMM/FI_H/run10/samples/nlive40000_eff0.1_noCONST_noMM_noIS_tol-1post_equal_weights.dat
Note that in the raw MultiNest files containing parameter vectors (i.e.,
samples), the mass is the column corresponding to the second parameter, and the
radius is the column corresponding to the third parameter.
For all models, the order of the parameter vector is:
    [mass, radius, distance, Earth inclination to rotation axis, ...],
and the parameter vector maps to columns contiguously. So once you identify the
column at which the parameter vector begins (which varies between MultiNest
file types), you should be able to extract the mass and radius, and other
parameters of interest, manually.

From the article Raaijmakers et al. (2021; see top of this README), we have
contour files available primarily to ease M-R posterior plotting. On paths:
    * S_U/NICERxXMM/FI_H/run10/
        > NICER_x_XMM_J0740__XPSI_STU_NSX_FIH__contour_radius_mass_68.txt
        > NICER_x_XMM_J0740__XPSI_STU_NSX_FIH__contour_radius_mass_95.txt
    * STU/NICER/FI_H/run12/
        > NICER_J0740__XPSI_STU_NSX_FIH__contour_radius_mass_68.txt
        > NICER_J0740__XPSI_STU_NSX_FIH__contour_radius_mass_95.txt
the iso-density contours enclosing 68% and 95% of the posterior mass,
respectively, are reported as sequences of points. The first column is the
equatorial radius in km, whilst the second column is the mass in solar masses.

The geometric configuration of surface hot regions is more involved to extract
and work with, often requiring spherical trigonometric transformations. We have
attempted to transform to a geometric parameter space that is straightforward
to understand: this space is not the most straightforward for sampling
processes, and if it were used, would result in severe efficiency loss and,
ultimately, higher posterior computation error.

The space can be characterized as follows:
    * for a given model there are two physical surface hot regions;
    * for all but one model the hot regions have their own location parameters;
    * for each such hot region the spherical coordinate colatitude of the
        center of the member regions is given, relative to the stellar
        rotation axis, together with the angular radii;
    * the azimuth of the hot regions is defined as the initial phase, in the
        sense of rotation, of the center of a member region relative to a
        reference meridian between rotational poles;
    * for a primary hot region the meridian is aligned with the Earth direction;
    * for a secondary hot region the meridian is aligned with the negative
        Earth direction from the stellar center;
    * if the hot region is constructed from two member regions, the phase
        specifies the azimuth of the center of the superseding member (which
        may or may not radiate);
    * finally a convention is that
        > for a radiating superseding member, the azimuth of the ceding region
            center is given w.r.t the superseding member center;
        > for a non-radiating superseding member, the azimuth of its center is
            given relative to the ceding member center (thus the negative of
            azimuth of ceding member center w.r.t. superseding member center).

Consider the following example for ST+PST:
    * the primary is designated as ST;
    * inspect "CustomPrior.__doc__" on path ST_PST/CustomPrior.py;
    * vector elements "p[4]", "p[5]", and "p[6]" are respectively
        > the ST initial center phase in cycles w.r.t Earth meridian;
        > the ST center colatitude in radians
        > the ST angular radius in radians
    * on path ST_PST/NICER/FI_H/run1/samples/nlive2000_eff0.1_noCONST_noMM_noIS_tol-1.txt,
        columns seven, eight, and nine correspond to these parameters;
    * vector elements "p[8]", "p[9]", and "p[10]" are respectively
        > the PST ceding member initial phase in cycles w.r.t to the opposite
            meridian to the Earth
        > the colatitude of the center of the PST ceding member in radians
        > the PST ceding member angular radius in radians;
    * on path ST_PST/NICER/FI_H/run1/samples/nlive2000_eff0.1_noCONST_noMM_noIS_tol-1.txt,
        columns eleven, twelve, and thirteen ccorrespond to these parameters;
    * vector elements "p[11]", "p[12]", and "p[13]" are respectively
        > the center colatitude of the PST superseding member region in radians
            (the "omit" in the docstring, which is a term for a non-radiating or
            zero temperature member that masks a ceding member)
        > the PST superseding member angular radius in radians
        > the PST superseding member azimuth - in radians and in the sense of
            rotation - relative to the PST ceding member center
    * on path ST_PST/NICER/FI_H/run1/samples/nlive2000_eff0.1_noCONST_noMM_noIS_tol-1.txt,
        columns forteen, fifteen, and sixteen correspond to these parameters.


Model & job scripts
-------------------

The likelihood and prior were implemented as callbacks using v0.7 of X-PSI.
X-PSI is hosted @ https://github.com/ThomasEdwardRiley/xpsi
Documentation is hosted at https://thomasedwardriley.github.io/xpsi/
If you need the specific version v0.7, you can checkout the tag on the master
branch of the repo.

In each hot-region model directory (e.g., STU) there is a set of Python modules
together with an empty "__init__.py".
These files contain the custom subclasses that constitute the model.
All but one of the files are associated with the likelihood function.
The remaining file defines the joint prior distribution together with
the transformation from the unit hypercube to a physical space for likelihood
evaluation.

In each hot-region model directory there are "main" Python scripts/modules that
are imported into an analysis notebook, e.g., for post-processing or model
setup checks. These scripts import the custom subclasses and perform
initialisation of the likelihood callable and prior callable. Both callables
are then passed as callbacks to MultiNest via PyMultiNest
(https://github.com/JohannesBuchner/PyMultiNest). These main files correspond
to the main files in the run subdirectories, with slight changes to enable
importing various model variants into an analysis notebook. The main files in
the run subdirectories are passed to an MPI executable such as mpiexec.

In each run directory there are one or more simple job bash scripts used with
the Slurm batch scheduler on either the SURFsara Cartesius supercomputer or
the Toulouse CALMIP cluster. Such a script can in principle be translated
straightforwardly to other systems.

For runs that were resumed, there are separate job scripts for the initial stage
of the process ("init" directory) and each resume (e.g., "resume1").


Summary
=======

The information archived here should in principle permit reproduction of
all or parts of the work. All files necessary to reproduce headline posteriors
with X-PSI are available in this repository. The NSX fully-ionized helium
atmosphere table used for exploratory analyses is not included in this
repository, and access requests need to be directed to Wynn C. G. Ho.

In practice we appreciate that the model scripts are not pristine and verbose,
and the API defined by v0.7 of X-PSI is still under development. Improvements
will continue to be made actively on GitHub.

Finally, if you think there are missing files that you need, you can contact
us using the details below.


Contact
=======

If you have questions, comments, or requests, please contact:
Thomas E. Riley, t.riley.phd@gmail.com
(in cc) Anna L. Watts, A.L.Watts@uva.nl
(in cc) Paul Ray, paul.ray@nrl.navy.mil
(in cc) Slavko Bogdanov, slavko@astro.columbia.edu

These contacts are the first four authors of the accompanying journal article.
This should better ensure awareness of problems and in principle a faster
response. You are of course free to contact other authors if you are unsure
about where to send your message, but you may be redirected.

Once we have established the nature of the message we can determine who is
best positioned to respond and continue correspondence as appropriate.


References
==========

Skilling J. 2006, doi:10.1214/06-BA127
Feroz, F., Hobson, M. P., & Bridges, M. 2009, MNRAS, 398, 1601
Feroz, F., Hobson, M. P., Cameron, E., & Pettitt, A. N. 2013,
    arXiv e-prints, arXiv:1306.2144
Handley, W. J., Hobson, M. P., & Lasenby, A. N. 2015, MNRAS, 453, 4384
Higson, E., Handley, W., Hobson, M., & Lasenby, A. 2018a,
    Statistics and Computing, doi:10.1007/s11222-018-9844-0
—. 2018b, Bayesian Analysis, 13, 873
Speagle, J. S. 2019, arXiv e-prints, arXiv:1904.02180
Handley, W. J. 2019, doi:10.21105/joss.01414

