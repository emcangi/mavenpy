==========================
MAVENPy
==========================

Python routines to download, read, and interact with MAVEN spacecraft data. Routines are present for following data products:
    - Spice: Downloading from NAIF, loading via SpiceyPy
    - ephemeris: spacecraft ephemeris (from SSL SPRG) and orbit ephemeris (from Spice) 
    - EUV: Level-3 FISM-fitted spectra, Level-2 bands
    - MAG: Level-2 and 1 magnetic fields (all cadences) (all file formats)
    - NGIMS: Level-2 neutral and ion densities (Level-3 disrecommended)
    - SEP: Level-1 (raw counts), Level-2 (counts and processed spectra), Level-3 (ancillary data, PADs)
    - SWEA: Level-2 (spectra, PAD, 3d), Level-3 (spacecraft potential, topology)
    - SWIA: Level-2 (coarse/fine/onboardsvyspec) angle-averaged spectra (CDF), Level-2 (coarse/fine/onboardsvyspec) differential spectra (CDF), Level-2 moments (CDF).

Depends on Python3.11 and above.

This is a package heavily based on IDL programs in the Berkeley SSL tplot SPEDAS distribution that read and analyze MAVEN spacecraft data. We thank the original authors of those routines, and have thanked them and referenced the source IDL program where derived.


Development
=============================

This package is still in development. Follow these instructions to install:

0. From terminal or conda, use git to clone the repository to local:

    ``git clone <>``

1. Install the package:
    
    ``python -m pip install -e .``

2. Run ``demo_solarwind.py`` to confirm your install is working

    ``python demo_solarwind.py -d <<MAVEN_CONTAINING_DIR_HERE>> --start_date 2015-03-01 --n_days 1``


Running
==========

To access routines in a file:

    import mavenpy

For example, to download L2 30sec MAG data in payload-coordinates to a directory containing MAVEN data (data_directory) in a given date range (between start_date and end_date):

    ```
    from mavenpy import retrieve

    retrieve.sdc_retrieve(
        'mag', destination_dir=data_directory,
        ext='sav', res='30sec', level='l2', coord='pl',
        start_date=start_date, end_date=end_date)
    ```

To retrieve the file names of those MAG data based on the SSL_SPRG directory tree:

    ```
    from mavenpy import file_path

    mag_file_names = file_path.local_file_names(
        data_directory, 'mag', start_date=start_date, end_date=end_date,
        ext='sav', res='30sec', level='l2', coord='pl', source='ssl_sprg')
    ```

To load the data given the file names:

    ```
    from mavenpy import load

    mag = load.load_data(
        mag_file_names, ext='sav', res='30sec', level='l2', coord='pl')
    ```


To rotate the data into MSO coordinates from payload coordinates:

    ```
    from mavenpy import spice

    mag_epoch = mag["epoch"][0]
    bx = mag["Bx"][0]
    by = mag["By"][0]
    bz = mag["Bz"][0]
    b_mso = spice.bpl_to_bmso(mag_epoch, bx, by, bz)
    bx, by, bz = b_mso[:, 0], b_mso[:, 1], b_mso[:, 2]
    ```



