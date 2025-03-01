=============================
How to get MAVEN data
=============================

Below are the scripts to update/download MAVEN data (using update.py in scripts/), from inside the mavenpy folder. WARNING! Need to replace 'MAVEN_DATA_DIR' with the directory you want to store MAVEN data. This will also by default save the data according to the SSL SPRG directory scheme. If you don't want this, add ``--no_mirror `` (will save to a folder called 'data' regardless of directory structure on remote) or add ``local_source_tree lasp_sdc_public`` (will save data according to a directory structure based on the LASP sdc).

* Spice: Can download/update MAVEN-related kernels. Defaults to mirroring Spice kernels as organized in the SSL SPRG server, downloading MAVEN position info (SPK, FK) and *not* downloading the MAVEN pointing files (CK, for MAVEN and the APP).

  * Assumes a directory structure based on SSL SPRG: ``python scripts/update.py -d MAVEN_DATA_DIR --spice --start_date 2015-03-15 --n_days 1``
  * If you don't want a directory structure based on SSL SPRG: ``python scripts/update.py -d KERNEL_DIR --spice --start_date 2015-03-15 --n_days 1 --no_mirror``
    * To download MAVEN pointing information, add ``--sc_pointing``.
    * To download MAVEN articulated payload platform (APP) pointing information, add ``--app_pointing``.
    * To ask before downloading new files, add ``--prompt_for_spice_download``.

* EUV:

  * Level 2 (bands): ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument euv --remote lasp_sdc_public``
  * Level 3 (minute, spectra): ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument euv --remote lasp_sdc_public --level l3``

* SWIA:

  * Moments: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument swia --remote lasp_sdc_public --dataset_name onboardsvymom``
  * 3D distributions:

    * Coarse (survey mode): ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument swia --remote lasp_sdc_public --dataset_name coarsesvy3d``
    * Coarse (archive mode): ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument swia --remote lasp_sdc_public --dataset_name coarsearc3d``
    * Fine (survey mode): ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument swia --remote lasp_sdc_public --dataset_name finesvy3d``
    * Fine (archive mode): ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument swia --remote lasp_sdc_public --dataset_name finearc3d``

* SWEA: 

  * Spectra: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument swea --dataset_name svyspec``
  * Survey mode:

    * Pitch angle distribution: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument swea --dataset_name svyspec``
    * 3D distribution: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument swea --dataset_name svy3d``

  * Archive mode:

    * Pitch angle distribution: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument swea --dataset_name arcpad``
    * 3D distribution: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument swea --dataset_name arc3d``

* IUVS: The level 2 files are generally filled with NaNs, so best to download the L1cs.
 * Limb: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2016-05-20 --end_date 2016-05-22 --remote lasp_sdc_public  --instrument iuvs --dataset limb --orbit_segment periapse --level l1c``

* LPW:

  * Density and temperatures: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument lpw``
  * IV curve: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument lpw --dataset_name lpiv``
  * Passive spec: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument lpw --dataset_name specpas``


* SEP:  ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-13 --instrument sep --remote lasp_sdc_public``

* MAG (only works for ssl_sprg remote):  ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2022-07-10 --end_date 2022-07-10 --instrument mag --remote ssl_sprg``

* NGIMS:

  * Level 2 neutral data: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument ngims --remote lasp_sdc_public``
  * Level 2 ion data: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument ngims --remote lasp_sdc_public --dataset_name 'ion-abund-(.*)'``
