=============================
How to get MAVEN data
=============================

Below are the scripts to update/download MAVEN data (using update.py in scripts/), from inside the mavenpy folder. Replace MAVEN_DATA_DIR with the directory where MAVEN data is stored.

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

* LPW:

  * Density and temperatures: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument lpw``
  * IV curve: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument lpw --dataset_name lpiv``
  * Passive spec: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --remote lasp_sdc_public --instrument lpw --dataset_name specpas``


* SEP:  ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-13 --instrument sep --remote lasp_sdc_public``

* MAG (only works for ssl_sprg remote):  ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2022-07-10 --end_date 2022-07-10 --instrument mag --remote ssl_sprg``

* NGIMS:

  * Level 2 neutral data: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument ngims --remote lasp_sdc_public``
  * Level 2 ion data: ``python scripts/update.py -d MAVEN_DATA_DIR --start_date 2023-07-10 --end_date 2023-07-10 --instrument ngims --remote lasp_sdc_public --dataset_name 'ion-abund-(.*)'``
