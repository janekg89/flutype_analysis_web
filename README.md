# Analysis scripts for FluType
## Overview

The FluTypeDB project is a web application for the data management of binding assays 
for the classification of influenza subtypes.
 
FluTypeDB consists of a database and web interface with focus on various binding assays 
for the classification of influenza viruses and contains experimental data based on

* ELISA antibody binding assays
* Peptid binding assays on microwell plates
* Peptid bindings assays on peptide microarrays

FluTypeDB is developed for the data management and data analysis within the FluType project
by <b>Janek Grzegorzewski</b> (Universität Potsdam) and
<b><a href="https://livermetabolism.com" target="_blank">Matthias König</a></b> (Humboldt Universität Berlin).

This repository contains the pipelines for pre-processing and analysis.


### License
* Source Code: [LGPLv3](http://opensource.org/licenses/LGPL-3.0)
* Documentation: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)
* Data: Copyright and ownership FluType.

# Technical Documentation
In this section technical information for setup is provided. For most analysis a local version of
FlutypeDB (https://github.com/janekg89/flutype_analysis_web) is required to run.

## Setup
```
...
(flutype_webapp) python manage.py runserver
```
Necessary to create core database with


To fill the database with test data run the following script.
This applies all migrations and writes the database content.
```
(flutype_webapp) ./create_db.sh
```

Create repository
```
git clone https://github.com/janekg89/flutype_webapp_analysis.git
mkvirtualenv flutype_webapp_analysis --python=python3
(flutype_webapp_analysis) pip install -r requirements.txt
(flutype_webapp_analysis) pip install -e .
cd path/to/flutype_webapp
(flutype_webapp_analysis) pip install -e .
```
Install kernel
```
(flutype_webapp_analysis) ipython kernel install --user --name=flutype_webapp_analysis

