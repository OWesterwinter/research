********************************************************************************
* Do file replication "The complex structure of commercial peace
* contrasting trade interdependence,asymmetry, and multipolarity", 
* WEIS Goldstein data
********************************************************************************

* Set working directory

cd "I:\Conflict and Trade Dependencies\Replication_Files_JPR_2016"

* Load WEIS data (Goldstein scale)

use "I:\Conflict and Trade Dependencies\Replication_Files_JPR_2016\Replication_Data_WEIS_GOLDSTEIN_JPR_2016.dta", clear

********************************************************************************
* Zero-inflated negative binomial models, dyadic asymmetric dependence and
* interdependence and their interaction with extra-dyadic asymmetric dependence
* and interdependence
********************************************************************************

zinb countgoldmil logtradelo logdifftrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig gldconyrs _gldcnflt*, inflate(distance loggdplo) robust

fitstat

zinb countgoldnomil logtradelo reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig gldconyrs _gldcnflt*, inflate(distance loggdplo) robust

fitstat

zinb countgoldmil logtradelo c.logdifftrade##c.countdephidy reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig gldconyrs _gldcnflt*, inflate (distance loggdplo) robust

fitstat

zinb countgoldmil logtradelo c.logdifftrade##c.countdeplody reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig gldconyrs _gldcnflt*, inflate(distance loggdplo) robust

fitstat

zinb countgoldmil logdifftrade c.logtradelo##c.sumcountlargebitrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig gldconyrs _gldcnflt*, inflate(distance loggdplo) robust

fitstat

* Conduct vuong and zip test

zinb countgoldmil logtradelo logdifftrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig gldconyrs _gldcnflt*, inflate(distance loggdplo) vuong zip

zinb countgoldnomil logtradelo reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig gldconyrs _gldcnflt*, inflate(distance loggdplo) vuong zip

zinb countgoldmil logtradelo c.logdifftrade##c.countdephidy reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig gldconyrs _gldcnflt*, inflate (distance loggdplo) vuong zip

zinb countgoldmil logtradelo c.logdifftrade##c.countdeplody reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig gldconyrs _gldcnflt*, inflate(distance loggdplo) vuong zip

zinb countgoldmil logdifftrade c.logtradelo##c.sumcountlargebitrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig gldconyrs _gldcnflt*, inflate(distance loggdplo) vuong zip

********************************************************************************
* End of do file
********************************************************************************
