*********************************************************************
* Do file replication "The complex structure of commercial peace
* contrasting trade interdependence,asymmetry, and multipolarity", 
* COPDAB data
*********************************************************************

*********************************************************************
* Set working directory and load data
*********************************************************************

* Set working directory

cd "E:\Conflict and Trade Dependencies\Replication_Files_JPR_2016"

* Load COPDAB data

use "E:\Conflict and Trade Dependencies\Replication_Files_JPR_2016\Replication_Data_COPDAB_JPR_2016.dta", clear

********************************************************************************
* Zero-inflated negative binomial models, dyadic asymmetric dependence and
* interdependence
********************************************************************************

zinb countcopconfl logtradelo logdifftrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robust

fitstat

zinb countcopmil logtradelo logdifftrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robust

fitstat

zinb countcopnomil logtradelo reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robust

fitstat

* Conduct vuong and zip test

zinb countcopconfl logtradelo logdifftrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) vuong zip

zinb countcopmil logtradelo logdifftrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) vuong zip

zinb countcopnomil logtradelo reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) vuong zip

* Compute substantive effects

zinb countcopconfl logtradelo logdifftrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robust

listcoef, help

zinb countcopmil logtradelo logdifftrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robust

listcoef, help

zinb countcopnomil logtradelo reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robust

listcoef, help

* Compute predicted events

* Asymmetric dyadic dependence and militarized conflict events

zinb countcopmil logtradelo logdifftrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robust

margins, at(logdifftrade=(0(0.5)4)) atmeans vsquish

marginsplot, xtitle("Dependence difference") recast(line) recastci(rline) ///
ytitle("Predicted number of" "militarized conflict events") yline(0)  ///
title("") scheme(s1manual)

* Dyadic interdependence and non-militarized conflict events

zinb countcopnomil logtradelo reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robust

margins, at(logtradelo=(0(0.01)0.5)) atmeans vsquish

marginsplot, xtitle("Dependence low") recast(line) recastci(rline) ///
ytitle("Predicted number of" "non-militarized conflict events") yline(0)  ///
title("") scheme(s1manual) xscale(range(0 0.5)) xlabel(0[0.1]0.5)

********************************************************************************
* Zero-inflated negative binomial models, dyadic asymmetric dependence and
* interdependence and their interaction with extra-dyadic asymmetric dependence
* and interdependence
********************************************************************************

zinb countcopmil logtradelo c.logdifftrade##c.countdephidy reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate (distance loggdplo) robust

fitstat

zinb countcopmil logtradelo c.logdifftrade##c.countdeplody reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robust

fitstat

zinb countcopmil logdifftrade c.logtradelo##c.sumcountlargebitrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robust

fitstat

* Conduct vuong and zip test

zinb countcopmil logtradelo c.logdifftrade##c.countdephidy reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate (distance loggdplo) vuong zip

zinb countcopmil logtradelo c.logdifftrade##c.countdeplody reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) vuong zip

zinb countcopmil logdifftrade c.logtradelo##c.sumcountlargebitrade reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) vuong zip

* Compute marginal effects for interation terms

* Interaction extra-dyadic dependencies more dependent state in dyad

zinb countcopmil logtradelo c.logdifftrade##c.countdephidy reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate (distance loggdplo) robust

margins, dydx(logdifftrade) at(countdephidy=(0(5)80)) atmeans post

mat t = e(at)
mat t = t[1...,"countdephidy"]
mat t = t,J(17,3,.)

forvalue i=1/17 {
  mat t[`i',2] = _b[`i'._at]
  mat t[`i',3] = _b[`i'._at] - 1.96*_se[`i'._at]
  mat t[`i',4] = _b[`i'._at] + 1.96*_se[`i'._at]
  }
  
mat colnames t = at dif ll ul
svmat t, names(col)

twoway (line dif at)(line ll at)(line ul at), ///
       yline(0) legend(off) ///
       xtitle("Number high dependence trade relationships" "of more dependent state in dyad") ///
	   ytitle("Effect on predicted number" "of militarized conflict events") ///
       scheme(s1manual)
 
drop at dif ll ul
 
* Interaction extra-dyadic dependencies less dependent state in dyad

zinb countcopmil logtradelo c.logdifftrade##c.countdeplody reglo diffreg logrgdppclo loggdplo sun3cati capratio allydmy distance contig conflyrs _cnflct*, inflate(distance loggdplo) robus

margins, dydx(logdifftrade) at(countdeplody=(0(5)80)) atmeans post

mat t = e(at)
mat t = t[1...,"countdeplody"]
mat t = t,J(17,3,.)

forvalue i=1/17 {
  mat t[`i',2] = _b[`i'._at]
  mat t[`i',3] = _b[`i'._at] - 1.96*_se[`i'._at]
  mat t[`i',4] = _b[`i'._at] + 1.96*_se[`i'._at]
  }
  
mat colnames t = at dif ll ul
svmat t, names(col)

twoway (line dif at)(line ll at)(line ul at), ///
       yline(0) legend(off) ///
       xtitle("Number high dependence trade relationships" "of less dependent state in dyad") ///
	   ytitle("Effect on predicted number" "of militarized conflict events") ///
       scheme(s1manual)

drop dif at ll ul
	   
********************************************************************************
* End of do file
********************************************************************************
