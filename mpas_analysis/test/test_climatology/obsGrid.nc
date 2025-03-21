CDF       
      lat    �   lon   h         creation_date         Mon Apr  7 14:45:47 MDT 2008       creator       Dennis Shea, CGD       story        �
Clara Deser and Jim Hurrell wanted a SST and SEAICE datasets that spanned
1870-2005. The Hadley data from http://www.hadobs.org/ were downloaded.  
These SST values were used to compute a Hadley 1971-2000 climatology.    
This Hadley climatology  was used to compute anomalies. The Hadley       
anomalies were added to the NCEP-OI2 1971-2000 climatology to create     
a new SST dataset spanning 1870-2005.                                    
                                                                         
The 187001-198110 data from these new datasets were merged with the      
the NCEP-OI2 198111-Present dataset.                                      
                                                                         
All SST values were subjected to the Hurrell SST/SEAICE consistency      
criteria.                                                                
                                                                         
The same software was used for land interpolation Hadley and NCEP-OI2.   
Linear interpolation was used when only a few pts needed to be filled.   
Near coasts of large land areas, a nearest neighbor was used for a few   
grid points inland. Then a iterative bilinear interpolator               
was used. Results are visually ok. Central Siberia has some issues       
because areas like the Aral Sea influence the interpolations.            
                                                                         
     title         ASST: Merged: HADLEY 187001-198110 with  NCEP OI2: 198111-present       history      �Fri Mar 31 00:48:17 2017: ncks -v lon,lat SST_annual_1870-1900.nc /home/xylar/code/mpas-work/analysis/fix_climatology_dates/mpas_analysis/test/test_climatology/obsGrid.nc
Wed Aug 19 17:28:31 2015: ncks -O -F -d time,1,1,1 /ccs/home/milena/ACME/observations/obsdir/SST/MODEL.SST.HAD187001-198110.OI198111-201203.nc /ccs/home/milena/ACME/observations/obsdir/SST/SST_annual_1870-1900.nc
Tue Apr 10 11:08:42 2012: ncrcat MODEL.SST.HAD187001-198110.OI198111-201112.nc SST.update.nc MODEL.SST.HAD187001-198110.OI198111-201203.nc
Mon Jan  9 09:04:34 2012: ncrcat MODEL.SST.HAD187001-198110.OI198111-201108.nc SST.update.nc MODEL.SST.HAD187001-198110.OI198111-201112.nc
Tue Oct  4 16:44:56 2011: ncrcat MODEL.SST.HAD187001-198110.OI198111-201103.nc SST.update.nc MODEL.SST.HAD187001-198110.OI198111-201108.nc
Mon Feb  7 14:12:59 2011: ncrcat MODEL.SST.HAD187001-198110.OI198111-201003.nc SST.update.nc MODEL.SST.HAD187001-198110.OI198111-201101.nc
Wed Apr  7 10:52:47 2010: ncrcat MODEL.SST.HAD187001-198110.OI198111-200903.nc SST.update.nc MODEL.SST.HAD187001-198110.OI198111-201003.nc
Fri Apr 17 09:39:47 2009: ncrcat MODEL.SST.HAD187001-198110.OI198111-200803.nc SST.update.nc MODEL.SST.HAD187001-198110.OI198111-200903.nc      nco_openmp_thread_number            NCO       "4.6.4"          lat                 	long_name         	latitude       units         degrees_north        �  �   lon                units         degrees_east       	long_name         
longitude        �  �³  ±  ¯  ­  «  ©  §  ¥  £  ¡                                  �~  �z  �v  �r  �n  �j  �f  �b  �^  �Z  �V  �R  �N  �J  �F  �B  �>  �:  �6  �2  �.  �*  �&  �"  �  �  �  �  �  �
  �  �  ��  ��  ��  ��  ��  ��  ��  ��  ��  ��  ��  ��  ��  ��  ��  ��  �x  �h  �X  �H  �8  �(  �  �  ��  ��  ��  ��  �`  �   ��  �   ?   ?�  @   @`  @�  @�  @�  @�  A  A  A(  A8  AH  AX  Ah  Ax  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  B  B  B
  B  B  B  B  B  B"  B&  B*  B.  B2  B6  B:  B>  BB  BF  BJ  BN  BR  BV  BZ  B^  Bb  Bf  Bj  Bn  Br  Bv  Bz  B~  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  ?   ?�  @   @`  @�  @�  @�  @�  A  A  A(  A8  AH  AX  Ah  Ax  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  B  B  B
  B  B  B  B  B  B"  B&  B*  B.  B2  B6  B:  B>  BB  BF  BJ  BN  BR  BV  BZ  B^  Bb  Bf  Bj  Bn  Br  Bv  Bz  B~  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C � C� C� C� C� C� C� C� C� C	� C
� C� C� C� C� C� C� C� C� C� C� C� C� C� C� C� C� C� C� C� C� C� C � C!� C"� C#� C$� C%� C&� C'� C(� C)� C*� C+� C,� C-� C.� C/� C0� C1� C2� C3� C4� C5� C6� C7� C8� C9� C:� C;� C<� C=� C>� C?� C@� CA� CB� CC� CD� CE� CF� CG� CH� CI� CJ� CK� CL� CM� CN� CO� CP� CQ� CR� CS� CT� CU� CV� CW� CX� CY� CZ� C[� C\� C]� C^� C_� C`� Ca� Cb� Cc� Cd� Ce� Cf� Cg� Ch� Ci� Cj� Ck� Cl� Cm� Cn� Co� Cp� Cq� Cr� Cs� Ct� Cu� Cv� Cw� Cx� Cy� Cz� C{� C|� C}� C~� C� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� C�@ C�� 