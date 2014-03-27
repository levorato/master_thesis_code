***files .g

Describe a directed weighted signed graph for each vote section. 
Each signed graph was construct as follows:

For each pair of counstries i < j:
   sum_Pos=0
   sum_Neg=0
   For each resolution:
      if (i votes YES) e (j votes YES) OR (i votes NO) e (j vota NO) then
         sum_Pos += 1
      if (i abst) AND (j abst) then
         sum_Pos += 0.5
      if (i votes YES) and (j votes NO) OR (i votes NO) and (j votes YES) then
         sum_Neg -= 1
      if (i votes YES) and (j abst) OR (i abst) and (j votes YES) then
         sum_Neg -= 0.5
      se (i votes NO) and (j abst) OR (i abst) and (j votes NO) then
         sum_Neg -= 0.5
   End-For
   Perc_Pos = sum_Pos/tot_resolutions 
   Perc_Neg = sum_Neg/tot_resolutions
 
   If (Perc_Pos + Perc_Neg != 0.0) then
      creates an arc (i,j) with weight Perc_Pos + Perc_Neg
   
End-For

***files .g.ccode

Describe, for each vertex i in the signed graph described by the associated file .g, which
country (code and name) is defined by i.
