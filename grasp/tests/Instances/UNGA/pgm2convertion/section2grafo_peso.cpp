// Votes: 1- yes
//        2- Abst
//        3- No
//	  	  8- Absent
//        9- Not a member

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CCODE_MAX 1000
#define COUNTRY_MAX 200
#define EPSILON 0.001

typedef struct country{
    int pos;
    char cname[5];
} id_country;

int main(int argc, char *argv[])
{ 
    
    int n=0, 
        m=0;

    int rcid, //resolution number?
        actual_rcid,
        last_rcid,
        section,
        vote,
        cnum,
        tot_resolutions=0;

    char date[15],
         cname[5];

    id_country ccode[CCODE_MAX];
    int vote_country[COUNTRY_MAX];
    double Sum_pos[COUNTRY_MAX][COUNTRY_MAX],
           Sum_neg[COUNTRY_MAX][COUNTRY_MAX];
    double SG[COUNTRY_MAX][COUNTRY_MAX];
  
    int i, j, k, sinal;

    char fileNameIn[100],
         fileNameOut[100];
    FILE *fileIn, *fileOut;

    strcpy(fileNameIn,argv[1]);
    fileIn = fopen(fileNameIn,"r");
    if(!fileIn)
    {
       printf("Fail to open file\n");
       exit(0);
    }

    //initialization
    for(i=0;i<CCODE_MAX;i++)
       ccode[i].pos = -1;

    for(i=0;i<COUNTRY_MAX;i++)
	   for(j=0;j<COUNTRY_MAX;j++){
	      Sum_pos[i][j] = 0.0;
	      Sum_neg[i][j] = 0.0;
       }

    //sum over resolutions -> Sum_pos and Sum_neg
    last_rcid=-1;
    while(!feof(fileIn)){   
       //read a resolution
       fscanf(fileIn,"%d %d %s %d %s %d",&rcid,&section,&date,&cnum,&cname,&vote);
       if(actual_rcid!=rcid){
          tot_resolutions++;
          //a new resolution starts       
          last_rcid=actual_rcid;
          actual_rcid=rcid;
          
          //calculate weights for last resolution
          if(last_rcid!=-1){
             for(i=0;i<n;i++)
            	 if(vote_country[i]!=8){
            		 for(j=i+1;j<n;j++){
            			 if(vote_country[j]!=8){
            				 if(vote_country[i]==vote_country[j]){
            					 if(vote_country[i]!= 2)   
            						 Sum_pos[i][j]=Sum_pos[i][j]+ 1.0;
            					 else 
            						 Sum_pos[i][j]=Sum_pos[i][j]+ 0.5;                         
            				 }
            				 else{
            					 if(vote_country[i]!= 2 && vote_country[j]!= 2)
            						 Sum_neg[i][j]=Sum_neg[i][j]- 1.0;
            					 else
            						 Sum_neg[i][j]=Sum_neg[i][j]- 0.5;
            			 	 }
            			 }
            		 }
                }
          }
          
          //initialize array vote_country
          for(i=0;i<COUNTRY_MAX;i++)
             vote_country[i]=0;
       }          
       //register this vote   
       if(vote!=9){
          //country belongs to UNGA
          if(ccode[cnum].pos == -1){
             //first data from this country
             ccode[cnum].pos = n++;
             strcpy(ccode[cnum].cname,cname);
          }
          vote_country[ccode[cnum].pos]= vote;
       }

    }
    //calculate weights for last resolution in the file
    last_rcid=actual_rcid;
    if(last_rcid!=-1){
       for(i=0;i<n;i++)
          for(j=i+1;j<n;j++){
             if(vote_country[i]==vote_country[j]){
                if(vote_country[i]!= 2)   
                   Sum_pos[i][j]=Sum_pos[i][j]+ 1.0;
                else 
                   Sum_pos[i][j]=Sum_pos[i][j]+ 0.5;                         
                }
             else{
                if(vote_country[i]!= 2 && vote_country[j]!= 2)
                   Sum_neg[i][j]=Sum_neg[i][j]- 1.0;
                else
                   Sum_neg[i][j]=Sum_neg[i][j]- 0.5;
             }
          }
    }


    fclose(fileIn);

    //generate output with the positive and negative sums

    printf("tot_resolutions: %d\n",tot_resolutions);
    
    printf("Ligacoes positivas: \n");
    for(i=0;i<CCODE_MAX;i++)
       if(ccode[i].pos != -1){
          for(j=i+1;j<CCODE_MAX;j++)
             if(ccode[j].pos != -1){
                double PercPos = Sum_pos[ccode[i].pos][ccode[j].pos]/tot_resolutions;
                double PercNeg = Sum_neg[ccode[i].pos][ccode[j].pos]/tot_resolutions;
                double PercSum = PercPos + PercNeg;
                if(PercSum >= 0.0){
                   printf("(%s,%s): +%6.3lf | %6.3lf *** +%5.3lf | %5.3lf = %5.3lf\n",
                        ccode[i].cname,ccode[j].cname,
                        Sum_pos[ccode[i].pos][ccode[j].pos],
                        Sum_neg[ccode[i].pos][ccode[j].pos],
                        PercPos,
                        PercNeg, 
                        PercSum);
                }
             }
       }
       
    printf("Ligacoes Negativas: \n");
    for(i=0;i<CCODE_MAX;i++)
       if(ccode[i].pos != -1){
          for(j=i+1;j<CCODE_MAX;j++)
             if(ccode[j].pos != -1){
                double PercPos = Sum_pos[ccode[i].pos][ccode[j].pos]/tot_resolutions;
                double PercNeg = Sum_neg[ccode[i].pos][ccode[j].pos]/tot_resolutions;
                double PercSum = PercPos + PercNeg;
                if(PercSum < 0.0){
                   printf("(%s,%s): +%6.3lf | %6.3lf *** +%5.3lf | %5.3lf = %5.3lf\n",
                        ccode[i].cname,ccode[j].cname,
                        Sum_pos[ccode[i].pos][ccode[j].pos],
                        Sum_neg[ccode[i].pos][ccode[j].pos],
                        PercPos,
                        PercNeg, 
                        PercSum);
                }
             }
       }

             
    //generate signed graph
    for(i=0;i<COUNTRY_MAX;i++)
	   for(j=0;j<COUNTRY_MAX;j++){
	      SG[i][j] = (double) tot_resolutions;
       }
    for(i=0;i<CCODE_MAX;i++)
       if(ccode[i].pos != -1){
          for(j=i+1;j<CCODE_MAX;j++)
             if(ccode[j].pos != -1){
                double PercPos = Sum_pos[ccode[i].pos][ccode[j].pos]/tot_resolutions;
                double PercNeg = Sum_neg[ccode[i].pos][ccode[j].pos]/tot_resolutions;
                double PercSum = PercPos + PercNeg;
                if(PercSum < -EPSILON || PercSum > EPSILON){
                   SG[ccode[i].pos][ccode[j].pos] = SG[ccode[j].pos][ccode[i].pos] = PercSum;
                   m++;
                }
             }
       } 
         
              
//    //Print output files
//    //Signed graph
    strncpy(fileNameOut,fileNameIn,(strlen(fileNameIn)-3));
    fileNameOut[strlen(fileNameIn)-3]='\0';
    strcat(fileNameOut,"g");
    fileOut = fopen(fileNameOut,"w");
    if(!fileOut)
    {
       printf("Fail to open file\n");
       exit(0);
    }
    fprintf(fileOut,"%d %d\n",n,m);
    for(i=0;i<COUNTRY_MAX;i++)
       for(j=i+1;j<COUNTRY_MAX;j++){
          if(SG[i][j] < (double) tot_resolutions)
             fprintf(fileOut,"%d %d %5.4lf\n",i,j,SG[i][j]);
       }
    fclose(fileOut);
    //country codes for vertices in the signed graph
    strcat(fileNameOut,".ccode");    
    fileOut = fopen(fileNameOut,"w");
    if(!fileOut)
    {
       printf("Fail to open file");
       exit(0);
    }    
    fprintf(fileOut,"<ccode> <cname> <vertex_label>\n");
    for(i=0;i<CCODE_MAX;i++)
       if(ccode[i].pos!=-1)
          fprintf(fileOut,"%d %s %d\n",i,ccode[i].cname,ccode[i].pos);
    fclose(fileOut);

    
    return(0);
}
