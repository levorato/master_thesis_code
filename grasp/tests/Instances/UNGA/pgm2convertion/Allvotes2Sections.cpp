#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main(int argc, char *argv[])
{
    char fileInName[100],
         fileOutName[100];
    FILE *fileIn, 
         *fileOut=NULL;

    double rcid; //resolution number?
    int section,
        vote,
        cnum,
        actual_section=-1,
        tot_sections=0;

    char date[15],
         cname[5],
         resolution[20],
         aux_section[5];
    
    strcpy(fileInName,argv[1]);
    fileIn = fopen(fileInName,"r");
    if(!fileIn)
    {
       printf("Fail to open file In");
       exit(1);
    }

    while(!feof(fileIn)){   
       //read line
       fscanf(fileIn,"%lf %d %s %d %s %d",&rcid,&section,&date,&cnum,&cname,&vote);
       if(section!=actual_section){
          tot_sections++;
          //a new section starts
          if(fileOut!=NULL) fclose(fileOut);
          
          itoa(section,aux_section,10);
          strcpy(fileOutName,"Section");
          strcat(fileOutName,aux_section);
          strcat(fileOutName,".txt");
               
          actual_section= section;
          
          fileOut = fopen(fileOutName,"w");
          if(!fileOut)
          {
             printf("Fail to open file Out");
             exit(1);
          }
       }          
       
       //write this line   
       fprintf(fileOut,"%d %d %s %d %s %d\n",(int) ceil(rcid),section,date,cnum,cname,vote);

    }

    fclose(fileIn);
    
    printf("Tot Sections: %d\n",tot_sections);
     
    return(0);
}
