#include <string.h>
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv[]){

    char fileInName[100],
         fileOutName[100];
    FILE *fileIn,
         *fileOut=NULL;

    strcpy(fileInName,argv[1]);
    fileIn = fopen(fileInName,"r");
    if(!fileIn)
    {
       printf("Fail to open file In");
       exit(1);
    }

    printf("abriu");
    strcpy(fileOutName,argv[1]);
    strcat(fileOutName,".datmos");
    fileOut = fopen(fileOutName,"w");
    if(!fileOut)
    {
       printf("Fail to open file In");
       exit(1);
    }

    int i, j, k;
    int n, m;
    double w;

    fscanf(fileIn,"%d %d",&n,&m);
    fprintf(fileOut,"people: %d\n",n);
    fprintf(fileOut,"VarErr: 0.5\n");
    fprintf(fileOut,"Names: [\n");
    for (i=1;i<=n;i++)
        fprintf(fileOut,"%d ",i);
    fprintf(fileOut,"]\n");

    fprintf(fileOut,"Mrel: [\n");
    for (k=0;k<m;k++){ 
        fscanf(fileIn,"%d %d %lf",&i,&j,&w);
        fprintf(fileOut,"( %d , %d ) %d\n",i+1,j+1,(int) (w*10000));
    }
    fprintf(fileOut,"]\n");
   
    fclose(fileOut);
    return 1;
}
