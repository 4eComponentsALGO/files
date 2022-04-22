#include <stdio.h>
#include <stdlib.h>

int toint(char*);

void main(int n_args, char** args)
{
   int n = toint(args[1]);
   int m = toint(args[2]);
   srand(toint(args[3]));
   int* edges = (int*)malloc(sizeof(int)*2*m);
   for(int i=0;i<m;i++)
   {
      int x=rand()%n; int y=(x+1+rand()%(n-1))%n;
      edges[2*i]=x; edges[2*i+1]=y;
   }
   FILE* fp = fopen(args[4],"w");
   fprintf(fp,"%d %d\n",n,m);
   for(int i=0;i<m;i++){fprintf(fp,"%d %d\n",edges[2*i],edges[2*i+1]);}
   fclose(fp);
   free(edges);
}

int toint(char* x)
{
   int l=0;
   while(x[l]!=0){l++;}
   int pow10=1;
   int ret=0;
   l--;
   while(l!=-1){ret+=pow10*(x[l]-'0');l--;pow10*=10;}
   return ret;
}
