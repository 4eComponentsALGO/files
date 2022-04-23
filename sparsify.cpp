//Runs with ./sparsify <input_graph> <output_graph>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>

int find(int*,int);
void unite(int*,int*,int,int);

using namespace std::chrono;

int main(int nArgs, char** args)
{
   srand(time(NULL));
   //source filename, target filename
   int n; int m;
   FILE* fp = fopen(args[1],"r");
   fscanf(fp,"%d %d",&n,&m);
   int* graph_edges = (int*)malloc(sizeof(int)*2*m);
   for(int i=0;i<m;i++){fscanf(fp,"%d %d",graph_edges+2*i,graph_edges+2*i+1);}
   fclose(fp);
high_resolution_clock::time_point t1 = high_resolution_clock::now();
   int* edges = (int*)malloc(sizeof(int)*8*n);
   char* is_available = (char*)malloc(sizeof(char)*m);
   for(int i=0;i<m;i++){is_available[i]=1;}
   int* ufparent = (int*)malloc(sizeof(int)*n);
   int* ufsize = (int*)malloc(sizeof(int)*n);
   int edgeIndx=0;
   for(int t=0;t<4;t++)
   {
      for(int x=0;x<n;x++){ufparent[x]=x; ufsize[x]=1;}
      for(int i=0;i<m;i++)
      {
         if(!is_available[i]){continue;}
         int x=graph_edges[2*i]; int y=graph_edges[2*i+1];
         if(find(ufparent,x)!=find(ufparent,y))
         {
            unite(ufparent,ufsize,x,y);
            edges[2*edgeIndx]=x; edges[2*edgeIndx+1]=y; edgeIndx++;
            is_available[i]=0;
         } 
      }
   }
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
   fp = fopen(args[2],"w");
   fprintf(fp,"%d %d\n",n,edgeIndx);
   for(int i=0;i<edgeIndx;i++){fprintf(fp,"%d %d\n",edges[2*i],edges[2*i+1]);}
   fclose(fp);

   free(edges); free(graph_edges); free(is_available);
   free(ufparent); free(ufsize);
   printf("time: %f\n",time_span.count());
   return 0;
}

void unite(int* p, int* size, int x, int y)
{
   int r1=find(p,x);
   int r2=find(p,y);
   if(r1==r2){return;}
   if(size[r1]<size[r2]){p[r1]=r2;size[r2]+=size[r1];}
   else{p[r2]=r1;size[r1]+=size[r2];}
}

int find(int* p, int x)
{
   int r=x;
   while(p[r]!=r){r=p[r];}
   while(x!=r){int next=p[x];p[x]=r;x=next;}
   return r;
}
