//Runs with ./get4eComponentsGIK <input_graph> <output_4components>

#include <stdio.h>
#include <stdlib.h>
#include <chrono>

void get_adj(int,int,int*,int**,int**);
void read_graph(char*,int*,int**,int**);

void DFS(int,int*,int*,int**,int**,int**);
void get_low(int,int*,int*,int*,int*,int*,int**);
void get_l_and_bcount(int,int*,int*,int*,int*,int*,int**,int**);
void get_lowChildren(int,int*,int*,int*,int*,int**,int**);
void get_M(int,int*,int*,int*,int*,int*,int*,int**,int**);

int get_4e_components_high(int,int*,int*,int**);
int get_4e_components_connected_high(int,int*,int*,int**);
int get_4e_components_2econnected_high(int,int*,int*,int**);
int get_4e_components_3econnected_high(int,int*,int*,int**,int*,int**);
void get_high(int,int*,int*,int*,int*,int*,int**,int**);
void get_allM(int,int*,int*,int*,int*,int*,int*,int*,int**,int**,int**);
void get_2low(int,int*,int*,int*,int*,int*,int**,int**,int**,int**);
int find(int*,int);
void unite(int*,int*,int,int);

using namespace std::chrono;

int numberOfComponents=0;
int numberOf2eComponents=0;
int numberOf3eComponents=0;
int numberOf4eComponents=0;
int numberOf1Cuts=0;
int numberOf2Cuts=0;
int numberOf3Cuts=0;
double timeForComponents=0;
double timeFor2eComponents=0;
double timeFor3eComponents=0;
double timeFor4eComponents=0;

int main(int n_args, char** args)
{
   int n; int* adj; int* firstOut;
   read_graph(args[1],&n,&adj,&firstOut);
   int* C;
high_resolution_clock::time_point t1 = high_resolution_clock::now();
   int k=get_4e_components_high(n,adj,firstOut,&C);
high_resolution_clock::time_point t2 = high_resolution_clock::now();
   FILE* fp = fopen(args[2],"w");
   fprintf(fp,"%d\n",n);
   for(int i=0;i<n;i++){fprintf(fp,"%d\n",C[i]);}
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
   fprintf(fp,"%f\n",time_span.count());
   fprintf(fp,"%f\n",timeFor4eComponents);
   fclose(fp);
   free(adj); free(firstOut); free(C);

printf("number of components: %d\n",numberOfComponents);
printf("time to compute components: %f\n\n",timeForComponents);

printf("number of 1-cuts: %d\n",numberOf1Cuts);
printf("number of 2e-components: %d\n",numberOf2eComponents);
printf("time to compute 2e-components: %f\n\n",timeFor2eComponents);

printf("number of 2-cuts: %d\n",numberOf2Cuts);
printf("number of 3e-components: %d\n",numberOf3eComponents);
printf("time to compute 3e-components: %f\n\n",timeFor3eComponents);

printf("number of 3-cuts: %d\n",numberOf3Cuts);
printf("number of 4e-components: %d\n",numberOf4eComponents);
printf("time to compute 4e-components: %f\n\n",timeFor4eComponents);

printf("total time: %f\n",time_span.count());

   return 0;
}

int get_4e_components_high(int n, int* adj, int* firstOut, int** C)
{
high_resolution_clock::time_point t1;
high_resolution_clock::time_point t2;
duration<double> time_span;
double temp_time=0;

t1 = high_resolution_clock::now();
   (*C) = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){(*C)[i]=-1;}
   int* Q = (int*)malloc(sizeof(int)*n);
   int* edges = (int*)malloc(sizeof(int)*firstOut[n]);
   char* found = (char*)malloc(sizeof(char)*n);
   int* map = (int*)malloc(sizeof(int)*n);
   int* imap = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){found[i]=0;}
   int k=0;
   for(int r=0;r<n;r++)
   {
      if(found[r]){continue;}
      int Nr=0;
      int first=0; int last=0;
      Q[last++]=r; found[r]=1; map[r]=Nr; imap[Nr++]=r;
      while(first!=last)
      {
         int x=Q[first++];
         for(int i=firstOut[x];i<firstOut[x+1];i++)
         {
            int y=adj[i];
            if(!found[y])
            {
               Q[last++]=y; found[y]=1; map[y]=Nr; imap[Nr++]=y;
            }
         } 
      }
      int nC = last;
      if(nC>1)
      {
high_resolution_clock::time_point temp_t1 = high_resolution_clock::now();
         int edgeIndx=0;
         for(int i=0;i<nC;i++)
         {
            int x=Q[i];
            for(int j=firstOut[x];j<firstOut[x+1];j++)
            {
               int y=adj[j];
               if(x<y){edges[2*edgeIndx]=map[x];edges[2*edgeIndx+1]=map[y];edgeIndx++;}
            }
         }
         int* adjC; int* firstOutC; 
         get_adj(nC,edgeIndx,edges,&adjC,&firstOutC);
         int* C2;
         int k2 = get_4e_components_connected_high(nC,adjC,firstOutC,&C2);
         for(int i=0;i<nC;i++){(*C)[imap[i]]=C2[i]+k;}
         k+=k2;
         free(adjC); free(firstOutC); free(C2);
high_resolution_clock::time_point temp_t2 = high_resolution_clock::now();
duration<double> temp_time_span =  duration_cast<duration<double>>(temp_t2 - temp_t1);
temp_time += temp_time_span.count();
      }
      else
      {
         (*C)[r]=k; 
         k++;
numberOf2eComponents++;
numberOf3eComponents++;
numberOf4eComponents++;
      }    
numberOfComponents++; 
   }
   free(Q); free(edges); free(found);  
   free(map); free(imap);
t2 = high_resolution_clock::now();

time_span = duration_cast<duration<double>>(t2 - t1);
timeForComponents = time_span.count()-temp_time;
   return k;
}

int get_4e_components_connected_high(int n, int* adj, int* firstOut, int** C)
{
high_resolution_clock::time_point t1;
high_resolution_clock::time_point t2;
duration<double> time_span;
double temp_time=0;

t1 = high_resolution_clock::now();
   (*C) = (int*)malloc(sizeof(int)*n);
   int* dfs; int* idfs; int* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   int* low;
   get_low(n,adj,firstOut,dfs,idfs,p,&low);
   int* Q = (int*)malloc(sizeof(int)*n);
   int* map = (int*)malloc(sizeof(int)*n);
   int* imap = (int*)malloc(sizeof(int)*n);   
   char* found = (char*)malloc(sizeof(char)*n);
   int* edges = (int*)malloc(sizeof(int)*firstOut[n]);
   for(int i=0;i<n;i++){found[i]=0;}
   int k=0;
for(int v=1;v<n;v++){numberOf1Cuts+=low[v]==v;}
   for(int r=0;r<n;r++)
   {
      if(found[r]){continue;}
      int Nr=0;
      int first=0; int last=0;
      Q[last++]=r; found[r]=1; map[r]=Nr; imap[Nr++]=r;
      while(first!=last)
      {
         int x=Q[first++];
         for(int i=firstOut[x];i<firstOut[x+1];i++)
         {
            int y=adj[i];
            if(found[y]){continue;}
            if((x==p[y]&&low[y]==y)||(y==p[x]&&low[x]==x)){continue;}
            Q[last++]=y; found[y]=1; map[y]=Nr; imap[Nr++]=y;
         }
      }
      int nC = last;
      if(nC>1)
      {
high_resolution_clock::time_point temp_t1 = high_resolution_clock::now();
         int edgeIndx=0;
         for(int i=0;i<last;i++)
         {
            int x=Q[i];
            for(int j=firstOut[x];j<firstOut[x+1];j++)
            {
               int y=adj[j];
               if(y<x){continue;}
               if((x==p[y]&&low[y]==y)||(y==p[x]&&low[x]==x)){continue;}
               edges[2*edgeIndx]=map[x]; edges[2*edgeIndx+1]=map[y]; edgeIndx++;
            }
         }
         int* adjC; int* firstOutC;
         get_adj(nC,edgeIndx,edges,&adjC,&firstOutC);
         int* C3;
         int k3=get_4e_components_2econnected_high(nC,adjC,firstOutC,&C3);
         for(int i=0;i<nC;i++){(*C)[imap[i]]=C3[i]+k;}
         k+=k3;
         free(adjC); free(firstOutC); free(C3);
high_resolution_clock::time_point temp_t2 = high_resolution_clock::now();
duration<double> temp_time_span =  duration_cast<duration<double>>(temp_t2 - temp_t1);
temp_time += temp_time_span.count();
      }
      else
      { 
         (*C)[r]=k;
         k++;
numberOf3eComponents++;
numberOf4eComponents++;
      }
numberOf2eComponents++;
   }
   free(dfs); free(idfs); free(p); free(low);
   free(Q); free(found); free(map); free(imap); free(edges);
t2 = high_resolution_clock::now();

time_span = duration_cast<duration<double>>(t2 - t1);
timeFor2eComponents += time_span.count()-temp_time;
   return k;
}

int get_4e_components_2econnected_high(int n, int* adj, int* firstOut, int** C3)
{
high_resolution_clock::time_point t1;
high_resolution_clock::time_point t2;
duration<double> time_span;
double temp_time=0;

t1 = high_resolution_clock::now();
   int* dfs; int* idfs; int* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   int* low;
   get_low(n,adj,firstOut,dfs,idfs,p,&low);
   int* l; int* bcount;
   get_l_and_bcount(n,adj,firstOut,dfs,idfs,p,&l,&bcount);
   int* low1C; int* low2C;
   get_lowChildren(n,dfs,idfs,p,low,&low1C,&low2C);
   int* M; int* nextM;
   get_M(n,dfs,idfs,l,low,low1C,low2C,&M,&nextM);
   int* prevM = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){prevM[i]=-1;}
   for(int i=1;i<n;i++){if(nextM[i]!=-1){prevM[nextM[i]]=i;}}

   int* vEdgeStack = (int*)malloc(sizeof(int)*4*n);
   int* vEdgeFirst = (int*)malloc(sizeof(int)*n);
   int* vEdgeNext = (int*)malloc(sizeof(int)*4*n);
   int* firstVertex = (int*)malloc(sizeof(int)*4*n);
   for(int i=0;i<n;i++){vEdgeFirst[i]=-1;}
   int SP=0;

   char* isCutEdge = (char*)malloc(sizeof(char)*n);
   char* isCutEdgeLow = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){isCutEdge[i]=0;isCutEdgeLow[i]=0;}

   int* cutEdgeStack = (int*)malloc(sizeof(int)*4*n);
   int* cutEdgeFirst = (int*)malloc(sizeof(int)*n);
   int* cutEdgeNext = (int*)malloc(sizeof(int)*2*n);
   for(int i=0;i<n;i++){cutEdgeFirst[i]=-1;}
   int cutEdgeSP=0;
int* cutEdgeCount = (int*)malloc(sizeof(int)*n);
for(int i=0;i<n;i++){cutEdgeCount[i]=0;}
   for(int m=1;m<n;m++)
   {
      if(M[m]!=m){continue;}
      int u=m;
      while(u!=-1)
      {
         int z=nextM[u];
         if(z!=-1 && bcount[z]==bcount[u])
         {
            int last;
            while(z!=-1 && bcount[z]==bcount[u]){last=z; z=nextM[z];}
            if(bcount[u]==1)
            {
               isCutEdgeLow[m]=1;
               if(u!=m)
               {
                  vEdgeNext[SP]=vEdgeFirst[u]; vEdgeFirst[u]=SP; vEdgeStack[SP]=m; firstVertex[SP++]=u;
                  vEdgeNext[SP]=vEdgeFirst[m]; vEdgeFirst[m]=SP; vEdgeStack[SP]=u; firstVertex[SP++]=u;
               }
               if(p[last]!=low[m])
               {
                  vEdgeNext[SP]=vEdgeFirst[p[last]]; vEdgeFirst[p[last]]=SP; vEdgeStack[SP]=low[m]; firstVertex[SP++]=u;
                  vEdgeNext[SP]=vEdgeFirst[low[m]]; vEdgeFirst[low[m]]=SP; vEdgeStack[SP]=p[last]; firstVertex[SP++]=u;
               }
               cutEdgeNext[cutEdgeSP]=cutEdgeFirst[u]; cutEdgeFirst[u]=cutEdgeSP;
               cutEdgeStack[2*cutEdgeSP]=m; cutEdgeStack[2*cutEdgeSP+1]=low[m]; cutEdgeSP++;
               cutEdgeCount[u]++;
            }
            else
            {
               vEdgeNext[SP]=vEdgeFirst[u]; vEdgeFirst[u]=SP; vEdgeStack[SP]=p[last]; firstVertex[SP++]=u;
               vEdgeNext[SP]=vEdgeFirst[p[last]]; vEdgeFirst[p[last]]=SP; vEdgeStack[SP]=u; firstVertex[SP++]=u;
            }
            int v=u;
            while(v!=nextM[last])
            {
               isCutEdge[v]=1;
               if(v!=last&&p[v]!=nextM[v])
               {
                  vEdgeNext[SP]=vEdgeFirst[p[v]]; vEdgeFirst[p[v]]=SP; vEdgeStack[SP]=nextM[v]; firstVertex[SP++]=u;
                  vEdgeNext[SP]=vEdgeFirst[nextM[v]]; vEdgeFirst[nextM[v]]=SP; vEdgeStack[SP]=p[v]; firstVertex[SP++]=u;
               }
               cutEdgeNext[cutEdgeSP]=cutEdgeFirst[u]; cutEdgeFirst[u]=cutEdgeSP;
               cutEdgeStack[2*cutEdgeSP]=v; cutEdgeStack[2*cutEdgeSP+1]=p[v]; cutEdgeSP++;
               cutEdgeCount[u]++;
               v=nextM[v];
            }
         }
         else if(bcount[u]==1)
         {
            isCutEdge[u]=1; isCutEdgeLow[m]=1;
            if(u!=m)
            {
               vEdgeNext[SP]=vEdgeFirst[u]; vEdgeFirst[u]=SP; vEdgeStack[SP]=m; firstVertex[SP++]=u;
               vEdgeNext[SP]=vEdgeFirst[m]; vEdgeFirst[m]=SP; vEdgeStack[SP]=u; firstVertex[SP++]=u;
            }
            if(p[u]!=low[m])
            {
               vEdgeNext[SP]=vEdgeFirst[p[u]]; vEdgeFirst[p[u]]=SP; vEdgeStack[SP]=low[m]; firstVertex[SP++]=u;
               vEdgeNext[SP]=vEdgeFirst[low[m]]; vEdgeFirst[low[m]]=SP; vEdgeStack[SP]=p[u]; firstVertex[SP++]=u;
            }
            cutEdgeNext[cutEdgeSP]=cutEdgeFirst[u]; cutEdgeFirst[u]=cutEdgeSP;
            cutEdgeStack[2*cutEdgeSP]=u; cutEdgeStack[2*cutEdgeSP+1]=p[u]; cutEdgeSP++;
            cutEdgeNext[cutEdgeSP]=cutEdgeFirst[u]; cutEdgeFirst[u]=cutEdgeSP;
            cutEdgeStack[2*cutEdgeSP]=m; cutEdgeStack[2*cutEdgeSP+1]=low[m]; cutEdgeSP++;
            cutEdgeCount[u]=2;
         }
         u=z;
      } 
   } 

   int* C = (int*)malloc(sizeof(int)*n);
   int* Q = (int*)malloc(sizeof(int)*n);
   int* map = (int*)malloc(sizeof(int)*n);
   int* imapStack = (int*)malloc(sizeof(int)*n);
   int* imapFirst = (int*)malloc(sizeof(int)*n);
   int* imapNext = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){C[i]=-1;imapFirst[i]=-1;}
   SP=0;
   int k=0;
   for(int r=0;r<n;r++)
   {
      if(C[r]!=-1){continue;}
      int first=0; int last=0;
      C[r]=k; Q[last++]=r;
      int Nr=0;
      map[r]=Nr++;
      imapNext[SP]=imapFirst[k]; imapFirst[k]=SP; imapStack[SP++]=r;
      while(first!=last)
      {
         int v=Q[first++];
         for(int i=firstOut[v];i<firstOut[v+1];i++)
         {
            int u=adj[i];
            if((u==p[v]&&isCutEdge[v])||(v==p[u]&&isCutEdge[u])){continue;}
            if((M[u]==u&&v==low[u]&&isCutEdgeLow[u])||(M[v]==v&&u==low[v]&&isCutEdgeLow[v])){continue;}
            if(C[u]==-1)
            {
               C[u]=k; Q[last++]=u; map[u]=Nr++;
               imapNext[SP]=imapFirst[k]; imapFirst[k]=SP; imapStack[SP++]=u;
            }
         }
         for(int i=vEdgeFirst[v];i!=-1;i=vEdgeNext[i])
         {
            int u=vEdgeStack[i];
            if(C[u]==-1)
            {
               C[u]=k; Q[last++]=u; map[u]=Nr++;
               imapNext[SP]=imapFirst[k]; imapFirst[k]=SP; imapStack[SP++]=u;
            }
         }
      }
      k++;
   }

numberOf3eComponents+=k;
for(int v=1;v<n;v++){numberOf2Cuts+=bcount[v]==1;}
for(int x=1;x<n;x++)
{
   if(M[x]!=x){continue;}
   int u=x;
   while(u!=-1)
   {
      int v=u;
      int temp_num=1;
      while(nextM[v]!=-1 && bcount[nextM[v]]==bcount[v])
      {
         temp_num++;
         v=nextM[v];
      }
      numberOf2Cuts+=(temp_num*(temp_num-1))/2;
      u=nextM[v];
   }
}
t2 = high_resolution_clock::now();

time_span = duration_cast<duration<double>>(t2 - t1);
timeFor3eComponents += time_span.count();

t1 = high_resolution_clock::now();
   int* nC = (int*)malloc(sizeof(int)*k);
   int** adjC = (int**)malloc(sizeof(int*)*k);
   int** firstOutC = (int**)malloc(sizeof(int*)*k);
   int* edges = (int*)malloc(sizeof(int)*(firstOut[n]+n));
   for(int c=0;c<k;c++)
   {
      nC[c]=0;
      int eIndx=0;
      for(int t=imapFirst[c];t!=-1;t=imapNext[t])
      {
         nC[c]++;
         int x=imapStack[t];
         for(int i=firstOut[x];i<firstOut[x+1];i++)
         {
            int y=adj[i];
            if(dfs[x]<dfs[y]){continue;}
            if((y==p[x]&&isCutEdge[x])||(x==M[x]&&y==low[x]&&isCutEdgeLow[x])){continue;}
            edges[2*eIndx]=map[x]; edges[2*eIndx+1]=map[y]; eIndx++;
         }
         for(int i=vEdgeFirst[x];i!=-1;i=vEdgeNext[i])
         {
            int y=vEdgeStack[i];
            if(y>x){continue;}
            edges[2*eIndx]=map[x]; edges[2*eIndx+1]=map[y]; eIndx++; 
         }
      }
      if(nC[c]==1){continue;}
      get_adj(nC[c],eIndx,edges,adjC+c,firstOutC+c);
   }

   int* nCcuts = (int*)malloc(sizeof(int)*k);
   int** Ccuts = (int**)malloc(sizeof(int*)*k);

   (*C3) = (int*)malloc(sizeof(int)*n);
   int c_num=0;
   for(int c=0;c<k;c++)
   {
      if(nC[c]>1)
      {
         int* C4;
         int k4=get_4e_components_3econnected_high(nC[c],adjC[c],firstOutC[c],&C4,nCcuts+c,Ccuts+c);
         for(int t=imapFirst[c];t!=-1;t=imapNext[t])
         {
            int x=imapStack[t];
            (*C3)[x]=C4[map[x]]+c_num;
         }
         c_num+=k4;
         free(C4);
numberOf4eComponents+=k4;
      }
      else
      {
         int x=imapStack[imapFirst[c]];
         (*C3)[x]=c_num;
         c_num++;    
numberOf4eComponents++;
      }
   }
t2 = high_resolution_clock::now();

time_span = duration_cast<duration<double>>(t2 - t1);
timeFor4eComponents += time_span.count();

   /*get the number of 3-edge cuts*/
   int* imap = (int*)malloc(sizeof(int)*n);
   for(int c=0;c<k;c++)
   {
      if(nC[c]==1){continue;}
      SP=1; int temp;
      for(int i=imapFirst[c];i!=-1;i=imapNext[i])
      {
         imap[nC[c]-SP]=imapStack[i]; SP++;
      }
      for(int i=0;i<nCcuts[c];i++)
      {
         for(int t=0;t<6;t++){Ccuts[c][6*i+t]=imap[Ccuts[c][6*i+t]];}
      }
   }

   int* cut3Stack_ = (int*)malloc(sizeof(int)*12*n);
   int* cut3First_ = (int*)malloc(sizeof(int)*n);
   int* cut3Next_ = (int*)malloc(sizeof(int)*12*n);
   int* CIndx_ = (int*)malloc(sizeof(int)*12*n);
   int* cutIndx_ = (int*)malloc(sizeof(int)*12*n);
   int* edgeIndx_ = (int*)malloc(sizeof(int)*12*n);
   for(int i=0;i<n;i++){cut3First_[i]=-1;}
   SP=0;
   for(int c=0;c<k;c++)
   {
      if(nC[c]==1){continue;}
      for(int i=0;i<nCcuts[c];i++)
      {
         for(int t=0;t<3;t++)
         {
            int x=Ccuts[c][6*i+2*t]; int y=Ccuts[c][6*i+2*t+1];
            cut3Next_[SP]=cut3First_[x]; cut3First_[x]=SP; cut3Stack_[SP]=y; CIndx_[SP]=c; cutIndx_[SP]=i; edgeIndx_[SP++]=t;
            cut3Next_[SP]=cut3First_[y]; cut3First_[y]=SP; cut3Stack_[SP]=x; CIndx_[SP]=c; cutIndx_[SP]=i; edgeIndx_[SP++]=t;
         }
      }
   }
   int* cut3Stack = (int*)malloc(sizeof(int)*12*n);
   int* cut3First = (int*)malloc(sizeof(int)*n);
   int* cut3Next = (int*)malloc(sizeof(int)*12*n);
   int* CIndx = (int*)malloc(sizeof(int)*12*n);
   int* cutIndx = (int*)malloc(sizeof(int)*12*n);
   int* edgeIndx = (int*)malloc(sizeof(int)*12*n);
   for(int i=0;i<n;i++){cut3First[i]=-1;}
   SP=0;
   for(int x=0;x<n;x++)
   {
      for(int i=cut3First_[x];i!=-1;i=cut3Next_[i])
      {
         int y=cut3Stack_[i];
         cut3Next[SP]=cut3First[y]; cut3First[y]=SP; cut3Stack[SP]=x; 
         CIndx[SP]=CIndx_[i]; cutIndx[SP]=cutIndx_[i]; edgeIndx[SP++]=edgeIndx_[i];
      }
   } 
   int* vEdgesStackS = (int*)malloc(sizeof(int)*4*n);
   int* vEdgesFirstS = (int*)malloc(sizeof(int)*n);
   int* vEdgesNextS = (int*)malloc(sizeof(int)*4*n);
   int* firstVertexS = (int*)malloc(sizeof(int)*4*n);
   for(int i=0;i<n;i++){vEdgesFirstS[i]=-1;}
   SP=0;
   for(int x=0;x<n;x++)
   {
      for(int i=vEdgeFirst[x];i!=-1;i=vEdgeNext[i])
      {
         int y=vEdgeStack[i];
         vEdgesNextS[SP]=vEdgesFirstS[y]; vEdgesFirstS[y]=SP; vEdgesStackS[SP]=x; firstVertexS[SP++]=firstVertex[i];
      }   
   }

   int** corVertex = (int**)malloc(sizeof(int*)*k);
   for(int c=0;c<k;c++){if(nC[c]!=1){corVertex[c]=(int*)malloc(sizeof(int)*nCcuts[c]*3);}}
   for(int c=0;c<k;c++){if(nC[c]!=1){for(int i=0;i<3*nCcuts[c];i++){corVertex[c][i]=-1;}}}

   for(int x=0;x<n;x++)
   {
      for(int i=cut3First[x];i!=-1;i=cut3Next[i])
      {
         int y=cut3Stack[i];
         if(y<x){continue;}
         int j=vEdgesFirstS[x];
         while(j!=-1)
         {
            int z=vEdgesStackS[j];
            if(z>y){j=vEdgesNextS[j]; vEdgesFirstS[x]=j; continue;}
            if(z<y){break;}
            int c=CIndx[i]; int t=cutIndx[i]; int e=edgeIndx[i];
            corVertex[c][3*t+e]=firstVertexS[j];
            vEdgesFirstS[x]=j; 

            int i1=cut3Next[i];
            if(i1==-1){break;}
            int y1=cut3Stack[i1];
            if(y1!=y){break;}
            if(cutIndx[i1]!=t){break;}
            i=i1;
            j=vEdgesNextS[j]; vEdgesFirstS[x]=j;
            if(j==-1){break;}
            z=vEdgesStackS[j];
            if(z!=y){break;}
            e=edgeIndx[i];
            corVertex[c][3*t+e]=firstVertexS[j];

            i1=cut3Next[i];
            if(i1==-1){break;}
            y1=cut3Stack[i1];
            if(y1!=y){break;}
            if(cutIndx[i1]!=t){break;}
            i=i1;
            j=vEdgesNextS[j]; vEdgesFirstS[x]=j;
            if(j==-1){break;}
            z=vEdgesStackS[j];
            if(z!=y){break;}
            e=edgeIndx[i];
            corVertex[c][3*t+e]=firstVertexS[j];
            break;
         } 
      }
   }

   int n3Cuts=0;
   for(int c=0;c<k;c++)
   {
      if(nC[c]==1){continue;}
      for(int i=0;i<nCcuts[c];i++)
      {
         int num=1;
         for(int t=0;t<3;t++)
         {
            if(corVertex[c][3*i+t]!=-1){num*=cutEdgeCount[corVertex[c][3*i+t]];}
         }
         n3Cuts+=num;
      } 
   }
numberOf3Cuts+=n3Cuts;
   /*end: get the number of 3-edge cuts*/

   free(vEdgeStack); free(vEdgeFirst); free(vEdgeNext); free(firstVertex);
   free(vEdgesStackS); free(vEdgesFirstS); free(vEdgesNextS); free(firstVertexS);

   free(cutEdgeStack); free(cutEdgeFirst); free(cutEdgeNext); free(cutEdgeCount);

   free(cut3Stack_); free(cut3First_); free(cut3Next_); free(CIndx_); free(cutIndx_); free(edgeIndx_);
   free(cut3Stack); free(cut3First); free(cut3Next); free(CIndx); free(cutIndx); free(edgeIndx);

   free(imap);
   for(int i=0;i<k;i++){if(nC[i]>1){free(Ccuts[i]); free(corVertex[i]);}}
   free(nCcuts); free(Ccuts); free(corVertex);

   for(int i=0;i<k;i++){if(nC[i]>1){free(adjC[i]);free(firstOutC[i]);}}
   free(adjC); free(firstOutC); free(edges); free(nC);

   free(dfs); free(idfs); free(p); free(low);
   free(l); free(bcount); free(low1C); free(low2C);
   free(M); free(nextM); free(prevM);
   free(isCutEdge); free(isCutEdgeLow);
   free(C); free(Q); free(map); free(imapStack); free(imapFirst); free(imapNext);

   return c_num;
}

int get_4e_components_3econnected_high(int n, int* adj, int* firstOut, int** C, int* n3cuts, int** cuts3)
{
   int* cuts = (int*)malloc(sizeof(int)*n*12);
   int* dfs; int* idfs; int* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   int* low; int* lowD; int* low2; int* low2D;
   get_2low(n,adj,firstOut,dfs,idfs,p,&low,&lowD,&low2,&low2D);
   int* low1C; int* low2C;
   get_lowChildren(n,dfs,idfs,p,low,&low1C,&low2C);
   int* l; int* bcount;
   get_l_and_bcount(n,adj,firstOut,dfs,idfs,p,&l,&bcount);
   int* M; int* nextM;
   get_M(n,dfs,idfs,l,low,low1C,low2C,&M,&nextM);
   int* high; int* highD;
   get_high(n,adj,firstOut,dfs,idfs,p,&high,&highD);
   int* Ml; int* Mlow1; int* Mlow2;
   get_allM(n,dfs,idfs,l,low,M,low1C,low2C,&Ml,&Mlow1,&Mlow2);
   int* currentVertex = (int*)malloc(sizeof(int)*n);

   int* ND = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){ND[i]=1;}
   for(int i=n-1;i>0;i--){int v=idfs[i];ND[p[v]]+=ND[v];}

   int* size = (int*)malloc(sizeof(int)*2*n);
   char* type = (char*)malloc(sizeof(char)*2*n);
   int k=0;

   //one tree-edge
   for(int v=1;v<n;v++)
   {
      if(bcount[v]==2)
      {
         cuts[6*k+0]=v; cuts[6*k+1]=p[v];
         cuts[6*k+2]=low[v]; cuts[6*k+3]=lowD[v];
         cuts[6*k+4]=low2[v]; cuts[6*k+5]=low2D[v];
         type[k]=0; size[k]=ND[v];
         k++;      
      }
   }

   //two tree-edges
   //upper case
   for(int i=0;i<n;i++){currentVertex[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      int m=Ml[v];
      if(m==-1){continue;} 
      int u=currentVertex[m];
      while(nextM[u]!=-1 && dfs[nextM[u]]>dfs[v]){u=nextM[u];}
      currentVertex[m]=u;
      if(dfs[high[u]]<dfs[v] && bcount[v]==bcount[u]+1)
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=M[v]; cuts[6*k+5]=l[M[v]];
         type[k]=1; size[k]=ND[v]-ND[u];
         k++;
      }
   }
   for(int i=0;i<n;i++){currentVertex[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      int m=Mlow1[v];
      if(m==-1){continue;} 
      int u=currentVertex[m];
      while(nextM[u]!=-1 && dfs[nextM[u]]>dfs[v]){u=nextM[u];}
      currentVertex[m]=u;
      if(dfs[high[u]]<dfs[v] && bcount[v]==bcount[u]+1)
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=Mlow2[v]; cuts[6*k+5]=l[Mlow2[v]];
         type[k]=1; size[k]=ND[v]-ND[u];
         k++;
      }
   }
   for(int i=0;i<n;i++){currentVertex[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      int m=Mlow2[v];
      if(m==-1){continue;} 
      int u=currentVertex[m];
      while(nextM[u]!=-1 && dfs[nextM[u]]>dfs[v]){u=nextM[u];}
      currentVertex[m]=u;
      if(dfs[high[u]]<dfs[v] && bcount[v]==bcount[u]+1)
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=Mlow1[v]; cuts[6*k+5]=l[Mlow1[v]];
         type[k]=1; size[k]=ND[v]-ND[u];
         k++;
      }
   }

   //lower case
   for(int u=1;u<n;u++)
   {
      if(nextM[u]!=-1 && bcount[nextM[u]]==bcount[u]-1)
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=nextM[u]; cuts[6*k+3]=p[nextM[u]];
         cuts[6*k+4]=highD[u]; cuts[6*k+5]=high[u];
         type[k]=1; size[k]=ND[nextM[u]]-ND[u];
         k++;
      }
   }
   for(int i=0;i<n;i++){currentVertex[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int u=idfs[i];
      int m=Ml[u];
      if(m==-1){continue;} 
      int v=currentVertex[m];
      while(v!=-1 && dfs[v]>=dfs[u]){v=nextM[v];}
      currentVertex[m]=v;
      if(bcount[u]==bcount[v]+1 && M[u]!=M[v])
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=highD[u]; cuts[6*k+5]=high[u];
         type[k]=1; size[k]=ND[v]-ND[u];
         k++;
      }
   }
   for(int i=0;i<n;i++){currentVertex[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int u=idfs[i];
      int m=Mlow1[u];
      if(m==-1){continue;} 
      int v=currentVertex[m];
      while(v!=-1 && dfs[v]>=dfs[u]){v=nextM[v];}
      currentVertex[m]=v;
      if(bcount[u]==bcount[v]+1)
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=highD[u]; cuts[6*k+5]=high[u];
         type[k]=1; size[k]=ND[v]-ND[u];
         k++;
      }
   }


   //three tree-edges
   //u and v are not related as ancestor and descendant
   for(int i=0;i<n;i++){currentVertex[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int w=idfs[i];
      int m1=Mlow1[w]; int m2=Mlow2[w];
      if(m1==-1||m2==-1){continue;}
      int u=currentVertex[m1];
      while(nextM[u]!=-1 && dfs[nextM[u]]>dfs[w]){u=nextM[u];}
      currentVertex[m1]=u;
      int v=currentVertex[m2];
      while(nextM[v]!=-1 && dfs[nextM[v]]>dfs[w]){v=nextM[v];}
      currentVertex[m2]=v;
      if(bcount[w]==bcount[u]+bcount[v] && dfs[high[u]]<dfs[w] && dfs[high[v]]<dfs[w])
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=w; cuts[6*k+5]=p[w];
         type[k]=2; size[k]=ND[w]-ND[u]-ND[v];
         k++;
      }
   }

   int* lastM = (int*)malloc(sizeof(int)*n);
   for(int v=1;v<n;v++){if(nextM[v]==-1){lastM[M[v]]=v;}}   

   //u,v,w M[v]!=M[w]
   for(int i=0;i<n;i++){currentVertex[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      int m1=Mlow1[v]; int m2=Mlow2[v];
      if(m1==-1||m2==-1){continue;}
      int w=currentVertex[m1];
      while(w!=-1 && dfs[w]>=dfs[v]){w=nextM[w];}
      currentVertex[m1]=w;
      int u=lastM[m2];
      if(w!=-1 && dfs[high[u]]<dfs[v] && bcount[v]==bcount[u]+bcount[w])
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=w; cuts[6*k+5]=p[w];
         type[k]=2; size[k]=ND[w]-ND[v]+ND[u];
         k++;
      }
   }

   //w=nextM[v]
   int* A = (int*)malloc(sizeof(int)*firstOut[n]);
   for(int i=0;i<firstOut[n];i++){A[i]=-1;}
   int* Hstack = (int*)malloc(sizeof(int)*n);
   int* firstH = (int*)malloc(sizeof(int)*n);
   int* nextH = (int*)malloc(sizeof(int)*n); 
   for(int i=0;i<n;i++){firstH[i]=-1; nextH[i]=-1;}
   int SP=0;
   for(int i=1;i<n;i++)
   {
      int v=idfs[i];
      int h=high[v];
      nextH[SP]=firstH[h];
      firstH[h]=SP;
      Hstack[SP++]=v;
   }  
   int* stack = (int*)malloc(sizeof(int)*n);
   for(int h=0;h<n;h++)
   {
      int uIndx = firstH[h];
      SP=0;
      while(uIndx!=-1)
      {
         int zIndx = nextH[uIndx];
         if(zIndx==-1){break;}
         int u=Hstack[uIndx]; int z=Hstack[zIndx];
         if(nextM[u]==-1)
         {
            stack[SP++]=u;
            A[bcount[u]]=u; 
         }
         if(!(dfs[z]<=dfs[u] && dfs[u]<dfs[z]+ND[z]))
         {
            while(SP!=0)
            {
               int uTemp = stack[SP-1];
               A[bcount[uTemp]]=-1;
               SP--;
            } 
         }
         if(nextM[z]!=-1)
         {
            int v=z; int w=nextM[v];
            if(A[bcount[v]-bcount[w]]!=-1)
            {
               u=A[bcount[v]-bcount[w]];
               if(dfs[low[u]]>=dfs[w])
               {
                  cuts[6*k+0]=u; cuts[6*k+1]=p[u];
                  cuts[6*k+2]=v; cuts[6*k+3]=p[v];
                  cuts[6*k+4]=w; cuts[6*k+5]=p[w];
                  type[k]=2; size[k]=ND[w]-ND[v]+ND[u];
                  k++;              
               }
            }
         }
         uIndx=zIndx;
      }
      while(SP!=0)
      {
         int uTemp = stack[SP-1];
         A[bcount[uTemp]]=-1; 
         SP--;
      } 
   }

   //w!=nextM[v]
   int* Ustack = (int*)malloc(sizeof(int)*n);
   int* firstU = (int*)malloc(sizeof(int)*n);
   int* nextU = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){firstU[i]=-1; nextU[i]=-1;}
   int USP=0;
   for(int h=0;h<n;h++)
   {
      int uIndx=firstH[h];
      SP=0;
      while(uIndx!=-1)
      {
         int zIndx=nextH[uIndx];
         if(zIndx==-1){break;}
         int u=Hstack[uIndx]; int z=Hstack[zIndx];
         if(nextM[u]==-1){stack[SP++]=u;}
         if(!(dfs[z]<=dfs[u] && dfs[u]<dfs[z]+ND[z])){SP=0;}
         if(nextM[z]!=-1)
         {
            int v=z;
            while(SP!=0 && dfs[low[stack[SP-1]]]<dfs[lastM[M[v]]]){SP--;}
            while(SP!=0 && dfs[low[stack[SP-1]]]<dfs[nextM[v]])
            {
               SP--;
               int u=stack[SP];
               nextU[USP]=firstU[v];
               firstU[v]=USP;
               Ustack[USP++]=u;  
            }  
         }
         uIndx=zIndx;
      }
   }

   int* lowestW = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){lowestW[i]=nextM[i];}
   for(int i=1;i<n;i++)
   {
      int v=idfs[i];
      int indx=firstU[v];
      while(indx!=-1)
      {
         int u=Ustack[indx];
         int w=lowestW[v];
         while(dfs[w]>dfs[low[u]]){w=lowestW[w];}
         lowestW[v]=w;
         if(bcount[v]==bcount[u]+bcount[w])
         {
            cuts[6*k+0]=u; cuts[6*k+1]=p[u];
            cuts[6*k+2]=v; cuts[6*k+3]=p[v];
            cuts[6*k+4]=w; cuts[6*k+5]=p[w];
            type[k]=2; size[k]=ND[w]-ND[v]+ND[u];
            k++;  
         }
         indx=nextU[indx];
      }
   }

   int* sizeStack = (int*)malloc(sizeof(int)*k);
   int* sizeFirst = (int*)malloc(sizeof(int)*n);
   int* sizeNext = (int*)malloc(sizeof(int)*k);
   for(int i=0;i<n;i++){sizeFirst[i]=-1;}
   SP=0;
   for(int i=0;i<k;i++)
   {
      int s=size[i];
      sizeNext[SP]=sizeFirst[s]; sizeFirst[s]=SP; sizeStack[SP++]=i;
   } 
   int* cutStack = (int*)malloc(sizeof(int)*3*k);
   int* cutFirst = (int*)malloc(sizeof(int)*n);
   int* cutNext = (int*)malloc(sizeof(int)*3*k);
   for(int i=0;i<n;i++){cutFirst[i]=-1;}
   SP=0;
   for(int s=0;s<n;s++)
   {
      for(int t=sizeFirst[s];t!=-1;t=sizeNext[t])
      {
         int i=sizeStack[t];
         int u=-1; int v=-1; int w=-1;
         u=cuts[6*i+0];
         if(type[i]>0){v=cuts[6*i+2];}
         if(type[i]>1){w=cuts[6*i+4];}
         cutNext[SP]=cutFirst[u]; cutFirst[u]=SP; cutStack[SP++]=i;
         if(v!=-1)
         {
            cutNext[SP]=cutFirst[v]; cutFirst[v]=SP; cutStack[SP++]=i;
         }
         if(w!=-1)
         {
            cutNext[SP]=cutFirst[w]; cutFirst[w]=SP; cutStack[SP++]=i;
         }
      }
   }

   int* cactusParent = (int*)malloc(sizeof(int)*3*n);
   int* parentCut = (int*)malloc(sizeof(int)*3*n);
   char* marked = (char*)malloc(sizeof(char)*k);
   for(int i=0;i<k;i++){marked[i]=0;}
   int* phi = (int*)malloc(sizeof(int)*k);
   for(int i=0;i<k;i++){phi[i]=-1;}

   int* Cstack = (int*)malloc(sizeof(int)*3*n);
   int* cactusNode = (int*)malloc(sizeof(int)*n);
   (*C) = (int*)malloc(sizeof(int)*n);
   (*C)[0]=0; cactusNode[0]=0;
   int n_comp=1; int Nr=1;
   for(int i=1;i<n;i++)
   {
      int v=idfs[i];
      int x=cactusNode[p[v]];
      SP=0;
      while(x!=0)
      {
         int c=parentCut[x];
         char found=0;
         if(type[c]==0&&v==cuts[6*c+0])
         {
            marked[c]=1; Cstack[SP++]=c; found=1;
         }
         else if(type[c]==1&&(v==cuts[6*c+0]||v==cuts[6*c+2]))
         {
            marked[c]=1; Cstack[SP++]=c; found=1;
         }
         else if(type[c]==2&&(v==cuts[6*c+0]||v==cuts[6*c+2]||v==cuts[6*c+4]))
         {
            marked[c]=1; Cstack[SP++]=c; found=1;
         }
         if(!found){break;}
         x=cactusParent[x];
      }
      char new_comp=0;
      for(int t=cutFirst[v];t!=-1;t=cutNext[t])
      {
         int c=cutStack[t];
         if(marked[c]){continue;}
         if(phi[c]==-1)
         {
            phi[c]=Nr; parentCut[Nr]=c; cactusParent[Nr++]=x; new_comp=1;
         }
         x=phi[c];
      } 
      cactusNode[v]=x;
      while(SP!=0){marked[Cstack[SP-1]]=0;SP--;}
      if(new_comp){x=n_comp++;}
      (*C)[v]=x;
   }

   free(sizeStack); free(sizeFirst); free(sizeNext);
   free(cutStack); free(cutFirst); free(cutNext); free(Cstack);
   free(cactusParent); free(parentCut); free(marked); free(phi); free(cactusNode);

   free(dfs); free(idfs); free(p);
   free(low); free(lowD); free(low2); free(low2D);
   free(low1C); free(low2C); free(l); free(bcount);
   free(M); free(nextM); free(high); free(highD);
   free(Ml); free(Mlow1); free(Mlow2); free(currentVertex); 
   free(lastM); free(A); free(ND); 
   free(Hstack); free(firstH); free(nextH); free(stack);
   free(Ustack); free(firstU); free(nextU);
   free(lowestW);
   free(size); free(type);
   //free(cuts);
   *n3cuts=k; *cuts3=cuts;
   return n_comp;
}

void get_2low(int n, int* adj, int* firstOut, int* dfs, int* idfs, int* p, int** low1, int** low1D, int** low2, int** low2D)
{
   (*low1) = (int*)malloc(sizeof(int)*n);
   (*low1D) = (int*)malloc(sizeof(int)*n);
   (*low2) = (int*)malloc(sizeof(int)*n);
   (*low2D) = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){(*low1)[i]=i; (*low2)[i]=i;}
   char* foundP = (char*)malloc(sizeof(char)*n);
   char* foundC = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){foundP[i]=0;foundC[i]=0;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      for(int j=firstOut[v];j<firstOut[v+1];j++)
      {
         int u=adj[j];
         if(u==p[v]&&!foundP[v]){foundP[v]=1;continue;}
         if(dfs[u]<dfs[v])
         {
            if(dfs[u]<=dfs[(*low1)[v]])
            { 
               (*low2)[v]=(*low1)[v]; (*low2D)[v]=(*low1D)[v];
               (*low1)[v]=u; (*low1D)[v]=v;
            }
            else if(dfs[u]<dfs[(*low2)[v]])
            {
               (*low2)[v]=u; (*low2D)[v]=v;
            }  
         }
         else if(v==p[u]&&!foundC[u])
         {
            if(dfs[(*low2)[u]]<=dfs[(*low1)[v]])
            {
               (*low1)[v]=(*low1)[u]; (*low1D)[v]=(*low1D)[u];    
               (*low2)[v]=(*low2)[u]; (*low2D)[v]=(*low2D)[u];    
            }
            else if(dfs[(*low1)[u]]<=dfs[(*low1)[v]])
            {
               if(dfs[(*low1)[v]]<dfs[(*low2)[v]])
               {
                  (*low2)[v]=(*low1)[v]; (*low2D)[v]=(*low1D)[v];
               }
               else if(dfs[(*low2)[u]]<dfs[(*low2)[v]])
               {
                  (*low2)[v]=(*low2)[u]; (*low2D)[u]=(*low2D)[u];
               }
               (*low1)[v]=(*low1)[u]; (*low1D)[v]=(*low1D)[u];
            }
            else if(dfs[(*low1)[u]]<dfs[(*low2)[v]])
            {
               (*low2)[v]=(*low1)[u]; (*low2D)[v]=(*low1D)[u];
            }
            foundC[u]=1;
         }
      }
   }
   free(foundP); free(foundC);
}

void get_high(int n, int* adj, int* firstOut, int* dfs, int* idfs, int* p, int** high, int** highD)
{
   (*high) = (int*)malloc(sizeof(int)*n);
   (*highD) = (int*)malloc(sizeof(int)*n);
   int* P = (int*)malloc(sizeof(int)*n);
   int* size = (int*)malloc(sizeof(int)*n);
   int* repr = (int*)malloc(sizeof(int)*n);
   char* found_p = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){(*high)[i]=-1; P[i]=i; size[i]=1; repr[i]=i; found_p[i]=0;}
   for(int i=n-1;i>=0;i--)
   {
      int y=idfs[i];
      for(int j=firstOut[y];j<firstOut[y+1];j++)
      {
         int x=adj[j];
         if(p[x]==y&&!found_p[x]){found_p[x]=1;continue;}
         if(dfs[x]<dfs[y]){continue;}
         int u=repr[find(P,x)];
         while(u!=y)
         {
            (*high)[u]=y; (*highD)[u]=x;
            int next=repr[find(P,p[u])];
            unite(P,size,u,p[u]);
            repr[find(P,u)]=y;
            u=next;
         } 
      }
   }
   free(P); free(size); free(repr); free(found_p);
}

int find(int* p, int x)
{
   int r=x;
   while(p[r]!=r){r=p[r];}
   while(x!=r){int next=p[x];p[x]=r;x=next;}
   return r;
}
void unite(int* p, int* size, int x, int y)
{
   int r1=find(p,x);
   int r2=find(p,y);
   if(size[r1]<size[r2]){p[r1]=r2;size[r2]+=size[r1];}
   else{p[r2]=r1;size[r1]+=size[r2];}
}

void get_allM(int n, int* dfs, int* idfs, int* l, int* low, int* M, int* low1C, int* low2C, int** Ml, int** Mlow1, int** Mlow2)
{
   (*Ml) = (int*)malloc(sizeof(int)*n);
   (*Mlow1) = (int*)malloc(sizeof(int)*n);
   (*Mlow2) = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++)
   {
      (*Ml)[i]=-1; (*Mlow1)[i]=-1; (*Mlow2)[i]=-1;
   }

   int* currentM = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){currentM[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      int m=M[v];
      if(dfs[l[m]]>=dfs[v]){continue;}
      if(low1C[m]==-1 || dfs[low[low1C[m]]]>=dfs[v]){continue;}
      if(low2C[m]!=-1 && dfs[low[low2C[m]]]<dfs[v])
      {
         (*Ml)[v]=m; 
         continue;
      }
      int tempM=low1C[m];
      m=currentM[low1C[m]];
      while(1)
      {
         if(dfs[l[m]]<dfs[v]){break;}
         if(low2C[m]!=-1 && dfs[low[low2C[m]]]<dfs[v]){break;}
         m=currentM[low1C[m]];
      }
      (*Ml)[v]=m; 
      currentM[tempM]=m;
   }

   for(int i=0;i<n;i++){currentM[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      int m=M[v];
      if(dfs[l[m]]<dfs[v]){continue;}
      if(low1C[m]==-1 || dfs[low[low1C[m]]]>=dfs[v]){continue;}
      int tempM=low1C[m];
      m=currentM[low1C[m]];
      while(1)
      {
         if(dfs[l[m]]<dfs[v]){break;}
         if(low2C[m]!=-1 && dfs[low[low2C[m]]]<dfs[v]){break;}
         m=currentM[low1C[m]];
      }
      (*Mlow1)[v]=m; 
      currentM[tempM]=m;
   }

   for(int i=0;i<n;i++){currentM[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      int m=M[v];
      if(dfs[l[m]]<dfs[v]){continue;}
      if(low2C[m]==-1 || dfs[low[low2C[m]]]>=dfs[v]){continue;}
      int tempM=low2C[m];
      m=currentM[low2C[m]];
      while(1)
      {
         if(dfs[l[m]]<dfs[v]){break;}
         if(low2C[m]!=-1 && dfs[low[low2C[m]]]<dfs[v]){break;}
         m=currentM[low1C[m]];
      }
      (*Mlow2)[v]=m; 
      currentM[tempM]=m;
   }

   free(currentM);
}

void get_M(int n, int* dfs, int* idfs, int* l, int* low, int* low1C, int* low2C, int** M, int** nextM)
{
   (*M) = (int*)malloc(sizeof(int)*n);
   (*nextM) = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){(*M)[i]=-1; (*nextM)[i]=-1;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      int c=v; int m=v;
      while(1)
      {
         if(dfs[l[m]]<i){(*M)[v]=m;break;}
         if(low2C[m]!=-1 && dfs[low[low2C[m]]]<i){(*M)[v]=m;break;}
         c=low1C[m]; m=(*M)[c]; 
      }
      if(c!=v)
      {
         (*nextM)[c]=v;
      }
   } 
}

void get_lowChildren(int n, int* dfs, int* idfs, int* p, int* low, int** low1C, int** low2C)
{
   (*low1C) = (int*)malloc(sizeof(int)*n);
   (*low2C) = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){(*low1C)[i]=-1; (*low2C)[i]=-1;}
   for(int i=1;i<n;i++)
   {
      int x=idfs[i];
      int y=p[x];
      if((*low1C)[y]==-1){(*low1C)[y]=x;}
      else if(dfs[low[x]]<dfs[low[(*low1C)[y]]]){(*low1C)[y]=x;}
   }
   for(int i=1;i<n;i++)
   {
      int x=idfs[i];
      int y=p[x];
      if(x!=(*low1C)[y])
      {
        if((*low2C)[y]==-1){(*low2C)[y]=x;}
        else if(dfs[low[x]]<dfs[low[(*low2C)[y]]]){(*low2C)[y]=x;}
      }
   }
}

void get_l_and_bcount(int n, int* adj, int* firstOut, int* dfs, int* idfs, int* p, int** l, int** bcount)
{
   (*l) = (int*)malloc(sizeof(int)*n);
   (*bcount) = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){(*l)[i]=i; (*bcount)[i]=0;}
   char* found_p = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){found_p[i]=0;}
   char* found_c = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){found_c[i]=0;}
   for(int i=n-1;i>0;i--)
   {
      int x=idfs[i];
      for(int t=firstOut[x];t<firstOut[x+1];t++)
      {
         int y=adj[t];
         if(dfs[y]<dfs[x])
         {
            if(y==p[x]&&!found_p[x]) 
            {
               found_p[x]=1;
            }  
            else
            {
               if(dfs[y]<dfs[(*l)[x]]){(*l)[x]=y;}
               (*bcount)[x]++; (*bcount)[y]--;
            }
         }
         else if(x==p[y]&&!found_c[y])
         {
            (*bcount)[x]+=(*bcount)[y];
            found_c[y]=1;     
         }
      }
   }
   free(found_p); free(found_c);
}

void get_low(int n, int* adj, int* firstOut, int* dfs, int* idfs, int* p, int** low)
{
   *low = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){(*low)[i]=i;}
   char* found_p = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){found_p[i]=0;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      for(int t=firstOut[v];t<firstOut[v+1];t++)
      {
         int u=adj[t];
         if(dfs[u]<dfs[v])
         {
            if(u==p[v] && !found_p[v])
            {
               found_p[v]=1;
            } 
            else
            {
               if(dfs[u]<dfs[(*low)[v]]){(*low)[v]=u;}
            }
         } 
         else if(v==p[u])
         {
            if(dfs[(*low)[u]]<dfs[(*low)[v]]){(*low)[v]=(*low)[u];}
         }
      }
   }
   free(found_p);
}

void DFS(int n, int* adj, int* firstOut, int** dfs, int** idfs, int** p)
{
   *dfs = (int*)malloc(sizeof(int)*n);
   *idfs = (int*)malloc(sizeof(int)*n);
   *p = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){(*dfs)[i]=-1; (*p)[i]=-1;}
   int* temp_vertex = (int*)malloc(sizeof(int)*n);
   int* temp_out = (int*)malloc(sizeof(int)*n);
   
   int nr=0;
   (*dfs)[0]=nr; (*idfs)[nr++]=0;
   temp_vertex[0]=0; temp_out[0]=firstOut[0];
   int SP=0;
   while(SP!=-1)
   {
      int v=temp_vertex[SP];
      char descend=0;
      for(int i=temp_out[SP];i<firstOut[v+1];i++)
      {
         int u=adj[i];
         if((*dfs)[u]==-1)
         {
            (*dfs)[u]=nr; (*idfs)[nr++]=u; (*p)[u]=v;
            temp_vertex[SP+1]=u; temp_out[SP+1]=firstOut[u]; temp_out[SP]=i+1;
            descend=1; break;
         }
      }
      if(descend){SP++;continue;}
      SP--;
   }

   free(temp_vertex); free(temp_out);
}

void get_adj(int n, int m, int* edges, int** adj, int** firstOut)
{
   *adj = (int*)malloc(sizeof(int)*2*m);
   *firstOut = (int*)malloc(sizeof(int)*(n+1));
   for(int i=0;i<=n;i++){(*firstOut)[i]=0;}
   for(int i=0;i<m;i++){(*firstOut)[edges[2*i]+1]++; (*firstOut)[edges[2*i+1]+1]++;}
   for(int i=1;i<=n;i++){(*firstOut)[i]+=(*firstOut)[i-1];}
   int* nextOut = (int*)malloc(sizeof(int)*(n+1));
   for(int i=0;i<=n;i++){nextOut[i]=(*firstOut)[i];}
   for(int i=0;i<m;i++)
   {
      int x=edges[2*i]; int y=edges[2*i+1];
      (*adj)[nextOut[x]++]=y; (*adj)[nextOut[y]++]=x;
   } 
   free(nextOut);
}

void read_graph(char* filename, int* n, int** adj, int** firstOut)
{
   int m;
   FILE* fp = fopen(filename,"r");
   fscanf(fp,"%d %d",n,&m);
   int* edges = (int*)malloc(sizeof(int)*2*m);
   for(int i=0;i<m;i++){fscanf(fp,"%d %d",edges+2*i,edges+2*i+1);}
   fclose(fp);
   get_adj(*n,m,edges,adj,firstOut);
   free(edges);
}
