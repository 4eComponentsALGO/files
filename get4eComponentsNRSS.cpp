//Runs with ./get4eComponentsNRSS <input_graph> <output_4components>

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

int get_4e_components_rand(int,int*,int*,int**);
int get_4e_components_connected_rand(int,int*,int*,int**);
int get_4e_components_2econnected_rand(int,int*,int*,int**);
int get_4e_components_3econnected_rand(int,int*,int*,int**,int*,int**);
int get_3cuts_rand(int,int*,int*,int**);
int get_3cuts_2tree_rand(int,int*,int*,int**);
void get_2low(int,int*,int*,int*,int*,int*,int**,int**,int**,int**);
void get_LandR(int,int*,int*,int*,int*,int*,int**,int**,int**,int**,int**,int**);
void get_high_indx(int,int*,int*,int*,int*,int*,int**,int**,int**);
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
   int k=get_4e_components_rand(n,adj,firstOut,&C);
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

int get_4e_components_rand(int n, int* adj, int* firstOut, int** C)
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
         int k2 = get_4e_components_connected_rand(nC,adjC,firstOutC,&C2);
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

int get_4e_components_connected_rand(int n, int* adj, int* firstOut, int** C)
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
         int k3=get_4e_components_2econnected_rand(nC,adjC,firstOutC,&C3);
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

int get_4e_components_2econnected_rand(int n, int* adj, int* firstOut, int** C3)
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
         int k4=get_4e_components_3econnected_rand(nC[c],adjC[c],firstOutC[c],&C4,nCcuts+c,Ccuts+c);
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


int get_4e_components_3econnected_rand(int n, int* adj, int* firstOut, int** C, int* n3cuts, int** cuts3)
{
   int* cuts;
   int k = get_3cuts_rand(n,adj,firstOut,&cuts);
(*cuts3) = (int*)malloc(sizeof(int)*6*k);
*n3cuts=k; for(int i=0;i<6*k;i++){(*cuts3)[i]=cuts[i];}
   int* dfs; int* idfs; int* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   int* ND = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){ND[i]=1;}
   for(int i=n-1;i>0;i--){int v=idfs[i];ND[p[v]]+=ND[v];}
   char* type = (char*)malloc(sizeof(char)*k);
   int* size = (int*)malloc(sizeof(int)*k);
   for(int i=0;i<k;i++)
   {
      int temp;
      if(dfs[cuts[6*i+0]]<dfs[cuts[6*i+1]]){temp=cuts[6*i+0];cuts[6*i+0]=cuts[6*i+1];cuts[6*i+1]=temp;}
      if(dfs[cuts[6*i+2]]<dfs[cuts[6*i+3]]){temp=cuts[6*i+2];cuts[6*i+2]=cuts[6*i+3];cuts[6*i+3]=temp;}
      if(dfs[cuts[6*i+4]]<dfs[cuts[6*i+5]]){temp=cuts[6*i+4];cuts[6*i+4]=cuts[6*i+5];cuts[6*i+5]=temp;}
      
      for(int t=0;t<3;t++)
      {
         for(int s=t+1;s<3;s++)
         {
            if(dfs[cuts[6*i+2*t]]<dfs[cuts[6*i+2*s]])
            {
               temp=cuts[6*i+2*t]; cuts[6*i+2*t]=cuts[6*i+2*s]; cuts[6*i+2*s]=temp;
               temp=cuts[6*i+2*t+1]; cuts[6*i+2*t+1]=cuts[6*i+2*s+1]; cuts[6*i+2*s+1]=temp;
            }  
         }
      }
   }
   for(int i=0;i<k;i++)
   {
      int x1=cuts[6*i+0]; int y1=cuts[6*i+1];  
      int x2=cuts[6*i+2]; int y2=cuts[6*i+3];  
      int x3=cuts[6*i+4]; int y3=cuts[6*i+5];
      type[i]=0;
      if(y1==p[x1]){type[i]++;}
      if(y2==p[x2]){type[i]++;}
      if(y3==p[x3]){type[i]++;}
      if(type[i]==3) 
      {
         if(x1==x2&&x1==x3){type[i]=0;}
         else if(x1==x2){type[i]=1;x2=x3;cuts[6*i+2]=cuts[6*i+4];}
         else if(x2==x3){type[i]=1;}
         else{type[i]=2;}
      }
      else if(type[i]==2)
      {
         if(y1==p[x1]&&y2==p[x2])
         {
            if(x1==x2){type[i]=0;}
            else{type[i]=1;}
         }
         else if(y2==p[x2]&&y3==p[x3])
         {
            if(x2==x3){type[i]=0;x1=x2;cuts[6*i]=cuts[6*i+2];}
            else{type[i]=1;x1=x2;x2=x3;cuts[6*i]=cuts[6*i+2];cuts[6*i+2]=cuts[6*i+4];}
         }
         else if(y1==p[x1]&&y3==p[x3])
         {
            type[i]=1; x2=x3; cuts[6*i+2]=cuts[6*i+4];  
         } 
      }
      else if(type[i]==1)
      {
         type[i]=0;
         if(y2==p[x2]){x1=x2;cuts[6*i]=cuts[6*i+2];}
         else if(y3==p[x3]){x1=x3;cuts[6*i]=cuts[6*i+4];}
      } 

      if(type[i]==0){size[i]=ND[x1];}
      else if(type[i]==1){size[i]=ND[x2]-ND[x1];}
      else if(type[i]==2)
      {
         if(dfs[x2]<dfs[x1] && dfs[x1]<dfs[x2]+ND[x2])
         {
            size[i]=ND[x3]-ND[x2]+ND[x1]; 
         }
         else
         {
            size[i]=ND[x3]-ND[x2]-ND[x1];
         }
      }
   }

/*print_tree(n,adj,firstOut,dfs,p);
for(int i=0;i<k;i++){printf("*(%d,%d) (%d,%d) (%d,%d) %d %d\n",cuts[6*i+0],cuts[6*i+1],cuts[6*i+2],cuts[6*i+3],cuts[6*i+4],cuts[6*i+5],type[i],size[i]);}
printf("\n");*/

   int* sizeStack = (int*)malloc(sizeof(int)*k);
   int* sizeFirst = (int*)malloc(sizeof(int)*n);
   int* sizeNext = (int*)malloc(sizeof(int)*k);
   for(int i=0;i<n;i++){sizeFirst[i]=-1;}
   int SP=0;
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

   int* Cstack = (int*)malloc(sizeof(int)*3*n);
   int* cactusParent = (int*)malloc(sizeof(int)*3*n);
   int* parentCut = (int*)malloc(sizeof(int)*3*n);
   char* marked = (char*)malloc(sizeof(char)*k);
   for(int i=0;i<k;i++){marked[i]=0;}
   int* phi = (int*)malloc(sizeof(int)*k);
   for(int i=0;i<k;i++){phi[i]=-1;}

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

/*sort_cuts(k,cuts);
for(int i=0;i<k;i++){printf("r(%d,%d) (%d,%d) (%d,%d)\n",cuts[6*i+0],cuts[6*i+1],cuts[6*i+2],cuts[6*i+3],cuts[6*i+4],cuts[6*i+5]);}
for(int i=0;i<n;i++){printf("%d ",(*C)[i]);}printf("\n");*/

/*for(int x=0;x<n;x++)
{
   printf("%d: ",x);for(int i=firstOut[x];i<firstOut[x+1];i++){printf("%d ",adj[i]);}printf("\n");
}
printf("%d\n",n_comp);
for(int i=0;i<n;i++){printf("%d ",(*C)[i]);}printf("\n");*/

   free(sizeStack); free(sizeFirst); free(sizeNext);
   free(cutStack); free(cutFirst); free(cutNext); free(Cstack);
   free(cactusParent); free(parentCut); free(marked); free(phi); free(cactusNode);

   free(dfs); free(idfs); free(p); free(ND);
   free(cuts); free(type); free(size);
   return n_comp;
}

int get_3cuts_rand(int n, int* adj, int* firstOut, int** cuts)
{
   (*cuts) = (int*)malloc(sizeof(int)*n*12);
   int k=0; int k2;
   int* cuts_2tree;

   k2 = get_3cuts_2tree_rand(n,adj,firstOut,&cuts_2tree);

   for(int i=0;i<k2;i++)
   {
      for(int t=0;t<6;t++){(*cuts)[6*k+t]=cuts_2tree[6*i+t];}k++;
   }
   free(cuts_2tree);
   
   int* dfs; int* idfs; int* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   char* foundP = (char*)malloc(sizeof(char)*n);
   char* foundC = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){foundP[i]=0;foundC[i]=0;}
   int* C = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){C[i]=-1;}
   int* Q = (int*)malloc(sizeof(int)*n);
   int cIndx=0;
   for(int r=0;r<n;r++)
   {
      if(C[r]!=-1){continue;}
      int first=0; int last=0;
      Q[last++]=r; C[r]=cIndx;
      while(first!=last)
      {
         int v=Q[first++];
         for(int i=firstOut[v];i<firstOut[v+1];i++)
         {
            int u=adj[i];
            if(v==p[u]&&!foundC[u]){foundC[u]=1;continue;}
            if(u==p[v]&&!foundP[v]){foundP[v]=1;continue;}
            if(C[u]==-1)
            {
               Q[last++]=u; C[u]=cIndx;
            }
         }
      }
      cIndx++;
   }

   if(cIndx==1){free(dfs); free(idfs); free(p); free(foundP); free(C); free(Q); return k;}

   int* stackT_ = (int*)malloc(sizeof(int)*2*n);
   int* firstT_ = (int*)malloc(sizeof(int)*n);
   int* nextT_ = (int*)malloc(sizeof(int)*2*n);
   int* tree_edge_ = (int*)malloc(sizeof(int)*4*n);
   int SP=0;
   for(int i=0;i<n;i++){firstT_[i]=-1;}
   for(int v=1;v<n;v++)
   {
      int u=p[v];
      if(C[u]==C[v]){continue;}
      nextT_[SP]=firstT_[C[u]];
      firstT_[C[u]]=SP;
      stackT_[SP]=C[v];
      tree_edge_[2*SP]=u; tree_edge_[2*SP+1]=v; SP++;
      nextT_[SP]=firstT_[C[v]];
      firstT_[C[v]]=SP;
      stackT_[SP]=C[u];
      tree_edge_[2*SP]=u; tree_edge_[2*SP+1]=v; SP++;
   }

   int* stackT = (int*)malloc(sizeof(int)*2*n);
   int* firstT = (int*)malloc(sizeof(int)*n);
   int* nextT = (int*)malloc(sizeof(int)*2*n);
   int* tree_edge = (int*)malloc(sizeof(int)*4*n);
   SP=0;
   for(int i=0;i<n;i++){firstT[i]=-1;}
   for(int c=cIndx-1;c>=0;c--)
   {
      for(int s=firstT_[c];s!=-1;s=nextT_[s])
      {
         int d=stackT_[s];
         nextT[SP]=firstT[d];
         firstT[d]=SP;
         stackT[SP]=c;
         tree_edge[2*SP]=tree_edge_[2*s]; tree_edge[2*SP+1]=tree_edge_[2*s+1]; SP++;
      }
   }  

   int* adj_new = (int*)malloc(sizeof(int)*2*n);
   int* firstOut_new = (int*)malloc(sizeof(int)*(cIndx+1));
   for(int i=0;i<=cIndx;i++){firstOut_new[i]=0;}
   for(int v=1;v<n;v++)
   {
      int u=p[v];
      if(C[u]==C[v]){continue;}
      firstOut_new[C[u]+1]++; firstOut_new[C[v]+1]++;
   }
   int* currentOut = (int*)malloc(sizeof(int)*(cIndx+1));
   currentOut[0]=0;
   for(int i=1;i<=cIndx;i++){firstOut_new[i]+=firstOut_new[i-1]; currentOut[i]=firstOut_new[i];}
   for(int v=1;v<n;v++)
   {
      int u=p[v];
      if(C[u]==C[v]){continue;}
      adj_new[currentOut[C[v]]++]=C[u]; adj_new[currentOut[C[u]]++]=C[v];
   }

   int* cuts_new;
   k2 = get_3cuts_rand(cIndx,adj_new,firstOut_new,&cuts_new);

   int* stackC_ = (int*)malloc(sizeof(int)*cIndx*12);
   int* firstC_ = (int*)malloc(sizeof(int)*cIndx);
   int* nextC_ = (int*)malloc(sizeof(int)*cIndx*12);
   for(int i=0;i<cIndx;i++){firstC_[i]=-1;}
   int* cutIndx_ = (int*)malloc(sizeof(int)*cIndx*12);
   SP=0;
   for(int i=0;i<k2;i++)
   {
      for(int t=0;t<3;t++)
      {
         int x=cuts_new[6*i+2*t]; int y=cuts_new[6*i+2*t+1];
    
         nextC_[SP]=firstC_[x];
         firstC_[x]=SP;
         stackC_[SP]=y;
         cutIndx_[SP++]=i; 

         nextC_[SP]=firstC_[y];
         firstC_[y]=SP;
         stackC_[SP]=x;
         cutIndx_[SP++]=i; 
      }
   }

   int* stackC = (int*)malloc(sizeof(int)*cIndx*12);
   int* firstC = (int*)malloc(sizeof(int)*cIndx);
   int* nextC = (int*)malloc(sizeof(int)*cIndx*12);
   for(int i=0;i<cIndx;i++){firstC[i]=-1;}
   int* cutIndx = (int*)malloc(sizeof(int)*cIndx*12);
   SP=0;   
   for(int x=cIndx-1;x>=0;x--)
   {
      for(int indx=firstC_[x];indx!=-1;indx=nextC_[indx])
      {
         int y=stackC_[indx];
         nextC[SP]=firstC[y];
         firstC[y]=SP;
         stackC[SP]=x;
         cutIndx[SP++]=cutIndx_[indx];
      } 
   }

   int* currentEdge = (int*)malloc(sizeof(int)*k2);
   for(int i=0;i<k2;i++){currentEdge[i]=0;}
   for(int x=0;x<cIndx;x++)
   {
      for(int indx=firstC[x];indx!=-1;indx=nextC[indx])
      {
         int y=stackC[indx];
         if(y<x){continue;}
         int t=firstT[x];
         while(stackT[t]!=y){t=nextT[t];}
         firstT[x]=t;
         if(nextT[t]!=-1 && stackT[nextT[t]]==y){firstT[x]=nextT[t];}
         int cut_indx=cutIndx[indx];
         (*cuts)[6*k+6*cut_indx+2*currentEdge[cut_indx]]=tree_edge[2*t];
         (*cuts)[6*k+6*cut_indx+2*currentEdge[cut_indx]+1]=tree_edge[2*t+1];
         currentEdge[cut_indx]++;
      }
   }

   k+=k2;

   free(dfs); free(idfs); free(p); free(foundP);
   free(C); free(Q);
   free(stackT_); free(firstT_); free(nextT_); free(tree_edge_);
   free(stackT); free(firstT); free(nextT); free(tree_edge);
   free(adj_new); free(firstOut_new); free(currentOut);
   free(stackC_); free(firstC_); free(nextC_); free(cutIndx_);
   free(stackC); free(firstC); free(nextC); free(cutIndx);
   free(currentEdge);
   free(cuts_new);
   return k;
}

int get_3cuts_2tree_rand(int n, int* adj, int* firstOut, int** cuts)
{

   (*cuts) = (int*)malloc(sizeof(int)*12*n);

   int* dfs; int* idfs; int* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);

   int* high; int* highD; int* highIndx;

   get_high_indx(n,adj,firstOut,dfs,idfs,p,&high,&highD,&highIndx); 
   int* L; int* R; int* Lhead; int* Rhead; int* Lindx; int* Rindx;
   get_LandR(n,adj,firstOut,dfs,idfs,p,&L,&R,&Lhead,&Rhead,&Lindx,&Rindx); 

   int log2n=1; int tempN=n; while(tempN!=1){log2n++;tempN/=2;}
   int n_bits=3*log2n;
   int w=1;
   while(8*sizeof(int)*w<n_bits){w++;}

   int* hash = (int*)malloc(sizeof(int)*w*firstOut[n]);
   char* foundP = (char*)malloc(sizeof(char)*n); 
   for(int i=0;i<n;i++){foundP[i]=0;}
   for(int d=n-1;d>0;d--)
   {
      int x=idfs[d];
      for(int i=firstOut[x];i<firstOut[x+1];i++)
      {
         int y=adj[i];
         if(dfs[y]>dfs[x]){continue;}
         if(y==p[x]&&!foundP[x]){foundP[x]=1;continue;}
         for(int t=0;t<w;t++){hash[i*w+t]=rand();}
      }
   }
   int* H = (int*)malloc(sizeof(int)*w*n);
   for(int i=0;i<w*n;i++){H[i]=0;}
   char* foundC = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){foundP[i]=0; foundC[i]=0;}
   for(int d=n-1;d>0;d--)
   {
      int x=idfs[d];
      for(int i=firstOut[x];i<firstOut[x+1];i++)
      {
         int y=adj[i];
         if(y==p[x]&&!foundP[x]){foundP[x]=1;continue;}
         if(dfs[y]<dfs[x])
         {
            for(int t=0;t<w;t++){H[x*w+t]^=hash[i*w+t]; H[y*w+t]^=hash[i*w+t];}
         }
         else if(x==p[y]&&!foundC[y])
         {
            foundC[y]=1;
            for(int t=0;t<w;t++){H[x*w+t]^=H[y*w+t];}
         } 
      }
   }

   int N=1; for(int i=1;i<=log2n;i++){N*=2;}
   int** table = (int**)malloc(sizeof(int*)*N);
   int* capacity = (int*)malloc(sizeof(int)*N);
   int* taken = (int*)malloc(sizeof(int)*N);
   for(int i=0;i<N;i++){capacity[i]=1;}
   for(int i=0;i<N;i++){table[i]=(int*)malloc(sizeof(int)*w*capacity[i]);}
   int** corresponds = (int**)malloc(sizeof(int*)*N);
   for(int i=0;i<N;i++){corresponds[i]=(int*)malloc(sizeof(int)*capacity[i]);}
    
   for(int i=0;i<N;i++){taken[i]=0;}
   for(int u=1;u<n;u++)
   {
      int h=(uint)(H[u*w+0])%N;
      if(taken[h]==capacity[h])
      {
         capacity[h]*=2;
         table[h]=(int*)realloc(table[h],sizeof(int)*w*capacity[h]);
         corresponds[h]=(int*)realloc(corresponds[h],sizeof(int)*capacity[h]);
      }
      for(int t=0;t<w;t++){table[h][taken[h]*w+t]=H[u*w+t];}
      corresponds[h][taken[h]++]=u;
   }

   int cutIndx=0;
   int* low1; int* low1D; int* low2; int* low2D;
   get_2low(n,adj,firstOut,dfs,idfs,p,&low1,&low1D,&low2,&low2D);
   int* l; int* bcount;
   get_l_and_bcount(n,adj,firstOut,dfs,idfs,p,&l,&bcount);
   for(int v=1;v<n;v++)
   {
      if(bcount[v]==2)
      {
         (*cuts)[6*cutIndx+0]=v; (*cuts)[6*cutIndx+1]=p[v];
         (*cuts)[6*cutIndx+2]=low1D[v]; (*cuts)[6*cutIndx+3]=low1[v];
         (*cuts)[6*cutIndx+4]=low2D[v]; (*cuts)[6*cutIndx+5]=low2[v];
         cutIndx++;
      }
   }

   free(l); free(bcount); free(low1); free(low1D); free(low2); free(low2D);

   int* temp_hash = (int*)malloc(sizeof(int)*w);
   //L[v]
   for(int v=1;v<n;v++)
   {
      for(int t=0;t<w;t++){temp_hash[t]=H[v*w+t]^hash[Lindx[v]*w+t];}
      int h=(uint)(temp_hash[0])%N;
      int indx=-1;
      for(int i=0;i<taken[h];i++)
      {
         char found=1;
         for(int t=0;t<w;t++){if(table[h][i*w+t]!=temp_hash[t]){found=0;break;}}
         if(found){indx=i;break;}
      }
      if(indx==-1){continue;}
      int u=corresponds[h][indx];
      if(dfs[v]<dfs[u])
      {
         (*cuts)[6*cutIndx+0]=v; (*cuts)[6*cutIndx+1]=p[v];
         (*cuts)[6*cutIndx+2]=L[v]; (*cuts)[6*cutIndx+3]=Lhead[v];
         (*cuts)[6*cutIndx+4]=u; (*cuts)[6*cutIndx+5]=p[u]; cutIndx++;
      }      
   }
   //R[v]
   for(int v=1;v<n;v++)
   {
      for(int t=0;t<w;t++){temp_hash[t]=H[v*w+t]^hash[Rindx[v]*w+t];}
      int h=(uint)(temp_hash[0])%N;
      int indx=-1;
      for(int i=0;i<taken[h];i++)
      {
         char found=1;
         for(int t=0;t<w;t++){if(table[h][i*w+t]!=temp_hash[t]){found=0;break;}}
         if(found){indx=i;break;}
      }
      if(indx==-1){continue;}
      int u=corresponds[h][indx];
      if(dfs[v]<dfs[u])
      {
         (*cuts)[6*cutIndx+0]=v; (*cuts)[6*cutIndx+1]=p[v];
         (*cuts)[6*cutIndx+2]=R[v]; (*cuts)[6*cutIndx+3]=Rhead[v];
         (*cuts)[6*cutIndx+4]=u; (*cuts)[6*cutIndx+5]=p[u]; cutIndx++;
      }      
   }
   //highD[v]
   for(int v=1;v<n;v++)
   {
      for(int t=0;t<w;t++){temp_hash[t]=H[v*w+t]^hash[highIndx[v]*w+t];}
      int h=(uint)(temp_hash[0])%N;
      int indx=-1;
      for(int i=0;i<taken[h];i++)
      {
         char found=1;
         for(int t=0;t<w;t++){if(table[h][i*w+t]!=temp_hash[t]){found=0;break;}}
         if(found){indx=i;break;}
      }
      if(indx==-1){continue;}
      int u=corresponds[h][indx];
      if(dfs[v]>dfs[u])
      {
         (*cuts)[6*cutIndx+0]=v; (*cuts)[6*cutIndx+1]=p[v];
         (*cuts)[6*cutIndx+2]=highD[v]; (*cuts)[6*cutIndx+3]=high[v];
         (*cuts)[6*cutIndx+4]=u; (*cuts)[6*cutIndx+5]=p[u]; cutIndx++;
      }      
   }

   for(int i=0;i<N;i++){free(table[i]); free(corresponds[i]);}
   free(table); free(capacity); free(taken); free(corresponds);
   free(hash); free(H); free(temp_hash);

   free(dfs); free(idfs); free(p);
   free(high); free(highD); free(highIndx);
   free(L); free(R); free(Lhead); free(Rhead); free(Lindx); free(Rindx);
   free(foundP); free(foundC);
   return cutIndx;
}

void get_LandR(int n, int* adj, int* firstOut, int* dfs, int* idfs, int* p, int** L, int** R, int** Lhead, int** Rhead, int** Lindx, int** Rindx)
{
   (*L) = (int*)malloc(sizeof(int)*n);
   (*R) = (int*)malloc(sizeof(int)*n);
   (*Lhead) = (int*)malloc(sizeof(int)*n);
   (*Rhead) = (int*)malloc(sizeof(int)*n);
   (*Lindx) = (int*)malloc(sizeof(int)*firstOut[n]);
   (*Rindx) = (int*)malloc(sizeof(int)*firstOut[n]);

   int* next = (int*)malloc(sizeof(int)*n);
   int* prev = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){next[i]=i+1; prev[i]=i-1;}
   int* lastChild = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){lastChild[i]=-1;}
   for(int i=n-1;i>0;i--)
   {
      int c=idfs[i];
      if(lastChild[p[c]]==-1){lastChild[p[c]]=c;}
   }
   char* foundP = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){foundP[i]=0;}
   int* currentOut = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){currentOut[i]=firstOut[i];}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      int x=v;
      char found=0;
      while(!found)
      {
         for(int t=currentOut[x];t<firstOut[x+1];t++)
         {
            int y=adj[t];
            if(dfs[y]>dfs[x]){currentOut[x]=t+1;continue;}
            if(y==p[x]&&!foundP[x]){foundP[x]=1;currentOut[x]=t+1;continue;}
            if(dfs[y]<dfs[v])
            {
               (*L)[v]=x; (*Lhead)[v]=y; (*Lindx)[v]=t;
               found=1; break;
            }
            currentOut[x]=t+1;
         }
         if(!found)
         {
            int d=dfs[x];
            int next_vertex=idfs[next[d]];
            if(prev[d]!=-1){next[prev[d]]=next[d];}
            if(next[d]!=n){prev[next[d]]=prev[d];}
            x=next_vertex;
         }
      }
      if(lastChild[v]==-1){x=v;}
      else{x=(*R)[lastChild[v]];}
      found=0;
      while(!found)
      {
         for(int t=currentOut[x];t<firstOut[x+1];t++)
         {
            int y=adj[t];
            if(dfs[y]>dfs[x]){currentOut[x]=t+1;continue;}
            if(y==p[x]&&!foundP[x]){foundP[x]=1;currentOut[x]=t+1;continue;}
            if(dfs[y]<dfs[v])
            {
               (*R)[v]=x; (*Rhead)[v]=y; (*Rindx)[v]=t;
               found=1; break;
            }
            currentOut[x]=t+1;
         }
         if(!found)
         {
            int d=dfs[x];
            int next_vertex=idfs[prev[d]];
            if(prev[d]!=-1){next[prev[d]]=next[d];}
            if(next[d]!=n){prev[next[d]]=prev[d];}
            x=next_vertex;
         }
      }
   }
   free(next); free(prev); free(lastChild); free(foundP); free(currentOut);
}

void get_high_indx(int n, int* adj, int* firstOut, int* dfs, int* idfs, int* p, int** high, int** highD, int** highIndx)
{
   (*high) = (int*)malloc(sizeof(int)*n);
   (*highD) = (int*)malloc(sizeof(int)*n);
   (*highIndx) = (int*)malloc(sizeof(int)*firstOut[n]);
   for(int i=0;i<n;i++){(*high)[i]=-1;}
   int* ufparent = (int*)malloc(sizeof(int)*n);
   int* ufsize = (int*)malloc(sizeof(int)*n);
   int* representative = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){ufparent[i]=i; ufsize[i]=1; representative[i]=i;}
   char* foundP = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){foundP[i]=0;}
   for(int i=n-1;i>=0;i--)
   {
      int y=idfs[i];
      for(int j=firstOut[y];j<firstOut[y+1];j++)
      {
         int x=adj[j];
         if(dfs[x]<dfs[y]){continue;}
         if(y==p[x]&&!foundP[x]){foundP[x]=1;continue;}
         int v=representative[find(ufparent,x)];
         while(v!=y)
         {
            (*high)[v]=y; (*highD)[v]=x; //(*highIndx)[v]=j;
            int next=representative[find(ufparent,p[v])];
            unite(ufparent,ufsize,v,p[v]);
            representative[find(ufparent,v)]=y;
            v=next;
         }
      }
   }

   int* Hstack_ = (int*)malloc(sizeof(int)*2*n);
   int* Hfirst_ = (int*)malloc(sizeof(int)*n);
   int* Hnext_ = (int*)malloc(sizeof(int)*2*n);
   int* vertex_ = (int*)malloc(sizeof(int)*2*n);
   for(int i=0;i<n;i++){Hfirst_[i]=-1;}
   int SP=0;
   for(int v=1;v<n;v++)
   {
      int x=(*highD)[v]; int y=(*high)[v];
      Hnext_[SP]=Hfirst_[x]; Hfirst_[x]=SP; Hstack_[SP]=y; vertex_[SP++]=v;
      Hnext_[SP]=Hfirst_[y]; Hfirst_[y]=SP; Hstack_[SP]=x; vertex_[SP++]=v;
   }
   int* Hstack = (int*)malloc(sizeof(int)*n);
   int* Hfirst = (int*)malloc(sizeof(int)*n);
   int* Hnext = (int*)malloc(sizeof(int)*n);
   int* vertex = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){Hfirst[i]=-1;}
   SP=0;
   for(int y=n-1;y>=0;y--)
   {
      for(int i=Hfirst_[y];i!=-1;i=Hnext_[i])
      {
         int x=Hstack_[i];
         if(dfs[x]<dfs[y]){continue;}
         Hnext[SP]=Hfirst[x]; Hfirst[x]=SP; Hstack[SP]=y; vertex[SP++]=vertex_[i];
      }
   }

   int* adj2 = (int*)malloc(sizeof(int)*firstOut[n]);
   int* currentOut = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){currentOut[i]=firstOut[i];}
   int* index = (int*)malloc(sizeof(int)*firstOut[n]);
   for(int x=0;x<n;x++)
   {
      for(int i=firstOut[x];i<firstOut[x+1];i++)
      {
         int y=adj[i];
         adj2[currentOut[y]]=x; index[currentOut[y]++]=i; 
      }
   }

   int* adj3 = (int*)malloc(sizeof(int)*firstOut[n]);
   for(int i=0;i<n;i++){currentOut[i]=firstOut[i];}
   int* index_ = (int*)malloc(sizeof(int)*firstOut[n]);
   for(int y=0;y<n;y++)
   {
      for(int i=firstOut[y];i<firstOut[y+1];i++)
      {
         int x=adj2[i];
         adj3[currentOut[x]]=y; index_[currentOut[x]++]=index[i]; 
      }
   }   

   for(int i=0;i<n;i++){foundP[i]=0;}
   for(int x=0;x<n;x++)
   {
      for(int i=firstOut[x];i<firstOut[x+1];i++)
      {
         int y=adj3[i];
         if(dfs[x]<dfs[y]){continue;}
         if(y==p[x]&&!foundP[x]){foundP[x]=1;continue;}
         while(Hfirst[x]!=-1 && y==Hstack[Hfirst[x]])
         {
            (*highIndx)[vertex[Hfirst[x]]]=index_[i];
            Hfirst[x]=Hnext[Hfirst[x]];
         }
      }
   }

/*//flawless computation of highIndx
   for(int i=0;i<n;i++){foundP[i]=0;}
   for(int v=1;v<n;v++)
   {
      int x=(*highD)[v];
      for(int i=firstOut[x];i<firstOut[x+1];i++)
      {
         int y=adj[i];
         if(y==p[x]&&!foundP[x]){foundP[x]=1;continue;}
         if(y==(*high)[v]){(*highIndx)[v]=i;break;} 
      }
   }*/

   free(Hstack_); free(Hfirst_); free(Hnext_); free(vertex_);
   free(Hstack); free(Hfirst); free(Hnext); free(vertex);
   free(adj2); free(adj3); free(currentOut); free(index); free(index_);
   free(ufparent); free(ufsize); free(representative); free(foundP);
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
