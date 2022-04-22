#include <stdio.h>
#include <stdlib.h>
#include <chrono>

void get_adj(int,int,int*,int**,int**);
void read_graph(char*,int*,int**,int**);

void DFS(int,int*,int*,int**,int**,int**);
void get_2low(int,int*,int*,int*,int*,int*,int**,int**,int**,int**);
void get_l_and_bcount(int,int*,int*,int*,int*,int*,int**,int**);
void get_lowChildren(int,int*,int*,int*,int*,int**,int**);
void get_M(int,int*,int*,int*,int*,int*,int*,int**,int**);
void get_allM(int,int*,int*,int*,int*,int*,int*,int*,int**,int**,int**);
void get_high(int,int*,int*,int*,int*,int*,int**,int**);
int find(int*,int);
void unite(int*,int*,int,int);

int get_3cuts_high(int,int*,int*,int**);

using namespace std::chrono;

int main(int n_args, char** args)
{
   int n; int* adj; int* firstOut;
   read_graph(args[1],&n,&adj,&firstOut);
high_resolution_clock::time_point t1 = high_resolution_clock::now();
   int* cuts;
   int k=get_3cuts_high(n,adj,firstOut,&cuts);
high_resolution_clock::time_point t2 = high_resolution_clock::now();
   FILE* fp = fopen(args[2],"w");
   fprintf(fp,"%d\n",k);
   for(int i=0;i<k;i++){fprintf(fp,"%d %d %d %d %d %d\n",cuts[6*i+0],cuts[6*i+1],cuts[6*i+2],cuts[6*i+3],cuts[6*i+4],cuts[6*i+5]);}
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
   fprintf(fp,"%f\n",time_span.count());
   fclose(fp);
   printf("%f\n",time_span.count());;
   free(adj); free(firstOut); free(cuts);
   return 0;
}

int get_3cuts_high(int n, int* adj, int* firstOut, int** cuts)
{
   (*cuts) = (int*)malloc(sizeof(int)*n*12);
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

   int k=0;

   //one tree-edge
   for(int v=1;v<n;v++)
   {
      if(bcount[v]==2)
      {
         (*cuts)[6*k+0]=v; (*cuts)[6*k+1]=p[v];
         (*cuts)[6*k+2]=low[v]; (*cuts)[6*k+3]=lowD[v];
         (*cuts)[6*k+4]=low2[v]; (*cuts)[6*k+5]=low2D[v];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=M[v]; (*cuts)[6*k+5]=l[M[v]];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=Mlow2[v]; (*cuts)[6*k+5]=l[Mlow2[v]];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=Mlow1[v]; (*cuts)[6*k+5]=l[Mlow1[v]];
         k++;
      }
   }

   //lower case
   for(int u=1;u<n;u++)
   {
      if(nextM[u]!=-1 && bcount[nextM[u]]==bcount[u]-1)
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=nextM[u]; (*cuts)[6*k+3]=p[nextM[u]];
         (*cuts)[6*k+4]=highD[u]; (*cuts)[6*k+5]=high[u];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=highD[u]; (*cuts)[6*k+5]=high[u];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=highD[u]; (*cuts)[6*k+5]=high[u];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=w; (*cuts)[6*k+5]=p[w];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=w; (*cuts)[6*k+5]=p[w];
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
                  (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
                  (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
                  (*cuts)[6*k+4]=w; (*cuts)[6*k+5]=p[w];
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
            (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
            (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
            (*cuts)[6*k+4]=w; (*cuts)[6*k+5]=p[w];
            k++;  
         }
         indx=nextU[indx];
      }
   }

   free(dfs); free(idfs); free(p);
   free(low); free(lowD); free(low2); free(low2D);
   free(low1C); free(low2C); free(l); free(bcount);
   free(M); free(nextM); free(high); free(highD);
   free(Ml); free(Mlow1); free(Mlow2); free(currentVertex); 
   free(lastM); free(A); free(ND); 
   free(Hstack); free(firstH); free(nextH); free(stack);
   free(Ustack); free(firstU); free(nextU);
   free(lowestW);
   return k;
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
