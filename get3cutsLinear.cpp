#include <stdio.h>
#include <stdlib.h>
#include <chrono>

void get_adj(int,int,int*,int**,int**);
void read_graph(char*,int*,int**,int**);

void DFS(int,int*,int*,int**,int**,int**);
void get_2low(int,int*,int*,int*,int*,int*,int**,int**,int**,int**);
void get_l1l2_and_bcount(int,int*,int*,int*,int*,int*,int**,int**,int**);
void get_lowChildren(int,int*,int*,int*,int*,int**,int**);
void get_M(int,int*,int*,int*,int*,int*,int*,int**,int**);
void get_allM(int,int*,int*,int*,int*,int*,int*,int*,int**,int**,int**);
void get_lowM(int,int*,int*,int*,int*,int*,int*,int*,int*,int*,int**,int**);
void sortAdjInc(int,int*,int*,int*);

int get_3cuts_linear(int,int*,int*,int**);
int get_3cuts_2tree(int,int*,int*,int**);

using namespace std::chrono;

int main(int n_args, char** args)
{
   int n; int* adj; int* firstOut;
   read_graph(args[1],&n,&adj,&firstOut);
high_resolution_clock::time_point t1 = high_resolution_clock::now();
   int* cuts;
   int k=get_3cuts_linear(n,adj,firstOut,&cuts);
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

int get_3cuts_linear(int n, int* adj, int* firstOut, int** cuts)
{
   (*cuts) = (int*)malloc(sizeof(int)*n*12);
   int k=0; int k2;
   int* cuts_2tree;

   k2 = get_3cuts_2tree(n,adj,firstOut,&cuts_2tree);
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
   k2 = get_3cuts_linear(cIndx,adj_new,firstOut_new,&cuts_new);

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


int get_3cuts_2tree(int n, int* adj, int* firstOut, int** cuts)
{
   (*cuts) = (int*)malloc(sizeof(int)*n*12);
   int* dfs; int* idfs; int* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   int* low1; int* low1D; int* low2; int* low2D;
   get_2low(n,adj,firstOut,dfs,idfs,p,&low1,&low1D,&low2,&low2D);
   int* l1; int* l2; int* bcount;
   get_l1l2_and_bcount(n,adj,firstOut,dfs,idfs,p,&l1,&l2,&bcount);
   int* low1C; int* low2C;
   get_lowChildren(n,dfs,idfs,p,low1,&low1C,&low2C);
   int* M; int* nextM; 
   get_M(n,dfs,idfs,l1,low1,low1C,low2C,&M,&nextM);
   int* Ml; int* Mlow1; int* Mlow2;
   get_allM(n,dfs,idfs,l1,low1,M,low1C,low2C,&Ml,&Mlow1,&Mlow2);
  
   int* ND = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){ND[i]=1;}
   for(int i=n-1;i>0;i--){int v=idfs[i];ND[p[v]]+=ND[v];}
   int* prevM = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){prevM[i]=-1;}
   for(int i=0;i<n;i++){if(nextM[i]!=-1){prevM[nextM[i]]=i;}}

   int* lowM; int* lowMD;
   get_lowM(n,adj,firstOut,dfs,idfs,p,ND,low1C,M,prevM,&lowM,&lowMD);

   int* low3C = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){low3C[i]=-1;}
   for(int c=1;c<n;c++)
   {
      int v=p[c];
      if(low1C[v]==c||low2C[v]==c){continue;}
      if(low3C[v]==-1||dfs[low1[c]]<dfs[low1[low3C[v]]]){low3C[v]=c;}
   }

   int k=0;
   int* currentVertex = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){currentVertex[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      int m=Ml[v];
      if(m==-1){continue;}
      int u=currentVertex[m];
      while(nextM[u]!=-1 && dfs[nextM[u]]>dfs[v]){u=nextM[u];}
      currentVertex[m]=u;
      if(bcount[v]==bcount[u]+1 && dfs[l2[M[v]]]>=dfs[v] && (low2C[M[v]]==-1 || dfs[low1[low2C[M[v]]]]>=dfs[v]))
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=M[v]; (*cuts)[6*k+5]=l1[M[v]];
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
      if(bcount[v]==bcount[u]+1 && dfs[low2[Mlow2[v]]]>=dfs[v] && (low3C[M[v]]==-1 || dfs[low1[low3C[M[v]]]]>=dfs[v]))
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=Mlow2[v]; (*cuts)[6*k+5]=l1[Mlow2[v]];
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
      if(bcount[v]==bcount[u]+1 && dfs[low2[Mlow1[v]]]>=dfs[v] && (low3C[M[v]]==-1 || dfs[low1[low3C[M[v]]]]>=dfs[v]))
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=Mlow1[v]; (*cuts)[6*k+5]=l1[Mlow1[v]];
         k++;
      }
   }
  
   for(int u=1;u<n;u++)
   {
      if(nextM[u]!=-1 && bcount[u]==bcount[nextM[u]]+1)
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=nextM[u]; (*cuts)[6*k+3]=p[nextM[u]];
         (*cuts)[6*k+4]=lowMD[u]; (*cuts)[6*k+5]=lowM[u];
         k++;
      }
   }
   for(int i=0;i<n;i++){currentVertex[i]=i;}
   for(int i=n-1;i>0;i--)
   {
      int u=idfs[i];
      int m=Ml[u];
      if(m==-1 || Ml[u]==M[u]){continue;}
      int v=currentVertex[m];
      while(v!=-1 && dfs[v]>=dfs[u]){v=nextM[v];}
      currentVertex[m]=v;
      if(v!=-1 && bcount[u]==bcount[v]+1)
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=M[u]; (*cuts)[6*k+5]=l1[M[u]];
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
      if(v!=-1 && bcount[u]==bcount[v]+1)
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=Mlow2[u]; (*cuts)[6*k+5]=l1[Mlow2[u]];
         k++;
      }
   }
  
   for(int v=1;v<n;v++)
   {
      if(bcount[v]==2)
      {
         (*cuts)[6*k+0]=v; (*cuts)[6*k+1]=p[v];
         (*cuts)[6*k+2]=low1D[v]; (*cuts)[6*k+3]=low1[v];
         (*cuts)[6*k+4]=low2D[v]; (*cuts)[6*k+5]=low2[v];
         k++;   
      }
   }
 
   free(dfs); free(idfs); free(p); free(ND);
   free(low1); free(low1D); free(low2); free(low2D);  free(lowM); free(lowMD);
   free(l1); free(l2); free(bcount);
   free(low1C); free(low2C); free(low3C);
   free(M); free(nextM); free(prevM);
   free(Ml); free(Mlow1); free(Mlow2); free(currentVertex);
   return k;
}

void get_lowM(int n, int* adj, int* firstOut, int* dfs, int* idfs, int* p, int* ND, int* low1C, int* M, int* prevM, int** lowM, int** lowMD)
{
   sortAdjInc(n,adj,firstOut,idfs);
   (*lowMD) = (int*)malloc(sizeof(int)*n);
   (*lowM) = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){(*lowMD)[i]=-1; (*lowM)[i]=-1;}
   int* currentOut = (int*)malloc(sizeof(int)*n);
   char* foundC = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){foundC[i]=0;}
   for(int i=0;i<n;i++){currentOut[i]=firstOut[i];}
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
      if(prevM[v]==-1){continue;}
      int u=prevM[v];
      int y=v;
      while((*lowM)[u]==-1)
      {
         while(currentOut[y]!=firstOut[y+1])
         {
            int x=adj[currentOut[y]];
            if(dfs[x]<dfs[y]){currentOut[y]++;continue;}
            if(y==p[x]&&foundC[x]==0){foundC[x]=1;currentOut[y]++;continue;}
            if(dfs[x]<dfs[M[u]]){currentOut[y]++;continue;}
            if(dfs[x]<dfs[M[u]]+ND[M[u]])
            {
               (*lowMD)[u]=x; (*lowM)[u]=y;
            }
            break;
         } 
         if((*lowM)[u]==-1)
         {
            if(prevM[low1C[y]]==-1){y=low1C[y];}
            else{y=(*lowM)[prevM[low1C[y]]];}
         }
      }
   }
   free(currentOut); free(foundC);
}

void get_l1l2_and_bcount(int n, int* adj, int* firstOut, int* dfs, int* idfs, int* p, int** l1, int** l2, int** bcount)
{
   (*l1) = (int*)malloc(sizeof(int)*n);
   (*l2) = (int*)malloc(sizeof(int)*n);
   (*bcount) = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){(*l1)[i]=i; (*l2)[i]=i; (*bcount)[i]=0;}
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
               if(dfs[y]<=dfs[(*l1)[x]]){(*l2)[x]=(*l1)[x]; (*l1)[x]=y;}
               else if(dfs[y]<dfs[(*l2)[x]]){(*l2)[x]=y;}
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

void sortAdjInc(int n, int* adj, int* firstOut, int* idfs)
{
   int* adj_copy = (int*)malloc(sizeof(int)*firstOut[n]);
   int* currentOut = (int*)malloc(sizeof(int)*(n+1));
   for(int i=0;i<firstOut[n];i++){adj_copy[i]=adj[i];}
   for(int i=0;i<=n;i++){currentOut[i]=firstOut[i];}
   for(int i=0;i<n;i++)
   {
      int x=idfs[i];
      for(int t=firstOut[x];t<firstOut[x+1];t++)
      {
         int y=adj_copy[t];
         adj[currentOut[y]++]=x;
      }
   }
   free(adj_copy); free(currentOut);
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
