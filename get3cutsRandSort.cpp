#include <stdio.h>
#include <stdlib.h>
#include <chrono>

void get_adj(int,int,int*,int**,int**);
void read_graph(char*,int*,int**,int**);

void DFS(int,int*,int*,int**,int**,int**);
void get_2low(int,int*,int*,int*,int*,int*,int**,int**,int**,int**);
void get_bcount(int,int*,int*,int*,int*,int*,int**);
void get_LandR(int,int*,int*,int*,int*,int*,int**,int**,int**,int**,int**,int**);
void get_high_sort(int,int*,int*,int*,int*,int*,int**,int**,int**);
int find(int*,int);
void unite(int*,int*,int,int);

int get_3cuts_rand(int,int*,int*,int**);
int get_3cuts_2tree_rand(int,int*,int*,int**);

using namespace std::chrono;

int main(int n_args, char** args)
{
   int n; int* adj; int* firstOut;
   read_graph(args[1],&n,&adj,&firstOut);
high_resolution_clock::time_point t1 = high_resolution_clock::now();
   int* cuts;
   int k=get_3cuts_rand(n,adj,firstOut,&cuts);
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

   get_high_sort(n,adj,firstOut,dfs,idfs,p,&high,&highD,&highIndx); 
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
   int* bcount;
   get_bcount(n,adj,firstOut,dfs,idfs,p,&bcount);
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

   free(bcount); free(low1); free(low1D); free(low2); free(low2D);

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

void get_high_sort(int n, int* adj, int* firstOut, int* dfs, int* idfs, int* p, int** high, int** highD, int** highIndx)
{
   int* backEdges = (int*)malloc(sizeof(int)*(firstOut[n]-2*n+2));
   int* index = (int*)malloc(sizeof(int)*(firstOut[n]/2-n+1));
   int bp=0;  
   char* foundP = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){foundP[i]=0;}
   for(int x=0;x<n;x++)
   {
      for(int i=firstOut[x];i<firstOut[x+1];i++)
      {
         int y=adj[i];
         if(y==p[x]&&!foundP[x]){foundP[x]=1;continue;}
         if(dfs[y]>dfs[x]){continue;}
         backEdges[2*bp]=x; backEdges[2*bp+1]=y; index[bp++]=i;
      }
   }
   
   //sort the back-edges decreasingly w.r.t. their lower end
   int* backEdgesStack = (int*)malloc(sizeof(int)*2*bp);
   int* backEdgesFirst = (int*)malloc(sizeof(int)*n);
   int* backEdgesNext = (int*)malloc(sizeof(int)*2*bp);
   for(int i=0;i<n;i++){backEdgesFirst[i]=-1;}
   int SP=0;
   for(int i=0;i<bp;i++)
   {
      int d=dfs[backEdges[2*i+1]];
      backEdgesNext[SP]=backEdgesFirst[d]; backEdgesFirst[d]=SP; backEdgesStack[SP++]=i;
   }

   int* ufparent = (int*)malloc(sizeof(int)*n);
   int* ufsize = (int*)malloc(sizeof(int)*n);
   int* lowest = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){ufparent[i]=i; ufsize[i]=1; lowest[i]=i;}
   (*highIndx) = (int*)malloc(sizeof(int)*n);
   (*high) = (int*)malloc(sizeof(int)*n);
   (*highD) = (int*)malloc(sizeof(int)*n);
   for(int d=n-1;d>=0;d--)
   {
      for(int i=backEdgesFirst[d];i!=-1;i=backEdgesNext[i])
      {
         int t=backEdgesStack[i];
         int x=backEdges[2*t]; int y=backEdges[2*t+1];
         int u=lowest[find(ufparent,x)];
         while(u!=y)
         {
            (*high)[u]=y; (*highD)[u]=x; (*highIndx)[u]=index[t];
            int next=lowest[find(ufparent,p[u])];
            unite(ufparent,ufsize,u,p[u]);
            lowest[find(ufparent,p[u])]=y;
            u=next;
         }
      }
   }

   free(backEdges); free(index); free(foundP);
   free(backEdgesStack); free(backEdgesFirst); free(backEdgesNext);
   free(ufparent); free(ufsize); free(lowest);
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

void get_bcount(int n, int* adj, int* firstOut, int* dfs, int* idfs, int* p, int** bcount)
{
   (*bcount) = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){(*bcount)[i]=0;}
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
