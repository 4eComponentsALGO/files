#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void get_adj(int,int,int*,int**,int**);

void make_2e_connected(int,int*,int*,int**,int**);
int make_3e_connected(int,int*,int*,int**);

void DFS(int,int*,int*,int**,int**,int**);
void get_low(int,int*,int*,int*,int*,int*,int**);
void get_l_and_bcount(int,int*,int*,int*,int*,int*,int**,int**);
void get_lowChildren(int,int*,int*,int*,int*,int**,int**);
void get_M(int,int*,int*,int*,int*,int*,int*,int**,int**);

void main(int n_args, char** args)
{
   srand(time(NULL));
   FILE* fp = fopen(args[1],"r");
   int n; int m;
   fscanf(fp,"%d %d",&n,&m);
   int* initial_edges = (int*)malloc(sizeof(int)*m*2);
   for(int i=0;i<m;i++){fscanf(fp,"%d %d",initial_edges+2*i,initial_edges+2*i+1);}
   fclose(fp);

   int* adj; int* firstOut;
   get_adj(n,m,initial_edges,&adj,&firstOut);
   free(initial_edges);

   int* adj2; int* firstOut2;
   make_2e_connected(n,adj,firstOut,&adj2,&firstOut2);
   free(adj); free(firstOut);

   int* edges;
   m = make_3e_connected(n,adj2,firstOut2,&edges);

   fp = fopen(args[2],"w");
   fprintf(fp,"%d %d\n",n,m);
   for(int i=0;i<m;i++){fprintf(fp,"%d %d\n",edges[2*i],edges[2*i+1]);}
   fclose(fp);

   free(adj2); free(firstOut2);
   free(edges);
}

int make_3e_connected(int n, int* adj, int* firstOut, int** edges)
{
   int* dfs; int* idfs; int* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   int* low; int* l; int* bcount;
   get_low(n,adj,firstOut,dfs,idfs,p,&low);
   get_l_and_bcount(n,adj,firstOut,dfs,idfs,p,&l,&bcount);
   int* low1C; int* low2C;
   get_lowChildren(n,dfs,idfs,p,low,&low1C,&low2C);
   int* M; int* nextM;
   get_M(n,dfs,idfs,l,low,low1C,low2C,&M,&nextM);

   int* vEdgeStack = (int*)malloc(sizeof(int)*2*n);
   int* vEdgeFirst = (int*)malloc(sizeof(int)*n);
   int* vEdgeNext = (int*)malloc(sizeof(int)*2*n);
   for(int i=0;i<n;i++){vEdgeFirst[i]=-1;}
   char* isCutEdge = (char*)malloc(sizeof(char)*n);
   char* isCutEdge2 = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){isCutEdge[i]=0; isCutEdge2[i]=0;}
   int SP=0;
   int* corresponding = (int*)malloc(sizeof(int)*n);
   int* corresponding2 = (int*)malloc(sizeof(int)*n);

   for(int m=1;m<n;m++)
   {
      if(M[m]!=m){continue;}
      int v=m;
      while(v!=-1)
      {
         int z=v;
         if(bcount[v]==1){isCutEdge2[m]=1; isCutEdge[v]=1; corresponding2[m]=v; corresponding[v]=v;}
         while(nextM[z]!=-1 && bcount[nextM[z]]==bcount[v])
         {
            z=nextM[z];  
         }
         if(z==v){v=nextM[v];continue;}
         for(int u=v;u!=nextM[z];u=nextM[u])
         {
            corresponding[u]=v;
            isCutEdge[u]=1;
         }
         if(bcount[v]!=1)
         {
            vEdgeNext[SP]=vEdgeFirst[v]; vEdgeFirst[v]=SP; vEdgeStack[SP++]=p[z];
            vEdgeNext[SP]=vEdgeFirst[p[z]]; vEdgeFirst[p[z]]=SP; vEdgeStack[SP++]=v;
         }
         v=nextM[z];
      } 
   }

   int* C = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){C[i]=-1;}
   int k=0;
   int* Q = (int*)malloc(sizeof(int)*n);
   for(int r=0;r<n;r++)
   {
      if(C[r]!=-1){continue;}
      int first=0; int last=0;
      Q[last++]=r; C[r]=k;
      while(first!=last)
      {
         int x=Q[first++];
         for(int i=firstOut[x];i<firstOut[x+1];i++)
         {
            int y=adj[i];
            if(C[y]!=-1){continue;}
            if((x==p[y]&&isCutEdge[y])||(y==p[x]&&isCutEdge[x])){continue;}
            if((M[x]==x&&y==l[x]&&isCutEdge2[x])||(M[y]==y&&x==l[y]&&isCutEdge2[y])){continue;}
            Q[last++]=y; C[y]=k;
         }
         for(int i=vEdgeFirst[x];i!=-1;i=vEdgeNext[i])
         {
            int y=vEdgeStack[i];
            if(C[y]!=-1){continue;}
            Q[last++]=y; C[y]=k;
         }
      }
      k++;
   }

   int* compStack = (int*)malloc(sizeof(int)*n);
   int* compFirst = (int*)malloc(sizeof(int)*k);
   int* compNext = (int*)malloc(sizeof(int)*n);
   int* compSize = (int*)malloc(sizeof(int)*k);
   for(int i=0;i<k;i++){compFirst[i]=-1; compSize[i]=0;}
   SP=0;
   for(int v=0;v<n;v++)
   {
      int c=C[v];
      compNext[SP]=compFirst[c]; compFirst[c]=SP; compStack[SP++]=v; compSize[c]++;
   }

   int** adjC = (int**)malloc(sizeof(int*)*k);
   int** cname = (int**)malloc(sizeof(int*)*k);
   int* capacity = (int*)malloc(sizeof(int)*k);
   int* size = (int*)malloc(sizeof(int)*k);
   for(int i=0;i<k;i++){capacity[i]=1; size[i]=0;}
   for(int i=0;i<k;i++){adjC[i]=(int*)malloc(sizeof(int)*capacity[i]);}
   for(int i=0;i<k;i++){cname[i]=(int*)malloc(sizeof(int)*capacity[i]);}
   for(int x=0;x<n;x++)
   {
      for(int i=firstOut[x];i<firstOut[x+1];i++)
      {
         int y=adj[i];
         if(C[x]==C[y]){continue;}
         int c=C[x];
         if(size[c]==capacity[c])
         {
            capacity[c]*=2;
            adjC[c]=(int*)realloc(adjC[c],sizeof(int)*capacity[c]);
            cname[c]=(int*)realloc(cname[c],sizeof(int)*capacity[c]);
         }          
         int corr;
         if(dfs[x]>dfs[y]&&y==p[x]){corr=corresponding[x];}
         if(dfs[x]>dfs[y]&&y!=p[x]){corr=corresponding2[x];}
         if(dfs[y]>dfs[x]&&x==p[y]){corr=corresponding[y];}
         if(dfs[y]>dfs[x]&&x!=p[y]){corr=corresponding2[y];}
         adjC[c][size[c]]=C[y]; cname[c][size[c]]=corr; size[c]++;
      }
   }

   int* dfsC = (int*)malloc(sizeof(int)*k);
   int* idfsC = (int*)malloc(sizeof(int)*k);
   for(int i=0;i<k;i++){dfsC[i]=-1;}
   int* temp_node = (int*)malloc(sizeof(int)*k);
   int* temp_indx = (int*)malloc(sizeof(int)*k);
   int Nr=0;
   dfsC[0]=Nr; idfsC[Nr++]=0;
   temp_node[0]=0; temp_indx[0]=0;
   SP=0;
   while(SP!=-1)
   {
      int v=temp_node[SP];
      char descend=0;
      for(int i=temp_indx[SP];i<size[v];i++)
      {
         int u=adjC[v][i];
         if(dfsC[u]==-1)
         {
            dfsC[u]=Nr; idfsC[Nr++]=u;
            temp_node[SP+1]=u; temp_indx[SP+1]=0; temp_indx[SP]=i+1;
            for(int t=0;t<size[u];t++)
            {
               int z=adjC[u][t];
               if(z!=v && cname[u][t]==cname[v][i])
               {
                  int temp;
                  temp=adjC[u][t]; adjC[u][t]=adjC[u][size[u]-1]; adjC[u][size[u]-1]=temp;
                  temp=cname[u][t]; cname[u][t]=cname[u][size[u]-1]; cname[u][size[u]-1]=temp;
                  break; 
               }
            }
            descend=1; break;
         }
      }
      if(descend){SP++;continue;}
      SP--;
   }

   int* leaf_stack = (int*)malloc(sizeof(int)*k);
   char* is_available = (char*)malloc(sizeof(char)*k);
   for(int i=0;i<k;i++){is_available[i]=1;}
   SP=0;
   for(int i=0;i<k;i++)
   {
      int c=idfsC[i];
      if(size[c]==2 && is_available[c])
      {
         leaf_stack[SP++]=c; is_available[c]=0;
      }
   }

   (*edges) = (int*)malloc(sizeof(int)*(firstOut[n]+SP+2*(SP%2)));
   int edgeIndx=0;

   for(int i=0;i<SP/2;i++)
   {
      int c=leaf_stack[i];
      int d=leaf_stack[i+SP/2];
      int i1=rand()%compSize[c];
      int j1=rand()%compSize[d];
      int indx1=compFirst[c];
      while(i1-->0){indx1=compNext[indx1];}
      int x=compStack[indx1];
      int indx2=compFirst[d];
      while(j1-->0){indx2=compNext[indx2];}
      int y=compStack[indx2];
      (*edges)[2*edgeIndx]=x; (*edges)[2*edgeIndx+1]=y; edgeIndx++; //printf("added (%d,%d)\n",x,y);
   }
   if(SP%2==1)
   {
      int c=leaf_stack[SP-1];
      int d=leaf_stack[rand()%(SP-1)];
      int i1=rand()%compSize[c];
      int j1=rand()%compSize[d];
      int indx1=compFirst[c];
      while(i1-->0){indx1=compNext[indx1];}
      int x=compStack[indx1];
      int indx2=compFirst[d];
      while(j1-->0){indx2=compNext[indx2];}
      int y=compStack[indx2];
      (*edges)[2*edgeIndx]=x; (*edges)[2*edgeIndx+1]=y; edgeIndx++; //printf("added (%d,%d)\n",x,y);
   }

   for(int x=0;x<n;x++)
   {
      for(int i=firstOut[x];i<firstOut[x+1];i++)
      {
         int y=adj[i];
         if(x>y){continue;}
         (*edges)[2*edgeIndx]=x; (*edges)[2*edgeIndx+1]=y; edgeIndx++;
      }
   }

   free(dfs); free(idfs); free(p); free(low); free(l); free(bcount);
   free(low1C); free(low2C); free(M); free(nextM);
   free(vEdgeStack); free(vEdgeFirst); free(vEdgeNext);
   free(isCutEdge); free(isCutEdge2); free(corresponding); free(corresponding2);
   free(C); free(Q);
   free(compStack); free(compFirst); free(compNext); free(compSize);
   for(int i=0;i<k;i++){free(adjC[i]); free(cname[i]);}
   free(adjC); free(cname); free(capacity); free(size);
   free(dfsC); free(idfsC); free(temp_node); free(temp_indx);;
   free(leaf_stack); free(is_available);

   return edgeIndx;
}

void make_2e_connected(int n, int* adj, int* firstOut, int** adj_new, int** firstOut_new)
{
   int* edgeStack = (int*)malloc(sizeof(int)*4*n);
   int* edgeFirst = (int*)malloc(sizeof(int)*n);
   int* edgeNext = (int*)malloc(sizeof(int)*4*n);
   for(int i=0;i<n;i++){edgeFirst[i]=-1;}
   int edgeIndx=0;
   
   char* found = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){found[i]=0;}
   int* Q = (int*)malloc(sizeof(int)*n);
   int prevRoot=-1;
   for(int r=0;r<n;r++)
   {
      if(found[r]){continue;}
      int first=0; int last=0;
      Q[last++]=r; found[r]=1;
      while(first!=last)
      {
         int x=Q[first++];
         for(int i=firstOut[x];i<firstOut[x+1];i++)
         {
            int y=adj[i];
            if(!found[y]){Q[last++]=y; found[y]=1;}
         }
      }
      if(prevRoot!=-1)
      {
         edgeNext[edgeIndx]=edgeFirst[prevRoot]; edgeFirst[prevRoot]=edgeIndx; edgeStack[edgeIndx++]=r;
         edgeNext[edgeIndx]=edgeFirst[r]; edgeFirst[r]=edgeIndx; edgeStack[edgeIndx++]=prevRoot;
      }
      prevRoot=r;
   }

   int* dfs = (int*)malloc(sizeof(int)*n);
   int* idfs = (int*)malloc(sizeof(int)*n);
   int* p = (int*)malloc(sizeof(int)*n);
   int* low = (int*)malloc(sizeof(int)*n); 
   for(int i=0;i<n;i++){dfs[i]=-1;}
   int* temp_vertex = (int*)malloc(sizeof(int)*n);
   int* temp_out = (int*)malloc(sizeof(int)*n);
   int* temp_next = (int*)malloc(sizeof(int)*n);
   char* foundP = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){foundP[i]=0;low[i]=i;}
   int Nr=0;
   dfs[0]=Nr; idfs[Nr++]=0; p[0]=-1;
   temp_vertex[0]=0; temp_out[0]=firstOut[0]; temp_next[0]=edgeFirst[0];
   int SP=0;
   while(SP!=-1)
   {
      int v=temp_vertex[SP];
      char descend=0; 
      for(int i=temp_out[SP];i<firstOut[v+1];i++)
      {
         int u=adj[i];
         if(dfs[u]==-1)
         {
            dfs[u]=Nr; idfs[Nr++]=u; p[u]=v;
            temp_vertex[SP+1]=u; temp_out[SP+1]=firstOut[u]; temp_out[SP]=i; temp_next[SP+1]=edgeFirst[u];
            descend=1; break; 
         }
         if(u==p[v]&&!foundP[v]){foundP[v]=1;continue;}
         if(dfs[u]<dfs[v])
         {
            if(dfs[u]<dfs[low[v]]){low[v]=u;}
         }
         else if(v==p[u])
         {
            if(dfs[low[u]]<dfs[low[v]]){low[v]=low[u];}
         } 
      }
      if(descend){SP++;continue;}
      for(int i=temp_next[SP];i!=-1;i=edgeNext[i])
      {
         int u=edgeStack[i]; 
         if(dfs[u]==-1)
         {
            dfs[u]=Nr; idfs[Nr++]=u; p[u]=v; low[u]=u;
            temp_vertex[SP+1]=u; temp_out[SP+1]=firstOut[u]; temp_next[SP+1]=edgeFirst[u]; temp_next[SP]=edgeNext[i];
            descend=1; break;
         }
      }
      if(descend){SP++;continue;}
      SP--;
   }

   int* nBridgesToVertex = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){nBridgesToVertex[i]=0;}
   for(int v=1;v<n;v++)
   {
      if(low[v]==v)
      {
         nBridgesToVertex[v]++; nBridgesToVertex[p[v]]++;;
      }
   }
   int* nBridgesToComp = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){nBridgesToComp[i]=0;}
   int* C = (int*)malloc(sizeof(int)*n);
   C[0]=0; nBridgesToComp[0]=nBridgesToVertex[0]; int k=1;
   for(int i=1;i<n;i++)
   {
      int v=idfs[i];
      if(low[v]!=v){C[v]=C[p[v]];}
      else{C[v]=k++;}
      nBridgesToComp[C[v]]+=nBridgesToVertex[v];
   }

   int* compStack = (int*)malloc(sizeof(int)*n);
   int* compFirst = (int*)malloc(sizeof(int)*k);
   int* compNext = (int*)malloc(sizeof(int)*n);
   int* compSize = (int*)malloc(sizeof(int)*k);
   for(int i=0;i<k;i++){compFirst[i]=-1; compSize[i]=0;}
   SP=0;
   for(int v=0;v<n;v++)
   {
      int c=C[v];
      compNext[SP]=compFirst[c]; compFirst[c]=SP; compStack[SP++]=v; compSize[c]++;
   }

   char* available = (char*)malloc(sizeof(char)*k);
   for(int i=0;i<k;i++){available[i]=1;}
   int* comp_stack = (int*)malloc(sizeof(int)*k);
   SP=0;
   for(int i=0;i<n;i++)
   {
      int c=C[idfs[i]];
      if(nBridgesToComp[c]==1&&available[c])
      {
         available[c]=0; comp_stack[SP++]=c;
      }
   }

   for(int t=0;t<SP/2;t++)
   {
      int i=comp_stack[t];
      int j=comp_stack[t+SP/2];
      int i1 = rand()%compSize[i];
      int j1 = rand()%compSize[j];
      int indx1=compFirst[i];
      for(int t=0;t<i1;t++){indx1=compNext[indx1];}
      int indx2=compFirst[j];
      for(int t=0;t<j1;t++){indx2=compNext[indx2];}
      int x=compStack[indx1]; int y=compStack[indx2];
      edgeNext[edgeIndx]=edgeFirst[x]; edgeFirst[x]=edgeIndx; edgeStack[edgeIndx++]=y;
      edgeNext[edgeIndx]=edgeFirst[y]; edgeFirst[y]=edgeIndx; edgeStack[edgeIndx++]=x;
   }
   if(SP%2==1)
   {
      int i=comp_stack[SP-1];
      int j=comp_stack[rand()%(SP-1)];
      int i1 = rand()%compSize[i];
      int j1 = rand()%compSize[j];
      int indx1=compFirst[i];
      for(int t=0;t<i1;t++){indx1=compNext[indx1];}
      int indx2=compFirst[j];
      for(int t=0;t<j1;t++){indx2=compNext[indx2];}
      int x=compStack[indx1]; int y=compStack[indx2];
      edgeNext[edgeIndx]=edgeFirst[x]; edgeFirst[x]=edgeIndx; edgeStack[edgeIndx++]=y;
      edgeNext[edgeIndx]=edgeFirst[y]; edgeFirst[y]=edgeIndx; edgeStack[edgeIndx++]=x;
   }

   int* edges = (int*)malloc(sizeof(int)*(firstOut[n]+edgeIndx));
   edgeIndx=0;
   for(int x=0;x<n;x++)
   {
      for(int i=firstOut[x];i<firstOut[x+1];i++)
      {
         int y=adj[i];
         if(x>y){continue;}
         edges[2*edgeIndx]=x; edges[2*edgeIndx+1]=y; edgeIndx++;
      }
   }
   for(int x=0;x<n;x++)
   {
      for(int i=edgeFirst[x];i!=-1;i=edgeNext[i])
      {
         int y=edgeStack[i];
         if(x>y){continue;}
         edges[2*edgeIndx]=x; edges[2*edgeIndx+1]=y; edgeIndx++;
      }
   }

   get_adj(n,edgeIndx,edges,adj_new,firstOut_new);

   free(edgeStack); free(edgeFirst); free(edgeNext);
   free(found); free(Q); free(edges);
   free(dfs); free(idfs); free(p); free(low);
   free(temp_vertex); free(temp_out); free(temp_next); free(foundP);
   free(compStack); free(compFirst); free(compNext); free(compSize);
   free(nBridgesToComp); free(nBridgesToVertex); free(C);
   free(available); free(comp_stack);
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
