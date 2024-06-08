/**************************************************************************************************
 * If you use this code, please cite:
 
 G. Bianconi and S. N. Dorogovstev 
 "Nature of hypergraph k-core percolation problems" 
 Physical Review E, 109, 014307 (2024).
 
***************************************************************************************************
 * Code that  generates random hypergraph and performs (2ng neighbour) hyperedge (k,s)-core percolation 
 * This code uses:
 * N  Number of nodes
 * M number of hyperedges
 * Kcore  k of the core
 * Ncore  s of the core
 *
 * Nrunmax  Number of MonteCarlo simulations
 *************************************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define N 5000
#define Q 10000
#define Factor_graph  0
#define Nrunmax 2
#define Kcore  3
#define Ncore  2
#define mav  4.0
#define err_tol 0.01
int *vis1,*size_cluster1,**knn1,*k1,c1,c2,c3,*occ,*dam1,*dam2;


/**************************************************************************
 Recurrence is a subrutine for calulating the giant component
 **************************************************************************/
int Recurrence( int i , int cluster_size, int ncluster){
	int j, n3,aus_cluster_size,n,j2,aus,is,iaus;
    
    vis1[i]=ncluster;
    if (i<N){cluster_size++;}
    vis1[i]=ncluster;
    for(n3=0;n3<k1[i];n3++){
			j=knn1[i][n3];
            if(j<N){
                iaus=0;
                if(Factor_graph==0){iaus=1;}
                if((Factor_graph==1)&&(dam1[j]==1)){iaus=1;}
                if((vis1[j]==0)&&(iaus==1)){
				aus_cluster_size=Recurrence(j, cluster_size, ncluster);
				cluster_size=aus_cluster_size;
                }
                   }
            if (j>=N){
                if(dam1[j]>0){
                aus=1;
                
                if((aus==1)&&(vis1[j]==0)){
                    aus_cluster_size=Recurrence(j, cluster_size, ncluster);
                    cluster_size=aus_cluster_size;
                }
                }
                }
               
                   }
	
		return cluster_size;
}

int main(int argc, char** argv){
	int i,j,n,it,ncluster1,ncluster2,ncluster3, GC,nn,mu,*occ2,cluster_size,m1,m2,m3,m1_aus,m2_aus,m3_aus,c1_aus,c2_aus,c3_aus,*sigma,nrun,nc,np,GCold;
	int s1,s2,Nc1,Nc2,Nc3,aus,aus3,**adj1,**adj2,**adj3,N0,i_aus,ik,nd,nd2,xaus,n2,*q;
	float p,x,f,*xd1,*Sm,**n1,MGC1,nsum1,nsumold1,nsum2,aus1;
    int *GCm;
	
	   char filename[60],string1[50],string2[50];
	
	FILE *gp3,*gp;
    


  /**************************************************************************
  Generate the hypergraph
   **************************************************************************/

		srand48(time(NULL));
	    N0=N;
	
	

	vis1=(int*)calloc(N+Q,sizeof(int));
    occ2=(int*)calloc(N+Q,sizeof(int));
	k1=(int*)calloc(N+Q,sizeof(int));
    q=(int*)calloc(N+Q,sizeof(int));
	xd1=(float*)calloc(N+Q,sizeof(float));
	dam1=(int*)calloc(N+Q,sizeof(int));
    dam2=(int*)calloc(N+Q,sizeof(int));
	knn1=(int**)calloc(N+Q,sizeof(int*));
    GCm=(int*)calloc(100,sizeof(int));
		for(i=0;i<(N+Q);i++){
			knn1[i]=(int*)calloc(N+Q,sizeof(int));
        }
	size_cluster1=(int*)calloc(N+Q,sizeof(int));
	
    
    for(i=0;i<N+Q;i++){
        k1[i]=0;
    }
    for(mu=N;mu<N+Q;mu++){
        q[mu]=2;
        for(i=0;i<5*N;i++){
            if(drand48()<(mav-2)/((float)5*N))
                q[mu]++;
        }
    }

  

	for(mu=0;mu<Q;mu++){
        for (n=0;n<q[N+mu];n++){
            aus=0;
            while(aus==0){
            i=(int)((float)N*drand48());
                if(mu+N>N+Q) printf("problem\n");
            aus=1;
            for (nn=0;nn<n;nn++){
                if(occ2[nn]==i) aus=0;
            }
            }
            occ2[n]=i;

				
            k1[i]++;
            k1[mu+N]++;
            
            knn1[i][k1[i]-1]=mu+N;
            knn1[mu+N][k1[mu+N]-1]=i;
            
    
        }
    }
		
 
    
				
/*%%%%%%%%%%%%SIMULATION%%%%%%%%%%%%%%%%%*/
    /*Pruning algorithm*/
	for (nrun=0;nrun<Nrunmax;nrun++){	
    
        for(i=N;i<(N+Q);i++){
            xd1[i]=drand48();
        }
        for(i=0;i<N;i++){
            dam1[i]=1;
        }
        
        
        for(nc=0;nc<100;nc++){
		f=nc*0.01;
            nd=0;
            for(i=N;i<(N+Q);i++){
                if(xd1[i]<f){
                    dam1[i]=0;}
                else {dam1[i]=1;nd=nd+1;}
            }
            aus=1;
            while(aus==1){
                nd2=0;
                
            for(i=N;i<(N+Q);i++){
                dam2[i]=dam1[i];
                if(dam1[i]==1){
                ik=0;
                for(n=0;n<k1[i];n++){
                    j=knn1[i][n];
                    xaus=0;
                        for(n2=0;n2<k1[j];n2++){
                        xaus=xaus+dam1[knn1[j][n2]];
                    }
                    if(xaus>(Kcore-1)){
                        ik=ik+1;
                    }
                }
                
                if(ik<Ncore){
                    dam2[i]=0;
                }
                else{dam2[i]=1;nd2=nd2+1;}
                }
            }
                if(nd2==nd)aus=0;
                nd=nd2;
                for (i=N;i<(N+Q);i++){
                    dam1[i]=dam2[i];
                }
                }
            
            if(Factor_graph==1){
            for(i=0;i<N;i++){
                xaus=0;
                    for(n2=0;n2<k1[i];n2++){
                    xaus=xaus+dam1[knn1[i][n2]];
                }
                
                if(xaus<(Kcore)){
                    dam1[i]=0;
                }
            }
            }
            
            /*Check nodes in the giant component*/
		for(i=0;i<(N+Q);i++){
			vis1[i]=0;
		}
 
        
        m1=0;
		ncluster1=0;
		for(n=0;n<(N+Q);n++){
			if((vis1[n]==0)&&(dam1[n]==1)){
				cluster_size=0;
				ncluster1++;
				cluster_size=Recurrence(n, cluster_size, ncluster1);
				size_cluster1[ncluster1]=cluster_size;
				if(cluster_size>m1){m1=cluster_size;
					c1=ncluster1;}				
			}			
        }
               
        
               GCold=GC;
               Nc1=c1;
               GC=0;
               for(i=0;i<(N);i++){
                   if((vis1[i]==Nc1)&&(dam1[i]==1)){
                       GC++;
                      
                   }
               }
            
        GCm[nc]+=GC;
        }
    }

    gp3=fopen("Kcore_Hypergraph_2nd_neighbor_Hyperedge.txt","w");
        for (nc=0;nc<100;nc++){
            f=nc*0.01;
        
        fprintf(gp3,"%f %f\n",1-f,(float)GCm[nc]/(float)(Nrunmax*N));
        }
		fclose(gp3);
	
	return 0;	
}

