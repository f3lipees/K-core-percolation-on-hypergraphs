/**************************************************************************************************
 * If you use this code, please cite:
 
 G. Bianconi and S. N. Dorogovstev 
 "Nature of hypergraph k-core percolation problems" 
 Physical Review E, 109, 014307 (2024).
 
***************************************************************************************************
 * Code that  generates random hypergraph and performs (2nd neighbour) node (k,s)-core percolation 
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
#define Q 20000
#define m  3
#define Hypergraph 1
#define Second_neighbour 0
#define Nrunmax 100
#define Kcore  2
#define Ncore  4
#define mav  4
#define err_tol 0.01
int *vis1,*size_cluster1,**knn1,*k1,c1,c2,c3,*occ,*dam1,*dam2;


/**************************************************************************
 Recurrence is a subrutine for calulating the giant component
 **************************************************************************/
int Recurrence( int i , int cluster_size, int ncluster){
	int j, n3,aus_cluster_size,n,j2,aus,is,iaus;
    
    vis1[i]=ncluster;
    iaus=Ncore-1;
    if(Second_neighbour==1)iaus=1;
    if(i<N)	{cluster_size++;}
		vis1[i]=ncluster;
    for(n3=0;n3<k1[i];n3++){
			j=knn1[i][n3];
            if(j<N){
                if((vis1[j]==0)&&(dam1[j]==1)){
				aus_cluster_size=Recurrence(j, cluster_size, ncluster);
				cluster_size=aus_cluster_size;
			}
            }
            if (j>=N){
                if(k1[j]>iaus){
                aus=1;
                if(Hypergraph==1){
              //      printf("ci sono\n");
                    for(n=0;n<k1[j];n++){
                        j2=knn1[j][n];
                        aus=aus*dam1[j2];
                    }
                }
                
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
	float p,x,f,**xd1,*Sm,**n1,MGC1,nsum1,nsumold1,nsum2,aus1;
    int *GCm;
	
	   char filename[60],string1[50],string2[50];
	
	FILE *gp3,*gp;
    


  /**************************************************************************
   Generate the hypergraph
   **************************************************************************/

		srand48(time(NULL));
	
	N0=N;
	
	

	vis1=(int*)calloc(N+Q,sizeof(int));
	occ=(int*)calloc(N+Q,sizeof(int));
    occ2=(int*)calloc(m,sizeof(int));
	k1=(int*)calloc(N+Q,sizeof(int));
    q=(int*)calloc(N+Q,sizeof(int));
	xd1=(float**)calloc(N,sizeof(float*));
	dam1=(int*)calloc(N,sizeof(int));
    dam2=(int*)calloc(N,sizeof(int));
	knn1=(int**)calloc(N+Q,sizeof(int*));
    n1=(float**)calloc(N+Q,sizeof(float*));
    GCm=(int*)calloc(100,sizeof(int));
		for(i=0;i<N+Q;i++){
			knn1[i]=(int*)calloc(N,sizeof(int));
            n1[i]=(float*)calloc(N+Q,sizeof(float));
        }
    for(i=0;i<N;i++){
			xd1[i]=(float*)calloc(N,sizeof(float));
		}
	size_cluster1=(int*)calloc(N+Q,sizeof(int));
	
    
    for(i=0;i<N+Q;i++){
        k1[i]=0;
    }
    for(mu=N;mu<N+Q;mu++){
        q[mu]=2;
        for(i=0;i<N;i++){
            if(drand48()<(mav-2)/(float)N)
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
    
        for(i=0;i<N;i++){
            xd1[i][0]=drand48();
        }
        GC=N;
        for(nc=0;nc<100;nc++){
		f=nc*0.01;
            nd=0;
            for(i=0;i<N;i++){
                if(xd1[i][0]<f){
                    dam1[i]=0;}
                else {dam1[i]=1;nd=nd+1;}
            }
            aus=1;
            while(aus==1){
                nd2=0;
            for(i=0;i<N;i++){
                dam2[i]=dam1[i];
                if(dam1[i]==1){
                ik=0;
                for(n=0;n<k1[i];n++){
                    j=knn1[i][n];
                    xaus=0;
                    if(k1[j]>(Ncore-1)){
                    if(Hypergraph==1){
                    xaus=1;
                        for(n2=0;n2<k1[j];n2++){
                        xaus=xaus*dam1[knn1[j][n2]];
                    }
                    }
                    if(Hypergraph==0){
                    xaus=-1;
                        for(n2=0;n2<k1[j];n2++){
                        xaus=xaus+dam1[knn1[j][n2]];
                    }
                    }
                    }
                    
                    if(xaus>0){
                        ik=ik+1;
                    }
                }
                
                if(ik<Kcore){
                    dam2[i]=0;
                }
                else{dam2[i]=1;nd2=nd2+1;}
                }
            }
            
                if(nd2==nd)aus=0;
                nd=nd2;
                for (i=0;i<N;i++){
                    dam1[i]=dam2[i];
                }
            }
        
        
            
		for(i=0;i<N+Q;i++){
			vis1[i]=0;
		}
    
		/*Check for the giant component*/
        m1=0;
		ncluster1=0;
		for(n=0;n<N+Q;n++){
			if((vis1[n]==0)&&(dam1[n]==1)){
				cluster_size=0;
				ncluster1++;
				cluster_size=Recurrence(n, cluster_size, ncluster1);
				size_cluster1[ncluster1]=cluster_size;
				if(cluster_size>m1){m1=cluster_size;
					c1=ncluster1;}				
			}			
        }
               
          /*      Nc1=c1;
                for(i=0;i<N;i++){
                    if((dam1[i]==1)&&(vis1[i]==c1)){
                    ik=0;
                    for(n=0;n<k1[i];n++){
                        j=knn1[i][n];
                        if(vis1[j]==c1){
                            ik=ik+1;
                        }
                    }
                    if(ik<Kcore){
                        dam1[i]=0;
                    }
                    else{dam1[i]=1;}
                    }
                }*/
               GCold=GC;
               Nc1=c1;
               GC=0;
               for(i=N;i<N+Q;i++){
                   if(vis1[i]==Nc1){
                       GC++;
                   }
               }

            
        GCm[nc]+=GC;
        }
    }

    gp3=fopen("Kcore_Hypergraph_Node_2nd.txt","w");
        for (nc=0;nc<100;nc++){
            f=nc*0.01;
        
        fprintf(gp3,"%f %f\n",1-f,(float)GCm[nc]/(float)(Nrunmax*Q));
        }
		fclose(gp3);
	
	return 0;	
}

