//
//  main.c
//  main8
//
//  Created by TokimotoJun on 2013/12/27.
//  Copyright (c) 2013年 Jun Tokimoto. All rights reserved.
//  fft部分も並列化させたもの。k.pの不具合も修正済み。読み込みの際の不具合も修正済み（これにより、完璧な一様系の再現に成功。)。さらにEc以下のみの要素のみ読み込み、計算することでメモリー及び計算時間の大幅な節約を施したもの。これが完全版。2014/2/4

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <mkl.h>
#define N 400
#define Nwide 800
#define P 200000

int main(int argc,char *argv[])
{
    
    FILE *input1,*input2,*output;
    int M;
    int i,j,z,l,alphayindex,alphazindex,nu;
    int Eccount=0;
    int alphay[60*N];
    int alphaz[60*N];
    double val,dt,dx;
    double x[2*Nwide],density[2*Nwide];
    double  ***E;
    double _Complex **u,**v,**udash,**vdash;
    double T=20.0;
    double Lcal=30.0;
    double Ly=29.176886;
    double Lz=29.176886;
    double _Complex delta0[2*Nwide],delta[2*Nwide],phase,ubar0,vbar0,f0,g0,ubar,vbar,ufft[2*Nwide],vfft[2*Nwide],deltcheck,plus;
    double Ec=15.0;
    double mu=0.1332031250000000;
    double g;
    double d=0.0;
    double beta=0.017262;
    double gbefore;
    double gafter;
    double scatterbefore=0.55;
    double scatterafter;
    double tin= 0.0;
    double tout ;
    double ky,kz;
    double absgap0,normu,normv,binorm,vpote,lambda,zeta,epsilon,cs,sn,dens;
    double normcheck[60*Nwide];
    double sum;
    double Ntotal=300.0;
    int alphamax=20;
    
    MKL_LONG status;
    DFTI_DESCRIPTOR_HANDLE my_desc_handle;
    char fname[500];
    
    
    scatterafter=atof(argv[1]);
    tout=fabs((scatterbefore-scatterafter))/0.1*1.0;
    
    /*printf("scatterafter=%f tout=%f",scatterafter,tout);*/
    
    sprintf(fname,"3DplaintdbdgharmonicEc%fNtotal%lfL%fscatterbefore%fscatterafter%fscattertime%fNdiv%domega%f.dat",Ec,Ntotal,2.0*Lcal,scatterbefore,scatterafter,tout-tin,N,beta);
    
    
    
    input1 = fopen("../3Dplainbdg/5test-New-E_3DplainbdgharmonicEc15.000000Ntotal300.000000L60.000000scatter0.550000Ndiv400omega0.185804error0.000050.dat","r");
    input2 = fopen("../3Dplainbdg/5test-New-uv_3DplainbdgharmonicEc15.000000Ntotal300.000000L60.000000scatter0.550000Ndiv400omega0.185804error0.000050.dat","r");
    output = fopen(fname,"w");
    
    M = 2*N;
    
    /*printf("確保開始\n");*/
    /******E(エネルギー固有値）の配列の確保*******/
    
    E=(double ***)malloc(sizeof(double **)*2*N*2);
    
    if (E==NULL) {
        printf("えらーーーー\n");
    }
    
    /*printf("確保開始\n");*/
    for (i=0; i<2*N*2; i++) {
        
        E[i]=(double **)malloc(sizeof(double *)*alphamax);
        if (E[i]==NULL) {
            printf("えらーーーー\n");
        }
        
        for (j=0; j<alphamax; j++) {
            
            E[i][j]=(double *)malloc(sizeof(double )*alphamax);
            if (E[i][j]==NULL) {
                printf("えらーーーー\n");
            }
        }
    }
    
    
    /*printf("確保終了\n");*/
    /**************Eの読み込み***************/
    
    printf("E読み込み開始\n");
    
    for(alphayindex=0; alphayindex < alphamax; alphayindex++){
        for(alphazindex=0;alphazindex <= alphayindex; alphazindex++){
            for(nu = 0; nu < 2*M; nu++){
                fscanf(input1,"%lf\n",&val);
                E[nu][alphayindex][alphazindex]= val ;
                
                if(val<0.0&&fabs(val)<Ec){
                    Eccount++;
                    alphay[Eccount-1]= alphayindex;
                    alphaz[Eccount-1]= alphazindex;
                    
                }
                
            }
        }
    }
    
    printf("E読み込み終了\n");
    
    printf("Eccount=%d\n",Eccount);
    
    
    /************u',v'用の配列の確保***************/
    
    udash=(double _Complex**)malloc(sizeof(double _Complex*)*Eccount);
    vdash=(double _Complex**)malloc(sizeof(double _Complex*)*Eccount);
    
    for(i=0;i<Eccount;i++){
        udash[i]=(double _Complex*)malloc(sizeof(double _Complex)*2*N);
        vdash[i]=(double _Complex*)malloc(sizeof(double _Complex)*2*N);
    }
    
    /************u,v用の配列の確保***************/
    
    u=(double _Complex**)malloc(sizeof(double _Complex*)*Eccount);
    v=(double _Complex**)malloc(sizeof(double _Complex*)*Eccount);
    
    for(i=0;i<Eccount;i++){
        u[i]=(double _Complex*)malloc(sizeof(double _Complex)*2*Nwide);
        v[i]=(double _Complex*)malloc(sizeof(double _Complex)*2*Nwide);
    }
    
    /************u,v初期化************************/
    
    for (nu=0; nu<Eccount; nu++) {
        for (j=0; j<2*Nwide; j++) {
            
            u[nu][j]=0.0;
            v[nu][j]=0.0;
            
        }
    }
    
    /**************u',v'の読み込み(Ec以下のものだけ読み込み)***************/
    
    i=0;
    
    printf("uv読み込み開始\n");
    
    for(alphayindex = 0; alphayindex < alphamax; alphayindex++){
        for(alphazindex = 0; alphazindex <= alphayindex; alphazindex++){
            for(nu=0; nu < 2*M ; nu++){
                
                if(nu<2*N){
                    
                    fscanf(input2,"%lf\n",&val);
                    
                    if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec) udash[i][N-1]=val;
                    
                    
                    for(j=1 ; j < N-1 ; j++){
                        fscanf(input2,"%lf\n",&val);
                        
                        
                        if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec){
                            udash[i][j+N-1]=val;
                            udash[i][N-1-j]=-1.0*val;
                        }
                        
                    }
                    
                    fscanf(input2,"%lf\n",&val);
                    
                    if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec) udash[i][0]=-1.0*val;
                    
                    
                    fscanf(input2,"%lf\n",&val);
                    
                    if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec) vdash[i][N-1]=val;
                    
                    
                    
                    for(j=1 ; j < N-1 ; j++){
                        fscanf(input2,"%lf\n",&val);
                        if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec){
                            vdash[i][j+N-1]=val;
                            vdash[i][N-1-j]=-1.0*val;
                        }
                    }
                    
                    fscanf(input2,"%lf\n",&val);
                    
                    if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec) vdash[i][0]=-1.0*val;
                    
                    
                    if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec) i++;
                    
                    /*printf("i=%d\n",i);*/
                    
                }else{
                    
                    fscanf(input2,"%lf\n",&val);
                    
                    if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec) udash[i][N-1]=val;
                    
                    
                    for(j=1 ; j < N-1 ; j++){
                        fscanf(input2,"%lf\n",&val);
                        
                        
                        if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec){
                            udash[i][j+N-1]=val;
                            udash[i][N-1-j]=val;
                        }
                        
                    }
                    
                    fscanf(input2,"%lf\n",&val);
                    
                    if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec) udash[i][0]=val;
                    
                    
                    fscanf(input2,"%lf\n",&val);
                    
                    if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec) vdash[i][N-1]=val;
                    
                    
                    
                    for(j=1 ; j < N-1 ; j++){
                        fscanf(input2,"%lf\n",&val);
                        if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec){
                            vdash[i][j+N-1]=val;
                            vdash[i][N-1-j]=val;
                        }
                    }
                    
                    fscanf(input2,"%lf\n",&val);
                    
                    if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec) vdash[i][0]=val;
                    
                    
                    if(E[nu][alphayindex][alphazindex]<0.0 && fabs(E[nu][alphayindex][alphazindex])<Ec) i++;
                }
            }
            
        }
        
    }
    
    printf("u'v'読み込み終了\n");
    
    /************u,vへの代入******************************/
    
    for(nu=0;nu<Eccount;nu++){
        
        for(j=0;j<2*N;j++){
            
            u[nu][j+N]=udash[nu][j];
            v[nu][j+N]=vdash[nu][j];
            
        }
        
    }
    
    
    /************時間発展用のパラメーターの用意***************/
    dt = T/P;
    dx = Lcal/(Nwide-1);
    
    for(j = 0; j<= 2*Nwide-1; j++) x[j] = -Nwide+dx*j;
    
    gbefore= 8.0*M_PI/(scatterbefore - sqrt(4*Ec/(M_PI*M_PI)));
    gafter =8.0*M_PI/(scatterafter - sqrt(4*Ec/(M_PI*M_PI)));
    
    g = gbefore;
    /*printf("g=%f\n",g);*/
    
    /******************normの計算(check用)*********************/
    
    for(i=0; i< Eccount; i++){
        
        normcheck[i] = 0.0;
        
        for(j=0; j <= Nwide-2; j++){
            
            normcheck[i] += (u[i][j]*conj(u[i][j]) + v[i][j]*conj(v[i][j])+u[i][j+1]*conj(u[i][j+1]) + v[i][j+1]*conj(v[i][j+1]))*dx/2.0*Ly*Lz*2.0;
            
        }
        
        
        /*printf(" check pre normloop[%d]=%.15f\n",i,normcheck[i]);*/
        
        /*for(j=0;j<2*N;j++){
         u[i][j] = u[i][j]/sqrt(normcheck[i]);
         v[i][j] = v[i][j]/sqrt(normcheck[i]);
         }*/
        
    }
    
    /*******************Time evo ***************************/
    for(l=0; l <= P-1; l++){
        
        if(l*dt<=tout){   g = g + (gafter - gbefore)/(tout - tin) * dt;}
        
        printf("g=%f l=%d\n",g,l);
        
        
        /************y,z方向のエネルギーパート1st**************/
        
        for (nu=0; nu<Eccount; nu++) {
            
            for (j=0; j<2*Nwide-2; j++) {
                
                ky= 2.0*M_PI*(alphay[nu])/Ly;
                kz= 2.0*M_PI*(alphaz[nu])/Lz;
                
                u[nu][j]= cexp(-I*0.5*0.5*dt*(ky*ky+kz*kz))*u[nu][j];
                v[nu][j]= cexp(I*0.5*0.5*dt*(ky*ky+kz*kz))*v[nu][j];
            }
            
        }
        /******************normの計算(check用)*********************/
        
        for(i=0; i< Eccount; i++){
            
            normcheck[i] = 0.0;
            
            for(j=0; j <= Nwide-2; j++){
                
                normcheck[i] += (u[i][j]*conj(u[i][j]) + v[i][j]*conj(v[i][j])+u[i][j+1]*conj(u[i][j+1]) + v[i][j+1]*conj(v[i][j+1]))*dx/2.0*Ly*Lz*2.0;
                
            }
            
            
            /*printf(" check 1stafter normloop[%d]=%.15f\n",i,normcheck[i]);*/
            
            /*for(j=0;j<2*N;j++){
             u[i][j] = u[i][j]/sqrt(normcheck[i]);
             v[i][j] = v[i][j]/sqrt(normcheck[i]);
             }*/
            
        }
        
#pragma omp parallel for private(i,absgap0,normu,normv,binorm,vpote,lambda,zeta,epsilon,phase,cs,sn,ubar0,vbar0,f0,g0,ubar,vbar)
        
        
        /**************ポテンシャルパート 1st*****************/
        for(j = 0; j <2*Nwide-2; j++){
            
            delta0[j] = 0.0;
            
            
            for(nu=0;nu<Eccount;nu++){
                
                if(alphaz[nu]==0&&alphay[nu]==0){
                    
                    delta0[j] = delta0[j] + u[nu][j]*conj(v[nu][j]);
                    
                }
                
                else if(alphaz[nu]==alphay[nu] || (alphaz[nu]==0 || alphay[nu]==0)){
                    
                    delta0[j] = delta0[j] + 4.0*u[nu][j]*conj(v[nu][j]);
                    
                }else {
                    
                    delta0[j] = delta0[j] + 8.0*(u[nu][j]*conj(v[nu][j]));
                    
                }
                
            }
            
            
            delta0[j] *= g;
            
            absgap0 =creal(delta0[j])*creal(delta0[j]) + cimag(delta0[j])*cimag(delta0[j]);
            
            normu = 0.0;
            normv = 0.0;
            
            for (nu=0; nu<Eccount; nu++) {
                
                if(alphay[nu]==0&&alphaz[nu]==0){
                    normu += u[nu][j]*conj(u[nu][j]);
                    normv += v[nu][j]*conj(v[nu][j]);
                    
                } else if(alphaz[nu]==alphay[nu] || (alphaz[nu]==0 || alphay[nu]==0)){
                    normu += 4.0*u[nu][j]*conj(u[nu][j]);
                    normv += 4.0*v[nu][j]*conj(v[nu][j]);
                }else{
                    normu += 8.0*u[nu][j]*conj(u[nu][j]);
                    normv += 8.0*v[nu][j]*conj(v[nu][j]);
                }
                
            }
            
            
            binorm = -1.0*g*(normu - normv);
            vpote = 0.5*beta*(dx*j-Lcal)*(dx*j-Lcal);
            lambda =  2.0*(vpote - mu) + binorm;
            zeta =  vpote - mu - 0.5*lambda;
            epsilon = sqrt(zeta*zeta + absgap0);
            if(fabs(creal(delta0[j]))>1.0E-10 && fabs(cimag(delta0[j])) > 1.0E-10)
                phase = delta0[j]/sqrt(absgap0);
            else phase = 1.0;
            if (fabs(epsilon)<1.0E-10) {
                cs = sqrt(0.5*(1.0 + 1.0));
                sn = sqrt(0.5*(1.0 - 1.0));
            }else{
                cs = sqrt(0.5*(1.0 + zeta/epsilon));
                sn = sqrt(0.5*(1.0 - zeta/epsilon));
            }
            
            
            delta[j] = cexp(-I*lambda*dt/2.0)*delta0[j];
            
            
            for (nu=0; nu<Eccount; nu++) {
                
                
                ubar0 = u[nu][j];
                vbar0 = v[nu][j];
                f0 = cs*ubar0/phase + sn*vbar0;
                g0 = -sn*ubar0/phase + cs*vbar0;
                /*printf("i=%d,j=%d,realf0=%f,realg0[i][j]=%f\n",i,j,creal(f0),creal(g0));*/
                ubar = phase*(cs*cexp(-I*epsilon*dt/2.0)*f0 - sn*cexp(I*epsilon*dt/2.0)*g0);
                vbar = sn*cexp(-I*epsilon*dt/2.0)*f0 + cs*cexp(I*epsilon*dt/2.0)*g0;
                /*printf("i=%d,j=%d,realubar[i][j]=%f,realvbar[i][j]=%f\n",i,j,creal(ubar),creal(vbar));*/
                u[nu][j] = cexp(-0.5*I*lambda*dt/2.0)*ubar;
                v[nu][j] = cexp(0.5*I*lambda*dt/2.0)*vbar;
                /* printf("i=%d,j=%d,realu[i][j]=%f,realv[i][j]=%f\n",i,j,creal(u[i][j]),creal(v[i][j]));*/
                
                
            }
            
        }
        
        /******************normの計算(check用)*********************/
        
        for(i=0; i< Eccount; i++){
            
            normcheck[i] = 0.0;
            
            for(j=0; j <= Nwide-2; j++){
                
                normcheck[i] += (u[i][j]*conj(u[i][j]) + v[i][j]*conj(v[i][j])+u[i][j+1]*conj(u[i][j+1]) + v[i][j+1]*conj(v[i][j+1]))*dx/2.0*Ly*Lz*2.0;
                
            }
            
            
            /*printf(" check p-after normloop[%d]=%.15f\n",i,normcheck[i]);*/
            
            /*for(j=0;j<2*N;j++){
             u[i][j] = u[i][j]/sqrt(normcheck[i]);
             v[i][j] = v[i][j]/sqrt(normcheck[i]);
             }*/
            
        }
        
        /************y,z方向のエネルギーパート2nd**************/
        
        for (nu=0; nu<Eccount; nu++) {
            
            for (j=0; j<2*Nwide-2; j++) {
                
                ky= 2.0*M_PI*(alphay[nu])/Ly;
                kz= 2.0*M_PI*(alphaz[nu])/Lz;
                
                u[nu][j]= cexp(-I*0.5*0.5*dt*(ky*ky+kz*kz))*u[nu][j];
                v[nu][j]= cexp(I*0.5*0.5*dt*(ky*ky+kz*kz))*v[nu][j];
            }
            
        }
        
        /*************運動量パートwith FFT******************/
#pragma omp parallel for private(j,ufft,vfft,status,my_desc_handle)
        for (nu=0; nu<Eccount; nu++) {
            
            
            for(j=0;j<2*Nwide-2;j++){
                
                ufft[j]=u[nu][j];
                vfft[j]=v[nu][j];
                
            }
            
            
            status = DftiCreateDescriptor(&my_desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1,2*Nwide-2);
            status = DftiCommitDescriptor(my_desc_handle);
            status = DftiComputeForward(my_desc_handle, ufft);
            status = DftiFreeDescriptor(&my_desc_handle);
            
            status = DftiCreateDescriptor(&my_desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1,2*Nwide-2);
            status = DftiCommitDescriptor(my_desc_handle);
            status = DftiComputeForward(my_desc_handle, vfft);
            status = DftiFreeDescriptor(&my_desc_handle);
            
            
            
            for(j=0; j<Nwide; j++){
                
                ufft[j] = cexp(-I*dt*(2.0*M_PI*j/(2.0*Lcal))*(2.0*M_PI*j/(2.0*Lcal)))*ufft[j];
                vfft[j] = cexp(I*dt*(2.0*M_PI*j/(2.0*Lcal))*(2.0*M_PI*j/(2.0*Lcal)))*vfft[j];
                
            }
            
            for(j=Nwide; j<2*Nwide-2; j++){
                
                ufft[j] = cexp(-I*dt*(2.0*M_PI*(j-(2.0*Nwide-2))/(2.0*Lcal))*(2.0*M_PI*(j-(2.0*Nwide-2))/(2.0*Lcal)))*ufft[j];
                vfft[j] = cexp(I*dt*(2.0*M_PI*(j-(2.0*Nwide-2))/(2.0*Lcal))*(2.0*M_PI*(j-(2.0*Nwide-2))/(2.0*Lcal)))*vfft[j];
                
            }
            
            
            status = DftiCreateDescriptor(&my_desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1,2*Nwide-2);
            status = DftiSetValue( my_desc_handle, DFTI_BACKWARD_SCALE, 1.0/(2.0*(double)(Nwide)-2));
            status = DftiCommitDescriptor(my_desc_handle);
            status = DftiComputeBackward(my_desc_handle, ufft);
            status = DftiFreeDescriptor(&my_desc_handle);
            
            status = DftiCreateDescriptor(&my_desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1,2*Nwide-2);
            status = DftiSetValue( my_desc_handle, DFTI_BACKWARD_SCALE, 1.0/(2.0*(double)(Nwide)-2));
            status = DftiCommitDescriptor(my_desc_handle);
            status = DftiComputeBackward(my_desc_handle, vfft);
            status = DftiFreeDescriptor(&my_desc_handle);
            
            for(j=0; j < 2*Nwide-2 ; j++){
                
                u[nu][j] =  ufft[j];
                v[nu][j] =  vfft[j];
                
            }
            
        }
        
        
        /**************normの計算(check用)*****************/
        
        /*for(i=0; i< Eccount; i++){
         
         normcheck[i] = 0.0;
         
         for(j=0; j < 2*N-2; j++){
         
         normcheck[i] += (u[i][j]*conj(u[i][j]) + v[i][j]*conj(v[i][j])+u[i][j+1]*conj(u[i][j+1]) + v[i][j+1]*conj(v[i][j+1]))*dx/2.0*Ly*Lz;
         
         }
         
         printf("normafterfft[i]=%f\n",creal(normcheck[i]));
         
         for(j=0;j<2*N;j++){
         u[i][j] = u[i][j]/sqrt(normcheck[i]);
         v[i][j] = v[i][j]/sqrt(normcheck[i]);
         }
         
         }*/
        
        
        /************y,z方向のエネルギーパート3rd**************/
        
        for (nu=0; nu<Eccount; nu++) {
            
            for (j=0; j<2*Nwide-2; j++) {
                
                ky= 2.0*M_PI*(alphay[nu])/Ly;
                kz= 2.0*M_PI*(alphaz[nu])/Lz;
                
                u[nu][j]= cexp(-I*0.5*0.5*dt*(ky*ky+kz*kz))*u[nu][j];
                v[nu][j]= cexp(I*0.5*0.5*dt*(ky*ky+kz*kz))*v[nu][j];
            }
            
        }
        
        
#pragma omp parallel for private(i,absgap0,normu,normv,binorm,vpote,lambda,zeta,epsilon,phase,cs,sn,ubar0,vbar0,f0,g0,ubar,vbar)
        /**************ポテンシャルパート 2nd****************/
        for(j = 0; j <2*Nwide-2; j++){
            
            delta0[j] = 0.0;
            
            
            for(nu=0;nu<Eccount;nu++){
                
                if(alphaz[nu]==0&&alphay[nu]==0){
                    
                    delta0[j] = delta0[j] + u[nu][j]*conj(v[nu][j]);
                    
                }
                
                else if(alphaz[nu]==alphay[nu] || (alphaz[nu]==0 || alphay[nu]==0)){
                    
                    delta0[j] = delta0[j] + 4.0*u[nu][j]*conj(v[nu][j]);
                    
                }else {
                    
                    delta0[j] = delta0[j] + 8.0*(u[nu][j]*conj(v[nu][j]));
                    
                }
                
                
            }
            
            
            delta0[j] *= g;
            
            
            
            absgap0 =creal(delta0[j])*creal(delta0[j]) + cimag(delta0[j])*cimag(delta0[j]);
            
            normu = 0.0;
            normv = 0.0;
            
            for (nu=0; nu<Eccount; nu++) {
                
                if(alphay[nu]==0&&alphaz[nu]==0){
                    normu += u[nu][j]*conj(u[nu][j]);
                    normv += v[nu][j]*conj(v[nu][j]);
                    
                } else if(alphaz[nu]==alphay[nu] || (alphaz[nu]==0 || alphay[nu]==0)){
                    normu += 4.0*u[nu][j]*conj(u[nu][j]);
                    normv += 4.0*v[nu][j]*conj(v[nu][j]);
                }else{
                    normu += 8.0*u[nu][j]*conj(u[nu][j]);
                    normv += 8.0*v[nu][j]*conj(v[nu][j]);
                }
                
            }
            
            
            
            binorm = -1.0*g*(normu - normv);
            vpote = 0.5*beta*(dx*j-Lcal)*(dx*j-Lcal);
            lambda =  2.0*(vpote - mu) + binorm;
            zeta =  vpote - mu - 0.5*lambda;
            epsilon = sqrt(zeta*zeta + absgap0);
            if(fabs(creal(delta0[j]))>1.0E-10 && fabs(cimag(delta0[j])) > 1.0E-10)
                phase = delta0[j]/sqrt(absgap0);
            else phase = 1.0;
            if (fabs(epsilon)<1.0E-10) {
                cs = sqrt(0.5*(1.0 + 1.0));
                sn = sqrt(0.5*(1.0 - 1.0));
            }else{
                cs = sqrt(0.5*(1.0 + zeta/epsilon));
                sn = sqrt(0.5*(1.0 - zeta/epsilon));
            }
            
            
            
            delta[j] = cexp(-I*lambda*dt/2.0)*delta0[j];
            
            
            for (nu=0; nu<Eccount; nu++) {
                
                
                ubar0 = u[nu][j];
                vbar0 = v[nu][j];
                /* printf("i=%d,j=%d,realubar0=%f,realvbar0=%f\n",i,j,creal(ubar0),creal(vbar0)); */
                f0 = cs*ubar0/phase + sn*vbar0;
                g0 = -sn*ubar0/phase + cs*vbar0;
                /*printf("i=%d,j=%d,realf0=%f,realg0[i][j]=%f\n",i,j,creal(f0),creal(g0));*/
                ubar = phase*(cs*cexp(-I*epsilon*dt/2.0)*f0 - sn*cexp(I*epsilon*dt/2.0)*g0);
                vbar = sn*cexp(-I*epsilon*dt/2.0)*f0 + cs*cexp(I*epsilon*dt/2.0)*g0;
                /*printf("i=%d,j=%d,realubar[i][j]=%f,realvbar[i][j]=%f\n",i,j,creal(ubar),creal(vbar));*/
                u[nu][j] = cexp(-0.5*I*lambda*dt/2.0)*ubar;
                v[nu][j] = cexp(0.5*I*lambda*dt/2.0)*vbar;
                /* printf("i=%d,j=%d,realu[i][j]=%f,realv[i][j]=%f\n",i,j,creal(u[i][j]),creal(v[i][j]));*/
                
                
            }
            
        }
        
        /************y,z方向のエネルギーパート4th**************/
        
        for (nu=0; nu<Eccount; nu++) {
            
            for (j=0; j<2*Nwide-2; j++) {
                
                ky= 2.0*M_PI*(alphay[nu])/Ly;
                kz= 2.0*M_PI*(alphaz[nu])/Lz;
                
                u[nu][j]= cexp(-I*0.5*0.5*dt*(ky*ky+kz*kz))*u[nu][j];
                v[nu][j]= cexp(I*0.5*0.5*dt*(ky*ky+kz*kz))*v[nu][j];
            }
            
        }
        
        /******************normの計算(check用)*********************/
        
        for(i=0; i< Eccount; i++){
            
            normcheck[i] = 0.0;
            
            for(j=0; j <= Nwide-2; j++){
                
                normcheck[i] += (u[i][j]*conj(u[i][j]) + v[i][j]*conj(v[i][j])+u[i][j+1]*conj(u[i][j+1]) + v[i][j+1]*conj(v[i][j+1]))*dx/2.0*Ly*Lz*2.0;
                
            }
            
            
            /*printf(" check after normloop[%d]=%.15f\n",i,normcheck[i]);*/
            
            /*for(j=0;j<2*N;j++){
             u[i][j] = u[i][j]/sqrt(normcheck[i]);
             v[i][j] = v[i][j]/sqrt(normcheck[i]);
             }*/
            
        }
        /*************オーダーパラメーター(delta)の計算********************/
        for(j=0;j<2*Nwide-2;j++){
            
            delta[j] = 0.0;
            
            for(nu=0;nu<Eccount;nu++){
                if(alphaz[nu]==0&&alphay[nu]==0){
                    
                    delta[j] = delta[j] + u[nu][j]*conj(v[nu][j]);
                    
                }
                
                else if(alphaz[nu]==alphay[nu] || (alphaz[nu]==0 || alphay[nu]==0)){
                    
                    delta[j] = delta[j] + 4.0*u[nu][j]*conj(v[nu][j]);
                    
                }else {
                    delta[j] = delta[j] + 8.0*(u[nu][j]*conj(v[nu][j]));
                }
                
                
            }
            
            
            delta[j] *= g;
            
            /*printf("all after  j=%d delt=%f\n",j,creal(sqrt(delta[j]*conj(delta[j]))));*/
            
        }
        
        printf("all after delt=%f\n",creal(sqrt(delta[800]*conj(delta[800]))));
        
        /*************粒子数密度の計算***********************/
        for(j=0;j<Nwide;j++){
            
            dens = 0.0;
            
            
            for(nu=0;nu<Eccount;nu++){
                
                if (alphay[nu]==0&&alphaz[nu]==0) {
                    dens += u[nu][j]*conj(u[nu][j]);
                }else if(alphaz[nu]==alphay[nu] || (alphaz[nu]==0 || alphay[nu]==0)){
                    dens += 4.0*u[nu][j]*conj(u[nu][j]);
                }else{
                    dens += 8.0*u[nu][j]*conj(u[nu][j]);
                }
            }
            
            density[j]= 2.0*dens;
            
            
        }
        
        /*************ファイルへの書き込み*******************/
        if (l % 50 == 0){
            
            for(j=0; j <Nwide; j++){
                
                fprintf(output,"%f\t%f\t%f\t%f\n",dx*j-Lcal,creal(delta[j]),cimag(delta[j]),density[j]);
                
            }
            
            fprintf(output,"\n");
        }
        
        
        
    }
    
    
    return 0;
    
}
