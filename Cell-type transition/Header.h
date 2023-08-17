//
//  Header.h
//  StemCell
//
//  Created by Rongsheng Huang on 2022/12/1.
//

#ifndef Header_h
#define Header_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "Random.h"
using namespace std;

//Definition of functions
#define _SQUARE(x) ((x)*(x))
#define _MAX(x,y) ((x)>=(y) ? (x):(y))
#define _MIN(x,y) ((x)<=(y) ? (x):(y))
#define _Heaviside(x) ((x)>(0) ? (1.0):(0.0))

//Definition of constants
#define COMLENGTH 100
#define StrLength 1024
#define UNITMAX 2147483647
#define MAXCELL 12000
#define NUMGENE 2
#define NUMEPI 2
#define NUMCELLINFOINT 5
#define NUMCELLINFODOUBLE 4
#define NUMSYSTEMINFOINT 4
#define NUMSYSTEMINFODOUBLE 10
#define NUMCELLTYPE 3
#define NUMRAND 2
#define PI 3.1415926
#define POSITIVE 0
#define CODESC 0
#define CODETA1 1
#define CODETA2 3
#define PATHOUTPUT "output/"


extern Random _Rand;
extern struct Parameter _Parameter;

//Model parameters
struct Parameter{
    
   // int seed;                // Seed of the random numbers.
    
    int index1,index2,index3;
    
    char FileOutputPath[COMLENGTH];
    int InitialCellNumber;
    double dt;
    double MaxTime;
    
    double sigma;
    
    double ODE_b1,ODE_b2,ODE_pou1,ODE_pou2,ODE_k1,ODE_k2,ODE_s1,ODE_s2,ODE_alpha1,ODE_alpha2,ODE_lambda1,ODE_lambda2,ODE_zeta1,ODE_zeta2,ODE_n;
    
    double mu0,kappa0,beta0,tau0,theta0;
    
    double m0,m1,m2,m3,m4;
  
    
};
struct DaughterCell{

    double Daughter_1_Gene[NUMGENE];//Gene expression variables of daughter cell 1.
    double Daughter_2_Gene[NUMGENE];//Gene expression variables of daughter cell 2.
    double Daughter_1_Epi[NUMEPI];//Epigenetic variables of daughter cell 1.
    double Daughter_2_Epi[NUMEPI];//Epigenetic variables of daughter cell 2.
    
    int Daughter_CellInfoInt[NUMCELLINFOINT];//Intracellular information variable (int type).
    double Daughter_CellInfoDouble[NUMCELLINFODOUBLE];//Intracellular information variable (double type).
};


class Cell{
public:
    
    double Gene[NUMGENE];//Gene expression variables.
    double Epi[NUMEPI];//Epigenetic variables.
    
    int CellInfoInt[NUMCELLINFOINT];//Intracellular information variable (int type).
    double CellInfoDouble[NUMCELLINFODOUBLE];//Intracellular information variable (double type).
    
    //CellInfoInt[0] denote Cell stage；
                    //-1:Cell death;
                    //0:The original resting period is still in the resting period.
                    //1:The original resting period enters the proliferative period,
                    //2:Original proliferating stage, still proliferating stage.
                    //3:Proproliferative phase, entering the resting phase (division).
    //CellInfoInt[1] denote the current cell type： 0 for SC；1 for TA1；2 for TA2
    //CellInfoInt[2] denote the parent cell type： 0 for SC；1 for TA1；2 for TA2
    //CellInfoInt[3] denote the number of the previous step of the cell
    //CellInfoInt[4] cell division number; 

    //CellInfoDouble[0] Cell proliferation stage age
    //CellInfoDouble[1] Cell age
    //CellInfoDouble[2] eta_1
    //CellInfoDouble[3] eta_2
    
    
    double ODE_a1,ODE_a2;
    
    double GeneTrend[NUMGENE];
    double eta_1,eta_2;
    
    double mu,kappa,beta,tau;
    double r1,r2;


    double Beta_a,Beta_b;
    double PXY_phi;

    
    Cell()
    {
        
    };
    ~Cell()
    {
    };
    
    void Initialized()
    {
    
        for(int i=0;i<NUMGENE;i++)
        {
            Gene[i]=3.0*_Rand.Value();
        };
      
    
        for(int i=0;i<NUMEPI;i++)
        {
            Epi[i]=_Rand.Value();
        };
      
        for(int i=0;i<NUMCELLINFOINT;i++)
        {
            CellInfoInt[i]=0;
        };
        for(int i=0;i<NUMCELLINFODOUBLE;i++)
        {
            CellInfoDouble[i]=0;
        }

        UpdateCellType();
    };


    void RunOneStep()
    {
      
        eta_1=CellInfoDouble[2];
        eta_2=CellInfoDouble[3];

        eta_1=eta_1+(-eta_1*1.0/_Parameter.ODE_zeta1)*_Parameter.dt+sqrt(2.0*_Parameter.dt/_Parameter.ODE_zeta1)*_Rand.NormalDistribution();
        eta_2=eta_2+(-eta_2*1.0/_Parameter.ODE_zeta2)*_Parameter.dt+sqrt(2.0*_Parameter.dt/_Parameter.ODE_zeta2)*_Rand.NormalDistribution();

        CellInfoDouble[2]=eta_1;
        CellInfoDouble[3]=eta_2;


        ODE_a1=_Parameter.ODE_alpha1*exp(_Parameter.ODE_lambda1*Epi[0])*exp(_Parameter.sigma*eta_1-_Parameter.sigma*_Parameter.sigma*0.5);
        ODE_a2=_Parameter.ODE_alpha2*exp(_Parameter.ODE_lambda2*Epi[1])*exp(_Parameter.sigma*eta_2-_Parameter.sigma*_Parameter.sigma*0.5);
        
        GeneTrend[0]=(ODE_a1*1+0.0)*(_Parameter.ODE_pou1+(1-_Parameter.ODE_pou1)*pow(Gene[0],_Parameter.ODE_n)/(pow(_Parameter.ODE_s1,_Parameter.ODE_n)+pow(Gene[0],_Parameter.ODE_n)))+_Parameter.ODE_b1*pow(_Parameter.ODE_s2,_Parameter.ODE_n)/(pow(_Parameter.ODE_s2,_Parameter.ODE_n)+pow(Gene[1],_Parameter.ODE_n))-_Parameter.ODE_k1*Gene[0];
        GeneTrend[1]=(ODE_a2*1+0.0)*(_Parameter.ODE_pou2+(1-_Parameter.ODE_pou2)*pow(Gene[1],_Parameter.ODE_n)/(pow(_Parameter.ODE_s2,_Parameter.ODE_n)+pow(Gene[1],_Parameter.ODE_n)))+_Parameter.ODE_b2*pow(_Parameter.ODE_s1,_Parameter.ODE_n)/(pow(_Parameter.ODE_s1,_Parameter.ODE_n)+pow(Gene[0],_Parameter.ODE_n))-_Parameter.ODE_k2*Gene[1];
        
        for(int i=0;i<NUMGENE;i++)
        {
            Gene[i]=Gene[i]+GeneTrend[i]*_Parameter.dt;
        };
        
      
        if(CellInfoInt[0]==0||CellInfoInt[0]==3)//CellInfoInt[0]==0 or CellInfoInt[0]==3 means that cell in resting period 
        {
            UpdateCellType();
            
            if(CellInfoInt[1]==0)
            {
                CellInfoInt[4]=0;
            }
        }
        
    };
    
    void  UpdateCellType()
    {
        
        if(Gene[0]>0.5&&Gene[1]>0.5)
        {
            CellInfoInt[1]=CODESC;
        };
        if(Gene[0]>0.5&&Gene[1]<0.5)
        {
            CellInfoInt[1]=CODETA1;
        };
        if(Gene[0]<0.5&&Gene[1]>0.5)
        {
            CellInfoInt[1]=CODETA2;
        };
    }
 

    DaughterCell CellFateDecision(int NumSC)
    {
      
        DaughterCell NextCell;
    
        for(int i=0;i<NUMGENE;i++)
        {
            NextCell.Daughter_1_Gene[i]=Gene[i];
            NextCell.Daughter_2_Gene[i]=Gene[i];
        };
     
        for(int i=0;i<NUMEPI;i++)
        {
            NextCell.Daughter_1_Epi[i]=Epi[i];
            NextCell.Daughter_2_Epi[i]=Epi[i];
        };
   
        for(int i=0;i<NUMCELLINFOINT;i++)
        {
            NextCell.Daughter_CellInfoInt[i]=CellInfoInt[i];
        };
        for(int i=0;i<NUMCELLINFODOUBLE;i++)
        {
            NextCell.Daughter_CellInfoDouble[i]=CellInfoDouble[i];
        }
     
        r1=_Rand.Value();
        r2=_Rand.Value();
       
        switch(CellInfoInt[1])//CellInfoInt[1] cell type;；
        {
            case CODESC:
                mu=_Parameter.mu0;
                kappa=_Parameter.kappa0;
                beta=_Parameter.beta0/(1+1.0*NumSC/_Parameter.theta0);
                tau=_Parameter.tau0;
                break;
            case CODETA1:
            case CODETA2:
                mu=_Parameter.mu0;
                kappa=2*_Parameter.kappa0;
                beta=_Parameter.beta0;
                tau=_Parameter.tau0*0.5;
                break;
                
        }
        mu=mu*_Parameter.dt;
        kappa=kappa*_Parameter.dt;
        beta=beta*_Parameter.dt;

        switch (CellInfoInt[0])//CellInfoInt[0]  Cell stage；
        {
            case 0://The original resting period is still in the resting period.
            case 3://Proproliferative phase, entering the resting phase (cell division).
                if(r1<kappa)//Cell death in the resting phase
                {
                    NextCell.Daughter_CellInfoInt[0]=-1;
                }else
                {
                    if (r1<kappa+beta)//Cells in the resting phase enter the proliferative phase
                    {
                        NextCell.Daughter_CellInfoInt[0]=1;//Cell stage set 1
                        NextCell.Daughter_CellInfoDouble[0]=0;//Cell cycle age set to 0
                        
                    }else//The original resting period is still in the resting period.
                    {
                        NextCell.Daughter_CellInfoInt[0]=0;
                    }
                }
                break;
            case 1://The original resting period enters the proliferative period,
            case 2://Original proliferating stage, still proliferating stage.
                if(r2<mu)//Cell death in the proliferative phase
                {
                    NextCell.Daughter_CellInfoInt[0]=-1;
                }else
                {
                    if(CellInfoDouble[0]<tau)//Cells that are proliferating are still proliferating; CellInfoInt[0] denote Cell stage;CellInfoDouble[0]denote Cell proliferation stage age
                    {
                        NextCell.Daughter_CellInfoInt[0]=2;
                        NextCell.Daughter_CellInfoDouble[0]=CellInfoDouble[0]+_Parameter.dt;
                    }else//From the original growth phase to the resting phase (cell division)
                    {
                        NextCell.Daughter_CellInfoInt[0]=3;//CellInfoInt[0] Cell stage
                        
                        NextCell.Daughter_CellInfoInt[2]=CellInfoInt[1];//CellInfoInt[1] Current cell type，CellInfoInt[2] Parent cell type
                        
                        NextCell.Daughter_CellInfoInt[4]++;//CellInfoInt[4] Cell division number
                        
                        if(NextCell.Daughter_CellInfoInt[4]==15)
                        {
                            if(CellInfoInt[1]==CODETA1||CellInfoInt[1]==CODETA2)
                            {
                                NextCell.Daughter_CellInfoInt[0]=-1;
                            }
                        }
        
                        for(int i=0;i<NUMEPI;i++)
                        {
                            
                            PXY_phi = _Parameter.m1 +_Parameter.m2* pow(_Parameter.m3* Epi[i], _Parameter.m4) /(1+pow(_Parameter.m3 * Epi[i], _Parameter.m4));
            
                            Beta_a =_Parameter.m0 * PXY_phi;
                            Beta_b =_Parameter.m0 * (1-PXY_phi);
                            
                            NextCell.Daughter_1_Epi[i]=_Rand.BetaDistribution(Beta_a, Beta_b);
                            NextCell.Daughter_2_Epi[i]=_Rand.BetaDistribution(Beta_a, Beta_b);

                        };
                        
                    }
                }
                break;
        }
        
        
        return NextCell;
    };
};

class System{
public:
    
    Cell CellPool[MAXCELL];
    int SystemInfoInt[NUMSYSTEMINFOINT];//The information variables of the system (int type)
    double SystemInfoDouble[NUMSYSTEMINFODOUBLE];//The information variables of the system (double type)
    DaughterCell NextCell;
    
    int NumTotal;
    

    double GenePool[MAXCELL][NUMGENE];
    double EpiPool[MAXCELL][NUMEPI];
    double CellInfoIntPool[MAXCELL][NUMCELLINFOINT];
    double CellInfoDoublePool[MAXCELL][NUMCELLINFODOUBLE];
    

    System()
    {
    };
    ~System()
    {
    };

    void Initialized()
    {
        SystemInfoInt[0]=_Parameter.InitialCellNumber;//SystemInfoInt[0] The total number of cells in the cell pool
        
        for(int i=0;i<SystemInfoInt[0];i++)
        {
            CellPool[i].Initialized();
            
            CellPool[i].CellInfoInt[3]=i+1;//The initial number of the cell
        };
    };

    void Run()
    {
        double Time=0;
        int Step=0;

        char FileOutName[COMLENGTH];//Output file location
        snprintf(FileOutName,COMLENGTH,"%s",_Parameter.FileOutputPath);//Read the output file location
        
       //ofstream FileOutSys(FileOutName+to_string(0)+"/sys/"+to_string(_Parameter.index)+".sys");//Output system information
        
        ofstream FileOutSys("output/sys/"+to_string(_Parameter.index1)+"_"+to_string(_Parameter.index2)+".sys");// Output system information
        ofstream FileOutTran("output/tran/"+to_string(_Parameter.index1)+"_"+to_string(_Parameter.index2)+".tran");// Output the transition probability information
        
        for(Time=0;Time<_Parameter.MaxTime;Time=Time+_Parameter.dt)
        {
            Step++;

            //ofstream FileOutStep("output/step/"+to_string(_Parameter.index1)+"_"+to_string(_Parameter.index2)+"_"+to_string(Step)+".step");//Output information for each step

            NumTotal=SystemInfoInt[0];
            SystemInfoInt[1]=0;
            SystemInfoInt[2]=0;
            SystemInfoInt[3]=0;
            
            for(int i=0;i<NumTotal;i++)
            {
                 /*
                for(int j1=0;j1<NUMGENE;j1++)
                {
                    FileOutStep<<CellPool[i].Gene[j1]<<" ";
                }
                for(int j2=0;j2<NUMEPI;j2++)
                {
                    FileOutStep<<CellPool[i].Epi[j2]<<" ";
                }
                for(int j3=0;j3<NUMCELLINFOINT;j3++)
                {
                    FileOutStep<<CellPool[i].CellInfoInt[j3]<<" ";
                }
                for(int j4=0;j4<NUMCELLINFODOUBLE;j4++)
                {
                    FileOutStep<<CellPool[i].CellInfoDouble[j4]<<" ";
                }
                
                FileOutStep<<endl;
                */


                CellPool[i].RunOneStep();
                
                if(CellPool[i].CellInfoInt[0]==1&&Step>4*2000)
                {
                    FileOutTran<<CellPool[i].CellInfoInt[1]<<" "<<CellPool[i].CellInfoInt[2]<<" "<<Step<<endl;
                }
 
                switch (CellPool[i].CellInfoInt[1])
                {
                    case CODESC:
                        SystemInfoInt[1]++;
                        break;
                    case CODETA1:
                        SystemInfoInt[2]++;
                        break;
                    case CODETA2:
                        SystemInfoInt[3]++;
                        break;
                }
               
            };

          // FileOutStep.close();

            if(NumTotal==0)
                break;
            
            //Output the number of cells at each step
            FileOutSys<<Time<<" ";
            for(int i=0;i<NUMSYSTEMINFOINT;i++)
            {
                FileOutSys<<SystemInfoInt[i]<<" ";
            }
            FileOutSys<<endl;
            
            SystemUpdate(Time,Step,NumTotal);
           
        };
        
        FileOutSys.close();
        FileOutTran.close();
        
    };
    void  SystemUpdate(double Time, int Step,int NumTotal)
    {
        int NumTotalTemp=0;
        
        for(int i=0;i<NumTotal;i++)
        {
            NextCell=CellPool[i].CellFateDecision(SystemInfoInt[1]);
            
            switch(NextCell.Daughter_CellInfoInt[0]) //CellInfoInt[0] cell stage
            {
                case -1://Cell deat
                    break;
                case 0://The original resting period is still in the resting period.
                case 1://The original resting period enters the proliferative period,
                case 2://Original proliferating stage, still proliferating stage.
                case 3://Proproliferative phase, entering the resting phase (cell division).
                    
                //继承第一个子细胞的信息
                    for(int j1=0;j1<NUMGENE;j1++)
                    {
                        GenePool[NumTotalTemp][j1]=NextCell.Daughter_1_Gene[j1];
                    }
                    for(int j2=0;j2<NUMEPI;j2++)
                    {
                        EpiPool[NumTotalTemp][j2]=NextCell.Daughter_1_Epi[j2];
                    }
                    for(int j3=0;j3<NUMCELLINFOINT;j3++)
                    {
                        CellInfoIntPool[NumTotalTemp][j3]=NextCell.Daughter_CellInfoInt[j3];
                    }
                    for(int j4=0;j4<NUMCELLINFODOUBLE;j4++)
                    {
                        CellInfoDoublePool[NumTotalTemp][j4]=NextCell.Daughter_CellInfoDouble[j4];
                    }
                    
                    CellInfoIntPool[NumTotalTemp][3]=i+1;
                    NumTotalTemp++;
                    
                    if(NextCell.Daughter_CellInfoInt[0]==3)//cell division
                    {
                    
                        for(int j1=0;j1<NUMGENE;j1++)
                        {
                            GenePool[NumTotalTemp][j1]=NextCell.Daughter_2_Gene[j1];
                        }
                        for(int j2=0;j2<NUMEPI;j2++)
                        {
                            EpiPool[NumTotalTemp][j2]=NextCell.Daughter_2_Epi[j2];
                        }
                        for(int j3=0;j3<NUMCELLINFOINT;j3++)
                        {
                            CellInfoIntPool[NumTotalTemp][j3]=NextCell.Daughter_CellInfoInt[j3];
                        }
                        for(int j4=0;j4<NUMCELLINFODOUBLE;j4++)
                        {
                            CellInfoDoublePool[NumTotalTemp][j4]=NextCell.Daughter_CellInfoDouble[j4];
                        }
                        CellInfoIntPool[NumTotalTemp][3]=i+1;
                        NumTotalTemp++;
                    }
                    break;
            }
        };
        SystemInfoInt[0]=0;
        for(int i=0;i<NumTotalTemp && i<MAXCELL;i++)
        {
            SystemInfoInt[0]++;
            for(int j1=0;j1<NUMGENE;j1++)
            {
                CellPool[i].Gene[j1]=GenePool[i][j1];
            }
            for(int j2=0;j2<NUMEPI;j2++)
            {
                CellPool[i].Epi[j2]=EpiPool[i][j2];
            }
            for(int j3=0;j3<NUMCELLINFOINT;j3++)
            {
                CellPool[i].CellInfoInt[j3]=CellInfoIntPool[i][j3];
            }
            for(int j4=0;j4<NUMCELLINFODOUBLE;j4++)
            {
                CellPool[i].CellInfoDouble[j4]=CellInfoDoublePool[i][j4];
            }
        };
        
       
    };
};

#endif /* Header_h */

