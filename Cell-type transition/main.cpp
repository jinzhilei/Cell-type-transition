//
//  main.cpp
//  StemCell
//
//  Created by Rongsheng Huang on 2022/12/1.
//

#include "Header.h"

Random _Rand;
struct Parameter _Parameter;

//Model parameters
void SetDefault(int i ,int j, int k)
{
    _Parameter.InitialCellNumber=100;
    
    _Parameter.dt=0.25;
    
    _Parameter.MaxTime=4000;
    
    _Parameter.index1=i;

    _Parameter.index2=j;

    _Parameter.index3=k;

    
    _Parameter.sigma=0.05;
    
   
    _Parameter.ODE_pou1=0.1;
    _Parameter.ODE_pou2=0.1;
    _Parameter.ODE_b1=1;
    _Parameter.ODE_b2=1;
    _Parameter.ODE_k1=1;
    _Parameter.ODE_k2=1;
    _Parameter.ODE_s1=0.5;
    _Parameter.ODE_s2=0.5;
    _Parameter.ODE_n=4;
    
    _Parameter.ODE_alpha1=0.4;
    _Parameter.ODE_alpha2=0.4;
    _Parameter.ODE_lambda1=1.9;
    _Parameter.ODE_lambda2=1.9;

    _Parameter.ODE_zeta1=1;
    _Parameter.ODE_zeta2=1;

   
    
    
    _Parameter.mu0=0.0002;
    _Parameter.kappa0=0.01;
    _Parameter.beta0=0.04;
    _Parameter.tau0=4;
    _Parameter.theta0=200;


    _Parameter.m0=60;
    _Parameter.m1=0.08;
    _Parameter.m2=0.8+0.025*(_Parameter.index2-1);
    _Parameter.m3=2;
    _Parameter.m4=2;
    
    
    strcpy(_Parameter.FileOutputPath,PATHOUTPUT);//Output file location.
    
};


void Test()
{
    int N=30;//样本轨道数量

    int M=5;//参数变化数目
    
    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=M;j++)
        {
            SetDefault(i,j,0);

            System *sys;
        
            sys = new System();
       
            sys->Initialized();
        
            sys->Run();
        }

        
    }

};


void Backtracking();

int main(int argc, const char * argv[]) {
    // insert code here...
    cout << "Hello, World!\n";
    
    cout<<_Rand.Value()<<endl;
    
    Test();

   // Backtracking();
    
  
};



void Backtracking()
{
    //SetDefault();
    char line[COMLENGTH];
    char FileInName[COMLENGTH];
    snprintf(FileInName,COMLENGTH,"%s",_Parameter.FileOutputPath);

    char FileOutName[COMLENGTH];
    snprintf(FileOutName,COMLENGTH,"%s",_Parameter.FileOutputPath);
  
    int NumCell=1;
    int NumStep=17880;
    int NumColumn=11;
    double temp;
    int LineCode=1;
    
    for(int i=1;i<=NumCell;i++)
    {
        LineCode=i;
        LineCode=3004;
        ofstream FileOutCell(FileOutName+to_string(i)+".cell");
        for(int step=0;step<NumStep;step++)
        {
            ifstream FileIn;
            FileIn.open(FileInName+to_string(NumStep-step)+".step");
            
            if(LineCode==1)
            {
                for(int k=1;k<=NumColumn;k++)
                {
                    FileIn>>temp;
                    FileOutCell<<temp<<" ";
                    if(k==8)
                    {
                        LineCode=temp;
                    }
                }
                FileOutCell<<endl;
                
            }else
            {
                for(int j=1;j<LineCode;j++)
                {
                    FileIn.getline(line,COMLENGTH,'\n');
                }
                
                for(int k=1;k<=NumColumn;k++)
                {
                    FileIn>>temp;
                    FileOutCell<<temp<<" ";
                    if(k==8)
                    {
                        LineCode=temp;
                    }
                }
                FileOutCell<<endl;
            }
            FileIn.close();
        }
        FileOutCell.close();
    }
    
    
    
    
};




