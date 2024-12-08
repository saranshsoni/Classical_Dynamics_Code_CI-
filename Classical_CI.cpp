#include<bits/stdc++.h>
using namespace std;

std::ofstream outFile("output.txt");
std::ofstream poptrans_10_01("poptrans_10_01.txt");
std::ofstream poptrans_00_01("poptrans_00_01.txt");
std::ofstream poptrans_10_00("poptrans_10_00.txt");
std::ofstream poptrans_01_00("poptrans_01_00.txt");
std::ofstream poptrans_00_10("poptrans_00_10.txt");
std::ofstream poptrans_01_10("poptrans_01_10.txt");
std::ofstream poptrans("pop_trans.txt");
std::ofstream tempout("temp_out.txt");

int n_traj=2000;
const int timesteps=10000;
int P01_avg[timesteps]={};
int P10_avg[timesteps]={};
int P00_avg[timesteps]={};
int Q2_avg[timesteps]={};
int pop_trans_00_01[timesteps]={};
int pop_trans_00_10[timesteps]={};
int pop_trans_10_00[timesteps]={};
int pop_trans_10_01[timesteps]={};
int pop_trans_01_00[timesteps]={};
int pop_trans_01_10[timesteps]={};
long double mass=1.67*pow(10,-27);
long double c=3*pow(10,8);
long double h=6.6264*pow(10,-34);
long double h_=h/(2*3.14);
long double tot_en[timesteps]={};

//int h=1;
long double PI = 3.14;

//angular frequencies
long double f1=29000*c*2*PI;            //Q1
long double f2=21000*c*2*PI;            //Q2
long double f3=280000*c*2*PI;           //q1
long double f4=165000*c*2*PI;           //q2

//position scaling factors

long double Q10=sqrt(h_/(mass*f1));
long double Q20=sqrt(h_/(mass*f2));
long double q10=sqrt(h_/(mass*f3));
long double q20=sqrt(h_/(mass*f4));

// anharmonic contribution
long double alpha=0.05;
long double beta=0.05;

//coupling constant
long double f112=16500*c*2*PI;
long double f222=-1*13500*c*2*PI;
long double f121=40000*c*2*PI;

//force calculation
long double force(string mode,double Q1,double Q2,double q1,double q2){

        long double f=0;

        if(mode=="Q1"){
                f=-((h_*f1*Q1/(Q10*Q10))+(beta*2*h_*f1*Q1*Q1*Q1/(Q10*Q10*Q10*Q10))+(h_*f121*q1*q2/(q10*q20*Q10)));
//              cout<<h_*f1*Q1/(Q10*Q10)<<" "<<beta*2*h_*f1*Q1*Q1*Q1/(Q10*Q10*Q10*Q10)<<" "<<h_*f121*q1*q2/(q10*q20*Q10)<<endl;
        }
        if(mode=="Q2"){
                 f=-((h_*f2*Q2/(Q20*Q20))+(beta*2*h_*f2*Q2*Q2*Q2/(Q20*Q20*Q20*Q20))+(h_*f112*q1*q1/(q10*q10*Q20))+(h_*f222*q2*q2/(q20*q20*Q20)));
//               cout<<h_*f2*Q2/(Q20*Q20)<<" "<<beta*2*h_*f2*Q2*Q2*Q2/(Q20*Q20*Q20*Q20)<<" "<<h_*f112*q1*q1/(q10*q10*Q20)<<" "<<h_*f222*q2*q2/(q20*q20*Q20)<<endl;
        }
        if(mode=="q1"){
                 f=-((h_*f3*q1/(q10*q10))+(alpha*2*h_*f3*q1*q1*q1/(q10*q10*q10*q10))+(2*h_*f112*q1*Q2/(q10*q10*Q20))+(h_*f121*q2*Q1/(q20*Q10*q10)));
//               cout<<h_*f3*q1/(q10*q10)<<" "<<alpha*2*h_*f3*q1*q1*q1/(q10*q10*q10*q10)<<" "<<2*h_*f112*q1*Q2/(q10*q10*Q20)<<" "<<h_*f121*q2*Q1/(q20*Q10*q10)<<endl;
        }
        if(mode=="q2"){
                f=-((h_*f4*q2/(q20*q20))+(alpha*2*h_*f4*q2*q2*q2/(q20*q20*q20*q20))+(2*h_*f222*q2*Q2/(q20*q20*Q20))+(h_*f121*q1*Q1/(q10*q20*Q10)));
//              cout<<h_*f4*q2/(q20*q20)<<" "<<alpha*2*h_*f4*q2*q2*q2/(q20*q20*q20*q20)<<" "<<2*h_*f222*q2*Q2/(q20*q20*Q20)<<" "<<(h_*f121*q1*Q1/(q10*q20*Q10))<<endl;
        }

        return f;
}

void velocity_verlet(long double q1,long double p1,long double q2,long double p2,long double Q1,long double P1,long double Q2,long double P2){

        long double time=pow(10,-16);
        long double acce;
        long double new_acce;
        int i=0;

        int prev_P01=0;
        int prev_P10=1;
        int prev_P00=0;
//      cout<<q1<<" "<<p1<<" "<<q2<<" "<<p2<<" "<<Q1<<" "<<P1<<" "<<Q2<<" "<<P2<<endl;
//      long double Ez=(0.5)*mass*f2*f2*Q2*Q2+((0.5)*P2*P2)/mass;

        while(i<timesteps){

                long double in_q1=q1,in_q2=q2,in_p1=p1,in_p2=p2,in_Q1=Q1,in_Q2=Q2,in_P1=P1,in_P2=P2;

                long double acce1,acce2,acce3,acce4,new_acce1,new_acce2,new_acce3,new_acce4;

//              long double E3=((0.5)*mass*f3*f3*q1*q1)+(((0.5)*p1*p1)/mass);
//              long double E4=((0.5)*mass*f4*f4*q2*q2)+(((0.5)*p2*p2)/mass);
//              long double E1=((0.5)*mass*f1*f1*Q1*Q1)+(((0.5)*P1*P1)/mass);
//              long double E2=((0.5)*mass*f2*f2*Q2*Q2)+(((0.5)*P2*P2)/mass);

                long double E3=((0.5)*mass*f3*f3*q1*q1)+(((0.5)*p1*p1)/mass)+(0.5*h_*alpha*f3*q1*q1*q1*q1/(q10*q10*q10*q10));
                long double E4=((0.5)*mass*f4*f4*q2*q2)+(((0.5)*p2*p2)/mass)+(0.5*h_*alpha*f4*q2*q2*q2*q2/(q20*q20*q20*q20));
                long double E1=((0.5)*mass*f1*f1*Q1*Q1)+(((0.5)*P1*P1)/mass)+(0.5*h_*beta*f1*Q1*Q1*Q1*Q1/(Q10*Q10*Q10*Q10));
                long double E2=((0.5)*mass*f2*f2*Q2*Q2)+(((0.5)*P2*P2)/mass)+(0.5*h_*beta*f2*Q2*Q2*Q2*Q2/(Q20*Q20*Q20*Q20));

                long double Ec=h_*f112*q1*q1*Q2/(q10*q10*Q20)+h_*f222*q2*q2*Q2/(q20*q20*Q20)+h_*f121*q1*q2*Q1/(q10*q20*Q10);
                long double Eh=0.5*h_*beta*f1*Q1*Q1*Q1*Q1/(Q10*Q10*Q10*Q10)+0.5*h_*beta*f2*Q2*Q2*Q2*Q2/(Q20*Q20*Q20*Q20)+0.5*h_*alpha*f3*q1*q1*q1*q1/(q10*q10*q10*q10)+0.5*h_*alpha*f4*q2*q2*q2*q2/(q20*q20*q20*q20);

                acce1=force("Q1",Q1,Q2,q1,q2)/mass;
                acce2=force("Q2",Q1,Q2,q1,q2)/mass;
                acce3=force("q1",Q1,Q2,q1,q2)/mass;
                acce4=force("q2",Q1,Q2,q1,q2)/mass;

                Q1=Q1+(P1/mass)*time+(0.5)*acce1*time*time;
                Q2=Q2+(P2/mass)*time+(0.5)*acce2*time*time;
                q1=q1+(p1/mass)*time+(0.5)*acce3*time*time;
                q2=q2+(p2/mass)*time+(0.5)*acce4*time*time;

                new_acce1=force("Q1",Q1,Q2,q1,q2)/mass;
                new_acce2=force("Q2",Q1,Q2,q1,q2)/mass;
                new_acce3=force("q1",Q1,Q2,q1,q2)/mass;
                new_acce4=force("q2",Q1,Q2,q1,q2)/mass;

                P1=P1+(0.5)*mass*(acce1+new_acce1)*time;
                P2=P2+(0.5)*mass*(acce2+new_acce2)*time;
                p1=p1+(0.5)*mass*(acce3+new_acce3)*time;
                p2=p2+(0.5)*mass*(acce4+new_acce4)*time;

                int P01;
                int P10;
                int P00;
                int P;

//              tot_en[i]=(E1+E2+E3+E4);

//              cout<<((E1+E2+E3+E4)/(h_*f4))<<endl;
//              cout<<(E1+E2+E3+E4)/(0.5*h_*(f1+f2+f3+f4))<<endl;

//              tempout<<i<<" "<<(E1+E2+E3+E4+Ec)<<endl;
//              cout<<E1<<"    "<<h_*f3<<endl;
//              cout<<E2<<"    "<<h_*f4<<endl;
//              cout<<E3/(h_*f3)<<" "<<E4/(h_*f4)<<endl;

                if((E3>0 && E3<h_*f3) && (E4>0 && E4<h_*f4))
                        P00=1;
                else
                        P00=0;

                P00_avg[i]=P00_avg[i]+P00;


                if((E3>0 && E3<h_*f3) && (E4>h_*f4 && E4<2*h_*f4))
                        P01=1;
                else
                        P01=0;

                P01_avg[i]=P01_avg[i]+P01;

                if((E3>h_*f3 && E3<2*h_*f3) && (E4>0 && E4<h_*f4))
                        P10=1;
                else
                        P10=0;

                P10_avg[i]=P10_avg[i]+P10;
                i=i+1;
        }

//      cout<<q1<<" "<<p1<<" "<<q2<<" "<<p2<<" "<<Q1<<" "<<P1<<" "<<Q2<<" "<<P2<<endl;

}

//std::srand(static_cast<unsigned>(std::time(NULL)));

long double* init_cond(long double freq,long double n){

//      std::srand(static_cast<unsigned>(std::time(NULL)));
        double r2 = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/6.28));

        double angle=r2;
        long double energy=(n+0.5)*h_*freq;
        long double k=mass*freq*freq;

        long double q_=sqrt(2*energy)*cos(angle);
        long double p_=sqrt(2*energy)*sin(angle);

        long double p=p_*sqrt(mass);
        long double q=q_/sqrt(k);

        long double calc_energy=(0.5)*mass*freq*freq*q*q+((0.5)*p*p)/mass;


//      cout<<sqrt(energy)<<" "<<energy<<endl;
//      cout<<angle<<endl;

//      cout<<energy/calc_energy<<endl;

        long double* qp = new long double[2];
        qp[0]=q;
        qp[1]=p;

        return qp;
}
int main(){

        long double qp[10][2];

//      cout<<Q10<<endl;
//      cout<<Q20<<endl;
//      cout<<q10<<endl;
//      cout<<q20<<endl;

        std::srand(static_cast<unsigned>(std::time(NULL)));

        for(int i=0;i<n_traj;i++){

                // first mode

                long double* result=init_cond(f3,1);
                qp[0][0] = result[0];
                qp[0][1] = result[1];

                //second mode

                result=init_cond(f4,0);
                qp[1][0] = result[0];
                qp[1][1] = result[1];

                //third mode

                result=init_cond(f1,0);
                qp[2][0] = result[0];
                qp[2][1] = result[1];

                //fourth mode

                result=init_cond(f2,13);
                qp[3][0] = result[0];
                qp[3][1] = result[1];


                velocity_verlet(qp[0][0],qp[0][1],qp[1][0],qp[1][1],qp[2][0],qp[2][1],qp[3][0],qp[3][1]);

        }
                          for(int i=0;i<timesteps;i++){

                outFile<<i<<" ";
                outFile<<P01_avg[i]<<" ";
                outFile<<P10_avg[i]<<" ";
                outFile<<P00_avg[i]<<" ";
                outFile<<Q2_avg[i]<<" ";
                outFile<<endl;

                poptrans<<i<<" ";
                poptrans<<pop_trans_00_01[i]<<" ";
                poptrans<<pop_trans_00_10[i]<<" ";
                poptrans<<pop_trans_01_10[i]<<" ";
                poptrans<<endl;

        }


        outFile.close();

        poptrans_10_01.close();
        poptrans_00_01.close();
        poptrans_10_00.close();
        poptrans_01_00.close();
        poptrans_00_10.close();
        poptrans_01_10.close();

        return 0;

}
