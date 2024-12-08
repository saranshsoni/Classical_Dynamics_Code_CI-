#include<bits/stdc++.h>
using namespace std;


std::ofstream outFile("output.txt");
std::ofstream t0("t0.txt");
std::ofstream t12("t12.txt");
std::ofstream t20("t20.txt");
std::ofstream t28("t28.txt");
std::ofstream t36("t36.txt");
std::ofstream t100("t100.txt");
std::ofstream t900("t900.txt");
std::ofstream poptrans_10_01("poptrans_10_01.txt");
std::ofstream poptrans_00_01("poptrans_00_01.txt");
std::ofstream poptrans_10_00("poptrans_10_00.txt");
std::ofstream poptrans_01_00("poptrans_01_00.txt");
std::ofstream poptrans_00_10("poptrans_00_10.txt");
std::ofstream poptrans_01_10("poptrans_01_10.txt");
std::ofstream poptrans("pop_trans.txt");
std::ofstream tempout("temp_out.txt");

int n_traj=2000;
const int timesteps=100000;
int P01_avg[timesteps]={};
int P10_avg[timesteps]={};
int P00_avg[timesteps]={};
int pop_trans_00_01[timesteps]={};
int pop_trans_00_10[timesteps]={};
int pop_trans_10_00[timesteps]={};
int pop_trans_10_01[timesteps]={};
int pop_trans_01_00[timesteps]={};
int pop_trans_01_10[timesteps]={};
long double E_avg[timesteps]={};
long double KE_avg[timesteps]={};
long double PE_avg[timesteps]={};

long double mass=1.67*pow(10,-27);
long double c=3*pow(10,8);
long double h=6.6264*pow(10,-34);
long double h_=h/(2*3.14);

long double tot_en[timesteps]={};
//int count=0;
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

long double P10=sqrt(mass*h_*f1);
long double P20=sqrt(mass*h_*f2);
long double p10=sqrt(mass*h_*f3);
long double p20=sqrt(mass*h_*f4);
// atomic units

long double F0=8.23872*pow(10,-8);
long double x0=5.29*pow(10,-11);
long double P0=1.998*pow(10,-24);

// anharmonic contribution
long double alpha=0.05;
long double beta=0.05;

//coupling constant
long double f112=16500*c*2*PI;
long double f222=-1*13500*c*2*PI;
long double f121=40000*c*2*PI;

// numerical coefficents for friction

long double c0;
long double c1;
long double c2;
long double kb=1.38*pow(10,-23);
long double temp=298;
long double diff=0.187*pow(10,-9);
long double eta=2900000*c*2*PI;                 //((kb*temp)/(mass*diff));        //2900000*c;

// generate gaussian distribution with mean 0 and sigma 1

long double gaussian_rand(){

//      long double stddev=1;
//        long double mean=0;

        long double u1 = static_cast<double>(rand()) / RAND_MAX; // Uniform random number in (0,1)
        long double u2 = static_cast<double>(rand()) / RAND_MAX; // Uniform random number in (0,1)

        long double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2); // Standard normal variate

        return z0;

//      return z0 * stddev + mean; // Scale and shift to match desired mean and stddev

}

long double* stochastic_force(long double time) {

//      long double u1 = static_cast<double>(rand()) / RAND_MAX; // Uniform random number in (0,1)
//      long double u2 = static_cast<double>(rand()) / RAND_MAX; // Uniform random number in (0,1)
//      long double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2); // Standard normal variate

//    long double gdt=eta*time

        long double u1=gaussian_rand();
        long double u2=gaussian_rand();

        long double sig_r=time*sqrt((kb*temp/mass)*(1/(eta*time))*(2-(1/(eta*time))*(3-4*exp(-eta*time)+exp(-2*eta*time))));
        long double sig_v=sqrt((kb*temp/mass)*(1-exp(-2*eta*time)));
        long double sig_rv=(time*(kb*temp/mass)*(1/(eta*time))*pow((1-exp(-eta*time)),2))/(sig_r*sig_v);

        long double del_r=sig_r*u1;
        long double del_v=sig_v*(sig_rv*u1+sqrt(1-(sig_rv*sig_rv))*u2);

        long double* rv = new long double[2];
        rv[0]=del_r;
        rv[1]=del_v;

//      cout<<del_r<<"  "<<mass*del_v<<endl;

        return rv;

//      return z0 * stddev + mean; // Scale and shift to match desired mean and stddev
}


//force calculation
long double force(string mode,double Q1,double Q2,double q1,double q2,double P1,double P2){

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
//      int count=0;
//      cout<<q1<<" "<<p1<<" "<<q2<<" "<<p2<<" "<<Q1<<" "<<P1<<" "<<Q2<<" "<<P2<<endl;
//      long double Ez=(0.5)*mass*f2*f2*Q2*Q2+((0.5)*P2*P2)/mass;

        while(i<timesteps){

                long double in_q1=q1,in_q2=q2,in_p1=p1,in_p2=p2,in_Q1=Q1,in_Q2=Q2,in_P1=P1,in_P2=P2;

                long double acce1,acce2,acce3,acce4,new_acce1,new_acce2,new_acce3,new_acce4;

                long double E3=((0.5*h_*f3*q1*q1)/(q10*q10))+((0.5*h_*f3*p1*p1)/(p10*p10))+(0.5*h_*alpha*f3*q1*q1*q1*q1/(q10*q10*q10*q10));
                long double E4=((0.5*h_*f4*q2*q2)/(q20*q20))+((0.5*h_*f4*p2*p2)/(p20*p20))+(0.5*h_*alpha*f4*q2*q2*q2*q2/(q20*q20*q20*q20));
                long double E1=((0.5*h_*f1*Q1*Q1)/(Q10*Q10))+((0.5*h_*f1*P1*P1)/(P10*P10))+(0.5*h_*beta*f1*Q1*Q1*Q1*Q1/(Q10*Q10*Q10*Q10));
                long double E2=((0.5*h_*f2*Q2*Q2)/(Q20*Q20))+((0.5*h_*f2*P2*P2)/(P20*P20))+(0.5*h_*beta*f2*Q2*Q2*Q2*Q2/(Q20*Q20*Q20*Q20));

//              long double E3=((0.5)*mass*f3*f3*q1*q1)+(((0.5)*p1*p1)/mass)+(0.5*h_*alpha*f3*q1*q1*q1*q1/(q10*q10*q10*q10));
//                long double E4=((0.5)*mass*f4*f4*q2*q2)+(((0.5)*p2*p2)/mass)+(0.5*h_*alpha*f4*q2*q2*q2*q2/(q20*q20*q20*q20));
//                long double E1=((0.5)*mass*f1*f1*Q1*Q1)+(((0.5)*P1*P1)/mass)+(0.5*h_*beta*f1*Q1*Q1*Q1*Q1/(Q10*Q10*Q10*Q10));
//                long double E2=((0.5)*mass*f2*f2*Q2*Q2)+(((0.5)*P2*P2)/mass)+(0.5*h_*beta*f2*Q2*Q2*Q2*Q2/(Q20*Q20*Q20*Q20));
                long double Ec=h_*f112*q1*q1*Q2/(q10*q10*Q20)+h_*f222*q2*q2*Q2/(q20*q20*Q20)+h_*f121*q1*q2*Q1/(q10*q20*Q10);
                long double Eh=0.5*h_*beta*f1*Q1*Q1*Q1*Q1/(Q10*Q10*Q10*Q10)+0.5*h_*beta*f2*Q2*Q2*Q2*Q2/(Q20*Q20*Q20*Q20)+0.5*h_*alpha*f3*q1*q1*q1*q1/(q10*q10*q10*q10)+0.5*h_*alpha*f4*q2*q2*q2*q2/(q20*q20*q20*q20);

                c0=exp(-eta*time);
                c1=(1-c0)/(eta*time);
                c2=(1-c1)/(eta*time);
                long double* QP_gauss1=stochastic_force(time);
                long double* QP_gauss2=stochastic_force(time);


                acce1=force("Q1",Q1,Q2,q1,q2,P1,P2)/mass;
                acce2=force("Q2",Q1,Q2,q1,q2,P1,P2)/mass;
                acce3=force("q1",Q1,Q2,q1,q2,P1,P2)/mass;
                acce4=force("q2",Q1,Q2,q1,q2,P1,P2)/mass;

                Q1=(Q1)+(c1*(P1/mass)*time)+(c2*acce1*time*time)+(QP_gauss1[0]);
                Q2=(Q2)+(c1*(P2/mass)*time)+(c2*acce2*time*time)+(QP_gauss2[0]);
                q1=(q1)+((p1/mass)*time)+((0.5)*acce3*time*time);
                q2=(q2)+((p2/mass)*time)+((0.5)*acce4*time*time);

                new_acce1=force("Q1",Q1,Q2,q1,q2,P1,P2)/mass;
                new_acce2=force("Q2",Q1,Q2,q1,q2,P1,P2)/mass;
                new_acce3=force("q1",Q1,Q2,q1,q2,P1,P2)/mass;
                new_acce4=force("q2",Q1,Q2,q1,q2,P1,P2)/mass;

//              P1=c0*P1+(c1-c2)*(2*0.5)*mass*(acce1+0)*time+c2*time*new_acce1+mass*QP_gauss1[1];
//              P2=c0*P2+(c1-c2)*(2*0.5)*mass*(acce2+0)*time+c2*time*new_acce2+mass*QP_gauss2[1];

                P1=(c0*P1)+(mass*((c1-c2)*acce1+(c2)*new_acce1)*time)+(mass*QP_gauss1[1]);
                P2=(c0*P2)+(mass*((c1-c2)*acce2+(c2)*new_acce2)*time)+(mass*QP_gauss2[1]);
                p1=(p1)+((0.5)*mass*(acce3+new_acce3)*time);
                p2=(p2)+((0.5)*mass*(acce4+new_acce4)*time);


                int P01;
                int P10;
                int P00;
                int P;
//              tot_en[i]=(E1+E2+E3+E4);

//              cout<<((E1+E2+E3+E4)/(h_*f4))<<endl;
//              cout<<(E1+E2+E3+E4)/(0.5*h_*(f1+f2+f3+f4))<<endl;
//              outFile<<i<<" "<<Q1/Q10<<" "<<Q2/Q20<<" "<<endl;
//              outFile<<i<<" "<<(E1+E2+E3+E4+Ec)<<endl;
//              cout<<E1<<"    "<<h_*f3<<endl;
//              cout<<E2<<"    "<<h_*f4<<endl;
//              cout<<E3/(h_*f3)<<" "<<E4/(h_*f4)<<endl;
//
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

                E_avg[i]=E_avg[i]+(E1+E2+E3+E4+Ec);

                KE_avg[i]=KE_avg[i]+((0.5*h_*f1*P1*P1)/(P10*P10))+((0.5*h_*f2*P2*P2)/(P20*P20))+((0.5*h_*f3*p1*p1)/(p10*p10))+((0.5*h_*f4*p2*p2)/(p20*p20));
                PE_avg[i]=PE_avg[i]+((0.5*h_*f3*q1*q1)/(q10*q10))+((0.5*h_*f4*q2*q2)/(q20*q20))+((0.5*h_*f1*Q1*Q1)/(Q10*Q10))+((0.5*h_*f2*Q2*Q2)/(Q20*Q20))+Ec+Eh;

                prev_P01=P01;
                prev_P10=P10;
                prev_P00=P00;

                poptrans<<i<<" "<<Q1/Q10<<" "<<Q2/Q20<<" "<<endl;

//              cout<<acce<<endl;
//              cout<<i<<" "<<P01_avg[i]<<endl;
//
                tempout<<i<<" "<<E3<<" "<<E4<<endl;

//              cout<<i<<" "<<n1<<endl;
//              cout<<i<<" "<<n2<<endl;
//              cout<<n1<<" "<<n2<<endl;
//              cout<<E3/(h_*f3)<<endl;
//              cout<<E4/(h_*f4)<<endl;
//              cout<<E1/(h_*f3)<<endl;
//              cout<<E2/(h_*f3)<<endl;


//              cout<<Q1<<" "<<P1<<endl;
                i=i+1;
        }
        //cout<<count<<endl;

//      cout<<q1<<" "<<p1<<" "<<q2<<" "<<p2<<" "<<Q1<<" "<<P1<<" "<<Q2<<" "<<P2<<endl;

}
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
//      cout<<count<<endl;

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


//              cout<<qp[2][0]<<"  "<<qp[2][1]<<endl;

//              cout<<eta*qp[2][1]<<"  "<<(h_*f1*qp[2][0]/(Q10*Q10))<<" "<<generateGaussian(0,0.01*F0)<<endl;
//              cout<<P10<<" "<<P20<<endl;
                c0=exp(-eta*pow(10,-16));
                c1=(1-c0)/(eta*pow(10,-16));
                c2=(1-c1)/(eta*pow(10,-16));

//              cout<<c0<<"  "<<c1<<"  "<<c2<<endl;
//              cout<<eta*pow(10,-16)<<endl;

//              cout<<(mass*diff*290000*c)/kb<<endl;
                velocity_verlet(qp[0][0],qp[0][1],qp[1][0],qp[1][1],qp[2][0],qp[2][1],qp[3][0],qp[3][1]);


        }



        for(int i=0;i<timesteps;i++){

                outFile<<i<<" ";
                outFile<<E_avg[i]<<" ";
                outFile<<KE_avg[i]<<" ";
                outFile<<PE_avg[i]<<" ";
                outFile<<P01_avg[i]<<" ";
                outFile<<P10_avg[i]<<" ";
                outFile<<P00_avg[i]<<" ";

                outFile<<endl;


//      outFile<<"kjlngjkbkbf4jjkbfejkwn "<<endl;
        outFile.close();
        t0.close();
        t12.close();
        t20.close();
        t28.close();
        t36.close();
        t100.close();
        t900.close();
        poptrans_10_01.close();
        poptrans_00_01.close();
        poptrans_10_00.close();
        poptrans_01_00.close();
        poptrans_00_10.close();
        poptrans_01_10.close();
        return 0;

}
                
