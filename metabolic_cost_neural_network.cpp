/* file metabolic_cost_neural_network.cpp         Daniel Aguilar Velazquez danielvelaguil@gmail.com
    This code belongs to the paper "Low metabolic cost of rich-club neural networks at criticality" 
    The manuscript print the spikes and synaptic cost of avalanche activity in a hierarchical neural network
*/

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;


double base(int neighbors[][125], int indegree[], int counter)  /// function for all the 1 step R.B. model
 {
 int j;
 int i;
 for (i=0;i<4;i++)
  for (j=i+1;j<5;j++)
  {
   if(i!=j)
    {
    if(((double)rand()/((double)RAND_MAX+1.0))<0.5)
     {
     neighbors[int(counter+i)][int(indegree[counter+i])]=counter+j;
     indegree[counter+i]+=1;
     }
    if(((double)rand()/((double)RAND_MAX+1.0))<0.5)
     {
     neighbors[int(counter+j)][int(indegree[counter+j])]=counter+i;
     indegree[counter+j]+=1;
     }
    }
  }
 counter+=5;
 return counter;
 }


double hub(int neighbors[][125], int indegree[], int counter, int evol) // function for the 1,..,n step R.B. model
 {
 int j;
 int i;
 for (j=counter-pow(5,evol);j<counter-pow(5,evol-1);j+=pow(5,evol-1))
  for(i=j;i<j+pow(5,evol-1)-pow(5,evol-2);i++)   
  {
   if(((double)rand()/((double)RAND_MAX+1.0))<0.5)
    {
     neighbors[int(counter-1)][int(indegree[counter-1])]=i;
     indegree[counter-1]+=1;
    }
   if(((double)rand()/((double)RAND_MAX+1.0))<0.5)
    {
     neighbors[int(i)][int(indegree[i])]=counter-1;
     indegree[i]+=1;
    }
  }
 return counter;

 }



int main(int argc, char ** argv){
  int step=3;
  int evol=2,e;
  int clu=atof(argv[1]); // number of replicas of the R.B. network
  float kappa=atof(argv[2]); // rich-club connectivity
  float eta=atof(argv[3]) ; // percentage of excitatory hubs
  int nsteps=atoi(argv[4]); // number of time steps
  long seed=atoi(argv[5]); // seed
  float weight=atof(argv[6]); // synaptic weight (mV)
  int tau=int(10*atof(argv[7])); // synapic time course (ms)
  float global_inhibition=0.15,I; //percentage of low degree node inhibiton
  float h=0.1; // size of the time step (ms)
  int N=int(clu*pow(5,step)); //total size of the network
  int sum,hubex,sum_ex,counter=0;
  int indegree[1000];
  int neighbors[1000][125]; // maximum number of indegree connections=150
  int i,j,k,n,m,l,z;
  float A1,B1,A2,B2; // variables for numeric integration
  int gene=4;
  float w[N];
  float u[N],v[N],a[N],b[N],c[N],d[N],s; // variables of the Izhikevich model
  float im[N]; // state of neurons for time evolution
  int syn[N][tau]; //indicates if synaptic time course is actived 
  int timer[N]; // remember the number of time steps of synaptic time course
  int hubs_list[int(clu*5)];
  int nhubs=0;  
  int flag;
  float SC; // synaptic metabolic cost
  srand(seed);
  int trans=8000; // transitory time steps


/*/////////////////////////////////////////////////////
///////////// Structure of neural network ///////////
///////////////////////////////////////////////////*/
  counter=0;
  // creating a hierarchical network //////////////
  for (i=0; i<int(N/5);i+=5)
   for (j=0;j<5;j++)
    counter=base(neighbors,indegree,counter);

  while(evol<=gene-1)
  {
   for(i=pow(5,evol)-1;i<N;i+=pow(5,evol))
    counter=hub(neighbors,indegree,i+1,evol);
   evol+=1;
  }
 //////// rich club connectivity /////////////////////////
for(i=24;i<N;i+=25)
 for (j=24;j<N;j+=25)
  if(((double)rand()/((double)RAND_MAX+1.0))<kappa && i!=j && abs(i-j)<pow(5,3)) 
   {
   neighbors[int(i)][int(indegree[i])]=j;
   indegree[i]+=1;
   }


//////////////////////////////////////////////////////
///////////// Inhibition/excitation behavior ///////////
////////////////////////////////////////////////////



//////// low degree excitatory neurons /////////////////////////
for(i=0;i<N;i++)
 w[i]=weight;
 
 //////// low degree inhibition neurons /////////////////////////
for(i=0;i<N;i+=25)
 for (j=i;j<i+24;j++)
  if(((double)rand()/((double)RAND_MAX+1.0))<global_inhibition)
    w[j]=-weight;


 //////// hub inhibition neurons /////////////////////////

for(i=24;i<N;i+=25)
 w[i]=-weight;

// global hub = 5, local hub =1, each replica exhibits 1 global hub and 4 local hubs 
sum_ex=int(eta*(clu*9));
j=0;
nhubs=0;
while(j<sum_ex)
 {
  flag=0;
  hubex=25*int(clu*5*((double)rand()/((double)RAND_MAX+1.0)))+24;
  if(nhubs==0) {hubs_list[nhubs]=hubex; nhubs++;}
  else
  {
  for(i=0;i<int(clu*5);i++)
   if(hubex==hubs_list[i]) flag=1;
  }
  if(flag==0)
  {
   hubs_list[nhubs]=hubex;
   nhubs++; 
   if(indegree[hubex]>50)
    {
    if(j+5<sum_ex)
     {
     w[hubex]=weight;
     j+=5;
     }
    }
   else  
   {
    w[hubex]=weight;
    j+=1;
   }
  }
 }


/* /////////////////////////////////////
 /////////////  izhikevich model //////////
   ///////////////////////////////////   */


/////////////////// initialization ///////////////
for(i=0;i<N;i++)
 if(w[i]<0.0)
  {
  a[i]=0.02+0.08*((double)rand()/((double)RAND_MAX+1.0));
  b[i]=0.25-0.05*((double)rand()/((double)RAND_MAX+1.0));
  c[i]=-65;
  d[i]=2;
  }
 else
  {
  a[i]=0.02;
  b[i]=0.2;
  c[i]=-65+15*pow(((double)rand()/((double)RAND_MAX+1.0)),2);
  d[i]=8-6*pow(((double)rand()/((double)RAND_MAX+1.0)),2);
  }

for(i=0;i<N;i++)
 {
 timer[i]=tau-1;
 v[i]=c[i];
 u[i]=2;
 for(j=0;j<tau;j++)
  syn[i][j]=0;
 }

l=0;
e=0;
////////////////////// evolution //////////////////////////////
for(j=0;j<nsteps+trans;j++)
 {
  SC=0.0;
  for(i=0;i<N;i++) 
  {
  s=0;
  for(k=0;k<indegree[i];k++) // sum of incoming potentitals
   if(syn[neighbors[i][k]][l]==1)
    {
    s+=w[neighbors[i][k]];
    SC+=abs(w[neighbors[i][k]]);
    }
  if(w[i]>0.0)
   I=5*((double)rand()/((double)RAND_MAX+1.0));
  else
   I=2*((double)rand()/((double)RAND_MAX+1.0));
  A1=h*(0.04*pow(v[i],2)+5*v[i]+140-u[i]+I+s);  //second order Runge-kutta integration
  B1=h*(a[i]*(b[i]*v[i]-u[i]));
  A2=h*(0.04*pow((v[i]+A1/2.0),2)+5*(v[i]+A1/2.0)+140-(u[i]+B1/2.0)+I+s);
  B2=h*(a[i]*(b[i]*(v[i]+A1/2.0)-(u[i]+B1/2.0)));
  v[i]+=A2;
  u[i]+=B2;
  if(v[i]>=30.0)
   {
    for(z=0;z<tau;z++)
     if(z!=l)
      syn[i][z]=1;
   timer[i]=tau-1;
   v[i]=c[i];
   u[i]=u[i]+d[i];
   im[i]=30;
   }
  else if(v[i]<c[i])
   {
   if(l==tau-1 && timer[i]<=0)
    syn[i][0]=0;
   else if(l!=tau-1 && timer[i]<=0)
    syn[i][l+1]=0;
   v[i]=c[i];
   im[i]=c[i];
   timer[i]+=-1;
   }
  else
   {
   if(l==tau-1 && timer[i]<=0)
    syn[i][0]=0;
   else if(l!=tau-1 && timer[i]<=0)
    syn[i][l+1]=0;
   im[i]=v[i];
   timer[i]+=-1;
   }
  }
  if(j>=trans)
   {   
   sum=0;                 // counts the number of neurons that send synaptic potentials
   for(k=0;k<N;k++)
    {
    sum+=syn[k][l];
    if(im[k]==30)
     cout<<k<<" "<<0.0<<endl;  // print the neurons that fire (action potentials)
    }                     // column 1: number of spikes, column 2: synaptic cost
   if(sum==0)
    cout<<-1<<" "<<0.0<<endl; // if no neurons fire and send synaptic potentials it is printed "-1" to indicate the end of the avalanche
   else
    cout<<0.0<<" "<<SC<<endl;  // if synaptic potentials are sent it is printed the synaptic cost  
   }                         // column 1: number of spikes, column 2: synaptic cost
  if(l==tau-1)
   l=0;
  else
   l+=1;
  }
}
