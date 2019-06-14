// Branching network model with Poissonian external input of mean rate h
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <list>
#include <random>
#include <algorithm>
#include <zlib.h>

double lambda(double h){
  double dt=1;
  return 1.0-exp(-h*dt);
}

int main(int argc, char* argv[]){
  std::string path;  // output path
  unsigned N = 1e0;  // system size 
  double   T = 1e0;  // time steps 
  double   h = 1.00; // external drive
  double   m = 0;    // synaptic strength
  unsigned seed=1000;

  std::stringstream help;
  help << "usage:\n";
  help << "     -N : number of neurons     (N=1e4         )\n";
  help << "     -T : number of time steps  (T=1e7 ms goal )\n";
  help << "     -h : external drive/neuron (h<1)\n";
  help << "     -m : synaptic strength\n";
  help << "     -s : seed\n";
  help << "     -o : output\n";
  int valid_args=0;
  for(unsigned i=0; i<argc; i++){
    if(std::string(argv[i])=="-N"){ N    = atof(argv[i+1]); valid_args+=1;}
    if(std::string(argv[i])=="-T"){ T    = atof(argv[i+1]); valid_args+=1;}
    if(std::string(argv[i])=="-h"){ h    = atof(argv[i+1]); valid_args+=1;}
    if(std::string(argv[i])=="-m"){ m    = atof(argv[i+1]); valid_args+=1;}
    if(std::string(argv[i])=="-s"){ seed = atof(argv[i+1]); valid_args+=1;}
    if(std::string(argv[i])=="-o"){ path =      argv[i+1] ; valid_args+=1;}
  }
  if(valid_args<6){std::cout << "not enough arguments\n" << help.str() << std::endl; exit(1);}

  std::vector<int> neuron(N,0);
  std::vector<int> stimulus_int(N,0);
  std::vector<int> stimulus_ext(N,0);

  //random numbers
  std::mt19937                           mt(seed);   
  std::uniform_real_distribution<double> uni_01(0,1);  // r=uni_01(mt)

  //precompute m_fsc = m(A)
  std::vector<double> m_fsc(N+1,0.0);
  for(unsigned A=0; A<m_fsc.size(); A++){ 
    if(1-m*A/N > 0){
      m_fsc[A] = N*(1-std::pow(1-m*A/N, 1.0/A)); 
    }
    else{ // should only happen for m>=1
      m_fsc[A] = log(N);
    }
  }
  //precompute the corresponding probability of activation p_fsc=m_fsc/N
  std::vector<double> p_fsc(N+1,0.0);
  for(unsigned A=0; A<p_fsc.size(); A++){ p_fsc[A] = m_fsc[A]/static_cast<double>(N);}
  //initialize the binomial distributions from which the neurons-to-be-activated are drawn
  std::vector<std::binomial_distribution<int> > binom;
  for(unsigned A=0; A<p_fsc.size(); A++){ binom.push_back(std::binomial_distribution<int>(N,p_fsc[A])); }

  std::cout << "simulation\n";
  std::stringstream filename;
  filename << path << "/BNcc_driven_binomial";
  //filename << "_N" << std::fixed << std::scientific << std::setprecision(2) << float(N);
  filename << "_N" << std::setw(7) << std::setfill('0') << N;
  filename << "_m" << std::fixed << std::scientific << std::setprecision(2) << m;
  filename << "_h" << std::fixed << std::scientific << std::setprecision(2) << h;
  filename << "_T" << std::fixed << std::scientific << std::setprecision(2) << float(T);
  filename << "_seed" << std::setw(4) << std::setfill('0') << seed;
  filename << "_time-series.gz";
  std::cout << filename.str() << "\n";

  gzFile zfile = gzopen(filename.str().c_str(), "wb");
  std::stringstream line;
  line << "#N_a = number of active sites at time step t\n";
  line << "#N_int = number of internally activated sites at time step t+1 (excluding external drive)\n";
  line << "# N_a\t N_int";
  gzprintf(zfile, line.str().c_str());
  //observable
  int num_active=0, num_active_int=0;
  // postsynaptic acceptance
  double p_accept = m/N;
  unsigned time=0;
  //loop over time
  while(true){
    // internal signal processing (branching process)
    for(unsigned n=0;n<neuron.size(); n++){
      //annealed average: select random postsynaptic neurons for each active neuron 
      if(neuron[n] > 0 ){
        //int k_n = distribution(generator);
        // finite-size correction of m directly in binomial distribution
        int k_n = binom[num_active](mt) ;
        std::list<double> nconnect;
        int j=0;
        while(j<k_n){
          int nn = uni_01(mt)*neuron.size();
          if(std::find(nconnect.begin(), nconnect.end(), nn) == nconnect.end()){
            nconnect.push_back(nn);
            stimulus_int[nn] += 1;
            j+=1;
          }
        }
      }
      //random external activation with probability h
      if(uni_01(mt) < lambda(h)) { 
        stimulus_ext[n] = 1; 
      } 
    }
    // integration over internal activation and external activation
    // num_active remains the number of initially active neurons that lead to this activation
    num_active = 0;
    num_active_int = 0;
    for(unsigned n=0; n<neuron.size(); n++){
      if(stimulus_int[n]>0 or stimulus_ext[n]>0){
         neuron[n]=1; 
         num_active +=1;
         if(stimulus_int[n]>0 && stimulus_ext[n]==0){
           num_active_int += 1;
         }
         stimulus_int[n]=0;
         stimulus_ext[n]=0;
      }
      else if(neuron[n]==1){ neuron[n]=0; }
    }
    time++;
    std::stringstream data;
    data << std::scientific;
    data << num_active << " " << num_active_int <<  std::endl;
    gzprintf(zfile, data.str().c_str());
    if(time>T) break; 
  }
  gzclose(zfile);
}
