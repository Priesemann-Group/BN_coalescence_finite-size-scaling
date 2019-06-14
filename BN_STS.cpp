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
  double   m = 0;    // synaptic strength
  unsigned seed=1000;
  double   T = 1e0;  // time steps 
  unsigned avalanches = 0;

  std::stringstream help;
  help << "usage:\n";
  help << "     -N : number of neurons     (N=1e4         )\n";
  help << "     -m : synaptic strength\n";
  help << "     -s : seed\n";
  help << "     -T : number of time steps  (T=1e7 ms goal )\n";
  help << "     -A : number of externally driven avalanches (T=1e6 ms goal )\n";
  help << "     -o : output\n";
  int valid_args=0;
  for(unsigned i=0; i<argc; i++){
    if(std::string(argv[i])=="-N"){ N    = atof(argv[i+1]); valid_args+=1;}
    if(std::string(argv[i])=="-m"){ m    = atof(argv[i+1]); valid_args+=1;}
    if(std::string(argv[i])=="-s"){ seed = atof(argv[i+1]); valid_args+=1;}
    if(std::string(argv[i])=="-T"){ T    = atof(argv[i+1]); valid_args+=1;}
    if(std::string(argv[i])=="-A"){ avalanches = atof(argv[i+1]); valid_args+=1;}
    if(std::string(argv[i])=="-o"){ path =      argv[i+1] ; valid_args+=1;}
  }
  if(valid_args<5){std::cout << "not enough arguments\n" << help.str() << std::endl; exit(1);}
  
  bool flag_avalanches = false;
  if(avalanches>0) flag_avalanches=true;

  std::vector<int> neuron(N,0);
  std::vector<int> stimulus_int(N,0);

  //random numbers
  std::mt19937                           mt(seed);   
  std::uniform_real_distribution<double> uni_01(0,1);  // r=uni_01(mt)

  //initialize the binomial distributions from which the neurons-to-be-activated are drawn
  double p_std = m/static_cast<double>(N);
  std::binomial_distribution<int> binom(N,p_std);

  //loop over time
  std::cout << "simulation\n";
  std::stringstream filename;
  if(flag_avalanches){
    filename << path << "/BN_STS_binomial";
    filename << "_N" << std::setw(7) << std::setfill('0') << N;
    filename << "_m" << std::fixed << std::scientific << std::setprecision(2) << m;
    filename << "_A" << std::fixed << std::scientific << std::setprecision(2) << float(avalanches);
    filename << "_seed" << std::setw(4) << std::setfill('0') << seed;
    filename << "_avalanches.gz";
    std::cout << filename.str() << "\n";
  }
  else{
    filename << path << "/BN_STS_binomial";
    filename << "_N" << std::setw(7) << std::setfill('0') << N;
    filename << "_m" << std::fixed << std::scientific << std::setprecision(2) << m;
    filename << "_T" << std::fixed << std::scientific << std::setprecision(2) << float(T);
    filename << "_seed" << std::setw(4) << std::setfill('0') << seed;
    filename << "_time-series.gz";
    std::cout << filename.str() << "\n";
  }

  gzFile zfile = gzopen(filename.str().c_str(), "wb");
  std::stringstream line;
  if(flag_avalanches){ 
    line << "# size avalanche\n";
  }
  else{
    line << "#N_a = number of active sites at time step t\n";
    line << "#N_int = number of internally activated sites at time step t+1 (excluding external drive)\n";
    line << "# N_a\t N_int";
  }
  gzprintf(zfile, line.str().c_str());
  //observable
  int num_active=0, num_active_int=0;
  unsigned time=0;
  unsigned num_avalanches=0;
  double avalanche_time = 0, avalanche_size = 0;
  //loop over avalanche onsets
  while(true){
    // STS activation
    if(num_active == 0){
      unsigned n = uni_01(mt)*neuron.size(); 
      neuron[n]      = 1; 
      num_active     = 1;
      num_active_int = 0;
      avalanche_time = 0;
      avalanche_size = 0;
    }
    // internal signal processing (branching process)
    else{
      for(unsigned n=0;n<neuron.size(); n++){
        //annealed average: select random postsynaptic neurons for each active neuron 
        if(neuron[n] > 0 ){
          // finite-size correction of m directly in binomial distribution
          int k_n = binom(mt) ;
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
      }
      // integration over internal activation and external activation
      // num_active remains the number of initially active neurons that lead to this activation
      num_active = 0;
      num_active_int = 0;
      for(unsigned n=0; n<neuron.size(); n++){
        if(stimulus_int[n]>0 ){
           neuron[n]=1; 
           num_active +=1;
           if(stimulus_int[n]>0){
             num_active_int += 1;
           }
           stimulus_int[n]=0;
        }
        else if(neuron[n]==1){ neuron[n]=0; }
      }
    }
    time++;
    avalanche_time++;
    avalanche_size += num_active;

    //write to file
    if(flag_avalanches){
      if(avalanche_size>0 and num_active==0){
        std::stringstream data;
        data << avalanche_time-1 << " " << avalanche_size << std::endl;
        gzprintf(zfile, data.str().c_str());
        num_avalanches++;
      }
      if(num_avalanches+1 > avalanches) break;
    }
    else{
      std::stringstream data;
      data << std::scientific;
      data << num_active << " " << num_active_int << std::endl;
      gzprintf(zfile, data.str().c_str());
      if(time+1>T) break; 
    }
  }
  gzclose(zfile);
}
