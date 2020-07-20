// odometer.C -- generate all combinations using odometer algorithm

#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;

int main(){
  const int nParticle = 3, nSPState = 5;
  int SPStateVec[nParticle];
  for(int i = 0; i < nParticle; i++) SPStateVec[i] = i;
  for(int t : SPStateVec) cout << t << " ";
  cout << endl;

  while(1){
    int i = nParticle - 1;
    int maxMBSD = i + nSPState - nParticle; // max MBSD of the i-th bit
    // get the bit that is not full
    while(i >= 0 && SPStateVec[i--] == maxMBSD--){}
    SPStateVec[++i]++; // increment the unfull bit
    // recover the following bits
    while(i < nParticle - 1){
      SPStateVec[i + 1] = SPStateVec[i] + 1; i++;
    } // end while
    getchar();
    for(int t : SPStateVec) cout << t << " ";
    cout << endl;
    if(SPStateVec[0] == nSPState - nParticle) break;
  } // end while
} // end the main function
