#include <iostream>
#include <cmath>

using namespace std;

int main(){
    const int n = 20;
    double I[n+1]{};
//    I[20] = (1./(6.*21.) + 1./(5.*21.)) / 2.;
    I[20] = 0.008;

    for(int i = n; i--;){
        I[i] = 1./(5.*(i+1.)) - I[i+1]/5.;
        cout << "I[" << i << "]: ";
        cout << I[i] << endl;
    }

}
