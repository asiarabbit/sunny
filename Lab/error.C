#include <iostream>
#include <cmath>

using namespace std;

int main(){
    const int n = 20;
    double I[n+1] = {0.182322};
    
    for(int i = 0; i < n; i++){
        I[i + 1] = 1./(i + 1.) - 5.*I[i];
        cout << "I[" << i + 1 << "]: ";
        cout << I[i+1] << endl;
    }

}
