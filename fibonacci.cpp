#include <iostream>
using namespace std;

int fibo (int n){

    if((n==1) || (n==0))
    return n;

    else
        return fibo(n-1) + fibo(n-2);
}



int main (){
    int n;

    cout << "Insira um valor de n p/ a sequencia de Fibonacci:  " << endl;
    cin >> n;
    cout << fibo(n) << endl;

    return 0;
}