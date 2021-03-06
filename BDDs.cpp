#include <bits/stdc++.h>
using namespace std;

const int Lim = 100;
/*
return 0 means: input data error;
return 1 means: unsatisfied;
return 2 means: tansform secceed;
the value of varibale[n - 1] 0 means less then; is 1 means greater then;
the value of coefficient[n - 1] means RHS.
*/

class MWVC{
    public:
    int gcd_vector(vector <int> &coefficient){
        int k = *lower_bound(coefficient.begin(), coefficient.end(), Lim);
        if(k < 0) return 0;

        int n = coefficient.size();
        while(k--)
            for(int i = 0; i < n && (coefficient[i] % k != 0); i++)
                if(i == n-1) return k;
        return 0;
    }

    int to_PBC(vector <int> &coefficient, vector<bool> &variable){ 
        std::ofstream outFile;
        outFile.open("/Users/haoji/Downloads/test1.txt");
        
        //Definsive data error.
        if(coefficient.size() != variable.size()) return 0;
        
        //"less then" are changed into "greater then" by negating all constants.
        int n = coefficient.size();

        if(!variable[n - 1]){
            for(int i = 0; i < n; i++)
                coefficient[i] = 0 - coefficient[i]; 
            variable[n - 1] = variable[n - 1] ^ 1;
            return to_PBC(coefficient, variable);
        }

        //Negative coefficients are eliminated by changing p into not p and updating the RHS.
        int sum = 0;
        for(int i = 0; i < n - 1; i++)
            if(coefficient[i] < 0){
                variable[i] = variable[i] ^ 1;
                coefficient[i] = abs(coefficient[i]);
                sum += coefficient[i];
            }
        coefficient[n - 1] += sum;

        //The coefficient are sorted in ascending order.
        vector<pair<int, int> > temp;
        
        for(int i = 0; i < n; i++)
            temp.push_back(make_pair(coefficient[i], variable[i]));
        sort(temp.begin(), temp.end() - 1);

        for(int i = 0; i < n; i++)
            outFile << temp[i].first << ";" << temp[i].second << endl;
        
        //Defensive the sum of LHS's coefficient is greater then RHS.
        int sum_coe = 0;
        for(int i = 0; i < n - 1; i++)
            sum_coe += coefficient[i];
        if(sum_coe < coefficient[n - 1]) return 1;

        //Coefficients greater than the RHS are trimmed to (replaced with) RHS.
        for(int i = 0; i < n - 1; i++)
            if(coefficient[i] > coefficient[n - 1])
                coefficient[i] = coefficient[n - 1];

        //The coefficients of the LHS are divided by their greatest common divisor(gcd).
        int gcd = gcd_vector(coefficient);
        if(!gcd)
            for(int i = 0; i < n; i++)
                coefficient[i] %= gcd;

        return 2;
    }
};


int main(){
    std::ofstream outFile;
    srand (time(NULL));
    MWVC ojb;
    outFile.open("/Users/haoji/Downloads/test.txt");

    vector<int> a(Lim);
    vector<bool> b(Lim);
    for(int i = 0; i < Lim; i++)
        a[i] = i - Lim/2;

    for(int i = 0; i < Lim; i++){
        int idx = rand() % (Lim - i) + i;
        int temp = a[idx];
        a[idx] = a[i];
        a[i] = temp;
        b[i] = rand() % 2;
    }

    for(int i = 0; i < Lim; i++)
        outFile << a[i] <<  endl;
    outFile << b[Lim - 1] << endl; 
    int res = ojb.to_PBC(a, b);
    outFile << b[Lim - 1] << endl; 

    return res;
}
 