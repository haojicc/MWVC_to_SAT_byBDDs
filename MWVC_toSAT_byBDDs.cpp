#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <fstream>
#include <numeric>

const int Lim = 100;
/*
Input are Graph G(V,E) and Vertex‘s weight weight_v.
At first, We transfrom this problem to PBO(Pseudo-Boolean Optimization) problem. But OPtimizaiton problem is also hard problem.
So, the next step is find a minimun value from [0, sum of w_i], to satisfied PBC(Pseduo-Boolean Constrains)-function.
In this step, we can use the binary search method to find the optimal value within time complexity of Log(n).
return 0 means: input data error;
return 1 means: unsatisfied;
return 2 means: tansform secceed/satisfied;
*/

class MWVC{
    public:
    /*
    int Graph_to_PBO(vector<int> &weight_v, vector<vector<int>> &edge, vector <int> &coefficient, vector<bool> &variable){}
    */

    static int gcd_vector(std::vector <int> &coefficient){
        int k = *lower_bound(coefficient.begin(), coefficient.end(), Lim);
        if(k < 0) return 0;

        int n;
        n = coefficient.size();
        while(k > 0) {
            int i;
            for(i = 0; i < n; i++){
                if(coefficient[i] % k != 0) break;
            }
            if(i == n) return k;
            k--;
        }
        return 0;
    }
    /*
    coefficient<int>:the weight of the variables of PBO(Pseudo-Boolean Optimization). The last value is RHS.
    variable<bool>:the variables of PBO(Pseudo-Boolean Optimization). The last bool value(0 or 1) means less then or greater then.
    PBP<int>: the coefficients of PB(pseudo-boolean)-part.
    PBP_V<bool>: the variables of PB(pseudo-boolean)-part.
    CP<int>: the numbers of the variables of Clauses-part.
    CP_V<bool>: the variables of Clauses-part. The last bool value(0 or 1) means less then or greater then.
    */
    int to_PBC(std::vector <int> &coefficient, std::vector<bool> &variable, std::vector<int> &PBP,
               std::vector<bool> &PBP_V, int* CP, std::vector<bool> &CP_V){
        std::ofstream outFile;
        outFile.open("/Users/haoji/Downloads/test1.txt");

        //Definsive data error.
        if(coefficient.size() != variable.size()) return 0;

        //"less then" are changed into "greater then" by negating all constants.
        int n = coefficient.size();

        if(!variable[n - 1]){
            for(int i = 0; i < n; i++)
                coefficient[i] = 0 - coefficient[i];
            variable[n - 1] = !variable[n - 1];
            return to_PBC(coefficient, variable, PBP, PBP_V, CP, CP_V);
        }

        //Negative coefficients are eliminated by changing p into not p and updating the RHS.
        int sum = 0;
        for(int i = 0; i < n - 1; i++)
            if(coefficient[i] < 0){
                variable[i] = !variable[i];
                coefficient[i] = abs(coefficient[i]);
                sum += coefficient[i];
            }
        coefficient[n - 1] += sum;

        //The coefficient are sorted in ascending order.
        std::vector<std::pair<int, int> > temp;
        for(int i = 0; i < n; i++)
            temp.emplace_back(coefficient[i], variable[i]);
        sort(temp.begin(), temp.end() - 1);

        for(int i = 0; i < n; i++){
//        outFile << temp[i].first << ";" << temp[i].second << std::endl;
            coefficient[i] = temp[i].first;
            variable[i] = temp[i].second;
        }

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
        if(gcd > 0)
            for(int i = 0; i < n; i++)
                coefficient[i] %= gcd;

        //Decompose to Clauses part and PB part.
        int j = 0;
        int k = 0;
        for(int i = 0; i < n - 1; i++){
            if(coefficient[i] >= coefficient[n - 1]){
                CP_V[j++] = variable[i];
            }
            else{
                PBP[k] = coefficient[i];
                PBP_V[k] = variable[i];
                k++;
            }
        }
        CP_V[j++] = false; // Set the linked variable with the Clauses part and PB part.
        CP_V[j++] = false; // Set Inequality Symbol.
        k++;
        PBP[k] = coefficient[n - 1]; //Set the linked variable with the Clauses part and PB part. The coefficient is equal to the RHS, because the constraint is satisfied as long as one of the variables in Clause part is true.
        PBP_V[k] = true; //Set the linked variable with the Clauses part and PB part.
        k++;
        PBP[k] = coefficient[n - 1]; //Set the RHS of PB part.
        PBP_V[k] = true; //Set the inequality symbol for the constraint.
        *CP = j + 1;
        for (int i = 0; i < n; i++)outFile << PBP[i] << ";" << PBP_V[i] << CP_V[i] <<std::endl;

        //complete the job.
        return 2;
    }

    bool bulidBDDs(std::vector<int> &PBP, std::vector<bool> &PBP_V, int i, int sum, int leftover, std::map<std::pair<int, int>, bool> &dp){
        int n = PBP.size();

        if(sum >= PBP[n - 1]) return true;
        else if(sum + leftover < PBP[n - 1]) return false;

        if(dp[std::make_pair(i,sum)]){
            i--;
            leftover -= PBP[i];
            int hi_sum = PBP_V[i] ? sum : sum + PBP[i];
            int lo_sum = PBP_V[i] ? sum + PBP[i] : sum;
            auto hi_result = bulidBDDs(PBP, PBP_V, i, hi_sum, leftover, dp);
            auto lo_result = bulidBDDs(PBP, PBP_V, i, lo_sum, leftover, dp);
            bool res = hi_result | lo_result;

            dp[std::make_pair(i,sum)] = res;
        }
        return dp[std::make_pair(i,sum)];
    }
    /*
    I need a consturct whose name is BDDNode has 4 satellite_dataes{int leftsum, int rightsum}.
    For save the dp_table, I need a map<pair<int, int>, BDDNode*>.
    When the inequality don't hold, we save the nullptr in the map_table.
    When the inequality hold, we save a particular node.(the satellite data have to some value?)
    */
};

int main(){
    std::ofstream outFile;
    srand (time(nullptr));
    MWVC ojb;
    outFile.open("/Users/haoji/Downloads/test.txt");

    std::vector<int> a(Lim);
    std::vector<bool> b(Lim);
    std::vector<int> PBP(Lim + 1);
    std::vector<bool> PBP_V(Lim + 1);
    int CP;
    std::vector<bool> CP_V(Lim + 1);

    for(int i = 0; i < Lim; i++)
        a[i] = i - Lim/2;
    for(int i = 0; i < Lim; i++){
        int idx = rand() % (Lim - i) + i;
        int temp = a[idx];
        a[idx] = a[i];
        a[i] = temp;
        b[i] = rand() % 2;
    }
    for(int i = 0; i <= Lim; i++)
        outFile << a[i] << ";"  << b[i] << ";" << std::endl;
    int res = ojb.to_PBC(a, b, PBP, PBP_V, &CP, CP_V);
    std::map<std::pair<int, int>, bool> dp;
    int j = PBP.size() - 1;
    outFile << "-----" << j << std::endl;
    int sum = 0;
    int leftover = std::accumulate(std::begin(PBP), std::end(PBP), 0);
    leftover -= PBP[j];
    bool res1;
    res1 = ojb.bulidBDDs(PBP, PBP_V, j, sum, leftover, dp);
    for(int i = 0; i <= Lim; i++)
        outFile << PBP[i] << ";" << PBP_V[i] << ";" << CP_V[i] << ";" << std::endl;
    for(auto &element: dp){
        auto fi = element.first;
        auto se = element.second;
        outFile << fi.first << ";" << fi.second << "---" << se << std::endl;
    }

    return 0;
}
