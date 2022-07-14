#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <time.h>

/*
Input are Graph G(V,E) and Vertex‘s weight weight_v.
At first, We transfrom this problem to PBO(Pseudo-Boolean Optimization) problem. But OPtimizaiton problem is also hard problem.
So, the next step is find a minimun value from [0, sum of w_i], to satisfied PBC(Pseduo-Boolean Constrains)-function.
In this step, we can use the binary search method to find the optimal value within time complexity of Log(n).
return 0 means: input data error;
return 1 means: unsatisfied;
return 2 means: tansform secceed/satisfied;
*/

const int Lim = 10;
typedef struct{
    int leftsum;
    int rightsum;
    int output_vari;

}BDDNode;
BDDNode TrueNode = {1,1};
BDDNode FalseNode = {0,0};

std::map<std::pair<int, int>, BDDNode> dp;
std::vector<int> corfficient(Lim + 1, 0);
std::vector<bool> variable(Lim + 1, 0);
std::vector<int> clause(10 * Lim);
int index_vari = Lim;
int clause_posion = 1;

/*
 a=corfficient and b=variable;
 a/b[0] = RHS / inequality sign.
 a/b[1] = the numbers of Pseudo-boolean part. /#
 a/b[2] = the connector of Pseudo-boolean part and Clause part. /#
text data:
a = {24,0,0,-6,-16,2,-6,-6,-16,-2,6,-14,16,-4,16,-4}
b = {0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1}
 */

class MWVC{
    public:
    /*
    int Graph_to_PBO(vector<int> &weight_v, vector<vector<int>> &edge, vector <int> &coefficient, vector<bool> &variable){}
    */

    static int gcd_vector(std::vector <int> &coefficient){
        int n = coefficient.size();
        int k = coefficient[n-1];

        while(k > 0) {
            int i = 3;
            while(i < n && coefficient[i] % k == 0) i++;
            if(i == n) return k;
            k--;
        }
        return k;
    }
    /*
    coefficient<int>:the weight of the variables of PBO(Pseudo-Boolean Optimization). The last value is RHS.
    variable<bool>:the variables of PBO(Pseudo-Boolean Optimization). The last bool value(0 or 1) means less then or greater then.
    PBP<int>: the coefficients of PB(pseudo-boolean)-part.
    PBP_V<bool>: the variables of PB(pseudo-boolean)-part.
    CP<int>: the numbers of the variables of Clauses-part.
    CP_V<bool>: the variables of Clauses-part. The last bool value(0 or 1) means less then or greater then.
    */
    int PBOtoPBC(std::vector<int> &coefficient, std::vector<bool> &variable){
        //text val.
        std::ofstream outFile;
        outFile.open("/Users/haoji/Downloads/test1.txt");

        //Definsive data error.
        if(coefficient.size() != variable.size()) return 0;

        //"less then" are changed into "greater then" by negating all constants.
        int n = coefficient.size();
        if(!variable[0]){
            for(int i = 2; i < n; i++)
                coefficient[i] = 0 - coefficient[i];
            variable[0] = 1 - variable[0];
            coefficient[0] = 0 - coefficient[0];
            return PBOtoPBC(coefficient, variable);
        }
        outFile << "--调整不等式符号--" << std::endl;
        for (int i = 0; i < n; i++)outFile << corfficient[i] << ";" << variable[i] <<std::endl;
        outFile << "--调整不等式符号--" << std::endl;

        //Negative coefficients are eliminated by changing p into not p and updating the RHS.
        int sum = 0;
        for(int i = 3; i < n; i++)
            if(coefficient[i] < 0){
                variable[i] = 1 - variable[i];
                coefficient[i] = abs(coefficient[i]);
                sum += coefficient[i];
            }
        coefficient[0] += sum;

        outFile << "--将负值系数调整为正数--" << std::endl;
        for (int i = 0; i < n; i++)outFile << corfficient[i] << ";" << variable[i] <<std::endl;
        outFile << "--将负值系数调整为正数--" << std::endl;

        //The coefficient are sorted in descending order.
        std::vector<std::pair<int, int> > temp;
        std::vector<bool> tempv = variable;
        for(int i = 3; i < n; i++)
            temp.emplace_back(coefficient[i], i);

        sort(temp.begin(), temp.end(), std::greater<>());

        outFile << "--系数排序--" << std::endl;
        for(int i = 3; i < n; i++){
            outFile << temp[i].first << ";" << temp[i].second << std::endl;
            coefficient[i] = temp[i - 3].first;
            variable[i] = tempv[temp[i-3].second];
        }
        outFile << "--系数排序--" << std::endl;

        //Defensive the sum of LHS's coefficient is greater then RHS.
        int sum_coe = 0;
        for(int i = 3; i < n; i++)
            sum_coe += coefficient[i];
        if(sum_coe <= coefficient[0]) return 1;

        //Coefficients greater than the RHS are trimmed to (replaced with) RHS.
        for(int i = 3; i < n; i++)
            if(coefficient[i] > coefficient[0])
                coefficient[i] = coefficient[0];

        outFile << "--剪掉多余的系数--" << std::endl;
        for (int i = 0; i < n; i++)outFile << corfficient[i] << ";" << variable[i] <<std::endl;
        outFile << "--剪掉多余的系数--" << std::endl;

        //The coefficients of the LHS are divided by their greatest common divisor(gcd).
        int gcd = gcd_vector(coefficient);
        outFile << "-----------------" << std::endl;
        outFile << "gcd：" << gcd <<std::endl;
        outFile << "-----------------" << std::endl;
        if(gcd > 0){
            coefficient[0] /= gcd;
            for(int i = 3; i < n; i++)
                coefficient[i] /= gcd;
        }

        //Decompose to Clauses part and PB part.
        int x = 3;
        while(x < n && coefficient[x] == coefficient[0]) x++;
        coefficient[1] = x;
        //complete the job.
        outFile << "--系数处理结束后的结果--" << std::endl;
        for (int i = 0; i < n; i++)outFile << corfficient[i] << ";" << variable[i] <<std::endl;
        outFile << "--系数处理结束后的结果--" << std::endl;
        return 2;

    }

    BDDNode bulidBDDs(std::vector<int> &PBP, std::vector<bool> &PBP_V, int i, int sum, int leftover, std::map<std::pair<int, int>, BDDNode> &dp){
        if(sum >= PBP[0]) return dp[std::make_pair(i,sum)] = TrueNode;
        else if(sum + leftover < PBP[0]) return dp[std::make_pair(i,sum)] = FalseNode;

        if(dp.find(std::make_pair(i,sum)) == dp.end()){
            leftover -= PBP[i];
            int hi_sum = PBP_V[i] ? sum + PBP[i] : sum;
            int lo_sum = PBP_V[i] ? sum : sum + PBP[i];
            auto hi_result = bulidBDDs(PBP, PBP_V, i + 1, hi_sum, leftover, dp);
            auto lo_result = bulidBDDs(PBP, PBP_V, i + 1, lo_sum, leftover, dp);
            dp[std::make_pair(i,sum)] = {hi_sum, lo_sum, index_vari++};
        }
        return dp[std::make_pair(i,sum)];
    }

    void BDDtoSAT(std::vector<int> &clause, int vari, BDDNode &node, int leftnode, int lefttype, int rightnode, int righttype){
        if(lefttype == 1 && righttype == 0){
            clause[3 * clause_posion] = 2 * vari;
            clause[3 * clause_posion + 1] = 2 * node.output_vari + 1;
            clause_posion++;
            clause[3 * clause_posion] = 2 * vari + 1;
            clause[3 * clause_posion + 1] = 2 * node.output_vari;
            clause_posion++;
        }else if(lefttype == 1 && righttype != 1 && righttype != 0){
            clause[3 * clause_posion] = 2 * node.output_vari + 1;
            clause[3 * clause_posion + 1] = 2 * righttype + 1;
            clause_posion++;
            clause[3 * clause_posion] = 2 * node.output_vari + 1;
            clause[3 * clause_posion + 1] = 2 * vari;
            clause_posion++;
            clause[3 * clause_posion] = 2 * node.output_vari;
            clause[3 * clause_posion + 1] = 2 * righttype;
            clause[3 * clause_posion + 2] = 2 * vari;
            clause_posion++;
        }else if(lefttype != 1 && lefttype != 0 && righttype == 0){
            clause[3 * clause_posion] = 2 * node.output_vari + 1;
            clause[3 * clause_posion + 1] = 2 * righttype + 1;
            clause_posion++;
            clause[3 * clause_posion] = 2 * node.output_vari + 1;
            clause[3 * clause_posion + 1] = 2 * vari;
            clause_posion++;
            clause[3 * clause_posion] = 2 * node.output_vari;
            clause[3 * clause_posion + 1] = 2 * righttype;
            clause[3 * clause_posion + 2] = 2 * vari;
            clause_posion++;
        }else if(lefttype != 1 && lefttype != 0 && righttype != 1 && righttype != 0){
            clause[3 * clause_posion] = 2 * node.output_vari + 1;
            clause[3 * clause_posion + 1] = 2 * vari;
            clause[3 * clause_posion + 2] = 2 * rightnode;
            clause_posion++;
            clause[3 * clause_posion] = 2 * node.output_vari + 1;
            clause[3 * clause_posion + 1] = 2 * leftnode;
            clause[3 * clause_posion + 2] = 2 * vari + 1;
            clause_posion++;
            clause[3 * clause_posion] = 2 * vari + 1;
            clause[3 * clause_posion + 1] = 2 * leftnode;
            clause[3 * clause_posion + 2] = 2 * node.output_vari;
            clause_posion++;
            clause[3 * clause_posion] = 2 * vari;
            clause[3 * clause_posion + 1] = 2 * node.output_vari;
            clause[3 * clause_posion + 2] = 2 * leftnode + 1;
            clause_posion++;
            clause[3 * clause_posion] = 2 * leftnode + 1;
            clause[3 * clause_posion + 1] = 2 * rightnode + 1;
            clause[3 * clause_posion + 2] = 2 * node.output_vari;
            clause_posion++;
            clause[3 * clause_posion] = 2 * node.output_vari + 1;
            clause[3 * clause_posion + 1] = 2 * rightnode;
            clause[3 * clause_posion + 2] = 2 * leftnode;
            clause_posion++;
        }
    }
};

int main(){
    clock_t start, finish;
    start = clock();
    std::ofstream outFile;
    srand (1);
    MWVC ojb;
    outFile.open("/Users/haoji/Downloads/test.txt");

    for(int i = 3; i < Lim; i++){
        corfficient[i] = i;
        variable[i] = 1;
    }

    for(int i = 3; i < Lim; i++){
        int idx = rand() % (Lim - i) + i;
        int temp = corfficient[idx];
        corfficient[idx] = corfficient[i];
        corfficient[i] = temp;
    }

//  corfficient[0] = rand() % ((Lim + 1)* Lim / 2);
//  corfficient[0] = rand() % (10*Lim);
//  corfficient[0] = rand() % (100*Lim);
    corfficient[0] = rand() % Lim;

    outFile << "--初始化数据--" << std::endl;
    for(int i = 0; i <= Lim; i++)
        outFile << corfficient[i] << ";"  << variable[i] << ";" << std::endl;
    outFile << "--初始化数据--" << std::endl;

    int res = ojb.PBOtoPBC(corfficient, variable);

    int j = corfficient.size() - corfficient[1];
    outFile << "-----" << std::endl;
    outFile << "Pseudo Boolean Part的长度：" << j << std::endl;
    int sum = 0;
    int index = corfficient[1];
    int leftover = std::accumulate(std::begin(corfficient) + index, std::end(corfficient), 0);
    outFile <<"初始的index："<< index  << "剩余参数和：" << leftover << std::endl;
    BDDNode res1;
    res1 = ojb.bulidBDDs(corfficient, variable, index, sum, leftover, dp);
    for(int i = 0; i <= corfficient.size() - 1; i++)
        outFile << corfficient[i] << ";" << variable[i] << ";" << std::endl;
    outFile << "-----" << std::endl;
    for(auto &element: dp){
        auto fi = element.first;
        auto se = element.second;
        int leftnode = dp[std::make_pair(fi.first + 1,se.leftsum)].output_vari;
        int lefttype = dp[std::make_pair(fi.first + 1,se.leftsum)].leftsum;
        int rightnode = dp[std::make_pair(fi.first + 1,se.leftsum)].output_vari;
        int righttype = dp[std::make_pair(fi.first + 1,se.leftsum)].rightsum;
        ojb.BDDtoSAT(clause, fi.first, reinterpret_cast<BDDNode &>(element), leftnode, lefttype, rightnode, righttype);
        outFile <<"变量"<< fi.first << ";" <<"总和：" << fi.second << "---" << "左节点：" << se.leftsum << ";" << "右节点：" << se.rightsum << "变量编号：" << se.output_vari<< std::endl;
    }
    finish = clock();
    outFile << "number of variable:" << index_vari << std::endl;
    outFile << "clock tick:" << finish - start << std::endl;
    outFile << "running time:" << (double)(finish - start)/ CLOCKS_PER_SEC<< std::endl;

    return 0;
}
