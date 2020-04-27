//
//  main.cpp
//  PageRank
//
//  Created by ray on 2020/4/3.
//  Copyright © 2020 ray. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <math.h>

using namespace std;

#define MAX 100000

// 节点-得分
struct nscore {
    int node;
    double score;
};

// 得分降序
bool cmp(const nscore &a, const nscore &b) { return a.score > b.score; }

bool iPageRank(string inf, double beta) {
    ifstream file(inf);  // fileName内容读取到file中
    // 确定文件打开
    if (!file) { cout<<"open \"test.txt\" failed! "<<endl; return false; }
    
    int degree[MAX] = {0};  // 出度
    int idnode[MAX] = {0};  // node得id
    int nodeid[MAX] = {0};  // id得node
    
    multimap<int, int> iedge;  // src, dst1, dst2 优化稀疏矩阵，去掉为0的边
    
    int N = 0;  // 节点数
    // 文件读入 节点->节点
    int node1, node2;  // from node1 to node2
    for (int cnt = 0; file>>node1>>node2; cnt++) {
        if (!idnode[node1]) {  // 节点ids对应，方便查询
            nodeid[N] = node1;
            N++;  // 避免id赋为0时被认为未赋值
            idnode[node1] = N;
        }
        if (!idnode[node2]) {
            nodeid[N] = node2;
            N++;
            idnode[node2] = N;
        }
        degree[node1]++;
        int row = idnode[node1]-1;
        int col = idnode[node2]-1;
        iedge.insert(pair<int, int>(col, row));  // M矩阵 M_ji
    }
    file.close();

    // cout<<"N: "<<N<<endl;
    double r[N], tmpr[N];  // r矩阵、临时r矩阵
    for (int i=0; i<N; i++) {
        r[i] = (double)1/N;
        tmpr[i] = 0.0;
    }
    // double beta = 0.85;
    while (1) {  // 迭代至收敛
        double e = 0.0;
        for (int i=0; i<N; i++) {  // 分块更新
            // 每次只对iedge一行的对应r[i]更新
            multimap<int, int>::iterator iter = iedge.find(i);
            for (int j=0; j<iedge.count(i); j++, iter++) {
                int node = nodeid[iter->second];
                // cout<<node<<":"<<degree[node]<<" ";
                // 临时r矩阵更新: r_j += M_ji·r_i
                // cout<<r[iter->second]<<" * "<<1.0/degree[node]<< " = "<<r[iter->second]/degree[node]<<endl;
                tmpr[i] += r[iter->second]/degree[node];
                // cout<<"r["<<i<<"] = "<<tmpr[i]<<endl;
            }
            // 考虑dead ends和spider trap节点: r_j = beta*(M_ji·r_i) + (1-beta)/N
            // cout<<"r["<<i<<"] = "<<beta*tmpr[i]<<" + "<<(1.0-beta)/N<< " = ";
            tmpr[i] = beta * tmpr[i] + (1.0-beta)/N;
            // cout<<tmpr[i]<<endl;
            // 计算epsilon，迭代至收敛
            e += pow((r[i] - tmpr[i]), 2);
            r[i] = tmpr[i];  // r矩阵更新
            tmpr[i] = 0.0;  // 临时r矩阵重置
        }
        e = sqrt(e);
        // cout<<"e: "<<e<<endl;
        if (e<0.01) { break; }
    }
    
    nscore ndsc[N];  // 得分
    memset(ndsc, 0, sizeof(ndsc));
    for (int i=0; i<N; i++) {  // 节点原始id与得分对应
        ndsc[i].node = nodeid[i];
        ndsc[i].score = r[i];
    }
    
    vector<nscore> vec(ndsc, ndsc+N);
    sort(vec.begin(), vec.end(), cmp);  // 按得分高到低排序

    double all = 0.0;
    for (vector<nscore>::iterator it=vec.begin(); it!=vec.end(); it++) {
        cout<<it->node<<" "<<it->score<<endl;
        all += it->score;
    }
    cout<<"all: "<<all<<endl;
    
//    ofstream result;
//    result.open("out/pagerank.txt");
//    if (!result) { cout<<"open \"pagerank.txt\" failed! "<<endl; return false; }
//    int num = 0;  // 记录前100个
//    for (vector<nscore>::iterator it=vec.begin(); it!=vec.end() && num<100; it++, num++) {
//        cout<<it->node<<" "<<it->score<<endl;
//        result<<it->node<<" "<<it->score<<endl;
//    }
//    cout<<all<<endl;
//    result.close();
    
    return true;
}

int main(int argc, const char * argv[]) {
    // 使用相对路径，在Xcode->product->scheme->edit->run->options->working directory修改工作目录
    string fileName = "PageRank_WikiData.txt";
    // string fileName = "test.txt";

    double bt;
    cout<<"beta please: "; cin>>bt;
    // bt = 0.8;  // 0.85
    if (bt>0.5 && bt<0.95) { iPageRank(fileName, bt); }
    else { cout<<"error! "<<endl; return -1; }
    
    return 0;
}

/*
 (4037, 0.0046127158911675485)
 (15, 0.0036812207295292792)
 (6634, 0.003524813657640259)
 (2625, 0.0032863743692309023)
 (2398, 0.0026053331717250192)
 */
