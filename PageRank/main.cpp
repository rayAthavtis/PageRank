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
#define NUM 105000

// 节点-得分
struct nscore {
    int node;
    double score;
};

struct edge {
    int src;
    // r被分为的k块，分为几块就几个dst(int dst1; int dst2; ...)
    int dst;  // r分为N块，每块1个dst
};

// 得分降序
bool cmp(const nscore &a, const nscore &b) { return a.score > b.score; }

bool basicPr(string inf, double beta) {
    ifstream file(inf);  // fileName内容读取到file中
    // 确定文件打开
    if (!file) { cout<<"open \""<<inf<<"\" failed! "<<endl; return false; }

    int idnode[MAX] = {0};  // node得id
    int nodeid[MAX] = {0};  // id得node
    
    multimap<int, int> iedge;  // 优化稀疏矩阵，去掉为0的边
    
    int N = 0;  // 节点数
    // 文件读入 节点->节点
    int node1, node2;  // from node1 to node2
    for (int cnt = 0; file>>node1>>node2; cnt++) {
        if (!idnode[node1]) {  // 节点id对应，方便查询
            nodeid[N] = node1;
            N++;  // 避免id赋为0时被认为未赋值
            idnode[node1] = N;
        }
        if (!idnode[node2]) {
            nodeid[N] = node2;
            N++;
            idnode[node2] = N;
        }
        int row = idnode[node1]-1;
        int col = idnode[node2]-1;
        iedge.insert(pair<int, int>(row, col));  // M矩阵
    }
    file.close();

    // cout<<"N: "<<N<<endl;
    double r_old[N], r_new[N];  // r矩阵、临时r矩阵
    for (int i=0; i<N; i++) {
        r_old[i] = (double)1/N;
        r_new[i] = 0.0;
    }
    // double beta = 0.85;
    while (1) {  // 迭代至收敛
        double eps = 0.0;
        for (int i=0; i<N; i++) {  // 扫描r矩阵
            multimap<int, int>::iterator iter = iedge.find(i);
            for (int j=0; j<iedge.count(i); j++, iter++) {  // 扫描M矩阵
                r_new[iter->second] += r_old[i]/iedge.count(i);
            }
        }
        
        double r_sum = 0.0;  // 记录pagerank泄露情况
        for (int i=0; i<N; i++) { r_sum += r_new[i]; }

        for (int i=0; i<N; i++) {
            // 考虑dead ends和spider trap节点
            r_new[i] = beta * (r_new[i] + (1.0-r_sum)/N) + (1.0-beta)/N;
            // cout<<r_new[i]<<endl;
            // 计算epsilon，迭代至收敛
            eps += fabs(r_new[i] - r_old[i]);
            r_old[i] = r_new[i];  // r矩阵更新
            r_new[i] = 0.0;  // 临时r矩阵重置
        }
        // cout<<"eps: "<<eps<<endl;
        if (eps<1e-8) { break; }
    }
    
    nscore ndsc[N];  // 得分
    memset(ndsc, 0, sizeof(ndsc));
    for (int i=0; i<N; i++) {  // 节点原始id与得分对应
        ndsc[i].node = nodeid[i];
        ndsc[i].score = r_old[i];
    }
    
    vector<nscore> vec(ndsc, ndsc+N);
    sort(vec.begin(), vec.end(), cmp);  // 按得分高到低排序
    
    ofstream result;
    result.open("out/pagerank.txt");
    if (!result) { cout<<"open \"pagerank.txt\" failed! "<<endl; return false; }
    int num = 0;  // 记录前100个
    for (vector<nscore>::iterator it=vec.begin(); it!=vec.end() && num<100; it++, num++) {
        cout<<it->node<<" "<<it->score<<endl;
        result<<it->node<<" "<<it->score<<endl;
    }
    result.close();
    
    return true;
}

bool blockStripePr(string inf, double beta) {
    ifstream file(inf);  // fileName内容读取到file中
    // 确定文件打开
    if (!file) { cout<<"open \""<<inf<<"\" failed! "<<endl; return false; }
    
    edge m_edge[NUM];
    memset(m_edge, 0, sizeof(m_edge));

    int degree[MAX] = {0};  // 出度
    int idnode[MAX] = {0};  // node得id
    int nodeid[MAX] = {0};  // id得node

    int cnt = 0, N = 0;  // 边数、节点数
    // 文件读入 节点->节点
    int node1, node2;  // from node1 to node2
    for (; file>>node1>>node2; cnt++) {
        if (!idnode[node1]) {  // 节点id对应，方便查询
            nodeid[N] = node1;
            N++;  // 避免id赋为0时被认为未赋值
            idnode[node1] = N;
        }
        if (!idnode[node2]) {
            nodeid[N] = node2;
            N++;
            idnode[node2] = N;
        }
        m_edge[cnt].src = idnode[node1]-1;
        m_edge[cnt].dst = idnode[node2]-1;  // M矩阵分块
        degree[idnode[node1]-1]++;
    }
    file.close();

    // cout<<"N: "<<N<<endl;
    double r_old[N], r_new[N];  // r矩阵、旧r矩阵
    for (int i=0; i<N; i++) {
        r_old[i] = (double)1/N;
        r_new[i] = 0.0;
    }
    // double beta = 0.85;
    while (1) {  // 迭代至收敛
        double eps = 0.0;
        for (int i=0; i<cnt; i++) {  // M分为的cnt块
            r_new[m_edge[i].dst] += r_old[m_edge[i].src]/degree[m_edge[i].src];
        }
        
        double r_sum = 0.0;  // 记录pagerank泄露情况
        for (int i=0; i<N; i++) { r_sum += r_new[i]; }

        for (int i=0; i<N; i++) {
            // 考虑dead ends和spider trap节点
            r_new[i] = beta * (r_new[i] + (1.0-r_sum)/N) + (1.0-beta)/N;
            // cout<<r_new[i]<<endl;
            // 计算epsilon，迭代至收敛
            eps += fabs(r_new[i] - r_old[i]);
            r_old[i] = r_new[i];  // r矩阵更新
            r_new[i] = 0.0;  // 临时r矩阵重置
        }
        // cout<<"eps: "<<eps<<endl;
        if (eps<1e-8) { break; }
    }
    
    nscore ndsc[N];  // 得分
    memset(ndsc, 0, sizeof(ndsc));
    for (int i=0; i<N; i++) {  // 节点原始id与得分对应
        ndsc[i].node = nodeid[i];
        ndsc[i].score = r_old[i];
    }
    
    vector<nscore> vec(ndsc, ndsc+N);
    sort(vec.begin(), vec.end(), cmp);  // 按得分高到低排序
    
    ofstream result;
    result.open("out/pagerank.txt");
    if (!result) { cout<<"open \"pagerank.txt\" failed! "<<endl; return false; }
    int num = 0;  // 记录前100个
    for (vector<nscore>::iterator it=vec.begin(); it!=vec.end() && num<100; it++, num++) {
        cout<<it->node<<" "<<it->score<<endl;
        result<<it->node<<" "<<it->score<<endl;
    }
    result.close();
    
    return true;
}

int main(int argc, const char * argv[]) {
    // 使用相对路径，在Xcode->product->scheme->edit->run->options->working directory修改工作目录
    string fileName = "PageRank_WikiData.txt";
    // string fileName = "test.txt";

    // cout<<"beta please: "; cin>>bt;
    double bt = 0.85;
    // basicPr(fileName, bt);
    blockStripePr(fileName, bt);
    
    return 0;
}
