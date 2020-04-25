//
//  main.cpp
//  PageRank
//
//  Created by Vodka_rl on 2020/4/3.
//  Copyright © 2020 Vodka_rl. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#define NUM 150000
#define MAX 100000

struct edge {
    int u;
    int v;
};

struct ndscore {
    int node;
    double scr;
};

edge m_edge[NUM];  // 边
int degree[MAX] = {0};  // 出度
int idnode[MAX] = {0};  // node-id
int nodeid[MAX] = {0};  // id-node
ndscore ndsc[MAX];  // 得分

bool cmp(const ndscore &a, const ndscore &b) { return a.scr > b.scr; }  // 降序

int main(int argc, const char * argv[]) {
    memset(m_edge, 0, sizeof(m_edge));
    memset(ndsc, 0, sizeof(ndsc));
    
    string fileName = "test.txt";
    ifstream file(fileName);  // fileName内容读取到file中
    if (!file.is_open()) {  // 确定文件打开了；
        cout<<"open file failed. "<<endl;
        return 0;
    }
    
    int cnt = 0, N = 0;  // 边数，节点数
    for (; file>>m_edge[cnt].u>>m_edge[cnt].v; cnt++) {
        int node1 = m_edge[cnt].u;
        int node2 = m_edge[cnt].v;
        if (!idnode[node1]) {  // 节点ids对应，方便查询
            nodeid[N] = node1;
            N++;  // 避免id0重复
            idnode[node1] = N;
        }
        if (!idnode[node2]) {
            nodeid[N] = node2;
            N++;
            idnode[node2] = N;
        }
        degree[node1]++;
    }
    
    double r[N], tmpr[N];  // r矩阵、临时r矩阵
    for (int i=0; i<N; i++) {
        r[i] = (double)1/N;
        tmpr[i] = 0.0;
    }
    double e = 1.0;  // epsilon
    double beta = 0.8;
    while (e>0.000000001) {  // 限
        e = 0.0;
        for (int i=0; i<cnt; i++) {
            int tmpi = idnode[m_edge[i].v] - 1;
            int tmpj = idnode[m_edge[i].u] - 1;
            tmpr[tmpi] += r[tmpj]/degree[m_edge[i].u];  // 临时r矩阵更新 r = M·r
        }
        for (int i=0; i<N; i++) {
            cout<<tmpr[i]<<endl;
            tmpr[i] = beta * tmpr[i] + (1-beta)/N;  // 矩阵分配
            cout<<tmpr[i]<<endl;
            e += r[i] > tmpr[i] ? (r[i] - tmpr[i]) : (tmpr[i] - r[i]);  // 计算epsilon
            r[i] = tmpr[i];  // r矩阵更新
            tmpr[i] = 0.0;  // 临时r矩阵z重置
        }
        // cout<<"e: "<<e<<endl;
    }
    for (int i=0; i<N; i++) {  // 节点原始id与得分对应
        ndsc[i].node = nodeid[i];
        ndsc[i].scr = r[i];
    }
    
    vector<ndscore> vec(ndsc, ndsc+N);
    sort(vec.begin(), vec.end(), cmp);
    for (vector<ndscore>::iterator it=vec.begin(); it!=vec.end(); it++) {
        cout<<it->node<<" "<<it->scr<<endl;
    }
    return 0;
}
