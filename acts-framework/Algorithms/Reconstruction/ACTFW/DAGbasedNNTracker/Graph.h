#ifndef _GRAPHh_H_
#define _GRAPHh_H_
// Implements weighted directed graph
// M.Kunze, Heidelberg University, 2018

#include <stdio.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <algorithm>
#include <functional>

class XMLP;
class ReadMLP1;
class ReadMLP2;
class ReadMLP3;
namespace TMVA { class Reader; }

template <typename T>
class Graph
{
private:
    std::set<T>                      fNodes;
    std::map<T, std::map<T, float> > fEdges;
    std::map<T, std::vector<int> >   fData;
    std::set<long long>              fHash[48];
    XMLP                             *n1,
                                     *n2,
                                     *n3;
    ReadMLP1 *mlp1;
    ReadMLP2 *mlp2;
    ReadMLP3 *mlp3;
    TMVA::Reader *reader1,*reader2,*reader3;
    float x1[12],x2[12],x3[12];
    
public:
    
    Graph() : n1(NULL), n2(NULL), n3(NULL), mlp1(NULL), mlp2(NULL), mlp3(NULL), reader1(NULL), reader2(NULL), reader3(NULL) {};
    ~Graph() {};

    void setNet1(XMLP *net) {n1 = net;}
    void setNet2(XMLP *net) {n2 = net;}
    void setNet3(XMLP *net) {n3 = net;}
    XMLP *net1() const { return n1;}
    XMLP *net2() const { return n2;}
    XMLP *net3() const { return n3;}
  
    void setNet1(ReadMLP1 *net) {mlp1 = net;}
    void setNet2(ReadMLP2 *net) {mlp2 = net;}
    void setNet3(ReadMLP3 *net) {mlp3 = net;}
    ReadMLP1 *getNet1() const { return mlp1;}
    ReadMLP2 *getNet2() const { return mlp2;}
    ReadMLP3 *getNet3() const { return mlp3;}
    
    void setReader1(TMVA::Reader *net) {reader1 = net;}
    void setReader2(TMVA::Reader *net) {reader2 = net;}
    void setReader3(TMVA::Reader *net) {reader3 = net;}
    TMVA::Reader *getReader1() const { return reader1;}
    TMVA::Reader *getReader2() const { return reader2;}
    TMVA::Reader *getReader3() const { return reader3;}
    float *getX1() { return x1;}
    float *getX2() { return x2;}
    float *getX3() { return x3;}

    void add(T n) // node
    {
        //if (fEdges.find(n)==fEdges.end()) return; // The node exists
        fNodes.insert(n);
        (void)fEdges[n];
        (void)fData[n];
    }
    
    void add(const T& n1, const T& n2, float d = 0.0, bool increment = false) // edge
    {
        add(n1);
        add(n2);
        auto& adj = fEdges[n1];
        auto  n   = adj.find(n2);
        if (n != adj.end()) {
            float& d1 = n->second;
            if (increment)
                d1++;
            else {
                if (d1==0) d1 = d;
                d1 = 0.25*(3.*d1+d);
            }
        } else {
            adj[n2] = d;
        }
    }

    const std::set<T>& nodes() const
    {
        return fNodes;
    }
    
    const std::map<T, float>& edges(const T& n) const
    {
        static const std::map<T, float> null;
        if (fEdges.find(n)==fEdges.end()) return null; // The node does not exist
        return fEdges.at(n);
    }

    std::vector<int>& data(const T& n)
    {
        static std::vector<int> null;
        null.clear();
        if (fData.find(n)==fData.end()) return null; // The node does not exist
        return fData.at(n);
    }
    
    std::set<long long>& hash(int l)
    {
        static std::set<long long> null;
        if (l<0||l>=48) return null; // The set does not exist
        return fHash[l];
    }
    
    void clear()
    {
        for (auto &n : fNodes) fData[n].clear();
        for (int i=0;i<48;i++) fHash[i].clear();
    }
    
    bool areConnected(const T& n1, const T& n2, float& d) const
    {
        if (fEdges.find(n1)==fEdges.end()) return false; // The node does not exist
        auto c = fEdges.at(n1);
        auto q = c.find(n2);
        if (q != c.end()) {
            d = q->second;
            return true;
        } else {
            return false;
        }
    }
    
    bool areConnected(const T& n1, const T& n2) const
    {
        float d;
        return areConnected(n1, n2, d);
    }
    
    void print(std::ostream& file=std::cout)
    {
        std::string sep      = "";
        bool   nok = false;
        
        file << "Graph {";
        for (const T& n : nodes()) {
            nok    = true;
            bool cok = false;
            if (edges(n).size()==0) continue;
            for (const auto& c : edges(n)) {
                cok = true;
                if (c.second == 0) {
                    file << sep << n << "->" << (c.first);
                } else {
                    file << sep << n << '-' << c.second << "->" << (c.first);
                }
                sep = ", ";
            }
            if (!cok) {
                file << sep << n;
            }
            sep = ", ";
        }
        file << "}" << std::endl;
    }

};

#define GCOUNT 0
template <typename T>
inline std::vector<std::vector<T> > serialize(const Graph<T>& G)
{
    typedef std::function<void(const Graph<T>& G, const T& N, std::set<T>& V, std::vector<T>& R)> Visitfun;
    Visitfun visit = [&visit](const Graph<T>& G, const T& N, std::set<T>& V, std::vector<T>& R) {
        if (V.find(N) == V.end()) {
            V.insert(N);
            //const auto &e = G.edges(N);
            //const auto &best = std::max_element(e.begin(), e.end(),e.value_comp());
            //if (best->second >= GCOUNT) visit(G, best->first, V, R);
            for (const auto& e : G.edges(N)) {
                if (e.second >= GCOUNT) visit(G, e.first, V, R);
            }
            R.push_back(N);
        }
    };
    
    std::vector<std::vector<T> > v;
    std::set<T>    V;
    for (const T& N : G.nodes()) {
        std::vector<T> R;
        visit(G, N, V, R);
        std::sort(R.begin(),R.end());
        if (R.size()>0) v.push_back(std::vector<T> (R));
    }
    return v;
}


template <typename T>
inline std::vector<T> serialize(const Graph<T>& G, const T& N)
{
    typedef std::function<void(const Graph<T>& G, const T& N, std::set<T>& V, std::vector<T>& R)> Visitfun;
    Visitfun visit = [&visit](const Graph<T>& G, const T& N, std::set<T>& V, std::vector<T>& R) {
        if (V.find(N) == V.end()) {
            V.insert(N);
            const auto &e = G.edges(N);
            const auto &best = std::max_element(e.begin(), e.end(),e.value_comp());
            if (best->second >= GCOUNT) visit(G, best->first, V, R);
//            for (const auto& e : G.edges(N)) {
//                if (e.second >= GCOUNT) visit(G, e.first, V, R);
//            }
            R.push_back(N);
        }
    };
    
    std::vector<T> R;
    std::set<T>    V;
    visit(G, N, V, R);
    return R;
}

template <typename T>
inline std::vector<std::vector<T> > parallelize(const Graph<T>& g)
{
    //-----------------------------------------------------------
    // Find the level of a node n -> {m1,m2,...} such that
    //        level(n -> {})            = 0
    //        level(n -> {m1,m2,...})    = 1 + max(level(mi))
    //-----------------------------------------------------------
    typedef std::function<int(const Graph<T>& g, const T& n1, std::map<T, float>&)> Levelfun;
    
    Levelfun level = [&level](const Graph<T>& g, const T& n1, std::map<T, float>& levelcache) -> int {
        auto p = levelcache.find(n1);
        if (p != levelcache.end()) {
            return p->second;
        } else {
            int l = -1;
            for (const auto& e : g.edges(n1)) {
                l = std::max(l, level(g, e.first, levelcache));
            }
            return levelcache[n1] = l + 1;
        }
    };
    
    std::map<T, float> levelcache;
    // compute the level of each node in the graph
    int l = -1;
    for (const T& n : g.nodes()) {
        l = std::max(l, level(g, n, levelcache));
    }
    // create a graph for each level and place
    // each node in the appropriate level
    std::vector<std::vector<T> > v;
    v.resize(l + 1);
    for (const T& n : g.nodes()) {
        v[levelcache[n]].push_back(n);
    }
    
    return v;
}


// print on stream

inline std::ostream& operator<<(std::ostream& file, const Graph<std::pair<int,int> >& g)
{
    for (const std::pair<int,int>& n : g.nodes()) {
        if (g.edges(n).size()==0) continue;
        for (const auto& c : g.edges(n)) {
            file << n.first << " " << n.second << " " << c.first.first << " " << c.first.second << " " << c.second << '\n';
        }
    }
    
    return file;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& file, const Graph<T>& g)
{
    for (const T& n : g.nodes()) {
        if (g.edges(n).size()==0) continue;
        for (const auto& c : g.edges(n)) {
                file << n << " " << c.second << " " << c.first << '\n';
        }
    }
    
    return file;
}

inline std::istream& operator>>(std::istream& file, Graph<std::pair<int,int> >& g)
{
    int n1,n2,n3,n4;
    float weight;
    while (file.good()) {
        file >> n1 >> n2 >> n3 >> n4 >> weight;
        std::pair<int,int> p1 = std::make_pair(n1,n2);
        std::pair<int,int> p2 = std::make_pair(n3,n4);
        g.add(p1,p2,weight);
    }
    
    return file;
}

template <typename T>
inline std::istream& operator>>(std::istream& file, Graph<T>& g)
{
    T n1,n2;
    float weight;
    while (file.good()) {
        file >> n1 >> weight >> n2;
        g.add(n1,n2,weight);
    }
    
    return file;
}

#endif

