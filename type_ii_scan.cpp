// Deterministic nonzero-syndrome probes and translated Gray blocks.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
using U64 = std::uint64_t;
struct Space {
    std::string id;
    int rank, dimension;
    U64 offset;
    std::vector<U64> rows, basis, seeds;
    std::vector<int> rhs, free;
    std::vector<std::array<int,5>> edges;
    std::unordered_set<U64> visited;
    std::map<std::string,U64> local;
};
struct UF {
    std::array<int,99> p, rank;
    int groups=99;
    UF() { for(int i=0;i<99;i++){p[i]=i;rank[i]=0;} }
    int find(int a) { while(p[a]!=a){p[a]=p[p[a]];a=p[a];}return a; }
    void join(int a,int b) {
        a=find(a);b=find(b);if(a==b)return;
        if(rank[a]<rank[b])std::swap(a,b);
        p[b]=a;if(rank[a]==rank[b])rank[a]++;groups--;
    }
};
struct Hit { std::string tri, stage; U64 mask; int syndrome; };
std::map<std::string,Hit> hits;
struct Phase { int radius, recovered; U64 attempted, distinct, reads; double elapsed; };
std::vector<Phase> phases;
U64 attempted=0, distinct=0, dividing=0, reads=0;
auto started=std::chrono::steady_clock::now();
std::string output;
double seconds(){return std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();}
int syndrome(const Space& s,U64 mask){
    int value=0;
    for(int i=0;i<s.rank;i++)value|=(__builtin_parityll(s.rows[i]&mask)^s.rhs[i])<<i;
    return value;
}
bool target(const std::string& key){
    if(key[0]=='U'){int r=std::stoi(key.substr(1));return r>=0&&r<=14;}
    if(key[0]!='N')return false;
    auto comma=key.find(',');int s=std::stoi(key.substr(1,comma-1)),t=std::stoi(key.substr(comma+1));
    return s>=0&&t>=1&&s+t<=13;
}
std::string read(const Space& s,U64 mask){
    UF u;std::array<int,308> cross;int nc=0;
    for(int k=0;k<(int)s.edges.size();k++){
        auto e=s.edges[k];
        if(((mask>>e[2])^(mask>>e[3])^e[4])&1)cross[nc++]=k;
        else u.join(e[0],e[1]);
    }
    std::array<int,99> roots,degree{};
    std::array<unsigned char,99*99> seen{};
    for(int i=0;i<99;i++)roots[i]=u.find(i);
    int root=-1,ne=0;
    for(int k=0;k<nc;k++){
        auto e=s.edges[cross[k]];int a=roots[e[0]],b=roots[e[1]];
        if(a==b){if(root>=0&&root!=a)throw std::runtime_error("multiple pseudoline regions");root=a;continue;}
        int key=std::min(a,b)*99+std::max(a,b);
        if(!seen[key]){seen[key]=1;degree[a]++;degree[b]++;ne++;}
    }
    if(root<0||ne!=u.groups-1)throw std::runtime_error("invalid region tree");
    int container=-1;
    for(int i=0;i<99;i++)if(i!=root&&degree[i]>1){if(container>=0)return "deep";container=i;}
    if(container<0)return "U"+std::to_string(degree[root]);
    int outer=degree[root]-1,inner=degree[container]-1;
    if(outer+inner+2!=u.groups)throw std::runtime_error("invalid single nest");
    return "N"+std::to_string(outer)+","+std::to_string(inner);
}
void save(){
    std::ofstream out(output);
    out<<"{\"algorithm\":\"deterministic-nonzero-syndrome-v1\",\"attempted\":"<<attempted
       <<",\"distinct_probes\":"<<distinct<<",\"type_i_skipped\":"<<dividing
       <<",\"unsigned_reads\":"<<reads<<",\"elapsed_seconds\":"<<seconds()<<",\"records\":[";
    bool first=true;
    for(const auto& [key,h]:hits){
        if(!first)out<<",";first=false;
        out<<"{\"key\":\""<<key<<"\",\"triangulation\":\""<<h.tri<<"\",\"mask\":"<<h.mask
           <<",\"syndrome\":"<<h.syndrome<<",\"stage\":\""<<h.stage<<"\"}";
    }
    out<<"],\"phases\":[";first=true;
    for(const auto& p:phases){
        if(!first)out<<",";first=false;
        out<<"{\"radius\":"<<p.radius<<",\"recovered\":"<<p.recovered<<",\"attempted\":"<<p.attempted
           <<",\"distinct\":"<<p.distinct<<",\"unsigned_reads\":"<<p.reads<<",\"elapsed_seconds\":"<<p.elapsed<<"}";
    }
    out<<"],\"target_complete\":"<<(hits.size()==106?"true":"false")<<"}\n";
}
void probe(Space& s,U64 mask,const std::string& stage){
    attempted++;
    if(!s.visited.insert(mask).second)return;
    distinct++;int syn=syndrome(s,mask);
    if(!syn){dividing++;return;}
    reads++;auto key=read(s,mask);
    if(!target(key))throw std::runtime_error("type II result outside target: "+key);
    s.local.emplace(key,mask);
    if(hits.emplace(key,Hit{s.id,stage,mask,syn}).second){
        std::cout<<"FOUND "<<hits.size()<<" "<<key<<" "<<s.id<<" "<<stage<<" seconds="<<seconds()<<std::endl;
    }
}
int main(int argc,char**argv){
    try{
        if(argc<3)throw std::runtime_error("usage: type_ii_scan input.scan results.json [radius] [gray-block]");
        output=argv[2];int radius=argc>3?std::stoi(argv[3]):2;
        U64 block=argc>4?std::stoull(argv[4]):65536;
        std::ifstream in(argv[1]);in.exceptions(std::ios::failbit|std::ios::badbit);
        std::string version;int count;in>>version>>count;
        if(version!="type-ii-v1")throw std::runtime_error("input version");
        std::vector<Space> spaces(count);
        for(auto& s:spaces){
            in>>s.id>>s.rank>>s.dimension>>s.offset;
            s.rows.resize(s.rank);s.rhs.resize(s.rank);
            for(int i=0;i<s.rank;i++)in>>s.rows[i]>>s.rhs[i];
            s.basis.resize(s.dimension);for(auto& b:s.basis)in>>b;
            int n;in>>n;s.seeds.resize(n);for(auto& b:s.seeds)in>>b;
            s.free.resize(33);for(auto& b:s.free)in>>b;
            in>>n;s.edges.resize(n);for(auto& e:s.edges)for(auto& v:e)in>>v;
            if(n!=308||syndrome(s,s.offset))throw std::runtime_error("invalid space");
        }
        for(int r=1;r<=radius&&hits.size()<106;r++){
            for(auto& s:spaces){
                std::vector<U64> seeds=s.seeds;
                std::sort(seeds.begin(),seeds.end());seeds.erase(std::unique(seeds.begin(),seeds.end()),seeds.end());
                for(auto mask:seeds){
                    for(int a=0;a<33&&hits.size()<106;a++){
                        if(r==1)probe(s,mask^(1ULL<<s.free[a]),"local-1");
                        else for(int b=a+1;b<33&&hits.size()<106;b++){
                            if(r==2)probe(s,mask^(1ULL<<s.free[a])^(1ULL<<s.free[b]),"local-2");
                            else for(int c=b+1;c<33&&hits.size()<106;c++)
                                probe(s,mask^(1ULL<<s.free[a])^(1ULL<<s.free[b])^(1ULL<<s.free[c]),"local-3");
                        }
                    }
                    if(hits.size()==106)break;
                }
                if(hits.size()==106)break;
            }
            phases.push_back({r,int(hits.size()),attempted,distinct,reads,seconds()});
            std::cout<<"PHASE radius="<<r<<" count="<<hits.size()<<" reads="<<reads<<" seconds="<<seconds()<<std::endl;save();
        }
        // Translations by coordinate vectors are cheap guaranteed-type-II blocks.
        for(auto& s:spaces){
            if(hits.size()==106)break;
            std::set<int> used;
            for(int bit:s.free){
                int syn=syndrome(s,s.offset^(1ULL<<bit));
                if(!syn||!used.insert(syn).second)continue;
                U64 mask=s.offset^(1ULL<<bit),end=std::min<U64>(block,1ULL<<s.dimension);
                for(U64 n=0;n<end&&hits.size()<106;n++){
                    if(n)mask^=s.basis[__builtin_ctzll(n)];
                    probe(s,mask,"translated-gray");
                }
                if(hits.size()==106)break;
            }
            std::cout<<"GRAY "<<s.id<<" count="<<hits.size()<<" reads="<<reads<<" seconds="<<seconds()<<std::endl;save();
            if(hits.size()==106)break;
        }
        save();std::cout<<"DONE "<<hits.size()<<"/106 reads="<<reads<<" seconds="<<seconds()<<std::endl;
    }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
