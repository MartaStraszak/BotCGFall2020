#include <bits/stdc++.h>


using namespace std;


#define sim template < class c
#define ris return * this
#define dor > debug & operator <<
#define eni(x) sim > typename   enable_if<sizeof dud<c>(0) x 1, debug&>::type operator<<(c i) {
sim > struct rge { c b, e; };
sim > rge<c> range(c i, c j) { return rge<c>{i, j}; }
sim > auto dud(c* x) -> decltype(cerr << *x, 0);
sim > char dud(...);
struct debug {
~debug() { cerr << endl; }
eni(!=) cerr << boolalpha << i; ris; }
eni(==) ris << range(begin(i), end(i)); }
sim, class b dor(pair < b, c > d) {
  ris << "(" << d.first << ", " << d.second << ")";
}
sim dor(rge<c> d) {
  *this << "[";
  for (auto it = d.b; it != d.e; ++it)
    *this << ", " + 2 * (it == d.b) << *it;
  ris << "]";
}
};
#define var(...) " [" << #__VA_ARGS__ ": " << (__VA_ARGS__) << "] "

using vpvii=vector<pair<vector<int>, int>>;

#define pb push_back
#define mp make_pair
#define st first
#define nd second
#define ALL(x) (x).begin(), (x).end()
#define SZ(x) ((int)(x).size())
#define un(x) x.erase(unique(ALL(x)),x.end())
#define REP(i,x,v)for(int i=x;i<=v;i++)
#define REPD(i,x,v)for(int i=x;i>=v;i--)
#define FOR(i,v)for(int i=0;i<v;i++)
#define TRAV(i,v)for(auto &i:(v))
#define REMIN(x,v) x=min(x,v);
#define REMAX(x,v) x=max(x,v);

mt19937_64 rng(134569);
int rd(int l, int r) {return uniform_int_distribution<int>(l, r)(rng);}

using ll=long long;
using pii=pair<int,int>;
using pll=pair<ll,ll>;
using vi=vector<int>;
using vd=vector<double>;
using vll=vector<ll>;






map<vector<int>,int> stats;
vector<pair<double,vector<int>>> usefulness;

void addNewStat(vector<int> delta, int val){
	if(stats.count(delta)==0) stats[delta]=0;
	stats[delta]+=1;
}

string load1(const string& path) {
    ifstream file(path);
    return string((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
}

void computeScoreFromStats(){
	vector<pair<int, vector<int>>> vectStats;
	int sum=0;
	vector<vector<int>> badSpells = {{2,0,0,0},{-1,1,0,0},{0,-1,1,0},{0,0,-1,1}};

	for(auto it: stats){
		if(it.first==badSpells[0] || it.first==badSpells[1] || it.first ==badSpells[2]  || it.first== badSpells[3]|| it.second<100) continue;
		sum+=it.second;
		vectStats.push_back({it.second,it.first});
	}
	assert(SZ(vectStats)==42);
	for(auto p: vectStats){
		usefulness.push_back({1.0*p.first/sum * 42, p.second});
	}
	sort(usefulness.begin(), usefulness.end());
	for(int i=0; i<SZ(usefulness);i++){
		vector<int> d = usefulness[i].second;
		cout<<i+1<<" "<< usefulness[i].first << " ["<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<"]"<<endl; ;
	}
}

void presentStats(){
	vector<pair<int, vector<int>>> vectStats;
	for(auto it: stats){
		vectStats.push_back({-it.second,it.first});
	}
	sort(vectStats.begin(), vectStats.end());
	for(int i=0; i<SZ(vectStats);i++){
		vector<int> d = vectStats[i].second;
		cout<<i+1<<" "<< -vectStats[i].first << " ["<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<"]"<<endl; ;
	}
}



void extractData(string s){
	int pos=0;
	while(1){
		auto occ=s.find("LOGSPELLB", pos);
		if(occ == string::npos) break;
		auto occEnd = s.find("LOGSPELLE", occ);
		assert(occEnd!=string::npos);
		string d = s.substr(occ+9, occEnd - occ-9);
		pos=occEnd+1;
		stringstream ss(d);
		vector<int> delta;
		for(int i=0; i<4; i++){
			int x;
			ss >> x;
			delta.push_back(x);
		}
		// debug() << var(delta);
		addNewStat(delta,1);

	}
}

void readLogsFromFile(string fileName){
	string s = load1(fileName);
	// debug() <<var(fileName)<< var(SZ(s));
	extractData(s);

}


void collectDataFromFiles(int fNo){
	for(int i=1; i<fNo+1;i++){
		stringstream ss;
		ss << "logs/game"<<i<<".json";
		string filePath = ss.str();
		readLogsFromFile(filePath);
	}
}

void transformUseToStrVect(){
	cout << "map<vector<int>, double> usefulness =  \n";
	cout << "{";
	for(int i=0; i<SZ(usefulness);i++){
		if(i!=0 && i%6==0) cout<< "\n";
		auto p = usefulness[i];
		cout << "{" <<"{"<< p.second[0]<<","<< p.second[1]<<","<< p.second[2]<<","<< p.second[3]<<"},"<<p.first <<"},";
	}
	cout<<"};\n";
}

int main(){
	collectDataFromFiles(1000);
	// readLogsFromFile("logs/game1.json");
	presentStats();
	computeScoreFromStats();
	transformUseToStrVect();

}