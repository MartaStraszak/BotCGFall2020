#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <queue>
//#include <bits/stdc++.h>
#include <chrono>
#include <set>
#include <cassert>
#include <iomanip>
#include <map>
#include <cmath>
#include <bitset>

using namespace std;

#define SZ(x) ((int)(x).size())
#define STATES 1001
#define INF (1<<30)
#define PLEN 14
#define NPOTIONS 5
#define MAXSPELLS 128
#define MAXSTATES 1600
#define BEAMMAX 1100
#define LEARNTURNS 9

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

mt19937_64 rng(134569);
int rd(int l, int r) {return uniform_int_distribution<int>(l, r)(rng);}

map<vector<int>, double> usefulness =  
{{{-3,0,0,1},0.367639},{{2,3,-2,0},0.424406},{{2,-3,2,0},0.466306},{{-4,0,2,0},0.469009},{{-2,0,1,0},0.623093},{{2,1,-2,1},0.637961},
{{2,2,0,-1},0.682564},{{0,2,-1,0},0.69608},{{3,-2,1,0},0.706893},{{0,3,2,-2},0.733925},{{0,2,-2,1},0.73798},{{0,-3,0,2},0.7542},
{{-1,-1,0,1},0.756903},{{3,-1,0,0},0.766364},{{3,0,1,-1},0.783935},{{0,3,0,-1},0.808264},{{-5,0,0,2},0.815022},{{0,0,-3,3},0.817725},
{{1,-3,1,1},0.831242},{{0,0,-2,2},0.846109},{{1,1,1,-1},0.856922},{{-5,0,3,0},0.8799},{{-3,1,1,0},0.902877},{{0,-3,3,0},0.924503},
{{2,-2,0,1},0.924503},{{-4,0,1,1},0.932612},{{1,1,3,-2},0.943425},{{0,-2,2,0},0.969106},{{0,0,2,-1},0.992083},{{1,2,-1,0},1.01776},
{{4,1,-1,0},1.19347},{{-2,0,-1,2},1.20429},{{-3,3,0,0},1.21104},{{1,1,0,0},1.24348},{{0,0,1,0},1.29349},{{-2,2,0,0},1.47461},
{{3,0,0,0},1.51245},{{4,0,0,0},1.79494},{{0,2,0,0},1.84225},{{2,1,0,0},1.94767},{{1,0,1,0},2.08148},{{0,0,0,1},2.13149},};

int invIdx[11][11][11][11];
vector<int> adjSpell[STATES]; // (spellIndex)
vector<int> adjPotion[STATES];
vector<int> adjLearnSpell[STATES]; 
int spellT[STATES][MAXSPELLS]; 
int potionT[STATES][NPOTIONS]; 
int learnSpellT[STATES][MAXSPELLS];
vector<vector<int>> allInv;
int turnNo = 0;
int myInv[4];
double invValue[STATES];
int totTierAbove0[STATES];
const double invWeights[4] = {1.49686278, 2.44179293, 3.44179293, 4.44113734};
bitset<MAXSPELLS> spell1Mask[MAXSPELLS]; // 1s on positions with equal id
bitset<MAXSPELLS> spell0Mask[MAXSPELLS];  // 0s on pos with equal id
bitset<MAXSPELLS> maskRepeatOne;
int distToState[STATES];


int playerScores[103][2];
int potionIds[103][5];
int bonusCnt[2]={4,4};
int bonusVal[2]={3,1};
int bonusThisTurn[5];

string actions[3]= {"REST", "BREW", "CAST"};
const double probREST = 0.1;
const double probBREW = 0.5;
const double probCAST = 0.4;

const double gammaCoeff = 0.96;
double gammaPow[200];

const double bonusSpellActive = 0.0;//1.0/10;
const double penaltyInv = 0.7;

int nPotionsBrewed;
const double bonusFinished=30;
const double bonusLearned = 0.0;
const double tomeIdxPenalty = 0.3;

double maxLearnTime=0;
double maxBeamTime =0;
double maxMoveTime =0;


int sampleType(bool restPossible, bool brewPossible, bool castPossible){
    assert(restPossible || brewPossible || castPossible);
    int startI[3];
    int endI[3];
    int P = 1000;
    int u =0;
    startI[0]=0;
    if(restPossible) u+=P*probREST;
    endI[0]=u;
    startI[1]=u;
    if(brewPossible) u+=P*probBREW;
    endI[1]=u;
    startI[2]=u;
    if(castPossible) u+=P*probCAST;
    endI[2]=u;
    int r = rd(0,endI[2]-1);
    for(int i=0; i<3;i++){
        if(r>= startI[i] && r<endI[i]){
            return i;
        }
    }
    assert(0);
    return 0;
}

int vecToVertex(const vector<int>& v){
    return invIdx[v[0]][v[1]][v[2]][v[3]];
}



vector<vector<int>> generateAllInventoryStates(){
    vector<vector<int>> allInventory;
    for(int a=0; a<11; a++){
        for(int b=0; b<11; b++){
            for(int c=0; c<11; c++){
                for(int d=0; d<11; d++){
                    if(a+b+c+d<11){
                        allInventory.push_back({a,b,c,d});
                        invIdx[a][b][c][d] = SZ(allInventory)-1;
                    }
                }
            }
        }
    }
    return allInventory;
}






struct Potion{
    int id;
    int formula[4];
    int price;

    Potion(){}

    Potion(int _id, int _d0, int _d1, int _d2, int _d3, int _price ){
        id = _id;
        formula[0] = _d0;
        formula[1] = _d1;
        formula[2] = _d2;
        formula[3] = _d3;
        price = _price;
    }

    bool operator < (const Potion& potion2) const
    {
        return price < potion2.price;
    }

};

vector<Potion> clientOrders;

struct Spell{
    int id;
    int formula[4];
    bool castable;
    bool repeatable;
    int tomeIdx;
    int tax;
    int repeats;//how many times was the spell repeated. Default = 1
    bool isLearned;

    Spell(){}

    Spell(int _id, int _d0, int _d1, int _d2, int _d3, bool _castable, bool _repeatable, int _tomeIdx, int _tax, int _repeats , bool _isLearned){
        id = _id;
        formula[0] = _d0;
        formula[1] = _d1;
        formula[2] = _d2;
        formula[3] = _d3;
        castable = _castable;
        repeatable = _repeatable;
        tomeIdx = _tomeIdx;
        tax = _tax;
        repeats = _repeats;
        isLearned = _isLearned;
    }

    int value(){
        int weights[4] = {1,2,3,4};
        int val =0;
        for(int i=0; i<4; i++){
            val += formula[i]*weights[i];
        }
        return val;
    }
};

vector<Spell> inputSpells;
vector<Spell> mySpells;
vector<Spell> tome;


bool spellActive[MAXSPELLS];
bool potionActive[NPOTIONS];
int possibleSpells[MAXSPELLS];
int possiblePotions[NPOTIONS];

struct Path{
    int moveId[PLEN];
    int moveType[PLEN];
    double score;

    Path(){
        for(int i=0; i< PLEN;i++){
            moveId[i]=0;
            moveType[i]=0;
        }
    };

    bool operator < (const Path& path2) const {
        return score > path2.score;
    }

};


struct State{
    int vertex;
    int level;
    bitset<MAXSPELLS> activeSpells;
    bitset<MAXSPELLS> learnedSpells;
    bitset<NPOTIONS> activePotions;
    int parent;
    int prevType;
    int prevId;
    double score;
    double partialScore;
    int brewedSoFar;

    State(){}

    bool operator < (const State& state2) const {
        return partialScore > state2.partialScore;
    }
};

int BEAMW=BEAMMAX;

State stateByLevel[PLEN+1][MAXSTATES];
int currentLen[PLEN+1];
double currentThreshold[PLEN+1];

void considerCandidate(int level){
    int lastIdx = currentLen[level]-1;
    if(stateByLevel[level][lastIdx].partialScore <= currentThreshold[level]){
        currentLen[level]--;
        return;
    }
    if(lastIdx==MAXSTATES-1){
        sort(stateByLevel[level], stateByLevel[level]+MAXSTATES);
        currentLen[level] = BEAMW;
        currentThreshold[level] = stateByLevel[level][BEAMW-1].partialScore;
    }
}

int cntExpands[PLEN+1];

void expandRest(int level, int idx){
    cntExpands[level]++;
    // int nSpells = SZ(mySpells);
    State& s = stateByLevel[level][idx];
    State& ns = stateByLevel[level+1][currentLen[level+1]++];
    ns.vertex = s.vertex;
    distToState[ns.vertex] = min(distToState[ns.vertex], level+1);
    ns.level = level+1;
    ns.activeSpells.set();
    ns.activePotions = s.activePotions;
    ns.learnedSpells = s.learnedSpells;
    ns.brewedSoFar = s.brewedSoFar;
    ns.parent = idx;
    ns.prevType = 0;
    ns.prevId = 0;
    ns.score = s.score;
    if(ns.brewedSoFar >=6) ns.score += bonusFinished* gammaPow[level]; 
    ns.partialScore = ns.score;
    // int cntLearned=(ns.learnedSpells & maskRepeatOne).count();
    // ns.partialScore += cntLearned*bonusLearned * gammaPow[level];
    if (ns.brewedSoFar < 6)
        ns.partialScore += penaltyInv * invValue[ns.vertex] * gammaPow[level];
    considerCandidate(level+1);
}


void expandCast(int level, int idx, int spellIdx){
    cntExpands[level]++;
    // int nSpells = SZ(mySpells);
    State& s = stateByLevel[level][idx];
    State& ns = stateByLevel[level+1][currentLen[level+1]++];
    ns.vertex = spellT[s.vertex][spellIdx];
    distToState[ns.vertex] = min(distToState[ns.vertex], level+1);
    ns.level = s.level+1;
    ns.activeSpells = s.activeSpells & spell0Mask[spellIdx];
    ns.activePotions = s.activePotions;
    ns.learnedSpells = s.learnedSpells;
    ns.parent = idx;
    ns.prevType = 2;
    ns.prevId = spellIdx;
    ns.score = s.score;
    ns.partialScore = ns.score;
    ns.brewedSoFar = s.brewedSoFar;
    // int cntLearned=(ns.learnedSpells & maskRepeatOne).count();
    // ns.partialScore += cntLearned*bonusLearned * gammaPow[level];
    ns.partialScore += penaltyInv * invValue[ns.vertex] * gammaPow[level];
    considerCandidate(level+1);
}


void expandBrew(int level, int idx, int potionIdx){
    cntExpands[level]++;
    // int nSpells = SZ(mySpells);
    State& s = stateByLevel[level][idx];
    State& ns = stateByLevel[level+1][currentLen[level+1]++];
    ns.vertex = potionT[s.vertex][potionIdx];
    ns.level = level+1;
    ns.activeSpells=s.activeSpells;
    ns.activePotions= s.activePotions;
    ns.activePotions[potionIdx]=0;
    ns.learnedSpells = s.learnedSpells;
    ns.parent = idx;
    ns.prevType = 1;
    ns.prevId = potionIdx;
    ns.brewedSoFar = s.brewedSoFar +1;
    ns.score = s.score + (clientOrders[potionIdx].price + bonusThisTurn[potionIdx]) * gammaPow[level];
    if (ns.brewedSoFar>=6)
        ns.score += totTierAbove0[ns.vertex]*gammaPow[level];
    ns.partialScore = ns.score;
    // int cntLearned=(ns.learnedSpells & maskRepeatOne).count();
    // ns.partialScore += cntLearned*bonusLearned * gammaPow[level];
    ns.partialScore += penaltyInv * invValue[ns.vertex] * gammaPow[level];
    considerCandidate(level+1);
}


void expandLearn(int level, int idx, int spellIdx){
    cntExpands[level]++;
    // int nSpells = SZ(mySpells);
    State& s = stateByLevel[level][idx];
    State& ns = stateByLevel[level+1][currentLen[level+1]++];
    ns.vertex = learnSpellT[s.vertex][spellIdx];
    assert(ns.vertex !=-1);
    ns.level = s.level+1;
    ns.activeSpells = s.activeSpells | spell1Mask[spellIdx];
    ns.activePotions = s.activePotions;
    ns.learnedSpells = s.learnedSpells | spell1Mask[spellIdx];
    ns.parent = idx;
    ns.prevType = 3;
    ns.prevId = spellIdx;
    ns.score = s.score;
    ns.partialScore = ns.score;
    ns.brewedSoFar = s.brewedSoFar;
    // int cntLearned=(ns.learnedSpells & maskRepeatOne).count();
    // ns.partialScore += cntLearned*bonusLearned * gammaPow[level];
    ns.partialScore += penaltyInv * invValue[ns.vertex] * gammaPow[level];
    considerCandidate(level+1);
}
 


void expandState(int level, int idx, bool castOnly){
    expandRest( level, idx);
    State& s = stateByLevel[level][idx];
    if(s.brewedSoFar>=6) return;
    bitset<MAXSPELLS> canCast = s.activeSpells & s.learnedSpells;
    int vertex = s.vertex;
    for(int j=0; j<SZ(adjSpell[vertex]);j++){
        if(canCast[adjSpell[vertex][j]]){
            expandCast(level, idx, adjSpell[vertex][j]);
        }
    }
    if(castOnly) return;
    for(int j=0; j<SZ(adjPotion[vertex]);j++){
        if(s.activePotions[adjPotion[vertex][j]]){
            expandBrew(level, idx, adjPotion[vertex][j]);
        }
    }
    if(level>0) return;
    for(int j=0; j<SZ(adjLearnSpell[vertex]);j++){
        if(s.learnedSpells[adjLearnSpell[vertex][j]]==0){
            expandLearn(level, idx, adjLearnSpell[vertex][j]);
        }
    }
}

Path getPathFromState(State s){
    assert(s.level==PLEN);
    Path p;
    p.score = s.partialScore;
    while(s.level>0){
        p.moveType[s.level-1] = s.prevType;
        p.moveId[s.level-1] = s.prevId;
        s = stateByLevel[s.level-1][s.parent];
    }
    return p;
}

int castSpellFromVertex(int v, const Spell& s){
    int inv[4];
    int suma = 0;
    for(int i=0; i<4;i++){
        inv[i] = allInv[v][i]+ s.formula[i];
        if(inv[i]<0) return -1;
        suma += inv[i];
    }
    if(suma>10) return -1;
    return invIdx[inv[0]][inv[1]][inv[2]][inv[3]];
}


int brewPotionFromVertex(int v, const Potion& p){
    int inv[4];
    int suma = 0;
    for(int i=0; i<4;i++){
        inv[i] = allInv[v][i]+ p.formula[i];
        if(inv[i]<0) return -1;
        suma += inv[i];
    }
    assert(suma<=10);
    return invIdx[inv[0]][inv[1]][inv[2]][inv[3]];
}


int learnSpellFromVertex(int v, const Spell& s){ 
    assert(s.tomeIdx>=0);
    int inv0 = allInv[v][0] - s.tomeIdx;
    if(inv0<0) return -1;
    int sum = inv0 + allInv[v][1] + allInv[v][2] + allInv[v][3];
    int bonus = s.tax;
    inv0 += min(bonus, 10-sum);
    return invIdx[inv0][allInv[v][1]][allInv[v][2]][allInv[v][3]];
}


void precomputeTransition(){
    for(int i=0; i<STATES; i++){
        adjSpell[i].clear();
        adjLearnSpell[i].clear();
    }
    for(int i =0; i< STATES;i++){
        for(int j=0; j<SZ(mySpells); j++){
            Spell& s=mySpells[j];
            int u = castSpellFromVertex(i, s);
            spellT[i][j]=u;
            if(u!=-1){
                adjSpell[i].push_back(j);
            }
        }
    }

    assert(SZ(clientOrders)==5);


    for(int i=0; i<STATES; i++){
        adjPotion[i].clear();
    }
    for(int i =0; i< STATES;i++){
        for(int j=0; j<NPOTIONS; j++){
            Potion& p=clientOrders[j];
            int u = brewPotionFromVertex(i, p);
            potionT[i][j]=u;
            if(u!=-1){
                adjPotion[i].push_back(j);
            }
        }
    }
    for(int i=0; i<STATES; i++){
        invValue[i]=0;
        for(int j=0;j<4;j++){
            invValue[i]+=invWeights[j]*allInv[i][j];
        }
    }

    int nSpells = SZ(mySpells);

    maskRepeatOne.reset();
    for(int i=0; i< nSpells; i++){
        if(mySpells[i].repeats == 1) maskRepeatOne[i]=1;
    }
    for(int i =0; i<nSpells;i++){
        spell1Mask[i].reset();
        spell0Mask[i].reset();
        for(int j=0; j<nSpells;j++){
            spell0Mask[i][j] = (mySpells[i].id != mySpells[j].id);
            spell1Mask[i][j] = (mySpells[i].id == mySpells[j].id);
        }
    }

    for(int i =0; i< STATES;i++){
        for(int j=0; j<SZ(mySpells); j++){
            Spell& s=mySpells[j];
            if(s.isLearned == 1 || s.repeats != 1) continue; 
            int u = learnSpellFromVertex(i, s);
            learnSpellT[i][j]=u;
            if(u!=-1){
                adjLearnSpell[i].push_back(j);
            }
        }
    }

    for(int i=0; i<STATES;i++){
        totTierAbove0[i] = allInv[i][1]+allInv[i][2]+allInv[i][3];
    }
}




bool repeatSpell(Spell spell, int n, Spell & result){
    result = spell;
    for(int i=0; i<4; i++){
        result.formula[i]*=n;
        if (abs(result.formula[i])>10) return 0;
    }
    result.repeats = n;
    return 1;
}


struct PotionInfo{
    double score;
    int vertex;
    int idx;

    PotionInfo(double _score, int _vertex, int _idx){
        score = _score;
        vertex = _vertex;
        idx = _idx;
    }

    bool operator < (const PotionInfo& potion2) const {
        return score > potion2.score;
    }
};

void addCandidate(vector<Path>& population, const Path& candidate, int populationSize){
    if(SZ(population)<populationSize){
        population.push_back(candidate);
        return;
    }
    if(candidate.score> population.back().score){
        population.push_back(candidate);
        int pos=populationSize;
        while(1){
            if(pos==0 || population[pos-1].score>=population[pos].score) break;
            swap(population[pos-1], population[pos]);
            pos--;
        }
        population.pop_back();
    }
}

double spellGoodness(int spellId){ //this is a spell id as given in input

    // chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    for(int i=0; i< STATES; i++){
        distToState[i] = INF;
    }
    for(int i=0;i<PLEN+1;i++){
        cntExpands[i]=0;
        currentThreshold[i] = -1.0;
        currentLen[i]=0;
    }
    State& root = stateByLevel[0][0];
    currentLen[0]=1;
    root.vertex = invIdx[0][0][0][0];
    root.level = 0;
    root.learnedSpells.reset();
    root.activeSpells.reset();
    root.activePotions.reset();
    int nSpells = SZ(mySpells);
    for(int i=0; i<nSpells; i++){
        root.learnedSpells[i]=(mySpells[i].isLearned || mySpells[i].id==spellId);
    }
    root.activeSpells = root.learnedSpells;
    root.parent=-1;
    root.prevType=-1;
    root.prevId=-1;
    root.score=0.0;
    root.partialScore=0.0;
    root.brewedSoFar=0;
    BEAMW=BEAMMAX/5;
    for(int level=0; level<PLEN; level++){
        for(int idx=0; idx<currentLen[level];idx++){
            expandState(level, idx, 1);
        }
        sort(stateByLevel[level+1], stateByLevel[level+1]+currentLen[level+1]);
        currentLen[level+1]=min(BEAMW, currentLen[level+1]);
    }
    double sum=0;
    for(int i=0; i<STATES;i++){
        assert(distToState[i]!=0);
        sum+= invValue[i]/(distToState[i]);
    }
    // chrono::steady_clock::time_point end = chrono::steady_clock::now();
    // cerr << "spellGoodness Time = " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
    return sum/STATES;
}

void chooseAndLearn(){
    vector<pair<double,int>> spellScores;
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    for(auto s: tome){
        if(myInv[0]<s.tomeIdx) continue;
        double goodness = spellGoodness(s.id);
        // cerr << "id: " << s.id << " goodness: " << goodness << endl;
        double score = goodness +tomeIdxPenalty*(-s.tomeIdx+0.5*s.tax);
        spellScores.push_back({score ,s.id});
    }
    sort(spellScores.begin(),spellScores.end());
    reverse(spellScores.begin(),spellScores.end());
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    double learnTime = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
    maxLearnTime = max(maxLearnTime, learnTime);
    cerr << "chooseAndLearn Time = " << learnTime << "[ms] max" <<maxLearnTime << endl;
    // debug() << var(tomeAvailableScore);
    // debug()<<var(spellScores);
    
    cout << "LEARN " <<  spellScores[0].second << endl;
}



Path bestBeam(){
    BEAMW=BEAMMAX;
    vector<int> startInv(myInv, myInv + 4);
    int startVertex = vecToVertex(startInv);
    int nSpells = SZ(mySpells);
    for(int i=0;i<PLEN+1;i++) cntExpands[i]=0;
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    for(int i=0; i< PLEN+1; i++){
        currentThreshold[i] = -1.0;
    }
    for(int i=0; i< PLEN+1; i++){
        currentLen[i]=0;
    }
    State& root = stateByLevel[0][0];
    currentLen[0]=1;
    root.learnedSpells.reset();
    root.activeSpells.reset();
    root.activePotions.reset();
    root.vertex = startVertex;
    root.level = 0;
    for(int i=0; i<nSpells; i++){
        root.activeSpells[i]=mySpells[i].castable;
    }
    for(int i=0; i<NPOTIONS; i++){
        root.activePotions[i]=1;
    }
    for(int i=0; i<nSpells; i++){
        root.learnedSpells[i]=mySpells[i].isLearned;
    }
    root.parent=-1;
    root.prevType=-1;
    root.prevId=-1;
    root.score=0.0;
    root.partialScore=0.0;
    root.brewedSoFar=nPotionsBrewed;
    for(int level=0; level<PLEN; level++){
        for(int idx=0; idx<currentLen[level];idx++){
            expandState(level, idx,0);
        }
        sort(stateByLevel[level+1], stateByLevel[level+1]+currentLen[level+1]);
        auto currTime = chrono::steady_clock::now();
        double timeSpent =chrono::duration_cast<chrono::milliseconds>(currTime - begin).count();
        if (timeSpent > 27.0) {
            cerr<<"panic "<<timeSpent<<endl;
            BEAMW = 200;
        }

        currentLen[level+1]=min(BEAMW, currentLen[level+1]);
    }
    // debug()<<var(vector<int>(cntExpands, cntExpands+PLEN));
    // debug()<<var(vector<int>(currentLen, currentLen+PLEN+1));
    // debug()<<var(vector<double>(currentThreshold, currentThreshold+PLEN+1));
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    double beamTime =chrono::duration_cast<chrono::milliseconds>(end - begin).count();
    maxBeamTime = max(maxBeamTime, beamTime);
    cerr << "BeamSearch Time = " << beamTime << "[ms] max" << maxBeamTime<< endl;
    return getPathFromState(stateByLevel[PLEN][0]);
}




void move(){
    mySpells.clear();
    for(auto s: inputSpells){
        mySpells.push_back(s);
        if(s.repeatable){
            Spell rs;
            for(int i=2; i<11; i++){
                if (repeatSpell(s,i,rs)) {
                    mySpells.push_back(rs);
                }
            }
        }
    }
    for(auto s: tome){
        mySpells.push_back(s);
        if(s.repeatable){
            Spell rs;
            for(int i=2; i<11; i++){
                if (repeatSpell(s,i,rs)) {
                    mySpells.push_back(rs);
                }
            }
        }
    }
    debug()<<var(SZ(mySpells));
    precomputeTransition();

    if(turnNo<LEARNTURNS){
        chooseAndLearn();
        return;
    } 



    Path best = bestBeam();
    //Path bestE=bestEvol(startVertex, 25);


    debug()<< var(best.score);
    //debug()<< var(bestE.score);
    int bestMoveType = best.moveType[0];
    int bestMoveId = best.moveId[0];
    if(bestMoveType==0){
        cout << "REST"<< endl;
    }
    if(bestMoveType==1){
        nPotionsBrewed +=1;
        cout << "BREW " << clientOrders[bestMoveId].id << endl;
    }
    if(bestMoveType==2){
        cout << "CAST " << mySpells[bestMoveId].id << " " << mySpells[bestMoveId].repeats << endl;
    }
    if(bestMoveType==3){
        cout << "LEARN " << mySpells[bestMoveId].id  << endl;
    }
}



void updateBonus(){
    if(turnNo==1) return;
    int nPotionsDelivered=0;
    for(int i=0;i<2;i++){
        if(playerScores[turnNo-1][i] != playerScores[turnNo][i]) nPotionsDelivered+=1;
    }
    int potionsDelivered[5]={0,0,0,0,0};
    set<int> potionsThisTurn;
    for(int i=0;i<5;i++){
        potionsThisTurn.insert(potionIds[turnNo][i]);
    }
    int cntDelivered=0;
    for(int i=0; i<5;i++){
        int p = potionIds[turnNo-1][i];
        if(potionsThisTurn.count(p)==0){
            potionsDelivered[i]=1;
            cntDelivered++;    
        } 
    }
    int decrease=1;
    assert(nPotionsDelivered>=cntDelivered);
    if(nPotionsDelivered != cntDelivered){
        decrease=2;
        assert(nPotionsDelivered == 2);
        assert(cntDelivered == 1);
    }
    for(int pos=0; pos<2; pos++){
        if(potionsDelivered[pos]){
            bonusCnt[pos]-=decrease;
        }
    }
    if(bonusCnt[0]<=0){
        bonusCnt[0]=bonusCnt[1];
        bonusCnt[1]=0;
        bonusVal[0] = bonusVal[1];
        bonusVal[1]=0;
    }
    for(int i=0;i<5;i++){
        bonusThisTurn[i]=0;
        if(i<2){
            if(bonusCnt[i]>0){
                bonusThisTurn[i] = bonusVal[i];
            }
        }
    }
    // debug()<< var(bonusThisTurn[0]) << var(bonusThisTurn[1]);
}


void globalInit(){
    gammaPow[0]=1.0;
    for(int i=1; i<200; i++){
        gammaPow[i] = gammaPow[i-1]*gammaCoeff;
    }
    allInv = generateAllInventoryStates();

    nPotionsBrewed=0;
}



int main()
{
    globalInit();

    
    // vector<Spell> oppSpells;

    // game loop
    while (1) {
        turnNo += 1;

        clientOrders.clear();
        inputSpells.clear();
        tome.clear();

        int actionCount; // the number of spells and recipes in play
        cin >> actionCount; cin.ignore();
        int potCnt=0;
        for (int i = 0; i < actionCount; i++) {
            int actionId; // the unique ID of this spell or recipe
            string actionType; // in the first league: BREW; later: CAST, OPPONENT_CAST, LEARN, BREW
            int delta0; // tier-0 ingredient change
            int delta1; // tier-1 ingredient change
            int delta2; // tier-2 ingredient change
            int delta3; // tier-3 ingredient change
            int price; // the price in rupees if this is a potion
            int tomeIndex; // in the first two leagues: always 0; later: the index in the tome if this is a tome spell, equal to the read-ahead tax
            int taxCount; // in the first two leagues: always 0; later: the amount of taxed tier-0 ingredients you gain from learning this spell
            bool castable; // in the first league: always 0; later: 1 if this is a castable player spell
            bool repeatable; // for the first two leagues: always 0; later: 1 if this is a repeatable player spell
            cin >> actionId >> actionType >> delta0 >> delta1 >> delta2 >> delta3 >> price >> tomeIndex >> taxCount >> castable >> repeatable; cin.ignore();




            if( actionType=="BREW" ){
                clientOrders.push_back(Potion(actionId, delta0, delta1, delta2, delta3, price));
                potionIds[turnNo][potCnt++] = actionId;
            }
            if(actionType == "CAST"){
                inputSpells.push_back(Spell(actionId, delta0, delta1, delta2, delta3, castable, repeatable, tomeIndex, taxCount, 1, 1 )); // isLearned=1
            }

            if(actionType == "LEARN"){
                tome.push_back(Spell(actionId, delta0, delta1, delta2, delta3, castable, repeatable, tomeIndex, taxCount, 1, 0));
            }

        }


        for (int i = 0; i < 2; i++) {
            int inv0; // tier-0 ingredients in inventory
            int inv1;
            int inv2;
            int inv3;
            int score; // amount of rupees
            cin >> inv0 >> inv1 >> inv2 >> inv3 >> score; cin.ignore();
            playerScores[turnNo][i] = score;
            if(i==0){
                myInv[0]=inv0;
                myInv[1]=inv1;
                myInv[2]=inv2;
                myInv[3]=inv3;
            }

        }


        updateBonus();

        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        move();
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        double moveTime = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        maxMoveTime = max(maxMoveTime, moveTime);
        cerr << "MOVE Time = " << moveTime << "[ms]  MAX "<< maxMoveTime << endl;
    }
}